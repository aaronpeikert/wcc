/* windcrosscum.c
 *
 * Batched windowed cross-correlation via cumulative (prefix) sums.
 *
 * windcrosscum_batch computes the WCC grid for a list of series pairs.
 * Per-series prefix sums (of the centered series and its squares) are
 * computed once per series and shared across all pairs; per pair, only
 * the lagged-product prefix sums are built.  Each cell of the
 * (window position x lag) grid then costs O(1), so total cost is
 * O(nPairs * n * nLags) independent of the window size.  The pair loop
 * is parallelized with OpenMP when available.
 *
 * Geometry matches the pure-R reference path in wccCalc.R:
 *   tStart = wMax + tMax
 *   row k (1-based) uses the base window [tStart + (k-1)*wInc - wMax + 1,
 *                                         tStart + (k-1)*wInc]
 *   column centerCol + i is cor(x[window], y[window - i*tInc])
 *   column centerCol - i is cor(y[window], x[window - i*tInc])
 *
 * The single-dyad entry point windcrosscum is a shim over the batch
 * routine with a single (1, 1) pair.
 *
 * Missing data is not supported; the R wrappers reject NA input.
 */

#include <R.h>
#include <Rinternals.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Compute one output column (signed lag) for one pair.
 * x, y: centered series (length n).
 * csA, csA2: prefix sums (length n+1, leading zero) of the lagged series
 *            and its squares; csB, csB2: same for the base series.
 * base = x and lagged = y for columns right of center, swapped left of it.
 * csP: scratch buffer (length n+1) for the lagged-product prefix sums. */
static void corrColumn(const double *base, const double *lagged,
                       const double *csB, const double *csB2,
                       const double *csA, const double *csA2,
                       double *csP, R_xlen_t n, R_xlen_t tStart,
                       R_xlen_t windowSize, R_xlen_t windowIncrement,
                       R_xlen_t L, R_xlen_t nRow, double *outCol) {
    double W = (double) windowSize;
    csP[0] = 0.0;
    for (R_xlen_t t = 0; t < n; t++) {
        csP[t + 1] = csP[t] + (t >= L ? base[t] * lagged[t - L] : 0.0);
    }
    for (R_xlen_t k = 0; k < nRow; k++) {
        R_xlen_t hi = tStart + k * windowIncrement;  /* prefix index of window end */
        R_xlen_t lo = hi - windowSize;
        double sxy = csP[hi] - csP[lo];
        double sb  = csB[hi] - csB[lo];
        double qb  = csB2[hi] - csB2[lo];
        double sa  = csA[hi - L] - csA[lo - L];
        double qa  = csA2[hi - L] - csA2[lo - L];
        double varb = W * qb - sb * sb;
        double vara = W * qa - sa * sa;
        double num = W * sxy - sb * sa;
        double den = varb * vara;
        outCol[k] = (den > 0.0) ? num / sqrt(den) : NA_REAL;
    }
}

/* seriesArray1/2: D x n numeric matrices (rows = series), pairs: P x 2
 * integer matrix of 1-based row indices into array1 / array2.
 * Returns an nRow x nCol x P numeric array. */
SEXP windcrosscum_batch(SEXP seriesArray1, SEXP seriesArray2, SEXP pairs,
                        SEXP wMax, SEXP tMax, SEXP wInc, SEXP tInc) {
    if (!isReal(seriesArray1) || !isReal(seriesArray2) || !isMatrix(seriesArray1) || !isMatrix(seriesArray2)) {
        error("seriesArray1 and seriesArray2 must be numeric matrices.");
    }
    if (!isInteger(pairs) || !isMatrix(pairs) || ncols(pairs) != 2) {
        error("pairs must be an integer matrix with two columns.");
    }
    if (!isReal(wMax) || !isReal(tMax) || !isReal(wInc) || !isReal(tInc)) {
        error("wMax, tMax, wInc, and tInc must be numeric.");
    }

    R_xlen_t D1 = nrows(seriesArray1);
    R_xlen_t D2 = nrows(seriesArray2);
    R_xlen_t n = ncols(seriesArray1);
    if ((R_xlen_t) ncols(seriesArray2) != n) {
        error("seriesArray1 and seriesArray2 must have the same number of columns.");
    }
    R_xlen_t nPairs = nrows(pairs);
    const int *pairIdx = INTEGER(pairs);
    for (R_xlen_t p = 0; p < nPairs; p++) {
        if (pairIdx[p] < 1 || pairIdx[p] > D1 || pairIdx[p + nPairs] < 1 || pairIdx[p + nPairs] > D2) {
            error("pairs contains an index outside the series arrays.");
        }
    }

    R_xlen_t windowSize = (R_xlen_t) REAL(wMax)[0];
    R_xlen_t windowIncrement = (R_xlen_t) REAL(wInc)[0];
    R_xlen_t maxLag = (R_xlen_t) REAL(tMax)[0];
    R_xlen_t lagIncrement = (R_xlen_t) REAL(tInc)[0];

    R_xlen_t tStart = windowSize + maxLag;
    R_xlen_t nRow = (n - tStart) / windowIncrement;
    R_xlen_t nLagSteps = maxLag / lagIncrement;
    R_xlen_t nCol = 2 * nLagSteps + 1;
    R_xlen_t centerCol0 = nLagSteps;  /* 0-based center column */
    if (nRow < 1) {
        error("Bad choice of parameters: the result matrix has %ld rows.", (long) nRow);
    }

    double resBytes = (double) nRow * (double) nCol * (double) nPairs * 8.0;
    if ((double) nRow * (double) nCol * (double) nPairs > 2147483647.0) {
        error("Result array would need %.1f GB (%ld x %ld x %ld doubles); "
              "split the pairs into smaller batches.",
              resBytes / 1073741824.0, (long) nRow, (long) nCol, (long) nPairs);
    }

    SEXP dims = PROTECT(allocVector(INTSXP, 3));
    INTEGER(dims)[0] = (int) nRow;
    INTEGER(dims)[1] = (int) nCol;
    INTEGER(dims)[2] = (int) nPairs;
    SEXP corResult = PROTECT(allocArray(REALSXP, dims));
    double *out = REAL(corResult);

    /* Stage 1: per-series centered copies and prefix sums, computed once.
     * Input matrices are D x n column-major, so series are strided;
     * transpose into contiguous per-series rows of length n, with the
     * prefix arrays as contiguous rows of length n + 1. */
    double *x1   = (double *) R_alloc(D1 * n, sizeof(double));
    double *x2   = (double *) R_alloc(D2 * n, sizeof(double));
    double *cs1  = (double *) R_alloc(D1 * (n + 1), sizeof(double));
    double *cs1q = (double *) R_alloc(D1 * (n + 1), sizeof(double));
    double *cs2  = (double *) R_alloc(D2 * (n + 1), sizeof(double));
    double *cs2q = (double *) R_alloc(D2 * (n + 1), sizeof(double));

    const double *in1 = REAL(seriesArray1);
    const double *in2 = REAL(seriesArray2);
    for (R_xlen_t d = 0; d < D1; d++) {
        double *xs = x1 + d * n;
        double m = 0.0;
        for (R_xlen_t t = 0; t < n; t++) m += in1[d + D1 * t];
        m /= (double) n;
        double *cs = cs1 + d * (n + 1);
        double *cq = cs1q + d * (n + 1);
        cs[0] = cq[0] = 0.0;
        for (R_xlen_t t = 0; t < n; t++) {
            double v = in1[d + D1 * t] - m;
            xs[t] = v;
            cs[t + 1] = cs[t] + v;
            cq[t + 1] = cq[t] + v * v;
        }
    }
    for (R_xlen_t d = 0; d < D2; d++) {
        double *ys = x2 + d * n;
        double m = 0.0;
        for (R_xlen_t t = 0; t < n; t++) m += in2[d + D2 * t];
        m /= (double) n;
        double *cs = cs2 + d * (n + 1);
        double *cq = cs2q + d * (n + 1);
        cs[0] = cq[0] = 0.0;
        for (R_xlen_t t = 0; t < n; t++) {
            double v = in2[d + D2 * t] - m;
            ys[t] = v;
            cs[t + 1] = cs[t] + v;
            cq[t + 1] = cq[t] + v * v;
        }
    }

#ifdef _OPENMP
    int maxThreads = omp_get_max_threads();
#else
    int maxThreads = 1;
#endif
    /* One product-prefix scratch buffer per thread. */
    double *csPAll = (double *) R_alloc((R_xlen_t) maxThreads * (n + 1), sizeof(double));

    /* Stage 2: independent pairs. */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (R_xlen_t p = 0; p < nPairs; p++) {
#ifdef _OPENMP
        double *csP = csPAll + (R_xlen_t) omp_get_thread_num() * (n + 1);
#else
        double *csP = csPAll;
#endif
        R_xlen_t d1 = pairIdx[p] - 1;
        R_xlen_t d2 = pairIdx[p + nPairs] - 1;
        const double *x  = x1 + d1 * n;
        const double *y  = x2 + d2 * n;
        const double *csX  = cs1 + d1 * (n + 1);
        const double *csX2 = cs1q + d1 * (n + 1);
        const double *csY  = cs2 + d2 * (n + 1);
        const double *csY2 = cs2q + d2 * (n + 1);
        double *slab = out + p * nRow * nCol;

        /* Center column: lag 0. */
        corrColumn(x, y, csX, csX2, csY, csY2, csP, n, tStart,
                   windowSize, windowIncrement, 0, nRow,
                   slab + nRow * centerCol0);
        for (R_xlen_t i = 1; i <= nLagSteps; i++) {
            R_xlen_t L = i * lagIncrement;
            /* Column centerCol + i: cor(x[window], y[window - L]). */
            corrColumn(x, y, csX, csX2, csY, csY2, csP, n, tStart,
                       windowSize, windowIncrement, L, nRow,
                       slab + nRow * (centerCol0 + i));
            /* Column centerCol - i: cor(y[window], x[window - L]). */
            corrColumn(y, x, csY, csY2, csX, csX2, csP, n, tStart,
                       windowSize, windowIncrement, L, nRow,
                       slab + nRow * (centerCol0 - i));
        }
    }

    UNPROTECT(2);
    return corResult;
}

/* Single-dyad entry point: shim over the batch routine with one (1, 1) pair. */
SEXP windcrosscum(SEXP inSeries1, SEXP inSeries2, SEXP wMax, SEXP tMax, SEXP wInc, SEXP tInc) {
    if (!isReal(inSeries1) || !isReal(inSeries2)) {
        error("inSeries1 and inSeries2 must be numeric.");
    }
    R_xlen_t n = XLENGTH(inSeries1);
    if (XLENGTH(inSeries2) != n) {
        error("inSeries1 and inSeries2 must have equal length.");
    }

    SEXP m1 = PROTECT(allocMatrix(REALSXP, 1, (int) n));
    SEXP m2 = PROTECT(allocMatrix(REALSXP, 1, (int) n));
    memcpy(REAL(m1), REAL(inSeries1), n * sizeof(double));
    memcpy(REAL(m2), REAL(inSeries2), n * sizeof(double));
    SEXP onePair = PROTECT(allocMatrix(INTSXP, 1, 2));
    INTEGER(onePair)[0] = 1;
    INTEGER(onePair)[1] = 1;

    SEXP res3d = PROTECT(windcrosscum_batch(m1, m2, onePair, wMax, tMax, wInc, tInc));

    int *d = INTEGER(getAttrib(res3d, R_DimSymbol));
    SEXP corResult = PROTECT(allocMatrix(REALSXP, d[0], d[1]));
    memcpy(REAL(corResult), REAL(res3d), (size_t) d[0] * d[1] * sizeof(double));

    UNPROTECT(5);
    return corResult;
}
