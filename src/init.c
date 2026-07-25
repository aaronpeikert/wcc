#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP windcross(SEXP inSeries1, SEXP inSeries2, SEXP wMax, SEXP tMax, SEXP wInc, SEXP tInc);
extern SEXP windcrosscum(SEXP inSeries1, SEXP inSeries2, SEXP wMax, SEXP tMax, SEXP wInc, SEXP tInc);
extern SEXP windcrosscum_batch(SEXP seriesArray1, SEXP seriesArray2, SEXP pairs, SEXP wMax, SEXP tMax, SEXP wInc, SEXP tInc);

static const R_CallMethodDef CallEntries[] = {
    {"windcross", (DL_FUNC) &windcross, 6},
    {"windcrosscum", (DL_FUNC) &windcrosscum, 6},
    {"windcrosscum_batch", (DL_FUNC) &windcrosscum_batch, 7},
    {NULL, NULL, 0}
};

void R_init_wcc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
