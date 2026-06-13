#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "windcross.h"

SEXP windcross(SEXP inSeries1, SEXP inSeries2, SEXP wMax, SEXP tMax, SEXP wInc, SEXP tInc) {
	int     debug = 0;	/* debug infomation */
	int     calcCov = 0;	/* 0=calculate correlations, 1=calculate covariances */
    // Check types
    if (!isReal(inSeries1) || !isReal(inSeries2) || !isReal(wMax) || !isReal(tMax) || !isReal(wInc) || !isReal(tInc)) {
        error("All arguments must be numeric.");
    }

    R_xlen_t maxRowLen = XLENGTH(inSeries1);
    R_xlen_t maxRowLen2 = XLENGTH(inSeries2);

    // Assume wMax, tMax, wInc, and tInc are scalar numeric parameters
    int windowSize = REAL(wMax)[0];
    int windowIncrement = REAL(wInc)[0];
    int maxLag = REAL(tMax)[0];
    int lagIncrement = REAL(tInc)[0];

    // Create an nRow-by-nCol numeric matrix
    int tStart = windowSize + maxLag - 1;
	int nRow = floor((maxRowLen - maxLag - windowSize) / windowIncrement);
	int nCol = (floor(maxLag / lagIncrement) * 2) + 1;
    int centerCol = (floor(maxLag / lagIncrement)) + 1;
    SEXP corResult = PROTECT(allocMatrix(REALSXP, nRow, nCol));
    double *tempResult = REAL(corResult);
    double *person1InputVector = REAL(inSeries1);
    double *person2InputVector = REAL(inSeries2);
//    double *r = REAL(corResult);

    // Fill matrix in column-major order, matching R
    int tRow = 0;
    int tCol = 0;
    int person1Index = 0;
    int person2Index = 0;
    double r = 0;
    /*
     * second person's window selecting so as to get negative lag for person1
     */
    tRow = 0;
    /* time t1 */
    for (person2Index = maxLag; person2Index <= (maxRowLen - windowSize - 1); person2Index += windowIncrement) {
    	tCol = 0;
    	/* first window selected */
    	/* negative lags */
    	for (person1Index = (person2Index - maxLag); person1Index <= person2Index; person1Index += lagIncrement) {
    		r = corr(person1InputVector, person2InputVector, person1Index, person2Index, windowSize, debug, calcCov);

    		/* write from left to right in a row */
    		tempResult[tRow + nRow * tCol] = r;

    		/* positive lags */
    		/* write from right to left in a row */
    		r = corr(person2InputVector, person1InputVector, person1Index, person2Index, windowSize, debug, calcCov);

    		tempResult[tRow + (nRow * ((nCol - 1) - tCol))] = r;
    		tCol++;
    	}
    	tRow++;
    } 



    UNPROTECT(1);
    return corResult;
}
