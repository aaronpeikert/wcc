

/* file name ---corr.c */
/* date---2/20/2001 */
/* author---Steven M. Boker, Minquan Xu */
/* calculate correlation beteen two variables */
/**input:
*  1). twoColData* pdata---data link list, this list include all the data read from a two
*  column ASCII file.
*  2). w1---the start position of a correlation window for person 1.
*  3). w2---the start position pf a correlation window for person 2.
*  4). size---window size for cross correlation.
*  5). debug---if dubug=1, print out debug information.
*  6). fout---output file
***how it works:
*  Based on w1, w2 and size, we know which window is selected for cross correlation.
*  fetch each value in the window from data, and using the formula:
*  r =(n*sum(xy)-sum(x)*sum(y))/sqrt((n*sum(x^2)-sum(x)^2)*(n*sum(y^2)-sum(y)^2))
*  to compute r.
*  two data points are used, one for the data for person 1, the oher for person2.
*  data[][0] for first person, data[][1] for second person.
*  to fetch the a specific data element, first find which memory block it is in,
*  then find the offset in that block. See readData.c for memory allocation scheme.
***output: a single pearson correlation.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "windcross.h"

/* calculate correlation for two windows of data, which start at w1 abd w2 */
double
corr(double *person1InputVector, double *person2InputVector, long int w1, long int w2, long int size, int debug, int calcCov)
{

	long int nonMissingPars = 0;	/* non missing pairs of data */
	long int i, j;		        /* loop indices */
	double  r = 0.0;	/* corr coeff */
	double  sumxy = 0.0;	/* sum of x*y */
	double  sumx = 0.0;	/* sum of x */
	double  sumy = 0.0;	/* sum of y */
	double  sumx2 = 0.0;	/* sum of x square */
	double  sumy2 = 0.0;	/* sum of y square */
	double  rdeno = 0.0;	/* denomitor of r */

/*	if (debug > 2) {
		printf("window1index=%ld, window2index=%ld, size=%ld\n", w1, w2, size);
		if (fout)
			fprintf(fout, "window1index=%ld, window2index=%ld, size=%ld\n", \
				w1, w2, size);
	}
*/


	/* use i control the loop */
	for (i = w1, j = w2; i < (w1 + size); i++, j++) {

		if (person1InputVector[i] != MISSINGVALUE && person2InputVector[j] != MISSINGVALUE) {
			sumxy += person1InputVector[i] * person2InputVector[j];

/*			if (debug > 3) {
				printf("sumxy = %.5f\n", sumxy);
				if (fout)
					fprintf(fout, "sumxy = %.5f\n", sumxy);
			}
*/

			sumx += person1InputVector[i];
			/* temp test--> real data used to compute r */


/*			if (debug > 3) {
				printf("sumx = %.5f\n", sumx);
				if (fout) {
					fprintf(fout, "sumx = %.5f\n", sumx);
				}
			}
*/

			sumy += person2InputVector[j];
                         

/*			if (debug > 3) {
				printf("sumy = %.5f\n", sumy);
				if (fout) {
					fprintf(fout, "sumy = %.5f\n", sumy);
				}
			}
*/


			sumx2 += person1InputVector[i] * person1InputVector[i];

/*			if (debug > 3) {
				printf("sum x square = %.5f\n", sumx2);
				if (fout)
					fprintf(fout, "sumx2 = %.5f\n", sumx2);
			}
*/

			sumy2 += person2InputVector[j] * person2InputVector[j];

/*			if (debug > 3) {
				printf("sum y square = %.5f\n", sumy2);
				if (fout)
					fprintf(fout, "sum y square = %.5f\n", sumy2);
			}
*/

			nonMissingPars++;
		}		/* if !NA */

	}			/* for i=w1... */


	if (nonMissingPars > 5) {

		rdeno = ((nonMissingPars * sumx2 - sumx * sumx) * (nonMissingPars * sumy2 - sumy * sumy));
		if (rdeno == 0.0)
			r = MISSINGVALUE;

		else if (rdeno != 0.0) {
		    if (calcCov==1) {
                r = (sumxy - sumx * sumy / nonMissingPars) / (nonMissingPars - 1);
		    }
		    else {
			    r = (nonMissingPars * sumxy - sumx * sumy) / sqrt((nonMissingPars * sumx2 - sumx * sumx) * (nonMissingPars * sumy2 - sumy * sumy));
		    }

/*			if (debug > 3) {
				printf("r = %.5f\n", r);
				if (fout)
					fprintf(fout, "r= %.5f\n", r);
			}
*/
		}		/* else rdeno!=0 */
	}
	else	/* too many missing values */
		r = MISSINGVALUE;

	return r;
}
