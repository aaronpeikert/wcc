
/* windCrossCorr.h */
/* head file for windcross */
/* Date: Feb 20/2001 */
/* Author: Steven M Boker, Minquan Xu */

//#define MAXFILENAMELEN 256 
//#define MAXOBS 100000
//#define NROW 50000
//#define NCOL 1000
//#define MAXDATALEN    20 	/* number of digits + sign + '.' */
#define MISSINGVALUE -999999.0  /* value as missing */
//#define DATALEN  (long int)(100000/(2*sizeof(double)))	
/* 2 columns --must cast this way, to match the data type for division */

/* for store data read from a file */
//struct twoColData {
//	double  data[DATALEN][2];  /* for direct access data */
//	long int fileRowLength;
//	struct twoColData *next;
//};

/* hold parameters */
/*struct cmdLineParams {
	char   *inputFile;
	char   *outputFile;
	int     inputFileFlag;
	int     outputFileFlag;
	int     maxRowLenFlag;
	int     debug;
	int     calcCov;
	long int windowSize;
	long int windowIncrement;
	long int maxLag;
	long int lagIncrement;
	long int maxRowLen;
};
*/

/* protype for corr */
double  corr(double *person1InputVector, double *person2InputVector, long int w1, long int w2, long int size, int debug, int calcCov);

/* ptotype for cmdLindParse */
// int cmdLineParse(int argc, char *argv[], struct cmdLineParams *cmdRet);

/* readData, read from fp */
//struct twoColData *readData(FILE * fp);

/* input number length check during cmd line parsing*/
//void    numLenCheck(char *argv[]);
