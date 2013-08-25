/*
 * clapack_ssyevd.c
 *
 *This program is a C interface to ssyevd.
 *
 * Written by Remi Delmas.
 *
 */

#include "mex.h"

#define min(a,b) ((a)>(b) ? (b) : (a))
#define max(a,b) ((a)<(b) ? (b) : (a))

#if defined(ADD_)
	#define f77_ssyevd ssyevd_
#elif defined(UPCASE)
	#define f77_ssyevd SSYEVD
#else
	#define f77_ssyevd ssyevd
#endif

char *select_name;

void mexFunction(int nlhs_m, mxArray *plhs_m[], int nrhs_m, const mxArray *prhs_m[])
{

	int i;

	char * jobz;
	char * uplo;
	int * n;
	float * a;
	int * lda;
	float * w;
	float * work;
	int * lwork;
	int * iwork;
	int * liwork;
	int * info;

	plhs_m[2]=mxDuplicateArray(prhs_m[2]);
	n=(int *)mxGetPr(plhs_m[2]);

	plhs_m[4]=mxDuplicateArray(prhs_m[4]);
	lda=(int *)mxGetPr(plhs_m[4]);

	plhs_m[7]=mxDuplicateArray(prhs_m[7]);
	lwork=(int *)mxGetPr(plhs_m[7]);

	plhs_m[9]=mxDuplicateArray(prhs_m[9]);
	liwork=(int *)mxGetPr(plhs_m[9]);

	plhs_m[10]=mxDuplicateArray(prhs_m[10]);
	info=(int *)mxGetPr(plhs_m[10]);

	plhs_m[0] = mxDuplicateArray(prhs_m[0]);
	jobz = (char*) mxArrayToString(plhs_m[0]);

	plhs_m[1] = mxDuplicateArray(prhs_m[1]);
	uplo = (char*) mxArrayToString(plhs_m[1]);



	plhs_m[3]=mxDuplicateArray(prhs_m[3]);
	a = (float*) mxGetPr(plhs_m[3]);



	plhs_m[5]=mxDuplicateArray(prhs_m[5]);
	w = (float*) mxGetPr(plhs_m[5]);

	plhs_m[6]=mxDuplicateArray(prhs_m[6]);
	work = (float*) mxGetPr(plhs_m[6]);



	plhs_m[8]=mxDuplicateArray(prhs_m[8]);
	iwork = (int *) mxGetPr(plhs_m[8]);





#ifdef F77_INT
	F77_INT* F77_n = n , F77_lda = lda , F77_lwork = lwork , F77_iwork = iwork , F77_liwork = liwork , F77_info = info ;
#else
	#define F77_n n 
	#define F77_lda lda 
	#define F77_lwork lwork 
	#define F77_iwork iwork 
	#define F77_liwork liwork 
	#define F77_info info 
#endif

#ifdef F77_CHAR
	F77_CHAR* F77_jobz = jobz , F77_uplo = uplo ;
#else
	#define F77_jobz jobz 
	#define F77_uplo uplo 
#endif

	f77_ssyevd(F77_jobz, F77_uplo, F77_n, a, F77_lda, w, work, F77_lwork, F77_iwork, F77_liwork, F77_info);

	return;
}
