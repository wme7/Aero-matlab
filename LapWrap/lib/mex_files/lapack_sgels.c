/*
 * clapack_sgels.c
 *
 *This program is a C interface to sgels.
 *
 * Written by Remi Delmas.
 *
 */

#include "mex.h"

#define min(a,b) ((a)>(b) ? (b) : (a))
#define max(a,b) ((a)<(b) ? (b) : (a))

#if defined(ADD_)
	#define f77_sgels sgels_
#elif defined(UPCASE)
	#define f77_sgels SGELS
#else
	#define f77_sgels sgels
#endif

char *select_name;

void mexFunction(int nlhs_m, mxArray *plhs_m[], int nrhs_m, const mxArray *prhs_m[])
{

	int i;

	char * trans;
	int * m;
	int * n;
	int * nrhs;
	float * a;
	int * lda;
	float * b;
	int * ldb;
	float * work;
	int * lwork;
	int * info;

	plhs_m[1]=mxDuplicateArray(prhs_m[1]);
	m=(int *)mxGetPr(plhs_m[1]);

	plhs_m[2]=mxDuplicateArray(prhs_m[2]);
	n=(int *)mxGetPr(plhs_m[2]);

	plhs_m[3]=mxDuplicateArray(prhs_m[3]);
	nrhs=(int *)mxGetPr(plhs_m[3]);

	plhs_m[5]=mxDuplicateArray(prhs_m[5]);
	lda=(int *)mxGetPr(plhs_m[5]);

	plhs_m[7]=mxDuplicateArray(prhs_m[7]);
	ldb=(int *)mxGetPr(plhs_m[7]);

	plhs_m[9]=mxDuplicateArray(prhs_m[9]);
	lwork=(int *)mxGetPr(plhs_m[9]);

	plhs_m[10]=mxDuplicateArray(prhs_m[10]);
	info=(int *)mxGetPr(plhs_m[10]);

	plhs_m[0] = mxDuplicateArray(prhs_m[0]);
	trans = (char*) mxArrayToString(plhs_m[0]);







	plhs_m[4]=mxDuplicateArray(prhs_m[4]);
	a = (float*) mxGetPr(plhs_m[4]);



	plhs_m[6]=mxDuplicateArray(prhs_m[6]);
	b = (float*) mxGetPr(plhs_m[6]);



	plhs_m[8]=mxDuplicateArray(prhs_m[8]);
	work = (float*) mxGetPr(plhs_m[8]);





#ifdef F77_INT
	F77_INT* F77_m = m , F77_n = n , F77_nrhs = nrhs , F77_lda = lda , F77_ldb = ldb , F77_lwork = lwork , F77_info = info ;
#else
	#define F77_m m 
	#define F77_n n 
	#define F77_nrhs nrhs 
	#define F77_lda lda 
	#define F77_ldb ldb 
	#define F77_lwork lwork 
	#define F77_info info 
#endif

#ifdef F77_CHAR
	F77_CHAR* F77_trans = trans ;
#else
	#define F77_trans trans 
#endif

	f77_sgels(F77_trans, F77_m, F77_n, F77_nrhs, a, F77_lda, b, F77_ldb, work, F77_lwork, F77_info);

	return;
}
