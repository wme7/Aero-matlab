/*
 * clapack_ssyevx.c
 *
 *This program is a C interface to ssyevx.
 *
 * Written by Remi Delmas.
 *
 */

#include "mex.h"

#define min(a,b) ((a)>(b) ? (b) : (a))
#define max(a,b) ((a)<(b) ? (b) : (a))

#if defined(ADD_)
	#define f77_ssyevx ssyevx_
#elif defined(UPCASE)
	#define f77_ssyevx SSYEVX
#else
	#define f77_ssyevx ssyevx
#endif

char *select_name;

void mexFunction(int nlhs_m, mxArray *plhs_m[], int nrhs_m, const mxArray *prhs_m[])
{

	int i;

	char * jobz;
	char * range;
	char * uplo;
	int * n;
	float * a;
	int * lda;
	float * vl;
	float * vu;
	int * il;
	int * iu;
	float * abstol;
	int * m;
	float * w;
	float * z;
	int * ldz;
	float * work;
	int * lwork;
	int * iwork;
	int * ifail;
	int * info;

	plhs_m[3]=mxDuplicateArray(prhs_m[3]);
	n=(int *)mxGetPr(plhs_m[3]);

	plhs_m[5]=mxDuplicateArray(prhs_m[5]);
	lda=(int *)mxGetPr(plhs_m[5]);

	plhs_m[8]=mxDuplicateArray(prhs_m[8]);
	il=(int *)mxGetPr(plhs_m[8]);

	plhs_m[9]=mxDuplicateArray(prhs_m[9]);
	iu=(int *)mxGetPr(plhs_m[9]);

	plhs_m[11]=mxDuplicateArray(prhs_m[11]);
	m=(int *)mxGetPr(plhs_m[11]);

	plhs_m[14]=mxDuplicateArray(prhs_m[14]);
	ldz=(int *)mxGetPr(plhs_m[14]);

	plhs_m[16]=mxDuplicateArray(prhs_m[16]);
	lwork=(int *)mxGetPr(plhs_m[16]);

	plhs_m[19]=mxDuplicateArray(prhs_m[19]);
	info=(int *)mxGetPr(plhs_m[19]);

	plhs_m[0] = mxDuplicateArray(prhs_m[0]);
	jobz = (char*) mxArrayToString(plhs_m[0]);

	plhs_m[1] = mxDuplicateArray(prhs_m[1]);
	range = (char*) mxArrayToString(plhs_m[1]);

	plhs_m[2] = mxDuplicateArray(prhs_m[2]);
	uplo = (char*) mxArrayToString(plhs_m[2]);



	plhs_m[4]=mxDuplicateArray(prhs_m[4]);
	a = (float*) mxGetPr(plhs_m[4]);



	plhs_m[6]=mxDuplicateArray(prhs_m[6]);
	vl = (float*) mxGetPr(plhs_m[6]);

	plhs_m[7]=mxDuplicateArray(prhs_m[7]);
	vu = (float*) mxGetPr(plhs_m[7]);





	plhs_m[10]=mxDuplicateArray(prhs_m[10]);
	abstol = (float*) mxGetPr(plhs_m[10]);



	plhs_m[12]=mxDuplicateArray(prhs_m[12]);
	w = (float*) mxGetPr(plhs_m[12]);

	plhs_m[13]=mxDuplicateArray(prhs_m[13]);
	z = (float*) mxGetPr(plhs_m[13]);



	plhs_m[15]=mxDuplicateArray(prhs_m[15]);
	work = (float*) mxGetPr(plhs_m[15]);



	plhs_m[17]=mxDuplicateArray(prhs_m[17]);
	iwork = (int *) mxGetPr(plhs_m[17]);

	plhs_m[18]=mxDuplicateArray(prhs_m[18]);
	ifail = (int *) mxGetPr(plhs_m[18]);



#ifdef F77_INT
	F77_INT* F77_n = n , F77_lda = lda , F77_il = il , F77_iu = iu , F77_m = m , F77_ldz = ldz , F77_lwork = lwork , F77_iwork = iwork , F77_ifail = ifail , F77_info = info ;
#else
	#define F77_n n 
	#define F77_lda lda 
	#define F77_il il 
	#define F77_iu iu 
	#define F77_m m 
	#define F77_ldz ldz 
	#define F77_lwork lwork 
	#define F77_iwork iwork 
	#define F77_ifail ifail 
	#define F77_info info 
#endif

#ifdef F77_CHAR
	F77_CHAR* F77_jobz = jobz , F77_range = range , F77_uplo = uplo ;
#else
	#define F77_jobz jobz 
	#define F77_range range 
	#define F77_uplo uplo 
#endif

	f77_ssyevx(F77_jobz, F77_range, F77_uplo, F77_n, a, F77_lda, vl, vu, F77_il, F77_iu, abstol, F77_m, w, z, F77_ldz, work, F77_lwork, F77_iwork, F77_ifail, F77_info);

	return;
}
