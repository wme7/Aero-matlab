/*
 * clapack_cgelsy.c
 *
 *This program is a C interface to cgelsy.
 *
 * Written by Remi Delmas.
 *
 */

#include "mex.h"

#define min(a,b) ((a)>(b) ? (b) : (a))
#define max(a,b) ((a)<(b) ? (b) : (a))

#if defined(ADD_)
	#define f77_cgelsy cgelsy_
#elif defined(UPCASE)
	#define f77_cgelsy CGELSY
#else
	#define f77_cgelsy cgelsy
#endif

char *select_name;

void mexFunction(int nlhs_m, mxArray *plhs_m[], int nrhs_m, const mxArray *prhs_m[])
{

	int i;

	int * m;
	int * n;
	int * nrhs;
	float * a;
	int * lda;
	float * b;
	int * ldb;
	int * jpvt;
	float * rcond;
	int * rank;
	float * work;
	int * lwork;
	float * rwork;
	int * info;

	plhs_m[0]=mxDuplicateArray(prhs_m[0]);
	m=(int *)mxGetPr(plhs_m[0]);

	plhs_m[1]=mxDuplicateArray(prhs_m[1]);
	n=(int *)mxGetPr(plhs_m[1]);

	plhs_m[2]=mxDuplicateArray(prhs_m[2]);
	nrhs=(int *)mxGetPr(plhs_m[2]);

	plhs_m[4]=mxDuplicateArray(prhs_m[4]);
	lda=(int *)mxGetPr(plhs_m[4]);

	plhs_m[6]=mxDuplicateArray(prhs_m[6]);
	ldb=(int *)mxGetPr(plhs_m[6]);

	plhs_m[9]=mxDuplicateArray(prhs_m[9]);
	rank=(int *)mxGetPr(plhs_m[9]);

	plhs_m[11]=mxDuplicateArray(prhs_m[11]);
	lwork=(int *)mxGetPr(plhs_m[11]);

	plhs_m[13]=mxDuplicateArray(prhs_m[13]);
	info=(int *)mxGetPr(plhs_m[13]);







	a=malloc((*lda)*(*n)*2*sizeof(float));
	for (i=0; i<(*lda)*(*n); i++)
	{
		a[i*2]=((float*)mxGetData(prhs_m[3]))[i];
		if (mxIsComplex(prhs_m[3]))
			a[i*2+1]=((float*)mxGetImagData(prhs_m[3]))[i];
		else
			a[i*2+1]=0;
	}



	b=malloc((*ldb)*(*nrhs)*2*sizeof(float));
	for (i=0; i<(*ldb)*(*nrhs); i++)
	{
		b[i*2]=((float*)mxGetData(prhs_m[5]))[i];
		if (mxIsComplex(prhs_m[5]))
			b[i*2+1]=((float*)mxGetImagData(prhs_m[5]))[i];
		else
			b[i*2+1]=0;
	}



	plhs_m[7]=mxDuplicateArray(prhs_m[7]);
	jpvt = (int *) mxGetPr(plhs_m[7]);

	plhs_m[8]=mxDuplicateArray(prhs_m[8]);
	rcond = (float*) mxGetPr(plhs_m[8]);



	work=malloc((max(1,*lwork))*(1)*2*sizeof(float));
	for (i=0; i<(max(1,*lwork))*(1); i++)
	{
		work[i*2]=((float*)mxGetData(prhs_m[10]))[i];
		if (mxIsComplex(prhs_m[10]))
			work[i*2+1]=((float*)mxGetImagData(prhs_m[10]))[i];
		else
			work[i*2+1]=0;
	}



	plhs_m[12]=mxDuplicateArray(prhs_m[12]);
	rwork = (float*) mxGetPr(plhs_m[12]);



#ifdef F77_INT
	F77_INT* F77_m = m , F77_n = n , F77_nrhs = nrhs , F77_lda = lda , F77_ldb = ldb , F77_jpvt = jpvt , F77_rank = rank , F77_lwork = lwork , F77_info = info ;
#else
	#define F77_m m 
	#define F77_n n 
	#define F77_nrhs nrhs 
	#define F77_lda lda 
	#define F77_ldb ldb 
	#define F77_jpvt jpvt 
	#define F77_rank rank 
	#define F77_lwork lwork 
	#define F77_info info 
#endif

	f77_cgelsy(F77_m, F77_n, F77_nrhs, a, F77_lda, b, F77_ldb, F77_jpvt, rcond, F77_rank, work, F77_lwork, rwork, F77_info);

	plhs_m[3] = mxCreateNumericMatrix((*lda),(*n),mxSINGLE_CLASS,mxCOMPLEX);
	for (i=0; i<(*lda)*(*n); i++)
	{
		((float*)mxGetData(plhs_m[3]))[i]=a[2*i];
		if (mxIsComplex(plhs_m[3]))
			((float*)mxGetImagData(plhs_m[3]))[i]=a[2*i+1];
	}

	plhs_m[5] = mxCreateNumericMatrix((*ldb),(*nrhs),mxSINGLE_CLASS,mxCOMPLEX);
	for (i=0; i<(*ldb)*(*nrhs); i++)
	{
		((float*)mxGetData(plhs_m[5]))[i]=b[2*i];
		if (mxIsComplex(plhs_m[5]))
			((float*)mxGetImagData(plhs_m[5]))[i]=b[2*i+1];
	}

	plhs_m[10] = mxCreateNumericMatrix((max(1,*lwork)),(1),mxSINGLE_CLASS,mxCOMPLEX);
	for (i=0; i<(max(1,*lwork))*(1); i++)
	{
		((float*)mxGetData(plhs_m[10]))[i]=work[2*i];
		if (mxIsComplex(plhs_m[10]))
			((float*)mxGetImagData(plhs_m[10]))[i]=work[2*i+1];
	}

	return;
}
