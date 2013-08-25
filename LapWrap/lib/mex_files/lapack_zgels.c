/*
 * clapack_zgels.c
 *
 *This program is a C interface to zgels.
 *
 * Written by Remi Delmas.
 *
 */

#include "mex.h"

#define min(a,b) ((a)>(b) ? (b) : (a))
#define max(a,b) ((a)<(b) ? (b) : (a))

#if defined(ADD_)
	#define f77_zgels zgels_
#elif defined(UPCASE)
	#define f77_zgels ZGELS
#else
	#define f77_zgels zgels
#endif

char *select_name;

void mexFunction(int nlhs_m, mxArray *plhs_m[], int nrhs_m, const mxArray *prhs_m[])
{

	int i;

	char * trans;
	int * m;
	int * n;
	int * nrhs;
	double * a;
	int * lda;
	double * b;
	int * ldb;
	double * work;
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







	a=malloc((*lda)*(*n)*2*sizeof(double));
	for (i=0; i<(*lda)*(*n); i++)
	{
		a[i*2]=(mxGetPr(prhs_m[4]))[i];
		if (mxIsComplex(prhs_m[4]))
			a[i*2+1]=(mxGetPi(prhs_m[4]))[i];
		else
			a[i*2+1]=0;
	}



	b=malloc((*ldb)*(*nrhs)*2*sizeof(double));
	for (i=0; i<(*ldb)*(*nrhs); i++)
	{
		b[i*2]=(mxGetPr(prhs_m[6]))[i];
		if (mxIsComplex(prhs_m[6]))
			b[i*2+1]=(mxGetPi(prhs_m[6]))[i];
		else
			b[i*2+1]=0;
	}



	work=malloc((max(1,*lwork))*(1)*2*sizeof(double));
	for (i=0; i<(max(1,*lwork))*(1); i++)
	{
		work[i*2]=(mxGetPr(prhs_m[8]))[i];
		if (mxIsComplex(prhs_m[8]))
			work[i*2+1]=(mxGetPi(prhs_m[8]))[i];
		else
			work[i*2+1]=0;
	}





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

	f77_zgels(F77_trans, F77_m, F77_n, F77_nrhs, a, F77_lda, b, F77_ldb, work, F77_lwork, F77_info);

	plhs_m[4] = mxCreateDoubleMatrix((*lda),(*n),mxCOMPLEX);
	for (i=0; i<(*lda)*(*n); i++)
	{
		(mxGetPr(plhs_m[4]))[i]=a[2*i];
		if (mxIsComplex(plhs_m[4]))
			(mxGetPi(plhs_m[4]))[i]=a[2*i+1];
	}

	plhs_m[6] = mxCreateDoubleMatrix((*ldb),(*nrhs),mxCOMPLEX);
	for (i=0; i<(*ldb)*(*nrhs); i++)
	{
		(mxGetPr(plhs_m[6]))[i]=b[2*i];
		if (mxIsComplex(plhs_m[6]))
			(mxGetPi(plhs_m[6]))[i]=b[2*i+1];
	}

	plhs_m[8] = mxCreateDoubleMatrix((max(1,*lwork)),(1),mxCOMPLEX);
	for (i=0; i<(max(1,*lwork))*(1); i++)
	{
		(mxGetPr(plhs_m[8]))[i]=work[2*i];
		if (mxIsComplex(plhs_m[8]))
			(mxGetPi(plhs_m[8]))[i]=work[2*i+1];
	}

	return;
}
