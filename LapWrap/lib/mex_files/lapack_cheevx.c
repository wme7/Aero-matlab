/*
 * clapack_cheevx.c
 *
 *This program is a C interface to cheevx.
 *
 * Written by Remi Delmas.
 *
 */

#include "mex.h"

#define min(a,b) ((a)>(b) ? (b) : (a))
#define max(a,b) ((a)<(b) ? (b) : (a))

#if defined(ADD_)
	#define f77_cheevx cheevx_
#elif defined(UPCASE)
	#define f77_cheevx CHEEVX
#else
	#define f77_cheevx cheevx
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
	float * rwork;
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

	plhs_m[20]=mxDuplicateArray(prhs_m[20]);
	info=(int *)mxGetPr(plhs_m[20]);

	plhs_m[0] = mxDuplicateArray(prhs_m[0]);
	jobz = (char*) mxArrayToString(plhs_m[0]);

	plhs_m[1] = mxDuplicateArray(prhs_m[1]);
	range = (char*) mxArrayToString(plhs_m[1]);

	plhs_m[2] = mxDuplicateArray(prhs_m[2]);
	uplo = (char*) mxArrayToString(plhs_m[2]);



	a=malloc((*lda)*(* n)*2*sizeof(float));
	for (i=0; i<(*lda)*(* n); i++)
	{
		a[i*2]=((float*)mxGetData(prhs_m[4]))[i];
		if (mxIsComplex(prhs_m[4]))
			a[i*2+1]=((float*)mxGetImagData(prhs_m[4]))[i];
		else
			a[i*2+1]=0;
	}



	plhs_m[6]=mxDuplicateArray(prhs_m[6]);
	vl = (float*) mxGetPr(plhs_m[6]);

	plhs_m[7]=mxDuplicateArray(prhs_m[7]);
	vu = (float*) mxGetPr(plhs_m[7]);





	plhs_m[10]=mxDuplicateArray(prhs_m[10]);
	abstol = (float*) mxGetPr(plhs_m[10]);



	plhs_m[12]=mxDuplicateArray(prhs_m[12]);
	w = (float*) mxGetPr(plhs_m[12]);

	z=malloc((*ldz)*( max(1,*n))*2*sizeof(float));
	for (i=0; i<(*ldz)*( max(1,*n)); i++)
	{
		z[i*2]=((float*)mxGetData(prhs_m[13]))[i];
		if (mxIsComplex(prhs_m[13]))
			z[i*2+1]=((float*)mxGetImagData(prhs_m[13]))[i];
		else
			z[i*2+1]=0;
	}



	work=malloc((max(1,*lwork))*(1)*2*sizeof(float));
	for (i=0; i<(max(1,*lwork))*(1); i++)
	{
		work[i*2]=((float*)mxGetData(prhs_m[15]))[i];
		if (mxIsComplex(prhs_m[15]))
			work[i*2+1]=((float*)mxGetImagData(prhs_m[15]))[i];
		else
			work[i*2+1]=0;
	}



	plhs_m[17]=mxDuplicateArray(prhs_m[17]);
	rwork = (float*) mxGetPr(plhs_m[17]);

	plhs_m[18]=mxDuplicateArray(prhs_m[18]);
	iwork = (int *) mxGetPr(plhs_m[18]);

	plhs_m[19]=mxDuplicateArray(prhs_m[19]);
	ifail = (int *) mxGetPr(plhs_m[19]);



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

	f77_cheevx(F77_jobz, F77_range, F77_uplo, F77_n, a, F77_lda, vl, vu, F77_il, F77_iu, abstol, F77_m, w, z, F77_ldz, work, F77_lwork, rwork, F77_iwork, F77_ifail, F77_info);

	plhs_m[4] = mxCreateNumericMatrix((*lda),(* n),mxSINGLE_CLASS,mxCOMPLEX);
	for (i=0; i<(*lda)*(* n); i++)
	{
		((float*)mxGetData(plhs_m[4]))[i]=a[2*i];
		if (mxIsComplex(plhs_m[4]))
			((float*)mxGetImagData(plhs_m[4]))[i]=a[2*i+1];
	}

	plhs_m[13] = mxCreateNumericMatrix((*ldz),( max(1,*n)),mxSINGLE_CLASS,mxCOMPLEX);
	for (i=0; i<(*ldz)*( max(1,*n)); i++)
	{
		((float*)mxGetData(plhs_m[13]))[i]=z[2*i];
		if (mxIsComplex(plhs_m[13]))
			((float*)mxGetImagData(plhs_m[13]))[i]=z[2*i+1];
	}

	plhs_m[15] = mxCreateNumericMatrix((max(1,*lwork)),(1),mxSINGLE_CLASS,mxCOMPLEX);
	for (i=0; i<(max(1,*lwork))*(1); i++)
	{
		((float*)mxGetData(plhs_m[15]))[i]=work[2*i];
		if (mxIsComplex(plhs_m[15]))
			((float*)mxGetImagData(plhs_m[15]))[i]=work[2*i+1];
	}

	return;
}
