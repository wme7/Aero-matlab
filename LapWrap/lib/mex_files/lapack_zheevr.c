/*
 * clapack_zheevr.c
 *
 *This program is a C interface to zheevr.
 *
 * Written by Remi Delmas.
 *
 */

#include "mex.h"

#define min(a,b) ((a)>(b) ? (b) : (a))
#define max(a,b) ((a)<(b) ? (b) : (a))

#if defined(ADD_)
	#define f77_zheevr zheevr_
#elif defined(UPCASE)
	#define f77_zheevr ZHEEVR
#else
	#define f77_zheevr zheevr
#endif

char *select_name;

void mexFunction(int nlhs_m, mxArray *plhs_m[], int nrhs_m, const mxArray *prhs_m[])
{

	int i;

	char * jobz;
	char * range;
	char * uplo;
	int * n;
	double * a;
	int * lda;
	double * vl;
	double * vu;
	int * il;
	int * iu;
	double * abstol;
	int * m;
	double * w;
	double * z;
	int * ldz;
	int * isuppz;
	double * work;
	int * lwork;
	double * rwork;
	int * lrwork;
	int * iwork;
	int * liwork;
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

	plhs_m[17]=mxDuplicateArray(prhs_m[17]);
	lwork=(int *)mxGetPr(plhs_m[17]);

	plhs_m[19]=mxDuplicateArray(prhs_m[19]);
	lrwork=(int *)mxGetPr(plhs_m[19]);

	plhs_m[21]=mxDuplicateArray(prhs_m[21]);
	liwork=(int *)mxGetPr(plhs_m[21]);

	plhs_m[22]=mxDuplicateArray(prhs_m[22]);
	info=(int *)mxGetPr(plhs_m[22]);

	plhs_m[0] = mxDuplicateArray(prhs_m[0]);
	jobz = (char*) mxArrayToString(plhs_m[0]);

	plhs_m[1] = mxDuplicateArray(prhs_m[1]);
	range = (char*) mxArrayToString(plhs_m[1]);

	plhs_m[2] = mxDuplicateArray(prhs_m[2]);
	uplo = (char*) mxArrayToString(plhs_m[2]);



	a=malloc((*lda)*(* n)*2*sizeof(double));
	for (i=0; i<(*lda)*(* n); i++)
	{
		a[i*2]=(mxGetPr(prhs_m[4]))[i];
		if (mxIsComplex(prhs_m[4]))
			a[i*2+1]=(mxGetPi(prhs_m[4]))[i];
		else
			a[i*2+1]=0;
	}



	plhs_m[6]=mxDuplicateArray(prhs_m[6]);
	vl = mxGetPr(plhs_m[6]);

	plhs_m[7]=mxDuplicateArray(prhs_m[7]);
	vu = mxGetPr(plhs_m[7]);





	plhs_m[10]=mxDuplicateArray(prhs_m[10]);
	abstol = mxGetPr(plhs_m[10]);



	plhs_m[12]=mxDuplicateArray(prhs_m[12]);
	w = mxGetPr(plhs_m[12]);

	z=malloc((*ldz)*( max(1,*n))*2*sizeof(double));
	for (i=0; i<(*ldz)*( max(1,*n)); i++)
	{
		z[i*2]=(mxGetPr(prhs_m[13]))[i];
		if (mxIsComplex(prhs_m[13]))
			z[i*2+1]=(mxGetPi(prhs_m[13]))[i];
		else
			z[i*2+1]=0;
	}



	plhs_m[15]=mxDuplicateArray(prhs_m[15]);
	isuppz = (int *) mxGetPr(plhs_m[15]);

	work=malloc((max(1,*lwork))*(1)*2*sizeof(double));
	for (i=0; i<(max(1,*lwork))*(1); i++)
	{
		work[i*2]=(mxGetPr(prhs_m[16]))[i];
		if (mxIsComplex(prhs_m[16]))
			work[i*2+1]=(mxGetPi(prhs_m[16]))[i];
		else
			work[i*2+1]=0;
	}



	plhs_m[18]=mxDuplicateArray(prhs_m[18]);
	rwork = mxGetPr(plhs_m[18]);



	plhs_m[20]=mxDuplicateArray(prhs_m[20]);
	iwork = (int *) mxGetPr(plhs_m[20]);





#ifdef F77_INT
	F77_INT* F77_n = n , F77_lda = lda , F77_il = il , F77_iu = iu , F77_m = m , F77_ldz = ldz , F77_isuppz = isuppz , F77_lwork = lwork , F77_lrwork = lrwork , F77_iwork = iwork , F77_liwork = liwork , F77_info = info ;
#else
	#define F77_n n 
	#define F77_lda lda 
	#define F77_il il 
	#define F77_iu iu 
	#define F77_m m 
	#define F77_ldz ldz 
	#define F77_isuppz isuppz 
	#define F77_lwork lwork 
	#define F77_lrwork lrwork 
	#define F77_iwork iwork 
	#define F77_liwork liwork 
	#define F77_info info 
#endif

#ifdef F77_CHAR
	F77_CHAR* F77_jobz = jobz , F77_range = range , F77_uplo = uplo ;
#else
	#define F77_jobz jobz 
	#define F77_range range 
	#define F77_uplo uplo 
#endif

	f77_zheevr(F77_jobz, F77_range, F77_uplo, F77_n, a, F77_lda, vl, vu, F77_il, F77_iu, abstol, F77_m, w, z, F77_ldz, F77_isuppz, work, F77_lwork, rwork, F77_lrwork, F77_iwork, F77_liwork, F77_info);

	plhs_m[4] = mxCreateDoubleMatrix((*lda),(* n),mxCOMPLEX);
	for (i=0; i<(*lda)*(* n); i++)
	{
		(mxGetPr(plhs_m[4]))[i]=a[2*i];
		if (mxIsComplex(plhs_m[4]))
			(mxGetPi(plhs_m[4]))[i]=a[2*i+1];
	}

	plhs_m[13] = mxCreateDoubleMatrix((*ldz),( max(1,*n)),mxCOMPLEX);
	for (i=0; i<(*ldz)*( max(1,*n)); i++)
	{
		(mxGetPr(plhs_m[13]))[i]=z[2*i];
		if (mxIsComplex(plhs_m[13]))
			(mxGetPi(plhs_m[13]))[i]=z[2*i+1];
	}

	plhs_m[16] = mxCreateDoubleMatrix((max(1,*lwork)),(1),mxCOMPLEX);
	for (i=0; i<(max(1,*lwork))*(1); i++)
	{
		(mxGetPr(plhs_m[16]))[i]=work[2*i];
		if (mxIsComplex(plhs_m[16]))
			(mxGetPi(plhs_m[16]))[i]=work[2*i+1];
	}

	return;
}
