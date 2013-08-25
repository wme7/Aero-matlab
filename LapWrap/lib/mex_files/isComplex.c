/*
 * clapack_dgetrf.c
 *
 *This program is a C interface to dgetrf.
 *
 * Written by Remi Delmas.
 *
 */

#include "mex.h"

#define min(a,b) ((a)>(b) ? (b) : (a))
#define max(a,b) ((a)<(b) ? (b) : (a))

char *select_name;

void mexFunction(int nlhs_m, mxArray *plhs_m[], int nrhs_m, const mxArray *prhs_m[])
{
		
	if (nlhs_m > 1)
	{
		printf("Error, only one ouput parameter allowed!\n");
		return;
	}

	if (nrhs_m != 1)
	{
		printf("Error, only one input parameter allowed!\n");
		return;
	}
	
	plhs_m[0] = mxCreateLogicalMatrix(1,1);		
	
	*(mxGetLogicals(plhs_m[0])) = mxIsComplex(prhs_m[0]);

	return;
}
