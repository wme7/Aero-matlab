/* 2D triangular mesh
   elem(:) comes from MATLAB
   index starts from 1
   Matlab version 1.0
   elem = int32(elem)
  function [neighbor,elem2edge,edge,bdEdge] = auxstructurec(elem)
 */

#include "mex.h"
#include "stdlib.h"

void auxstructurec ( int N,
					 int NT,
					 int *elem, 
					 int *bdEdgeSizePtr,
					 int *edgeSizePtr,
					 int **neighbor,
					 int **elem2edge,
					 int **edge,
					 int **bdEdge
				   )
{
	int i, j, k, p, t, lji, tt, pt, nnz, globalidx, edgeidx, bdSize;
	int *valence, *rowptr, *columnidx, *e2tidx, *localidx;
	int *tempEdge, *edge2elem, *isbdEdge;

	valence = (int *) calloc (   N, sizeof(int) );
	rowptr  = (int *) calloc ( N+1, sizeof(int) );
	
	/* valence = accumarray(elem(:),ones(3*NT,1),[N 1]);*/
	for ( i = 0; i < 3*NT; i++ )
	{
		valence[elem[i]-1]++;
	}
	
	/* sparse pattern of edge2elem */
	rowptr[0] = 1;
	for ( i = 1; i <= N; i++ )
	{
		rowptr[i] = rowptr[i-1] + valence[i-1];
	}
	
	nnz = rowptr[N] - 1;
	
	/* form edge2elem and e2tidx sparse matrix */
	columnidx = (int *) calloc ( nnz, sizeof(int) );
	e2tidx    = (int *) calloc ( nnz, sizeof(int) );
	localidx  = (int *) calloc (   N, sizeof(int) );
	edge2elem = (int *) calloc ( nnz, sizeof(int) );
	
	for ( t = 0; t < NT; t++ )
	{
		for ( p = 0; p < 3; p++ )
		{
			switch (p)
			{
				case 0:
					i = elem[t+NT];
					j = elem[t+NT*2];
					break;
				case 1:
					i = elem[t+NT*2];
					j = elem[t];
					break;
				case 2:
					i = elem[t];
					j = elem[t+NT];
					break;
			}
			globalidx = rowptr[i-1] + localidx[i-1];
			columnidx[globalidx-1] = j;
			edge2elem[globalidx-1] = t + 1;
			e2tidx[globalidx-1] = p + 1;
			localidx[i-1]++;
		}
	}
	
	/* form edge and neighbor matrix */
	(*elem2edge) = (int *) calloc ( NT*3, sizeof(int) );
	(*neighbor)  = (int *) calloc ( NT*3, sizeof(int) ); 
	tempEdge     = (int *) calloc (  N*8, sizeof(int) );
	isbdEdge     = (int *) calloc (  N*4, sizeof(int) );

	edgeidx = 0;	
	bdSize = 0;
	for ( t = 0; t < NT; t++ )
	{
		for ( p = 0; p < 3; p++ )
		{
			if ( (*elem2edge)[t+NT*p] == 0 )
			{
				switch (p)
				{
					case 0:
						i = elem[t+NT];
						j = elem[t+NT*2];
						break;
					case 1:
						i = elem[t+NT*2];
						j = elem[t];
						break;
					case 2:
						i = elem[t];
						j = elem[t+NT];
						break;
				}

				/* (*edge)[edgeidx-1,:] = sort(i,j);*/			
				if ( i < j )
				{
					tempEdge[edgeidx] = i;
					tempEdge[edgeidx+4*N] = j;
				}
				else
				{
					tempEdge[edgeidx] = j;
					tempEdge[edgeidx+4*N] = i;
				}
			
				(*elem2edge)[t+NT*p] = edgeidx+1;
			
				/* search the index of (j,i) in the sparse matrix*/
				lji = 0;
				for ( k = rowptr[j-1]; k <= rowptr[j]-1; k++ )
				{
					if ( columnidx[k-1] == i )
					{
						lji = k;
						break;
					}
				}
			
				if ( lji == 0 ) /* boundary edges */
				{
					(*neighbor)[t+NT*p] = t + 1;
					isbdEdge[edgeidx] = 1;
					bdSize++;
				}
				
				else /* lji ~= 0 interior edges*/
				{
					tt = edge2elem[lji-1];
					pt = e2tidx[lji-1];
					(*elem2edge)[(tt-1)+NT*(pt-1)] = edgeidx + 1;
					(*neighbor)[t+NT*p] = tt;
					(*neighbor)[(tt-1)+NT*(pt-1)] = t + 1;
				}
				
				edgeidx++;
			}
		}
	}
	
	(*bdEdge) = (int *) calloc (  bdSize*2, sizeof(int) );
	(*edge)   = (int *) calloc ( edgeidx*2, sizeof(int) );
	
	j = 0;
	
	for ( i = 0; i < edgeidx; i++ )
	{
		if ( isbdEdge[i] == 1 )
		{
			tempEdge[i];
			
			(*bdEdge)[j] = tempEdge[i];
			(*bdEdge)[j+bdSize] = tempEdge[i+4*N];
			j++;
		}

		(*edge)[i] = tempEdge[i];
		(*edge)[i+edgeidx] = tempEdge[i+4*N];
	}
	
	(*bdEdgeSizePtr) = bdSize;
	(*edgeSizePtr)   = edgeidx;
	
	free(valence);
	free(rowptr);
	free(columnidx);
	free(e2tidx);
	free(edge2elem);
	free(localidx);
	free(isbdEdge);
	free(tempEdge);
}

void mexFunction( int nlhs,       mxArray *plhs[], 
				  int nrhs, const mxArray *prhs[] 
				) 
{
	int i, NT, bdEdgeSize, edgeSize, N = 0;
	int *elem, *neighbor, *elem2edge, *edge, *bdEdge;
	int *output0, *output1, *output2, *output3;
	
	/* Check for proper number of arguments.*/
	if ( nrhs != 1 ) 
	{
		mexErrMsgTxt("Need one input argument.");
	} 
	else if ( nlhs > 4 ) 
	{
		mexErrMsgTxt("Too many output arguments.");
	}
	else ;

    NT   = mxGetM  ( prhs[0] );
	elem = (int *) mxGetData ( prhs[0] );
	
	/* Find largest index for nodes: n = max(elem(:))*/
	for ( i = 0; i < 3*NT; i++ )
	{
		if ( elem[i] > N )
		{
			N = elem[i];
		}
	}

	auxstructurec ( N, NT, elem, &bdEdgeSize, &edgeSize,
					&neighbor, &elem2edge, &edge, &bdEdge );

	/* Create matrix for the return argument. */
	plhs[0] = mxCreateNumericMatrix(NT, 3, mxINT32_CLASS, mxREAL);	
	plhs[1] = mxCreateNumericMatrix(NT, 3, mxINT32_CLASS, mxREAL);
	plhs[2] = mxCreateNumericMatrix(edgeSize,   2, mxINT32_CLASS, mxREAL);	
	plhs[3] = mxCreateNumericMatrix(bdEdgeSize, 2, mxINT32_CLASS, mxREAL);
	output0 = (int *) mxGetData (plhs[0]);
	output1 = (int *) mxGetData (plhs[1]);
	output2 = (int *) mxGetData (plhs[2]);
	output3 = (int *) mxGetData (plhs[3]);
	
	for ( i = 0; i < NT*3; i++ )
	{
		output0[i] = neighbor[i];
	}
	
	for ( i = 0; i < NT*3; i++ )
	{
		output1[i] = elem2edge[i];
	}
	
	for ( i = 0; i < edgeSize*2; i++ )
	{
		output2[i] = edge[i];
	}
	
	for ( i = 0; i < bdEdgeSize*2; i++ )
	{
		output3[i] = bdEdge[i];
	}

	free(neighbor); 
	free(elem2edge);
	free(edge);
	free(bdEdge);
}

