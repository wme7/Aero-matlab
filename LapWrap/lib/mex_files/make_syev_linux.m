%
mex isComplex.c -DADD_    
mex isDouble.c  -DADD_   
mex isinteger.c  -DADD_   
%
% SYEV / HEEV
%
mex lapack_cheev.c  -DADD_
mex lapack_cheevd.c -DADD_
mex lapack_cheevr.c -DADD_
mex lapack_cheevx.c -DADD_
%
mex lapack_dsyev.c  -DADD_
mex lapack_dsyevd.c -DADD_
mex lapack_dsyevr.c -DADD_
mex lapack_dsyevx.c -DADD_
%
mex lapack_ssyev.c  -DADD_
mex lapack_ssyevd.c -DADD_
mex lapack_ssyevr.c -DADD_
mex lapack_ssyevx.c -DADD_
%
mex lapack_zheev.c  -DADD_
mex lapack_zheevd.c -DADD_
mex lapack_zheevr.c -DADD_
mex lapack_zheevx.c -DADD_
%
