%
mex isComplex.c -DADD_    
mex isinteger.c -DADD_    
mex isDouble.c  -DADD_   
%
% GELS
%
mex lapack_cgels.c  -DADD_
mex lapack_cgelsd.c -DADD_
mex lapack_cgelss.c -DADD_
mex lapack_cgelsy.c -DADD_
%
mex lapack_zgels.c  -DADD_
mex lapack_zgelsd.c -DADD_
mex lapack_zgelss.c -DADD_
mex lapack_zgelsy.c -DADD_
%
mex lapack_dgels.c  -DADD_
mex lapack_dgelsd.c -DADD_
mex lapack_dgelss.c -DADD_
mex lapack_dgelsy.c -DADD_
%
mex lapack_sgels.c  -DADD_
mex lapack_sgelsd.c -DADD_
mex lapack_sgelss.c -DADD_
mex lapack_sgelsy.c -DADD_
