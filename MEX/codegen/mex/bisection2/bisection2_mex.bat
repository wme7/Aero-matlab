@echo off
set MATLAB=C:\PROGRA~1\MATLAB\R2012a
set MATLAB_ARCH=win64
set MATLAB_BIN="C:\Program Files\MATLAB\R2012a\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=.\
set LIB_NAME=bisection2_mex
set MEX_NAME=bisection2_mex
set MEX_EXT=.mexw64
call mexopts.bat
echo # Make settings for bisection2 > bisection2_mex.mki
echo COMPILER=%COMPILER%>> bisection2_mex.mki
echo COMPFLAGS=%COMPFLAGS%>> bisection2_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> bisection2_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> bisection2_mex.mki
echo LINKER=%LINKER%>> bisection2_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> bisection2_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> bisection2_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> bisection2_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> bisection2_mex.mki
echo BORLAND=%BORLAND%>> bisection2_mex.mki
echo OMPFLAGS= >> bisection2_mex.mki
echo OMPLINKFLAGS= >> bisection2_mex.mki
echo EMC_COMPILER=msvc100>> bisection2_mex.mki
echo EMC_CONFIG=optim>> bisection2_mex.mki
"C:\Program Files\MATLAB\R2012a\bin\win64\gmake" -B -f bisection2_mex.mk
