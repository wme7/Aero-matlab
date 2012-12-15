/*
 * bisection2.c
 *
 * Code generation for function 'bisection2'
 *
 * C source code generated on: Fri Dec 14 01:42:57 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "bisection2.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 11, "bisection2",
  "E:/Documents/Matlab/aero-matlab/MEX/bisection2.m" };

static emlrtRSInfo b_emlrtRSI = { 16, "error",
  "C:/Program Files/MATLAB/R2012a/toolbox/eml/lib/matlab/lang/error.m" };

static emlrtMCInfo emlrtMCI = { 16, 7, "error",
  "C:/Program Files/MATLAB/R2012a/toolbox/eml/lib/matlab/lang/error.m" };

/* Function Declarations */
static void error(const mxArray *b, emlrtMCInfo *location);

/* Function Definitions */
static void error(const mxArray *b, emlrtMCInfo *location)
{
  const mxArray *pArray;
  pArray = b;
  emlrtCallMATLAB(0, NULL, 1, &pArray, "error", TRUE, location);
}

real_T bisection2(real_T a, real_T b, real_T tol)
{
  real_T p;
  const mxArray *y;
  static const int32_T iv0[2] = { 1, 16 };

  const mxArray *m0;
  static const char_T cv0[16] = { 'W', 'r', 'o', 'n', 'g', ' ', 'c', 'h', 'o',
    'i', 'c', 'e', ' ', 'b', 'r', 'o' };

  real_T err;

  /*  provide the equation you want to solve with R.H.S = 0 form.  */
  /*  Write the L.H.S by using inline function */
  /*  Give initial guesses. */
  /*  Solves it by method of bisection. */
  /*  A very simple code and very handy! */
  if ((muDoubleScalarPower(a - 1.5, 3.0) + 10.0) * (muDoubleScalarPower(b - 1.5,
        3.0) + 10.0) > 0.0) {
    EMLRTPUSHRTSTACK(&emlrtRSI);
    EMLRTPUSHRTSTACK(&b_emlrtRSI);
    y = NULL;
    m0 = mxCreateCharArray(2, iv0);
    emlrtInitCharArray(16, m0, cv0);
    emlrtAssign(&y, m0);
    error(y, &emlrtMCI);
    EMLRTPOPRTSTACK(&b_emlrtRSI);
    EMLRTPOPRTSTACK(&emlrtRSI);
  } else {
    p = (a + b) / 2.0;
    err = muDoubleScalarAbs(b - a);
    while (err >= tol) {
      if ((muDoubleScalarPower(a - 1.5, 3.0) + 10.0) * (muDoubleScalarPower(p -
            1.5, 3.0) + 10.0) < 0.0) {
        b = p;
      } else {
        a = p;
      }

      p = (a + b) / 2.0;
      err = muDoubleScalarAbs(b - a);
      emlrtBreakCheck();
    }
  }

  return p;
}

/* End of code generation (bisection2.c) */
