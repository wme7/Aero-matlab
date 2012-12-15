/*
 * bisection2_api.c
 *
 * Code generation for function 'bisection2_api'
 *
 * C source code generated on: Fri Dec 14 01:42:57 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "bisection2.h"
#include "bisection2_api.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static real_T b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId);
static real_T c_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId);
static real_T emlrt_marshallIn(const mxArray *a, const char_T *identifier);
static const mxArray *emlrt_marshallOut(real_T u);

/* Function Definitions */
static real_T b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId)
{
  real_T y;
  y = c_emlrt_marshallIn(emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T c_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId)
{
  real_T ret;
  emlrtCheckBuiltInCtxR2011b(&emlrtContextGlobal, msgId, src, "double", FALSE,
    0U, 0);
  ret = *(real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T emlrt_marshallIn(const mxArray *a, const char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = b_emlrt_marshallIn(emlrtAlias(a), &thisId);
  emlrtDestroyArray(&a);
  return y;
}

static const mxArray *emlrt_marshallOut(real_T u)
{
  const mxArray *y;
  const mxArray *m1;
  y = NULL;
  m1 = mxCreateDoubleScalar(u);
  emlrtAssign(&y, m1);
  return y;
}

void bisection2_api(const mxArray * const prhs[3], const mxArray *plhs[1])
{
  real_T a;
  real_T b;
  real_T tol;

  /* Marshall function inputs */
  a = emlrt_marshallIn(emlrtAliasP(prhs[0]), "a");
  b = emlrt_marshallIn(emlrtAliasP(prhs[1]), "b");
  tol = emlrt_marshallIn(emlrtAliasP(prhs[2]), "tol");

  /* Invoke the target function */
  a = bisection2(a, b, tol);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(a);
}

/* End of code generation (bisection2_api.c) */
