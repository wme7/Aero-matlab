/*
 * bisection2_terminate.c
 *
 * Code generation for function 'bisection2_terminate'
 *
 * C source code generated on: Fri Dec 14 01:42:57 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "bisection2.h"
#include "bisection2_terminate.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void bisection2_atexit(void)
{
  emlrtEnterRtStack(&emlrtContextGlobal);
  emlrtLeaveRtStack(&emlrtContextGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void bisection2_terminate(void)
{
  emlrtLeaveRtStack(&emlrtContextGlobal);
}

/* End of code generation (bisection2_terminate.c) */
