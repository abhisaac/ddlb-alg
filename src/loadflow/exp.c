/*
 * exp.c
 *
 * Code generation for function 'exp'
 *
 * C source code generated on: Sun Dec 14 15:48:09 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "matlabrun.h"
#include "exp.h"
#include <stdio.h>

/* Function Definitions */
void b_exp(creal_T *x)
{
  real_T r;
  real_T x_im;
  real_T b_x_im;
  r = exp(x->re / 2.0);
  x_im = x->im;
  b_x_im = x->im;
  x->re = r * (r * cos(x_im));
  x->im = r * (r * sin(b_x_im));
}

/* End of code generation (exp.c) */
