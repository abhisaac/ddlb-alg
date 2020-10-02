/*
 * cos.c
 *
 * Code generation for function 'cos'
 *
 * C source code generated on: Sun Dec 14 15:48:09 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "matlabrun.h"
#include "cos.h"
#include <stdio.h>

/* Function Definitions */
void b_cos(real_T *x, uint_T sz)
{
  uint_T k;
  for (k = 0; k < sz; k++) {
    x[k] = cos(x[k]);
  }
}

/* End of code generation (cos.c) */
