/*
 * sin.c
 *
 * Code generation for function 'sin'
 *
 * C source code generated on: Sun Dec 14 15:48:09 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "matlabrun.h"
#include "sin.h"
#include <stdio.h>

/* Function Definitions */
void b_sin(real_T *x, uint_T sz)
{
  uint_T k;
  for (k = 0; k < sz; k++) {
    x[k] = sin(x[k]);
  }
}

/* End of code generation (sin.c) */
