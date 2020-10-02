/*
 * diag.c
 *
 * Code generation for function 'diag'
 *
 * C source code generated on: Sun Dec 14 15:48:09 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "matlabrun.h"
#include "diag.h"
#include <stdio.h>

/* Function Definitions */
void diag(const creal_T *v, creal_T *d, uint_T sz)
{
  uint_T j;
  for (j = 0; j < sz*sz; j++) {
    d[j].re = 0.0;
    d[j].im = 0.0;
  }

  for (j = 0; j < sz; j++) {
    d[j + sz * j] = v[j];
  }
}

/* End of code generation (diag.c) */
