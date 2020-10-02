/*
 * norm.c
 *
 * Code generation for function 'norm'
 *
 * C source code generated on: Sun Dec 14 15:48:09 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "matlabrun.h"
#include "norm.h"
#include <stdio.h>

/* Function Definitions */
real_T norm(const real_T *x, uint_T sz)
{
  real_T y;
  uint_T k;
  boolean_T exitg1;
  real_T absx;
  y = 0.0;
  k = 0;
  exitg1 = FALSE;
  while ((exitg1 == FALSE) && (k < sz)) {
    absx = fabs(x[k]);
    if (rtIsNaN(absx)) {
      y = rtNaN;
      exitg1 = TRUE;
    } else {
      if (absx > y) {
        y = absx;
      }

      k++;
    }
  }

  return y;
}

/* End of code generation (norm.c) */
