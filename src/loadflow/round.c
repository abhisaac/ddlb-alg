/*
 * round.c
 *
 * Code generation for function 'round'
 *
 * C source code generated on: Sun Dec 14 15:48:09 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "matlabrun.h"
#include "round.h"
#include "loadflow.h"
#include "matlabrun_rtwutil.h"
#include <stdio.h>

/* Function Definitions */
void b_round(real_T *x, uint_T sz)
{
  uint_T k;
  for (k = 0; k < sz; k++) {
    x[k] = rt_roundd_snf(x[k]);
  }
}

/* End of code generation (round.c) */
