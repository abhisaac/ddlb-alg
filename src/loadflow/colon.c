/*
 * colon.c
 *
 * Code generation for function 'colon'
 *
 * C source code generated on: Sun Dec 14 15:48:09 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "matlabrun.h"
#include "colon.h"
#include "matlabrun_emxutil.h"
#include <stdio.h>

/* Function Definitions */
void eml_signed_integer_colon(int32_T b, emxArray_int32_T *y)
{
  int32_T n;
  int32_T yk;
  int32_T k;
  if (b < 1) {
    n = 0;
  } else {
    n = b;
  }

  yk = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = n;
  emxEnsureCapacity((emxArray__common *)y, yk, (int32_T)sizeof(int32_T));
  if (n > 0) {
    y->data[0] = 1;
    yk = 1;
    for (k = 2; k <= n; k++) {
      yk++;
      y->data[k - 1] = yk;
    }
  }
}

/* End of code generation (colon.c) */
