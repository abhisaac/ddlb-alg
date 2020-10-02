/*
 * rdivide.c
 *
 * Code generation for function 'rdivide'
 *
 * C source code generated on: Mon Dec 15 20:15:42 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "matlabrun.h"
#include "rdivide.h"
#include "matlabrun_emxutil.h"
#include <stdio.h>

/* Function Definitions */
void rdivide(const emxArray_real_T *x, const creal_T *y, emxArray_creal_T *z, uint_T sz)
{
  uint_T i2;
  real_T x_re;
  real_T brm;
  real_T bim;
  real_T d;
  i2 = z->size[0];
  z->size[0] = sz;
  emxEnsureCapacity((emxArray__common *)z, i2, (int32_T)sizeof(creal_T));
  for (i2 = 0; i2 < sz; i2++) {
    x_re = x->data[i2];
    if (y[i2].im == 0.0) {
      z->data[i2].re = x_re / y[i2].re;
      z->data[i2].im = 0.0;
    } else if (y[i2].re == 0.0) {
      if (x_re == 0.0) {
        z->data[i2].re = 0.0 / y[i2].im;
        z->data[i2].im = 0.0;
      } else {
        z->data[i2].re = 0.0;
        z->data[i2].im = -(x_re / y[i2].im);
      }
    } else {
      brm = fabs(y[i2].re);
      bim = fabs(y[i2].im);
      if (brm > bim) {
        bim = y[i2].im / y[i2].re;
        d = y[i2].re + bim * y[i2].im;
        z->data[i2].re = (x_re + bim * 0.0) / d;
        z->data[i2].im = (0.0 - bim * x_re) / d;
      } else if (bim == brm) {
        if (y[i2].re > 0.0) {
          bim = 0.5;
        } else {
          bim = -0.5;
        }

        if (y[i2].im > 0.0) {
          d = 0.5;
        } else {
          d = -0.5;
        }

        z->data[i2].re = x_re * bim / brm;
        z->data[i2].im = (-0.0 - x_re * d) / brm;
      } else {
        bim = y[i2].re / y[i2].im;
        d = y[i2].im + bim * y[i2].re;
        z->data[i2].re = bim * x_re / d;
        z->data[i2].im = (bim * 0.0 - x_re) / d;
      }
    }
  }
}

/* End of code generation (rdivide.c) */
