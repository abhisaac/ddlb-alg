/*
 * mldivide.c
 *
 * Code generation for function 'mldivide'
 *
 * C source code generated on: Sun Dec 14 15:48:09 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "matlabrun.h"
#include "mldivide.h"
#include "matlabrun_emxutil.h"
#include "colon.h"
#include <stdio.h>

/* Function Declarations */
static real_T eml_matlab_zlarfg(int32_T n, real_T *alpha1, emxArray_real_T *x,
  int32_T ix0);
static void eml_qrsolve(const emxArray_real_T *A, emxArray_real_T *B,
  emxArray_real_T *Y);
static real_T eml_xnrm2(int32_T n, const emxArray_real_T *x, int32_T ix0);
static void eml_xscal(int32_T n, real_T a, emxArray_real_T *x, int32_T ix0);
static real_T rt_hypotd_snf(real_T u0, real_T u1);

/* Function Definitions */
static real_T eml_matlab_zlarfg(int32_T n, real_T *alpha1, emxArray_real_T *x,
  int32_T ix0)
{
  real_T tau;
  real_T xnorm;
  int32_T knt;
  real_T d0;
  int32_T k;
  tau = 0.0;
  if (n <= 0) {
  } else {
    xnorm = eml_xnrm2(n - 1, x, ix0);
    if (xnorm != 0.0) {
      xnorm = rt_hypotd_snf(*alpha1, xnorm);
      if (*alpha1 >= 0.0) {
        xnorm = -xnorm;
      }

      if (fabs(xnorm) < 1.0020841800044864E-292) {
        knt = 0;
        do {
          knt++;
          eml_xscal(n - 1, 9.9792015476736E+291, x, ix0);
          xnorm *= 9.9792015476736E+291;
          *alpha1 *= 9.9792015476736E+291;
        } while (!(fabs(xnorm) >= 1.0020841800044864E-292));

        xnorm = eml_xnrm2(n - 1, x, ix0);
        xnorm = rt_hypotd_snf(*alpha1, xnorm);
        if (*alpha1 >= 0.0) {
          xnorm = -xnorm;
        }

        tau = (xnorm - *alpha1) / xnorm;
        d0 = 1.0 / (*alpha1 - xnorm);
        eml_xscal(n - 1, d0, x, ix0);
        for (k = 1; k <= knt; k++) {
          xnorm *= 1.0020841800044864E-292;
        }

        *alpha1 = xnorm;
      } else {
        tau = (xnorm - *alpha1) / xnorm;
        d0 = 1.0 / (*alpha1 - xnorm);
        eml_xscal(n - 1, d0, x, ix0);
        *alpha1 = xnorm;
      }
    }
  }

  return tau;
}

static void eml_qrsolve(const emxArray_real_T *A, emxArray_real_T *B,
  emxArray_real_T *Y)
{
  emxArray_real_T *b_A;
  real_T smax;
  real_T s;
  real_T y;
  int32_T mn;
  int32_T i1;
  int32_T k;
  emxArray_real_T *tau;
  emxArray_int32_T *jpvt;
  int32_T pvt;
  int32_T b_mn;
  emxArray_real_T *work;
  emxArray_real_T *vn1;
  emxArray_real_T *vn2;
  int32_T nmi;
  int32_T i;
  int32_T i_i;
  int32_T mmi;
  int32_T ix;
  int32_T iy;
  int32_T i_ip1;
  int32_T lastv;
  int32_T lastc;
  boolean_T exitg3;
  int32_T exitg2;
  real_T t;
  boolean_T exitg1;
  uint32_T unnamed_idx_0;
  b_emxInit_real_T(&b_A, 2);
  smax = A->size[0];
  s = A->size[1];
  if (smax <= s) {
    y = smax;
  } else {
    y = s;
  }

  mn = (int32_T)y;
  i1 = b_A->size[0] * b_A->size[1];
  b_A->size[0] = A->size[0];
  b_A->size[1] = A->size[1];
  emxEnsureCapacity((emxArray__common *)b_A, i1, (int32_T)sizeof(real_T));
  k = A->size[0] * A->size[1];
  for (i1 = 0; i1 < k; i1++) {
    b_A->data[i1] = A->data[i1];
  }

  emxInit_real_T(&tau, 1);
  emxInit_int32_T(&jpvt, 2);
  k = A->size[0];
  pvt = A->size[1];
  if (k <= pvt) {
    b_mn = k;
  } else {
    b_mn = pvt;
  }

  i1 = tau->size[0];
  tau->size[0] = b_mn;
  emxEnsureCapacity((emxArray__common *)tau, i1, (int32_T)sizeof(real_T));
  eml_signed_integer_colon(A->size[1], jpvt);
  if ((A->size[0] == 0) || (A->size[1] == 0)) {
  } else {
    emxInit_real_T(&work, 1);
    k = A->size[1];
    i1 = work->size[0];
    work->size[0] = k;
    emxEnsureCapacity((emxArray__common *)work, i1, (int32_T)sizeof(real_T));
    for (i1 = 0; i1 < k; i1++) {
      work->data[i1] = 0.0;
    }

    emxInit_real_T(&vn1, 1);
    emxInit_real_T(&vn2, 1);
    k = A->size[1];
    i1 = vn1->size[0];
    vn1->size[0] = k;
    emxEnsureCapacity((emxArray__common *)vn1, i1, (int32_T)sizeof(real_T));
    i1 = vn2->size[0];
    vn2->size[0] = k;
    emxEnsureCapacity((emxArray__common *)vn2, i1, (int32_T)sizeof(real_T));
    k = 1;
    for (nmi = 0; nmi + 1 <= A->size[1]; nmi++) {
      vn1->data[nmi] = eml_xnrm2(A->size[0], A, k);
      vn2->data[nmi] = vn1->data[nmi];
      k += A->size[0];
    }

    for (i = 0; i + 1 <= b_mn; i++) {
      i_i = i + i * A->size[0];
      nmi = A->size[1] - i;
      mmi = (A->size[0] - i) - 1;
      if (nmi < 1) {
        pvt = -1;
      } else {
        pvt = 0;
        if (nmi > 1) {
          ix = i;
          smax = fabs(vn1->data[i]);
          for (k = 0; k + 2 <= nmi; k++) {
            ix++;
            s = fabs(vn1->data[ix]);
            if (s > smax) {
              pvt = k + 1;
              smax = s;
            }
          }
        }
      }

      pvt += i;
      if (pvt + 1 != i + 1) {
        ix = A->size[0] * pvt;
        iy = A->size[0] * i;
        for (k = 1; k <= A->size[0]; k++) {
          smax = b_A->data[ix];
          b_A->data[ix] = b_A->data[iy];
          b_A->data[iy] = smax;
          ix++;
          iy++;
        }

        k = jpvt->data[pvt];
        jpvt->data[pvt] = jpvt->data[i];
        jpvt->data[i] = k;
        vn1->data[pvt] = vn1->data[i];
        vn2->data[pvt] = vn2->data[i];
      }

      if (i + 1 < A->size[0]) {
        s = b_A->data[i_i];
        tau->data[i] = eml_matlab_zlarfg(mmi + 1, &s, b_A, i_i + 2);
      } else {
        smax = b_A->data[i_i];
        s = b_A->data[i_i];
        b_A->data[i_i] = smax;
        tau->data[i] = 0.0;
      }

      b_A->data[i_i] = s;
      if (i + 1 < A->size[1]) {
        s = b_A->data[i_i];
        b_A->data[i_i] = 1.0;
        i_ip1 = (i + (i + 1) * A->size[0]) + 1;
        if (tau->data[i] != 0.0) {
          lastv = mmi;
          pvt = i_i + mmi;
          while ((lastv + 1 > 0) && (b_A->data[pvt] == 0.0)) {
            lastv--;
            pvt--;
          }

          lastc = nmi - 1;
          exitg3 = FALSE;
          while ((exitg3 == FALSE) && (lastc > 0)) {
            k = i_ip1 + (lastc - 1) * A->size[0];
            nmi = k;
            do {
              exitg2 = 0;
              if (nmi <= k + lastv) {
                if (b_A->data[nmi - 1] != 0.0) {
                  exitg2 = 1;
                } else {
                  nmi++;
                }
              } else {
                lastc--;
                exitg2 = 2;
              }
            } while (exitg2 == 0);

            if (exitg2 == 1) {
              exitg3 = TRUE;
            }
          }
        } else {
          lastv = -1;
          lastc = 0;
        }

        if (lastv + 1 > 0) {
          if (lastc == 0) {
          } else {
            for (iy = 1; iy <= lastc; iy++) {
              work->data[iy - 1] = 0.0;
            }

            iy = 0;
            i1 = i_ip1 + A->size[0] * (lastc - 1);
            for (pvt = i_ip1; pvt <= i1; pvt += A->size[0]) {
              ix = i_i;
              smax = 0.0;
              k = pvt + lastv;
              for (nmi = pvt; nmi <= k; nmi++) {
                smax += b_A->data[nmi - 1] * b_A->data[ix];
                ix++;
              }

              work->data[iy] += smax;
              iy++;
            }
          }

          if (-tau->data[i] == 0.0) {
          } else {
            pvt = 0;
            for (nmi = 1; nmi <= lastc; nmi++) {
              if (work->data[pvt] != 0.0) {
                smax = work->data[pvt] * -tau->data[i];
                ix = i_i;
                i1 = lastv + i_ip1;
                for (k = i_ip1; k <= i1; k++) {
                  b_A->data[k - 1] += b_A->data[ix] * smax;
                  ix++;
                }
              }

              pvt++;
              i_ip1 += A->size[0];
            }
          }
        }

        b_A->data[i_i] = s;
      }

      for (nmi = i + 1; nmi + 1 <= A->size[1]; nmi++) {
        k = (i + A->size[0] * nmi) + 1;
        if (vn1->data[nmi] != 0.0) {
          s = fabs(b_A->data[i + b_A->size[0] * nmi]) / vn1->data[nmi];
          y = s * s;
          s = 1.0 - s * s;
          if (1.0 - y < 0.0) {
            s = 0.0;
          }

          smax = vn1->data[nmi] / vn2->data[nmi];
          if (s * (smax * smax) <= 1.4901161193847656E-8) {
            if (i + 1 < A->size[0]) {
              y = 0.0;
              if (mmi < 1) {
              } else if (mmi == 1) {
                y = fabs(b_A->data[k]);
              } else {
                smax = 2.2250738585072014E-308;
                pvt = k + mmi;
                while (k + 1 <= pvt) {
                  s = fabs(b_A->data[k]);
                  if (s > smax) {
                    t = smax / s;
                    y = 1.0 + y * t * t;
                    smax = s;
                  } else {
                    t = s / smax;
                    y += t * t;
                  }

                  k++;
                }

                y = smax * sqrt(y);
              }

              vn1->data[nmi] = y;
              vn2->data[nmi] = vn1->data[nmi];
            } else {
              vn1->data[nmi] = 0.0;
              vn2->data[nmi] = 0.0;
            }
          } else {
            vn1->data[nmi] *= sqrt(s);
          }
        }
      }
    }

    emxFree_real_T(&vn2);
    emxFree_real_T(&vn1);
    emxFree_real_T(&work);
  }

  t = 0.0;
  if (mn > 0) {
    k = 0;
    exitg1 = FALSE;
    while ((exitg1 == FALSE) && (k <= mn - 1)) {
      smax = A->size[0];
      s = A->size[1];
      if (smax >= s) {
        y = smax;
      } else {
        y = s;
      }

      if (fabs(b_A->data[k + b_A->size[0] * k]) <= y * fabs(b_A->data[0]) *
          2.2204460492503131E-16) {
        exitg1 = TRUE;
      } else {
        t++;
        k++;
      }
    }
  }

  unnamed_idx_0 = (uint32_T)A->size[1];
  i1 = Y->size[0];
  Y->size[0] = (int32_T)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)Y, i1, (int32_T)sizeof(real_T));
  k = (int32_T)unnamed_idx_0;
  for (i1 = 0; i1 < k; i1++) {
    Y->data[i1] = 0.0;
  }

  for (nmi = 0; nmi < mn; nmi++) {
    if (tau->data[nmi] != 0.0) {
      smax = B->data[nmi];
      i1 = A->size[0] + (int32_T)(1.0 - ((1.0 + (real_T)nmi) + 1.0));
      for (i = 0; i < i1; i++) {
        unnamed_idx_0 = ((uint32_T)nmi + i) + 2U;
        smax += b_A->data[((int32_T)unnamed_idx_0 + b_A->size[0] * nmi) - 1] *
          B->data[(int32_T)unnamed_idx_0 - 1];
      }

      smax *= tau->data[nmi];
      if (smax != 0.0) {
        B->data[nmi] -= smax;
        i1 = A->size[0] + (int32_T)(1.0 - ((1.0 + (real_T)nmi) + 1.0));
        for (i = 0; i < i1; i++) {
          unnamed_idx_0 = ((uint32_T)nmi + i) + 2U;
          B->data[(int32_T)unnamed_idx_0 - 1] -= b_A->data[((int32_T)
            unnamed_idx_0 + b_A->size[0] * nmi) - 1] * smax;
        }
      }
    }
  }

  emxFree_real_T(&tau);
  for (i = 0; i < (int32_T)t; i++) {
    Y->data[jpvt->data[(int32_T)(1.0 + (real_T)i) - 1] - 1] = B->data[(int32_T)
      (1.0 + (real_T)i) - 1];
  }

  for (nmi = 0; nmi < (int32_T)-(1.0 + (-1.0 - t)); nmi++) {
    smax = t + -(real_T)nmi;
    Y->data[jpvt->data[(int32_T)smax - 1] - 1] /= b_A->data[((int32_T)smax +
      b_A->size[0] * ((int32_T)smax - 1)) - 1];
    for (i = 0; i <= (int32_T)smax - 2; i++) {
      Y->data[jpvt->data[(int32_T)(1.0 + (real_T)i) - 1] - 1] -= Y->data
        [jpvt->data[(int32_T)smax - 1] - 1] * b_A->data[((int32_T)(1.0 + (real_T)
        i) + b_A->size[0] * ((int32_T)smax - 1)) - 1];
    }
  }

  emxFree_int32_T(&jpvt);
  emxFree_real_T(&b_A);
}

static real_T eml_xnrm2(int32_T n, const emxArray_real_T *x, int32_T ix0)
{
  real_T y;
  real_T scale;
  int32_T kend;
  int32_T k;
  real_T absxk;
  real_T t;
  y = 0.0;
  if (n < 1) {
  } else if (n == 1) {
    y = fabs(x->data[ix0 - 1]);
  } else {
    scale = 2.2250738585072014E-308;
    kend = (ix0 + n) - 1;
    for (k = ix0; k <= kend; k++) {
      absxk = fabs(x->data[k - 1]);
      if (absxk > scale) {
        t = scale / absxk;
        y = 1.0 + y * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }
    }

    y = scale * sqrt(y);
  }

  return y;
}

static void eml_xscal(int32_T n, real_T a, emxArray_real_T *x, int32_T ix0)
{
  int32_T i4;
  int32_T k;
  i4 = (ix0 + n) - 1;
  for (k = ix0; k <= i4; k++) {
    x->data[k - 1] *= a;
  }
}

static real_T rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T a;
  real_T b;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = b * sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * sqrt(b * b + 1.0);
  } else if (rtIsNaN(b)) {
    y = b;
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

void mldivide(const emxArray_real_T *A, emxArray_real_T *B)
{
  emxArray_real_T *b_A;
  emxArray_int32_T *ipiv;
  emxArray_real_T *b_B;
  emxArray_real_T *r0;
  uint32_T unnamed_idx_0;
  int32_T i2;
  int32_T iy;
  int32_T u1;
  int32_T j;
  int32_T mmj;
  int32_T c;
  int32_T ix;
  real_T temp;
  real_T s;
  int32_T i3;
  int32_T jy;
  int32_T b_j;
  int32_T ijA;
  b_emxInit_real_T(&b_A, 2);
  emxInit_int32_T(&ipiv, 2);
  emxInit_real_T(&b_B, 1);
  emxInit_real_T(&r0, 1);
  if ((A->size[0] == 0) || (A->size[1] == 0) || (B->size[0] == 0)) {
    unnamed_idx_0 = (uint32_T)A->size[1];
    i2 = B->size[0];
    B->size[0] = (int32_T)unnamed_idx_0;
    emxEnsureCapacity((emxArray__common *)B, i2, (int32_T)sizeof(real_T));
    iy = (int32_T)unnamed_idx_0;
    for (i2 = 0; i2 < iy; i2++) {
      B->data[i2] = 0.0;
    }
  } else if (A->size[0] == A->size[1]) {
    i2 = b_A->size[0] * b_A->size[1];
    b_A->size[0] = A->size[0];
    b_A->size[1] = A->size[1];
    emxEnsureCapacity((emxArray__common *)b_A, i2, (int32_T)sizeof(real_T));
    iy = A->size[0] * A->size[1];
    for (i2 = 0; i2 < iy; i2++) {
      b_A->data[i2] = A->data[i2];
    }

    iy = A->size[1];
    u1 = A->size[1];
    if (iy <= u1) {
    } else {
      iy = u1;
    }

    eml_signed_integer_colon(iy, ipiv);
    iy = A->size[1] - 1;
    u1 = A->size[1];
    if (iy <= u1) {
      i2 = iy;
    } else {
      i2 = u1;
    }

    for (j = 0; j + 1 <= i2; j++) {
      mmj = A->size[1] - j;
      c = j * (A->size[1] + 1);
      if (mmj < 1) {
        iy = -1;
      } else {
        iy = 0;
        if (mmj > 1) {
          ix = c;
          temp = fabs(b_A->data[c]);
          for (u1 = 1; u1 + 1 <= mmj; u1++) {
            ix++;
            s = fabs(b_A->data[ix]);
            if (s > temp) {
              iy = u1;
              temp = s;
            }
          }
        }
      }

      if (b_A->data[c + iy] != 0.0) {
        if (iy != 0) {
          ipiv->data[j] = (j + iy) + 1;
          ix = j;
          iy += j;
          for (u1 = 1; u1 <= A->size[1]; u1++) {
            temp = b_A->data[ix];
            b_A->data[ix] = b_A->data[iy];
            b_A->data[iy] = temp;
            ix += A->size[1];
            iy += A->size[1];
          }
        }

        i3 = c + mmj;
        for (iy = c + 1; iy + 1 <= i3; iy++) {
          b_A->data[iy] /= b_A->data[c];
        }
      }

      u1 = (A->size[1] - j) - 1;
      iy = c + A->size[1];
      jy = c + A->size[1];
      for (b_j = 1; b_j <= u1; b_j++) {
        temp = -b_A->data[jy];
        if (b_A->data[jy] != 0.0) {
          ix = c + 1;
          i3 = mmj + iy;
          for (ijA = 1 + iy; ijA + 1 <= i3; ijA++) {
            b_A->data[ijA] += b_A->data[ix] * temp;
            ix++;
          }
        }

        jy += A->size[1];
        iy += A->size[1];
      }
    }

    for (iy = 0; iy + 1 <= A->size[1]; iy++) {
      if (ipiv->data[iy] != iy + 1) {
        temp = B->data[iy];
        B->data[iy] = B->data[ipiv->data[iy] - 1];
        B->data[ipiv->data[iy] - 1] = temp;
      }
    }

    for (u1 = 0; u1 + 1 <= A->size[1]; u1++) {
      c = A->size[1] * u1;
      if (B->data[u1] != 0.0) {
        for (iy = u1 + 1; iy + 1 <= A->size[1]; iy++) {
          B->data[iy] -= B->data[u1] * b_A->data[iy + c];
        }
      }
    }

    for (u1 = A->size[1] - 1; u1 + 1 > 0; u1--) {
      c = A->size[1] * u1;
      if (B->data[u1] != 0.0) {
        B->data[u1] /= b_A->data[u1 + c];
        for (iy = 0; iy + 1 <= u1; iy++) {
          B->data[iy] -= B->data[u1] * b_A->data[iy + c];
        }
      }
    }
  } else {
    i2 = b_B->size[0];
    b_B->size[0] = B->size[0];
    emxEnsureCapacity((emxArray__common *)b_B, i2, (int32_T)sizeof(real_T));
    iy = B->size[0];
    for (i2 = 0; i2 < iy; i2++) {
      b_B->data[i2] = B->data[i2];
    }

    eml_qrsolve(A, b_B, r0);
    i2 = B->size[0];
    B->size[0] = r0->size[0];
    emxEnsureCapacity((emxArray__common *)B, i2, (int32_T)sizeof(real_T));
    iy = r0->size[0];
    for (i2 = 0; i2 < iy; i2++) {
      B->data[i2] = r0->data[i2];
    }
  }

  emxFree_real_T(&r0);
  emxFree_real_T(&b_B);
  emxFree_int32_T(&ipiv);
  emxFree_real_T(&b_A);
}

/* End of code generation (mldivide.c) */
