/*
 * loadflow.c
 *
 * Code generation for function 'loadflow'
 *
 * C source code generated on: Mon Dec 15 20:15:41 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "matlabrun.h"
#include "loadflow.h"
#include "exp.h"
#include "matlabrun_emxutil.h"
#include "norm.h"
#include "sin.h"
#include "cos.h"
#include "mldivide.h"
#include "diag.h"
#include "rdivide.h"
#include <stdio.h>
#include "../constants.h"

static real_T rt_roundd_snf(real_T u);

static real_T rt_roundd_snf(real_T u)
{
  real_T y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = -0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

uint_T loadflow(const real_T *bus, const real_T *line, const uint_T noofgens,
    const uint_T SB, const uint_T noofslack,
    real_T *bus_sol)
{
  emxArray_creal_T *Y;
  int32_T i0 = 0;
  int32_T loop_ub;
  int32_T ixstart;
  real_T tap_ratio;
  boolean_T exitg1;
  emxArray_real_T *bus_int;
  emxArray_real_T *r0;
  creal_T b_line[LINENUM];
  emxArray_creal_T *y;
  real_T x;
  real_T b_x;
  real_T temp_im;
  real_T ai;
  creal_T temp;
  real_T y_re;
  real_T brm;
  real_T V_rect_im;
  real_T A_re;
  real_T V[BUSNUM];
  real_T ang[BUSNUM];
  real_T Pg[BUSNUM];
  real_T Qg[BUSNUM];
  real_T Pl[BUSNUM];
  real_T Ql[BUSNUM];
  real_T bus_type[BUSNUM];
  emxArray_real_T *PQV_no;
  emxArray_real_T *PQ_no;
  uint32_T PQVptr;
  uint32_T PQptr;
  emxArray_real_T *ang_red;
  emxArray_real_T *volt_red;
  real_T P[BUSNUM];
  real_T Q[BUSNUM];
  creal_T V_rect[BUSNUM];
  int32_T i1;
  int32_T ic;
  int32_T br;
  int32_T ar;
  int32_T ib;
  int32_T ia;
  real_T delP[BUSNUM];
  real_T delQ[BUSNUM];
  int32_T conv_flag;
  emxArray_creal_T *V_2;
  emxArray_creal_T *XX_1;
  emxArray_real_T *J11;
  emxArray_real_T *J12;
  emxArray_real_T *J21;
  emxArray_real_T *J22;
  emxArray_real_T *red_delQ;
  emxArray_real_T *b_temp;
  emxArray_real_T *delAng;
  emxArray_real_T *delV;
  emxArray_creal_T *b;
  emxArray_real_T *a;
  emxArray_real_T *b_J11;
  creal_T ang_pert[BUSNUM];
  creal_T V_pert[BUSNUM];
  creal_T X_1[BUSNUM*BUSNUM];
  boolean_T b_b;
  creal_T b_y[BUSNUM*BUSNUM];
  emxInit_creal_T(&Y, 2);

  Y->size[0] = (int32_T)BUSNUM;
  Y->size[1] = (int32_T)BUSNUM;
  emxEnsureCapacity((emxArray__common *)Y, i0, (int32_T)sizeof(creal_T));
  loop_ub = (int32_T)BUSNUM * (int32_T)BUSNUM;
  for (i0 = 0; i0 < loop_ub; i0++) {
    Y->data[i0].re = 0.0;
    Y->data[i0].im = 0.0;
  }

  ixstart = 1;
  tap_ratio = bus[0];
  if (rtIsNaN(bus[0])) {
    loop_ub = 2;
    exitg1 = FALSE;
    while ((exitg1 == FALSE) && (loop_ub < (BUSNUM+1))) {
      ixstart = loop_ub;
      if (!rtIsNaN(bus[loop_ub - 1])) {
        tap_ratio = bus[loop_ub - 1];
        exitg1 = TRUE;
      } else {
        loop_ub++;
      }
    }
  }

  if (ixstart < BUSNUM) {
    while (ixstart + 1 < (BUSNUM+1)) {
      if (bus[ixstart] > tap_ratio) {
        tap_ratio = bus[ixstart];
      }

      ixstart++;
    }
  }

  emxInit_real_T(&bus_int, 1);
  i0 = bus_int->size[0];
  bus_int->size[0] = (int32_T)tap_ratio;
  emxEnsureCapacity((emxArray__common *)bus_int, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)tap_ratio;
  for (i0 = 0; i0 < loop_ub; i0++) {
    bus_int->data[i0] = 0.0;
  }

  /*  ibus = [1:1:BUSNUM]'; */
  for (ixstart = 0; ixstart < (int32_T)BUSNUM; ixstart++) {
    bus_int->data[(int32_T)rt_roundd_snf(bus[(int32_T)(1.0 + (real_T)ixstart) -
        1]) - 1] = 1.0 + (real_T)ixstart;
  }

  emxInit_real_T(&r0, 1);

  /*  process line data and build admittance matrix Y */
  /*  line impedance */
  i0 = r0->size[0];
  r0->size[0] = (int32_T)LINENUM;
  emxEnsureCapacity((emxArray__common *)r0, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)LINENUM;
  for (i0 = 0; i0 < loop_ub; i0++) {
    r0->data[i0] = 1.0;
  }

  for (i0 = 0; i0 < LINENUM; i0++) {
    b_line[i0].re = line[LINENUM*2 + i0] + 0.0 * line[LINENUM*3 + i0];
    b_line[i0].im = line[LINENUM*3 + i0];
  }

  b_emxInit_creal_T(&y, 1);
  rdivide(r0, b_line, y, LINENUM);
  ixstart = 0;
  emxFree_real_T(&r0);
  while (ixstart <= (int32_T)LINENUM - 1) {
    x = rt_roundd_snf(line[(int32_T)(1.0 + (real_T)ixstart) - 1]);
    b_x = rt_roundd_snf(line[(int32_T)(1.0 + (real_T)ixstart) + (LINENUM-1)]);
    tap_ratio = line[(int32_T)(1.0 + (real_T)ixstart) + (LINENUM*5-1)];
    if (line[(int32_T)(1.0 + (real_T)ixstart) + (LINENUM*5-1)] == 0.0) {
      /*  this line has no transformer */
      tap_ratio = 1.0;
    }

    temp_im = line[(int32_T)(1.0 + (real_T)ixstart) + (LINENUM*6-1)] * 0.0 *
      3.1415926535897931;
    ai = line[(int32_T)(1.0 + (real_T)ixstart) + (LINENUM*6-1)] * 3.1415926535897931;
    if (ai == 0.0) {
      temp.re = temp_im / 180.0;
      temp.im = 0.0;
    } else if (temp_im == 0.0) {
      temp.re = 0.0;
      temp.im = ai / 180.0;
    } else {
      temp.re = temp_im / 180.0;
      temp.im = ai / 180.0;
    }

    b_exp(&temp);
    temp.re *= tap_ratio;
    temp.im *= tap_ratio;
    temp_im = y->data[(int32_T)(1.0 + (real_T)ixstart) - 1].re;
    ai = y->data[(int32_T)(1.0 + (real_T)ixstart) - 1].im;
    if (-temp.im == 0.0) {
      if (ai == 0.0) {
        y_re = temp_im / temp.re;
        temp_im = 0.0;
      } else if (temp_im == 0.0) {
        y_re = 0.0;
        temp_im = ai / temp.re;
      } else {
        y_re = temp_im / temp.re;
        temp_im = ai / temp.re;
      }
    } else if (temp.re == 0.0) {
      if (temp_im == 0.0) {
        y_re = ai / -temp.im;
        temp_im = 0.0;
      } else if (ai == 0.0) {
        y_re = 0.0;
        temp_im = -(temp_im / -temp.im);
      } else {
        y_re = ai / -temp.im;
        temp_im = -(temp_im / -temp.im);
      }
    } else {
      brm = fabs(temp.re);
      V_rect_im = fabs(-temp.im);
      if (brm > V_rect_im) {
        brm = -temp.im / temp.re;
        V_rect_im = temp.re + brm * -temp.im;
        y_re = (temp_im + brm * ai) / V_rect_im;
        temp_im = (ai - brm * temp_im) / V_rect_im;
      } else if (V_rect_im == brm) {
        if (temp.re > 0.0) {
          tap_ratio = 0.5;
        } else {
          tap_ratio = -0.5;
        }

        if (-temp.im > 0.0) {
          V_rect_im = 0.5;
        } else {
          V_rect_im = -0.5;
        }

        y_re = (temp_im * tap_ratio + ai * V_rect_im) / brm;
        temp_im = (ai * tap_ratio - temp_im * V_rect_im) / brm;
      } else {
        brm = temp.re / -temp.im;
        V_rect_im = -temp.im + brm * temp.re;
        y_re = (brm * temp_im + ai) / V_rect_im;
        temp_im = (brm * ai - temp_im) / V_rect_im;
      }
    }

    Y->data[((int32_T)bus_int->data[(int32_T)x - 1] + Y->size[0] * ((int32_T)
          bus_int->data[(int32_T)b_x - 1] - 1)) - 1].re -= y_re;
    Y->data[((int32_T)bus_int->data[(int32_T)x - 1] + Y->size[0] * ((int32_T)
          bus_int->data[(int32_T)b_x - 1] - 1)) - 1].im -= temp_im;
    temp_im = y->data[(int32_T)(1.0 + (real_T)ixstart) - 1].re;
    ai = y->data[(int32_T)(1.0 + (real_T)ixstart) - 1].im;
    if (temp.im == 0.0) {
      if (ai == 0.0) {
        y_re = temp_im / temp.re;
        temp_im = 0.0;
      } else if (temp_im == 0.0) {
        y_re = 0.0;
        temp_im = ai / temp.re;
      } else {
        y_re = temp_im / temp.re;
        temp_im = ai / temp.re;
      }
    } else if (temp.re == 0.0) {
      if (temp_im == 0.0) {
        y_re = ai / temp.im;
        temp_im = 0.0;
      } else if (ai == 0.0) {
        y_re = 0.0;
        temp_im = -(temp_im / temp.im);
      } else {
        y_re = ai / temp.im;
        temp_im = -(temp_im / temp.im);
      }
    } else {
      brm = fabs(temp.re);
      V_rect_im = fabs(temp.im);
      if (brm > V_rect_im) {
        brm = temp.im / temp.re;
        V_rect_im = temp.re + brm * temp.im;
        y_re = (temp_im + brm * ai) / V_rect_im;
        temp_im = (ai - brm * temp_im) / V_rect_im;
      } else if (V_rect_im == brm) {
        if (temp.re > 0.0) {
          tap_ratio = 0.5;
        } else {
          tap_ratio = -0.5;
        }

        if (temp.im > 0.0) {
          V_rect_im = 0.5;
        } else {
          V_rect_im = -0.5;
        }

        y_re = (temp_im * tap_ratio + ai * V_rect_im) / brm;
        temp_im = (ai * tap_ratio - temp_im * V_rect_im) / brm;
      } else {
        brm = temp.re / temp.im;
        V_rect_im = temp.im + brm * temp.re;
        y_re = (brm * temp_im + ai) / V_rect_im;
        temp_im = (brm * ai - temp_im) / V_rect_im;
      }
    }

    Y->data[((int32_T)bus_int->data[(int32_T)b_x - 1] + Y->size[0] * ((int32_T)
          bus_int->data[(int32_T)x - 1] - 1)) - 1].re -= y_re;
    Y->data[((int32_T)bus_int->data[(int32_T)b_x - 1] + Y->size[0] * ((int32_T)
          bus_int->data[(int32_T)x - 1] - 1)) - 1].im -= temp_im;
    temp_im = 0.0 * line[(int32_T)(1.0 + (real_T)ixstart) + (LINENUM*4-1)];
    ai = line[(int32_T)(1.0 + (real_T)ixstart) + (LINENUM*4-1)];
    if (ai == 0.0) {
      V_rect_im = temp_im / 2.0;
      tap_ratio = 0.0;
    } else if (temp_im == 0.0) {
      V_rect_im = 0.0;
      tap_ratio = ai / 2.0;
    } else {
      V_rect_im = temp_im / 2.0;
      tap_ratio = ai / 2.0;
    }

    A_re = y->data[(int32_T)(1.0 + (real_T)ixstart) - 1].re + V_rect_im;
    y_re = y->data[(int32_T)(1.0 + (real_T)ixstart) - 1].im + tap_ratio;
    tap_ratio = temp.re * temp.re - temp.im * -temp.im;
    temp_im = temp.re * -temp.im + temp.im * temp.re;
    if (temp_im == 0.0) {
      if (y_re == 0.0) {
        ai = A_re / tap_ratio;
        y_re = 0.0;
      } else if (A_re == 0.0) {
        ai = 0.0;
        y_re /= tap_ratio;
      } else {
        ai = A_re / tap_ratio;
        y_re /= tap_ratio;
      }
    } else if (tap_ratio == 0.0) {
      if (A_re == 0.0) {
        ai = y_re / temp_im;
        y_re = 0.0;
      } else if (y_re == 0.0) {
        ai = 0.0;
        y_re = -(A_re / temp_im);
      } else {
        ai = y_re / temp_im;
        y_re = -(A_re / temp_im);
      }
    } else {
      brm = fabs(tap_ratio);
      V_rect_im = fabs(temp_im);
      if (brm > V_rect_im) {
        brm = temp_im / tap_ratio;
        V_rect_im = tap_ratio + brm * temp_im;
        ai = (A_re + brm * y_re) / V_rect_im;
        y_re = (y_re - brm * A_re) / V_rect_im;
      } else if (V_rect_im == brm) {
        if (tap_ratio > 0.0) {
          tap_ratio = 0.5;
        } else {
          tap_ratio = -0.5;
        }

        if (temp_im > 0.0) {
          V_rect_im = 0.5;
        } else {
          V_rect_im = -0.5;
        }

        ai = (A_re * tap_ratio + y_re * V_rect_im) / brm;
        y_re = (y_re * tap_ratio - A_re * V_rect_im) / brm;
      } else {
        brm = tap_ratio / temp_im;
        V_rect_im = temp_im + brm * tap_ratio;
        ai = (brm * A_re + y_re) / V_rect_im;
        y_re = (brm * y_re - A_re) / V_rect_im;
      }
    }

    Y->data[((int32_T)bus_int->data[(int32_T)x - 1] + Y->size[0] * ((int32_T)
          bus_int->data[(int32_T)x - 1] - 1)) - 1].re += ai;
    Y->data[((int32_T)bus_int->data[(int32_T)x - 1] + Y->size[0] * ((int32_T)
          bus_int->data[(int32_T)x - 1] - 1)) - 1].im += y_re;
    temp_im = 0.0 * line[(int32_T)(1.0 + (real_T)ixstart) + (LINENUM*4-1)];
    ai = line[(int32_T)(1.0 + (real_T)ixstart) + (LINENUM*4-1)];
    if (ai == 0.0) {
      V_rect_im = temp_im / 2.0;
      tap_ratio = 0.0;
    } else if (temp_im == 0.0) {
      V_rect_im = 0.0;
      tap_ratio = ai / 2.0;
    } else {
      V_rect_im = temp_im / 2.0;
      tap_ratio = ai / 2.0;
    }

    Y->data[((int32_T)bus_int->data[(int32_T)b_x - 1] + Y->size[0] * ((int32_T)
          bus_int->data[(int32_T)b_x - 1] - 1)) - 1].re = (Y->data[((int32_T)
          bus_int->data[(int32_T)b_x - 1] + Y->size[0] * ((int32_T)bus_int->data
            [(int32_T)b_x - 1] - 1)) - 1].re + y->data[(int32_T)(1.0 + (real_T)ixstart)
        - 1].re) + V_rect_im;
    Y->data[((int32_T)bus_int->data[(int32_T)b_x - 1] + Y->size[0] * ((int32_T)
          bus_int->data[(int32_T)b_x - 1] - 1)) - 1].im = (Y->data[((int32_T)
          bus_int->data[(int32_T)b_x - 1] + Y->size[0] * ((int32_T)bus_int->data
            [(int32_T)b_x - 1] - 1)) - 1].im + y->data[(int32_T)(1.0 + (real_T)ixstart)
        - 1].im) + tap_ratio;
    ixstart++;
  }

  /* %%%%%ybus.m ends */
  /*  process bus data */
  for (i0 = 0; i0 < BUSNUM; i0++) {
    V[i0] = bus[BUSNUM + i0];
    ang[i0] = bus[BUSNUM*2 + i0] * 3.1415926535897931 / 180.0;
    Pg[i0] = bus[BUSNUM*3 + i0];
    Qg[i0] = bus[BUSNUM*4 + i0];
    Pl[i0] = bus[BUSNUM*5 + i0];
    Ql[i0] = bus[BUSNUM*6 + i0];

    /* cyb = (Gb + jay*Bb)'; */
    bus_type[i0] = bus[BUSNUM*9 + i0];
  }

  b_emxInit_real_T(&PQV_no, 2);

  i0 = PQV_no->size[0] * PQV_no->size[1];
  PQV_no->size[0] = 1;
  PQV_no->size[1] = (int32_T)(BUSNUM - noofslack);
  emxEnsureCapacity((emxArray__common *)PQV_no, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)(BUSNUM - noofslack);
  for (i0 = 0; i0 < loop_ub; i0++) {
    PQV_no->data[i0] = 0.0;
  }

  b_emxInit_real_T(&PQ_no, 2);
  i0 = PQ_no->size[0] * PQ_no->size[1];
  PQ_no->size[0] = 1;
  PQ_no->size[1] = (int32_T)(BUSNUM - (noofgens + noofslack));
  emxEnsureCapacity((emxArray__common *)PQ_no, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)(BUSNUM - (noofgens + noofslack));
  for (i0 = 0; i0 < loop_ub; i0++) {
    PQ_no->data[i0] = 0.0;
  }

  PQVptr = 0U;

  /*  PQV_no pointer */
  PQptr = 0U;

  /*  PQ_no pointer */
  for (ixstart = 0; ixstart < (int32_T)BUSNUM; ixstart++) {
    if (bus_type[ixstart] == 3.0) {
      PQptr++;
      PQVptr++;
      PQV_no->data[(int32_T)PQVptr - 1] = 1.0 + (real_T)ixstart;
      PQ_no->data[(int32_T)PQptr - 1] = 1.0 + (real_T)ixstart;
    } else {
      if (bus_type[(int32_T)(1.0 + (real_T)ixstart) - 1] == 2.0) {
        PQVptr++;
        PQV_no->data[(int32_T)PQVptr - 1] = 1.0 + (real_T)ixstart;
      }
    }
  }

  b_emxInit_real_T(&ang_red, 2);

  /* % */
  /*  construct angle reduction matrix */
  ixstart = PQV_no->size[1];
  i0 = ang_red->size[0] * ang_red->size[1];
  ang_red->size[0] = ixstart;
  emxEnsureCapacity((emxArray__common *)ang_red, i0, (int32_T)sizeof(real_T));
  i0 = ang_red->size[0] * ang_red->size[1];
  ang_red->size[1] = (int32_T)BUSNUM;
  emxEnsureCapacity((emxArray__common *)ang_red, i0, (int32_T)sizeof(real_T));
  loop_ub = PQV_no->size[1] * (int32_T)BUSNUM;
  for (i0 = 0; i0 < loop_ub; i0++) {
    ang_red->data[i0] = 0.0;
  }

  for (ixstart = 0; ixstart < PQV_no->size[1]; ixstart++) {
    ang_red->data[ixstart + ang_red->size[0] * ((int32_T)PQV_no->data[ixstart] -
        1)] = 1.0;
  }

  b_emxInit_real_T(&volt_red, 2);

  /*  */
  /* il = length(PQV_no); */
  /* ii = 1:1:il; */
  /* ang_red = sparse(ii,PQV_no,ones(il,1),il,BUSNUM); */
  /*  construct voltage reduction matrix */
  ixstart = PQ_no->size[1];
  i0 = volt_red->size[0] * volt_red->size[1];
  volt_red->size[0] = ixstart;
  emxEnsureCapacity((emxArray__common *)volt_red, i0, (int32_T)sizeof(real_T));
  i0 = volt_red->size[0] * volt_red->size[1];
  volt_red->size[1] = (int32_T)BUSNUM;
  emxEnsureCapacity((emxArray__common *)volt_red, i0, (int32_T)sizeof(real_T));
  loop_ub = PQ_no->size[1] * (int32_T)BUSNUM;
  for (i0 = 0; i0 < loop_ub; i0++) {
    volt_red->data[i0] = 0.0;
  }

  for (ixstart = 0; ixstart < PQ_no->size[1]; ixstart++) {
    volt_red->data[ixstart + volt_red->size[0] * ((int32_T)PQ_no->data[ixstart]
        - 1)] = 1.0;
  }

  /*  */
  /* il = length(PQ_no); */
  /* ii = 1:1:il; */
  /* volt_red = sparse(ii,PQ_no,ones(il,1),il,BUSNUM); */
  uint_T iter = 0;

  /*  initialize iteration counter */
  /*  calculate the power mismatch and check convergence */
  /*   [delP,delQ,P,Q,conv_flag] =... */
  /*             calc(BUSNUM,bus_type,V,ang,Y,Pg,Qg,Pl,Ql,tol); */
  /* %%%%% calc.m */
  /*  voltage in rectangular coordinate */
  memcpy(&P[0], &ang[0], BUSNUM * sizeof(real_T));
  b_cos(P, BUSNUM);
  memcpy(&Q[0], &ang[0], BUSNUM * sizeof(real_T));
  b_sin(Q, BUSNUM);
  for (i0 = 0; i0 < BUSNUM; i0++) {
    V_rect[i0].re = V[i0] * (P[i0] + 0.0 * Q[i0]);
    V_rect[i0].im = V[i0] * Q[i0];
  }

  /*  bus current injection */
  if (Y->size[1] == 1) {
    i0 = y->size[0];
    y->size[0] = Y->size[0];
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(creal_T));
    loop_ub = Y->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      y->data[i0].re = 0.0;
      y->data[i0].im = 0.0;
      ixstart = Y->size[1];
      for (i1 = 0; i1 < ixstart; i1++) {
        tap_ratio = Y->data[i0 + Y->size[0] * i1].re * V_rect[i1].re - Y->
          data[i0 + Y->size[0] * i1].im * V_rect[i1].im;
        V_rect_im = Y->data[i0 + Y->size[0] * i1].re * V_rect[i1].im + Y->
          data[i0 + Y->size[0] * i1].im * V_rect[i1].re;
        y->data[i0].re += tap_ratio;
        y->data[i0].im += V_rect_im;
      }
    }
  } else {
    PQVptr = (uint32_T)Y->size[0];
    i0 = y->size[0];
    y->size[0] = (int32_T)PQVptr;
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(creal_T));
    loop_ub = (int32_T)PQVptr;
    for (i0 = 0; i0 < loop_ub; i0++) {
      y->data[i0].re = 0.0;
      y->data[i0].im = 0.0;
    }

    if (Y->size[0] == 0) {
    } else {
      loop_ub = 0;
      while ((Y->size[0] > 0) && (loop_ub <= 0)) {
        i0 = Y->size[0];
        for (ic = 1; ic <= i0; ic++) {
          y->data[ic - 1].re = 0.0;
          y->data[ic - 1].im = 0.0;
        }

        loop_ub = Y->size[0];
      }

      br = 0;
      loop_ub = 0;
      while ((Y->size[0] > 0) && (loop_ub <= 0)) {
        ar = 0;
        i0 = br + Y->size[1];
        for (ib = br; ib + 1 <= i0; ib++) {
          if ((V_rect[ib].re != 0.0) || (V_rect[ib].im != 0.0)) {
            temp.re = V_rect[ib].re - 0.0 * V_rect[ib].im;
            temp.im = V_rect[ib].im + 0.0 * V_rect[ib].re;
            ia = ar;
            i1 = Y->size[0];
            for (ic = 0; ic + 1 <= i1; ic++) {
              ia++;
              tap_ratio = temp.re * Y->data[ia - 1].re - temp.im * Y->data[ia -
                1].im;
              temp_im = temp.re * Y->data[ia - 1].im + temp.im * Y->data[ia - 1]
                .re;
              y->data[ic].re += tap_ratio;
              y->data[ic].im += temp_im;
            }
          }

          ar += Y->size[0];
        }

        br += Y->size[1];
        loop_ub = Y->size[0];
      }
    }
  }

  /*  power output  */
  for (i0 = 0; i0 < BUSNUM; i0++) {
    y_re = y->data[i0].re;
    temp_im = -y->data[i0].im;
    tap_ratio = V_rect[i0].re;
    V_rect_im = V_rect[i0].im;
    V_rect[i0].re = V_rect[i0].re * y_re - V_rect[i0].im * temp_im;
    V_rect[i0].im = tap_ratio * temp_im + V_rect_im * y_re;
  }

  for (ixstart = 0; ixstart < BUSNUM; ixstart++) {
    delP[ixstart] = (Pg[ixstart] - Pl[ixstart]) - V_rect[ixstart].re;
    delQ[ixstart] = (Qg[ixstart] - Ql[ixstart]) - V_rect[ixstart].im;
    P[ixstart] = V_rect[ixstart].re;
    Q[ixstart] = V_rect[ixstart].im;
  }

  /*  zero out mismatches on swing bus and generation bus */
  for (ixstart = 0; ixstart < (int32_T)BUSNUM; ixstart++) {
    if (bus_type[ixstart] == 1.0) {
      delP[(int32_T)(1.0 + (real_T)ixstart) - 1] = 0.0;
      delQ[(int32_T)(1.0 + (real_T)ixstart) - 1] = 0.0;
    } else {
      if (bus_type[(int32_T)(1.0 + (real_T)ixstart) - 1] == 2.0) {
        delQ[(int32_T)(1.0 + (real_T)ixstart) - 1] = 0.0;
      }
    }
  }

  /*   total mismatch */
  if (norm(delQ, BUSNUM) + norm(delP, BUSNUM) > 1.0E-11) {
    conv_flag = 1;
  } else {
    conv_flag = 0;
  }

  /* %%%%% calc.m ends */
  /* % start iteration process */
  emxInit_creal_T(&V_2, 2);
  emxInit_creal_T(&XX_1, 2);
  b_emxInit_real_T(&J11, 2);
  b_emxInit_real_T(&J12, 2);
  b_emxInit_real_T(&J21, 2);
  b_emxInit_real_T(&J22, 2);
  emxInit_real_T(&red_delQ, 1);
  emxInit_real_T(&b_temp, 1);
  b_emxInit_real_T(&delAng, 2);
  b_emxInit_real_T(&delV, 2);
  emxInit_creal_T(&b, 2);
  b_emxInit_real_T(&a, 2);
  b_emxInit_real_T(&b_J11, 2);
  while ((conv_flag == 1) && (iter < 100.0)) {
    (iter)++;

    /* %%% form_jac.m */
    /* [k,dum] = size(Y); */
    memcpy(&P[0], &ang[0], BUSNUM * sizeof(real_T));
    b_cos(P, BUSNUM);
    memcpy(&Q[0], &ang[0], BUSNUM * sizeof(real_T));
    b_sin(Q, BUSNUM);

    /*  voltage perturbation rectangular coordinates */
    for (ixstart = 0; ixstart < BUSNUM; ixstart++) {
      V_pert[ixstart].re = P[ixstart] + 0.0 * Q[ixstart];
      V_pert[ixstart].im = Q[ixstart];

      /*  Voltage rectangular coordinates */
      V_rect[ixstart].re = V[ixstart] * V_pert[ixstart].re;
      V_rect[ixstart].im = V[ixstart] * V_pert[ixstart].im;

      /*  angle and voltage perturbation rectangular coordinates */
      ang_pert[ixstart].re = -V[ixstart] * (Q[ixstart] - 0.0 * P[ixstart]);
      ang_pert[ixstart].im = -V[ixstart] * (0.0 - P[ixstart]);
    }

    if (Y->size[1] == 1) {
      i0 = y->size[0];
      y->size[0] = Y->size[0];
      emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(creal_T));
      loop_ub = Y->size[0];
      for (i0 = 0; i0 < loop_ub; i0++) {
        y->data[i0].re = 0.0;
        y->data[i0].im = 0.0;
        ixstart = Y->size[1];
        for (i1 = 0; i1 < ixstart; i1++) {
          tap_ratio = Y->data[i0 + Y->size[0] * i1].re * V_rect[i1].re - Y->
            data[i0 + Y->size[0] * i1].im * V_rect[i1].im;
          V_rect_im = Y->data[i0 + Y->size[0] * i1].re * V_rect[i1].im + Y->
            data[i0 + Y->size[0] * i1].im * V_rect[i1].re;
          y->data[i0].re += tap_ratio;
          y->data[i0].im += V_rect_im;
        }
      }
    } else {
      PQVptr = (uint32_T)Y->size[0];
      i0 = y->size[0];
      y->size[0] = (int32_T)PQVptr;
      emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(creal_T));
      loop_ub = (int32_T)PQVptr;
      for (i0 = 0; i0 < loop_ub; i0++) {
        y->data[i0].re = 0.0;
        y->data[i0].im = 0.0;
      }

      if (Y->size[0] == 0) {
      } else {
        loop_ub = 0;
        while ((Y->size[0] > 0) && (loop_ub <= 0)) {
          i0 = Y->size[0];
          for (ic = 1; ic <= i0; ic++) {
            y->data[ic - 1].re = 0.0;
            y->data[ic - 1].im = 0.0;
          }

          loop_ub = Y->size[0];
        }

        br = 0;
        loop_ub = 0;
        while ((Y->size[0] > 0) && (loop_ub <= 0)) {
          ar = 0;
          i0 = br + Y->size[1];
          for (ib = br; ib + 1 <= i0; ib++) {
            if ((V_rect[ib].re != 0.0) || (V_rect[ib].im != 0.0)) {
              temp.re = V_rect[ib].re - 0.0 * V_rect[ib].im;
              temp.im = V_rect[ib].im + 0.0 * V_rect[ib].re;
              ia = ar;
              i1 = Y->size[0];
              for (ic = 0; ic + 1 <= i1; ic++) {
                ia++;
                tap_ratio = temp.re * Y->data[ia - 1].re - temp.im * Y->data[ia
                  - 1].im;
                temp_im = temp.re * Y->data[ia - 1].im + temp.im * Y->data[ia -
                  1].re;
                y->data[ic].re += tap_ratio;
                y->data[ic].im += temp_im;
              }
            }

            ar += Y->size[0];
          }

          br += Y->size[1];
          loop_ub = Y->size[0];
        }
      }
    }

    i0 = y->size[0];
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(creal_T));
    loop_ub = y->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      y->data[i0].im = -y->data[i0].im;
    }

    diag(V_rect, X_1, BUSNUM);
    i0 = b->size[0] * b->size[1];
    b->size[0] = Y->size[0];
    b->size[1] = Y->size[1];
    emxEnsureCapacity((emxArray__common *)b, i0, (int32_T)sizeof(creal_T));
    loop_ub = Y->size[0] * Y->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b->data[i0].re = Y->data[i0].re;
      b->data[i0].im = -Y->data[i0].im;
    }

    if (b->size[0] == 1) {
      i0 = V_2->size[0] * V_2->size[1];
      V_2->size[0] = BUSNUM;
      V_2->size[1] = b->size[1];
      emxEnsureCapacity((emxArray__common *)V_2, i0, (int32_T)sizeof(creal_T));
      for (i0 = 0; i0 < BUSNUM; i0++) {
        loop_ub = b->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          V_2->data[i0 + V_2->size[0] * i1].re = 0.0;
          V_2->data[i0 + V_2->size[0] * i1].im = 0.0;
          for (ixstart = 0; ixstart < BUSNUM; ixstart++) {
            tap_ratio = X_1[i0 + BUSNUM * ixstart].re * b->data[ixstart + b->size[0]
              * i1].re - X_1[i0 + BUSNUM * ixstart].im * b->data[ixstart + b->size[0]
              * i1].im;
            V_rect_im = X_1[i0 + BUSNUM * ixstart].re * b->data[ixstart + b->size[0]
              * i1].im + X_1[i0 + BUSNUM * ixstart].im * b->data[ixstart + b->size[0]
              * i1].re;
            V_2->data[i0 + V_2->size[0] * i1].re += tap_ratio;
            V_2->data[i0 + V_2->size[0] * i1].im += V_rect_im;
          }
        }
      }
    } else {
      PQVptr = (uint32_T)b->size[1];
      i0 = V_2->size[0] * V_2->size[1];
      V_2->size[0] = BUSNUM;
      emxEnsureCapacity((emxArray__common *)V_2, i0, (int32_T)sizeof(creal_T));
      i0 = V_2->size[0] * V_2->size[1];
      V_2->size[1] = (int32_T)PQVptr;
      emxEnsureCapacity((emxArray__common *)V_2, i0, (int32_T)sizeof(creal_T));
      loop_ub = BUSNUM * (int32_T)PQVptr;
      for (i0 = 0; i0 < loop_ub; i0++) {
        V_2->data[i0].re = 0.0;
        V_2->data[i0].im = 0.0;
      }

      if (b->size[1] == 0) {
      } else {
        ixstart = BUSNUM * (b->size[1] - 1);
        for (loop_ub = 0; loop_ub <= ixstart; loop_ub += BUSNUM) {
          for (ic = loop_ub; ic + 1 <= loop_ub + BUSNUM; ic++) {
            V_2->data[ic].re = 0.0;
            V_2->data[ic].im = 0.0;
          }
        }

        br = 0;
        for (loop_ub = 0; loop_ub <= ixstart; loop_ub += BUSNUM) {
          ar = 0;
          for (ib = br; ib + 1 <= br + BUSNUM; ib++) {
            b_b = ((b->data[ib].re != 0.0) || (b->data[ib].im != 0.0));
            if (b_b) {
              temp.re = b->data[ib].re - 0.0 * b->data[ib].im;
              temp.im = b->data[ib].im + 0.0 * b->data[ib].re;
              ia = ar;
              for (ic = loop_ub; ic + 1 <= loop_ub + BUSNUM; ic++) {
                ia++;
                V_2->data[ic].re += temp.re * X_1[ia - 1].re - temp.im * X_1[ia
                  - 1].im;
                V_2->data[ic].im += temp.re * X_1[ia - 1].im + temp.im * X_1[ia
                  - 1].re;
              }
            }

            ar += BUSNUM;
          }

          br += BUSNUM;
        }
      }
    }

    for (ixstart = 0; ixstart < BUSNUM; ixstart++) {
      V_rect[ixstart].re = ang_pert[ixstart].re;
      V_rect[ixstart].im = -ang_pert[ixstart].im;
    }

    diag(V_rect, X_1, BUSNUM);
    if (V_2->size[1] == 1) {
      for (i0 = 0; i0 < BUSNUM; i0++) {
        for (i1 = 0; i1 < BUSNUM; i1++) {
          b_y[i0 + BUSNUM * i1].re = 0.0;
          b_y[i0 + BUSNUM * i1].im = 0.0;
          for (ixstart = 0; ixstart < BUSNUM; ixstart++) {
            tap_ratio = V_2->data[i0 + BUSNUM * ixstart].re * X_1[ixstart + BUSNUM * i1]
              .re - V_2->data[i0 + BUSNUM * ixstart].im * X_1[ixstart + BUSNUM * i1].im;
            V_rect_im = V_2->data[i0 + BUSNUM * ixstart].re * X_1[ixstart + BUSNUM * i1]
              .im + V_2->data[i0 + BUSNUM * ixstart].im * X_1[ixstart + BUSNUM * i1].re;
            b_y[i0 + BUSNUM * i1].re += tap_ratio;
            b_y[i0 + BUSNUM * i1].im += V_rect_im;
          }
        }
      }
    } else {
      for (i0 = 0; i0 < BUSNUM*BUSNUM; i0++) {
        b_y[i0].re = 0.0;
        b_y[i0].im = 0.0;
      }

      for (loop_ub = 0; loop_ub < (BUSNUM*BUSNUM - BUSNUM + 2); loop_ub += BUSNUM) {
        for (ic = loop_ub; ic + 1 <= loop_ub + BUSNUM; ic++) {
          b_y[ic].re = 0.0;
          b_y[ic].im = 0.0;
        }
      }

      br = 0;
      for (loop_ub = 0; loop_ub < (BUSNUM*BUSNUM - BUSNUM + 2); loop_ub += BUSNUM) {
        ar = 0;
        i0 = br + V_2->size[1];
        for (ib = br; ib + 1 <= i0; ib++) {
          if ((X_1[ib].re != 0.0) || (X_1[ib].im != 0.0)) {
            temp.re = X_1[ib].re - 0.0 * X_1[ib].im;
            temp.im = X_1[ib].im + 0.0 * X_1[ib].re;
            ia = ar;
            for (ic = loop_ub; ic + 1 <= loop_ub + BUSNUM; ic++) {
              ia++;
              tap_ratio = temp.re * V_2->data[ia - 1].re - temp.im * V_2->
                data[ia - 1].im;
              temp_im = temp.re * V_2->data[ia - 1].im + temp.im * V_2->data[ia
                - 1].re;
              b_y[ic].re += tap_ratio;
              b_y[ic].im += temp_im;
            }
          }

          ar += BUSNUM;
        }

        br += V_2->size[1];
      }
    }

    for (i0 = 0; i0 < BUSNUM; i0++) {
      V_rect[i0].re = y->data[i0].re * ang_pert[i0].re - y->data[i0].im *
        ang_pert[i0].im;
      V_rect[i0].im = y->data[i0].re * ang_pert[i0].im + y->data[i0].im *
        ang_pert[i0].re;
    }

    diag(V_rect, X_1, BUSNUM);
    for (i0 = 0; i0 < BUSNUM*BUSNUM; i0++) {
      X_1[i0].re += b_y[i0].re;
      X_1[i0].im += b_y[i0].im;
    }

    i0 = XX_1->size[0] * XX_1->size[1];
    XX_1->size[0] = BUSNUM;
    emxEnsureCapacity((emxArray__common *)XX_1, i0, (int32_T)sizeof(creal_T));
    ixstart = PQV_no->size[1];
    i0 = XX_1->size[0] * XX_1->size[1];
    XX_1->size[1] = ixstart;
    emxEnsureCapacity((emxArray__common *)XX_1, i0, (int32_T)sizeof(creal_T));
    loop_ub = BUSNUM * PQV_no->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      XX_1->data[i0].re = 0.0;
      XX_1->data[i0].im = 0.0;
    }

    for (ixstart = 0; ixstart < PQV_no->size[1]; ixstart++) {
      i0 = (int32_T)rt_roundd_snf(PQV_no->data[ixstart]);
      for (i1 = 0; i1 < BUSNUM; i1++) {
        XX_1->data[i1 + XX_1->size[0] * ixstart] = X_1[i1 + BUSNUM * (i0 - 1)];
      }
    }

    for (ixstart = 0; ixstart < BUSNUM; ixstart++) {
      V_rect[ixstart].re = V_pert[ixstart].re;
      V_rect[ixstart].im = -V_pert[ixstart].im;
    }

    diag(V_rect, X_1, BUSNUM);
    if (V_2->size[1] == 1) {
      for (i0 = 0; i0 < BUSNUM; i0++) {
        for (i1 = 0; i1 < BUSNUM; i1++) {
          b_y[i0 + BUSNUM * i1].re = 0.0;
          b_y[i0 + BUSNUM * i1].im = 0.0;
          for (ixstart = 0; ixstart < BUSNUM; ixstart++) {
            tap_ratio = V_2->data[i0 + BUSNUM * ixstart].re * X_1[ixstart + BUSNUM * i1]
              .re - V_2->data[i0 + BUSNUM * ixstart].im * X_1[ixstart + BUSNUM * i1].im;
            V_rect_im = V_2->data[i0 + BUSNUM * ixstart].re * X_1[ixstart + BUSNUM * i1]
              .im + V_2->data[i0 + BUSNUM * ixstart].im * X_1[ixstart + BUSNUM * i1].re;
            b_y[i0 + BUSNUM * i1].re += tap_ratio;
            b_y[i0 + BUSNUM * i1].im += V_rect_im;
          }
        }
      }
    } else {
      for (i0 = 0; i0 < BUSNUM*BUSNUM; i0++) {
        b_y[i0].re = 0.0;
        b_y[i0].im = 0.0;
      }

      for (loop_ub = 0; loop_ub < (BUSNUM*BUSNUM - BUSNUM + 2); loop_ub += BUSNUM) {
        for (ic = loop_ub; ic + 1 <= loop_ub + BUSNUM; ic++) {
          b_y[ic].re = 0.0;
          b_y[ic].im = 0.0;
        }
      }

      br = 0;
      for (loop_ub = 0; loop_ub < (BUSNUM*BUSNUM - BUSNUM + 2); loop_ub += BUSNUM) {
        ar = 0;
        i0 = br + V_2->size[1];
        for (ib = br; ib + 1 <= i0; ib++) {
          if ((X_1[ib].re != 0.0) || (X_1[ib].im != 0.0)) {
            temp.re = X_1[ib].re - 0.0 * X_1[ib].im;
            temp.im = X_1[ib].im + 0.0 * X_1[ib].re;
            ia = ar;
            for (ic = loop_ub; ic + 1 <= loop_ub + BUSNUM; ic++) {
              ia++;
              tap_ratio = temp.re * V_2->data[ia - 1].re - temp.im * V_2->
                data[ia - 1].im;
              temp_im = temp.re * V_2->data[ia - 1].im + temp.im * V_2->data[ia
                - 1].re;
              b_y[ic].re += tap_ratio;
              b_y[ic].im += temp_im;
            }
          }

          ar += BUSNUM;
        }

        br += V_2->size[1];
      }
    }

    for (i0 = 0; i0 < BUSNUM; i0++) {
      V_rect[i0].re = y->data[i0].re * V_pert[i0].re - y->data[i0].im *
        V_pert[i0].im;
      V_rect[i0].im = y->data[i0].re * V_pert[i0].im + y->data[i0].im *
        V_pert[i0].re;
    }

    diag(V_rect, X_1, BUSNUM);
    for (i0 = 0; i0 < BUSNUM*BUSNUM; i0++) {
      X_1[i0].re += b_y[i0].re;
      X_1[i0].im += b_y[i0].im;
    }

    i0 = V_2->size[0] * V_2->size[1];
    V_2->size[0] = BUSNUM;
    emxEnsureCapacity((emxArray__common *)V_2, i0, (int32_T)sizeof(creal_T));
    ixstart = PQ_no->size[1];
    i0 = V_2->size[0] * V_2->size[1];
    V_2->size[1] = ixstart;
    emxEnsureCapacity((emxArray__common *)V_2, i0, (int32_T)sizeof(creal_T));
    loop_ub = BUSNUM * PQ_no->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      V_2->data[i0].re = 0.0;
      V_2->data[i0].im = 0.0;
    }

    for (ixstart = 0; ixstart < PQ_no->size[1]; ixstart++) {
      i0 = (int32_T)rt_roundd_snf(PQ_no->data[ixstart]);
      for (i1 = 0; i1 < BUSNUM; i1++) {
        V_2->data[i1 + V_2->size[0] * ixstart] = X_1[i1 + BUSNUM * (i0 - 1)];
      }
    }

    /* % form submatrices of the Jacobian matrix */
    ixstart = PQV_no->size[1];
    i0 = J11->size[0] * J11->size[1];
    J11->size[0] = ixstart;
    emxEnsureCapacity((emxArray__common *)J11, i0, (int32_T)sizeof(real_T));
    ixstart = PQV_no->size[1];
    i0 = J11->size[0] * J11->size[1];
    J11->size[1] = ixstart;
    emxEnsureCapacity((emxArray__common *)J11, i0, (int32_T)sizeof(real_T));
    loop_ub = PQV_no->size[1] * PQV_no->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      J11->data[i0] = 0.0;
    }

    ixstart = PQV_no->size[1];
    i0 = J12->size[0] * J12->size[1];
    J12->size[0] = ixstart;
    emxEnsureCapacity((emxArray__common *)J12, i0, (int32_T)sizeof(real_T));
    ixstart = PQ_no->size[1];
    i0 = J12->size[0] * J12->size[1];
    J12->size[1] = ixstart;
    emxEnsureCapacity((emxArray__common *)J12, i0, (int32_T)sizeof(real_T));
    loop_ub = PQV_no->size[1] * PQ_no->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      J12->data[i0] = 0.0;
    }

    ixstart = PQ_no->size[1];
    i0 = J21->size[0] * J21->size[1];
    J21->size[0] = ixstart;
    emxEnsureCapacity((emxArray__common *)J21, i0, (int32_T)sizeof(real_T));
    ixstart = PQV_no->size[1];
    i0 = J21->size[0] * J21->size[1];
    J21->size[1] = ixstart;
    emxEnsureCapacity((emxArray__common *)J21, i0, (int32_T)sizeof(real_T));
    loop_ub = PQ_no->size[1] * PQV_no->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      J21->data[i0] = 0.0;
    }

    ixstart = PQ_no->size[1];
    i0 = J22->size[0] * J22->size[1];
    J22->size[0] = ixstart;
    emxEnsureCapacity((emxArray__common *)J22, i0, (int32_T)sizeof(real_T));
    ixstart = PQ_no->size[1];
    i0 = J22->size[0] * J22->size[1];
    J22->size[1] = ixstart;
    emxEnsureCapacity((emxArray__common *)J22, i0, (int32_T)sizeof(real_T));
    loop_ub = PQ_no->size[1] * PQ_no->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      J22->data[i0] = 0.0;
    }

    for (ixstart = 0; ixstart < PQV_no->size[1]; ixstart++) {
      loop_ub = XX_1->size[1];
      i0 = (int32_T)rt_roundd_snf(PQV_no->data[ixstart]);
      i1 = a->size[0] * a->size[1];
      a->size[0] = 1;
      a->size[1] = loop_ub;
      emxEnsureCapacity((emxArray__common *)a, i1, (int32_T)sizeof(real_T));
      for (i1 = 0; i1 < loop_ub; i1++) {
        a->data[a->size[0] * i1] = XX_1->data[(i0 + XX_1->size[0] * i1) - 1].re;
      }

      loop_ub = a->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        J11->data[ixstart + J11->size[0] * i0] = a->data[a->size[0] * i0];
      }

      loop_ub = V_2->size[1];
      i0 = (int32_T)rt_roundd_snf(PQV_no->data[ixstart]);
      i1 = a->size[0] * a->size[1];
      a->size[0] = 1;
      a->size[1] = loop_ub;
      emxEnsureCapacity((emxArray__common *)a, i1, (int32_T)sizeof(real_T));
      for (i1 = 0; i1 < loop_ub; i1++) {
        a->data[a->size[0] * i1] = V_2->data[(i0 + V_2->size[0] * i1) - 1].re;
      }

      loop_ub = a->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        J12->data[ixstart + J12->size[0] * i0] = a->data[a->size[0] * i0];
      }
    }

    for (ixstart = 0; ixstart < PQ_no->size[1]; ixstart++) {
      loop_ub = XX_1->size[1];
      i0 = (int32_T)rt_roundd_snf(PQ_no->data[ixstart]);
      i1 = a->size[0] * a->size[1];
      a->size[0] = 1;
      a->size[1] = loop_ub;
      emxEnsureCapacity((emxArray__common *)a, i1, (int32_T)sizeof(real_T));
      for (i1 = 0; i1 < loop_ub; i1++) {
        a->data[a->size[0] * i1] = XX_1->data[(i0 + XX_1->size[0] * i1) - 1].im;
      }

      loop_ub = a->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        J21->data[ixstart + J21->size[0] * i0] = a->data[a->size[0] * i0];
      }

      loop_ub = V_2->size[1];
      i0 = (int32_T)rt_roundd_snf(PQ_no->data[ixstart]);
      i1 = a->size[0] * a->size[1];
      a->size[0] = 1;
      a->size[1] = loop_ub;
      emxEnsureCapacity((emxArray__common *)a, i1, (int32_T)sizeof(real_T));
      for (i1 = 0; i1 < loop_ub; i1++) {
        a->data[a->size[0] * i1] = V_2->data[(i0 + V_2->size[0] * i1) - 1].im;
      }

      loop_ub = a->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        J22->data[ixstart + J22->size[0] * i0] = a->data[a->size[0] * i0];
      }
    }

    /* %%% form_jac.m ends */
    /*  reduced mismatch real and reactive power vectors */
    if (ang_red->size[1] == 1) {
      i0 = bus_int->size[0];
      bus_int->size[0] = ang_red->size[0];
      emxEnsureCapacity((emxArray__common *)bus_int, i0, (int32_T)sizeof(real_T));
      loop_ub = ang_red->size[0];
      for (i0 = 0; i0 < loop_ub; i0++) {
        bus_int->data[i0] = 0.0;
        ixstart = ang_red->size[1];
        for (i1 = 0; i1 < ixstart; i1++) {
          bus_int->data[i0] += ang_red->data[i0 + ang_red->size[0] * i1] *
            delP[i1];
        }
      }
    } else {
      PQVptr = (uint32_T)ang_red->size[0];
      i0 = bus_int->size[0];
      bus_int->size[0] = (int32_T)PQVptr;
      emxEnsureCapacity((emxArray__common *)bus_int, i0, (int32_T)sizeof(real_T));
      loop_ub = (int32_T)PQVptr;
      for (i0 = 0; i0 < loop_ub; i0++) {
        bus_int->data[i0] = 0.0;
      }

      if (ang_red->size[0] == 0) {
      } else {
        loop_ub = 0;
        while ((ang_red->size[0] > 0) && (loop_ub <= 0)) {
          i0 = ang_red->size[0];
          for (ic = 1; ic <= i0; ic++) {
            bus_int->data[ic - 1] = 0.0;
          }

          loop_ub = ang_red->size[0];
        }

        br = 0;
        loop_ub = 0;
        while ((ang_red->size[0] > 0) && (loop_ub <= 0)) {
          ar = 0;
          i0 = br + ang_red->size[1];
          for (ib = br; ib + 1 <= i0; ib++) {
            if (delP[ib] != 0.0) {
              ia = ar;
              i1 = ang_red->size[0];
              for (ic = 0; ic + 1 <= i1; ic++) {
                ia++;
                bus_int->data[ic] += delP[ib] * ang_red->data[ia - 1];
              }
            }

            ar += ang_red->size[0];
          }

          br += ang_red->size[1];
          loop_ub = ang_red->size[0];
        }
      }
    }

    if (volt_red->size[1] == 1) {
      i0 = red_delQ->size[0];
      red_delQ->size[0] = volt_red->size[0];
      emxEnsureCapacity((emxArray__common *)red_delQ, i0, (int32_T)sizeof(real_T));
      loop_ub = volt_red->size[0];
      for (i0 = 0; i0 < loop_ub; i0++) {
        red_delQ->data[i0] = 0.0;
        ixstart = volt_red->size[1];
        for (i1 = 0; i1 < ixstart; i1++) {
          red_delQ->data[i0] += volt_red->data[i0 + volt_red->size[0] * i1] *
            delQ[i1];
        }
      }
    } else {
      PQVptr = (uint32_T)volt_red->size[0];
      i0 = red_delQ->size[0];
      red_delQ->size[0] = (int32_T)PQVptr;
      emxEnsureCapacity((emxArray__common *)red_delQ, i0, (int32_T)sizeof(real_T));
      loop_ub = (int32_T)PQVptr;
      for (i0 = 0; i0 < loop_ub; i0++) {
        red_delQ->data[i0] = 0.0;
      }

      if (volt_red->size[0] == 0) {
      } else {
        loop_ub = 0;
        while ((volt_red->size[0] > 0) && (loop_ub <= 0)) {
          i0 = volt_red->size[0];
          for (ic = 1; ic <= i0; ic++) {
            red_delQ->data[ic - 1] = 0.0;
          }

          loop_ub = volt_red->size[0];
        }

        br = 0;
        loop_ub = 0;
        while ((volt_red->size[0] > 0) && (loop_ub <= 0)) {
          ar = 0;
          i0 = br + volt_red->size[1];
          for (ib = br; ib + 1 <= i0; ib++) {
            if (delQ[ib] != 0.0) {
              ia = ar;
              i1 = volt_red->size[0];
              for (ic = 0; ic + 1 <= i1; ic++) {
                ia++;
                red_delQ->data[ic] += delQ[ib] * volt_red->data[ia - 1];
              }
            }

            ar += volt_red->size[0];
          }

          br += volt_red->size[1];
          loop_ub = volt_red->size[0];
        }
      }
    }

    /* clear delP delQ */
    i0 = b_temp->size[0];
    b_temp->size[0] = bus_int->size[0] + red_delQ->size[0];
    emxEnsureCapacity((emxArray__common *)b_temp, i0, (int32_T)sizeof(real_T));
    loop_ub = bus_int->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_temp->data[i0] = bus_int->data[i0];
    }

    loop_ub = red_delQ->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_temp->data[i0 + bus_int->size[0]] = red_delQ->data[i0];
    }

    i0 = b_J11->size[0] * b_J11->size[1];
    b_J11->size[0] = J11->size[0] + J21->size[0];
    b_J11->size[1] = J11->size[1] + J12->size[1];
    emxEnsureCapacity((emxArray__common *)b_J11, i0, (int32_T)sizeof(real_T));
    loop_ub = J11->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      ixstart = J11->size[0];
      for (i1 = 0; i1 < ixstart; i1++) {
        b_J11->data[i1 + b_J11->size[0] * i0] = J11->data[i1 + J11->size[0] * i0];
      }
    }

    loop_ub = J12->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      ixstart = J12->size[0];
      for (i1 = 0; i1 < ixstart; i1++) {
        b_J11->data[i1 + b_J11->size[0] * (i0 + J11->size[1])] = J12->data[i1 +
          J12->size[0] * i0];
      }
    }

    loop_ub = J21->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      ixstart = J21->size[0];
      for (i1 = 0; i1 < ixstart; i1++) {
        b_J11->data[(i1 + J11->size[0]) + b_J11->size[0] * i0] = J21->data[i1 +
          J21->size[0] * i0];
      }
    }

    loop_ub = J22->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      ixstart = J22->size[0];
      for (i1 = 0; i1 < ixstart; i1++) {
        b_J11->data[(i1 + J11->size[0]) + b_J11->size[0] * (i0 + J21->size[1])] =
          J22->data[i1 + J22->size[0] * i0];
      }
    }

    mldivide(b_J11, b_temp);

    /*  expand solution vectors to all buses */
    if (1 > PQV_no->size[1]) {
      loop_ub = 0;
    } else {
      loop_ub = PQV_no->size[1];
    }

    i0 = a->size[0] * a->size[1];
    a->size[0] = 1;
    a->size[1] = loop_ub;
    emxEnsureCapacity((emxArray__common *)a, i0, (int32_T)sizeof(real_T));
    for (i0 = 0; i0 < loop_ub; i0++) {
      a->data[a->size[0] * i0] = b_temp->data[i0];
    }

    if ((a->size[1] == 1) || (ang_red->size[0] == 1)) {
      i0 = delAng->size[0] * delAng->size[1];
      delAng->size[0] = 1;
      delAng->size[1] = ang_red->size[1];
      emxEnsureCapacity((emxArray__common *)delAng, i0, (int32_T)sizeof(real_T));
      loop_ub = ang_red->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        delAng->data[delAng->size[0] * i0] = 0.0;
        ixstart = a->size[1];
        for (i1 = 0; i1 < ixstart; i1++) {
          delAng->data[delAng->size[0] * i0] += a->data[a->size[0] * i1] *
            ang_red->data[i1 + ang_red->size[0] * i0];
        }
      }
    } else {
      PQVptr = (uint32_T)ang_red->size[1];
      i0 = delAng->size[0] * delAng->size[1];
      delAng->size[0] = 1;
      emxEnsureCapacity((emxArray__common *)delAng, i0, (int32_T)sizeof(real_T));
      i0 = delAng->size[0] * delAng->size[1];
      delAng->size[1] = (int32_T)PQVptr;
      emxEnsureCapacity((emxArray__common *)delAng, i0, (int32_T)sizeof(real_T));
      loop_ub = (int32_T)PQVptr;
      for (i0 = 0; i0 < loop_ub; i0++) {
        delAng->data[i0] = 0.0;
      }

      if (ang_red->size[1] == 0) {
      } else {
        for (loop_ub = 0; loop_ub < ang_red->size[1]; loop_ub++) {
          for (ic = loop_ub; ic + 1 <= loop_ub + 1; ic++) {
            delAng->data[ic] = 0.0;
          }
        }

        br = 0;
        for (loop_ub = 0; loop_ub < ang_red->size[1]; loop_ub++) {
          ar = 0;
          i0 = br + a->size[1];
          for (ib = br; ib + 1 <= i0; ib++) {
            if (ang_red->data[ib] != 0.0) {
              ia = ar;
              for (ic = loop_ub; ic + 1 <= loop_ub + 1; ic++) {
                ia++;
                delAng->data[ic] += ang_red->data[ib] * a->data[ia - 1];
              }
            }

            ar++;
          }

          br += a->size[1];
        }
      }
    }

    ixstart = PQV_no->size[1];
    i0 = PQ_no->size[1];
    i1 = a->size[0] * a->size[1];
    a->size[0] = 1;
    a->size[1] = (int32_T)((real_T)i0 - 1.0) + 1;
    emxEnsureCapacity((emxArray__common *)a, i1, (int32_T)sizeof(real_T));
    loop_ub = (int32_T)((real_T)i0 - 1.0);
    for (i0 = 0; i0 <= loop_ub; i0++) {
      a->data[a->size[0] * i0] = b_temp->data[(int32_T)((real_T)ixstart + (1.0 +
            (real_T)i0)) - 1];
    }

    if ((a->size[1] == 1) || (volt_red->size[0] == 1)) {
      i0 = delV->size[0] * delV->size[1];
      delV->size[0] = 1;
      delV->size[1] = volt_red->size[1];
      emxEnsureCapacity((emxArray__common *)delV, i0, (int32_T)sizeof(real_T));
      loop_ub = volt_red->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        delV->data[delV->size[0] * i0] = 0.0;
        ixstart = a->size[1];
        for (i1 = 0; i1 < ixstart; i1++) {
          delV->data[delV->size[0] * i0] += a->data[a->size[0] * i1] *
            volt_red->data[i1 + volt_red->size[0] * i0];
        }
      }
    } else {
      PQVptr = (uint32_T)volt_red->size[1];
      i0 = delV->size[0] * delV->size[1];
      delV->size[0] = 1;
      emxEnsureCapacity((emxArray__common *)delV, i0, (int32_T)sizeof(real_T));
      i0 = delV->size[0] * delV->size[1];
      delV->size[1] = (int32_T)PQVptr;
      emxEnsureCapacity((emxArray__common *)delV, i0, (int32_T)sizeof(real_T));
      loop_ub = (int32_T)PQVptr;
      for (i0 = 0; i0 < loop_ub; i0++) {
        delV->data[i0] = 0.0;
      }

      if (volt_red->size[1] == 0) {
      } else {
        for (loop_ub = 0; loop_ub < volt_red->size[1]; loop_ub++) {
          for (ic = loop_ub; ic + 1 <= loop_ub + 1; ic++) {
            delV->data[ic] = 0.0;
          }
        }

        br = 0;
        for (loop_ub = 0; loop_ub < volt_red->size[1]; loop_ub++) {
          ar = 0;
          i0 = br + a->size[1];
          for (ib = br; ib + 1 <= i0; ib++) {
            if (volt_red->data[ib] != 0.0) {
              ia = ar;
              for (ic = loop_ub; ic + 1 <= loop_ub + 1; ic++) {
                ia++;
                delV->data[ic] += volt_red->data[ib] * a->data[ia - 1];
              }
            }

            ar++;
          }

          br += a->size[1];
        }
      }
    }

    /*  update voltage magnitude and phase angle */
    if (delV->size[1] == 0) {
      i0 = delV->size[0] * delV->size[1];
      delV->size[0] = 1;
      delV->size[1] = BUSNUM;
      emxEnsureCapacity((emxArray__common *)delV, i0, (int32_T)sizeof(real_T));
      for (i0 = 0; i0 < BUSNUM; i0++) {
        delV->data[i0] = 0.0;
      }
    }

    for (i0 = 0; i0 < BUSNUM; i0++) {
      tap_ratio = V[i0] + delV->data[i0];
      V[i0] = tap_ratio;
    }

    for (i0 = 0; i0 < BUSNUM; i0++) {
      /*  voltage higher than minimum */
      tap_ratio = V[i0];
      if (tap_ratio >= 0.0) {
      } else {
        tap_ratio = 0.0;
      }

      if (tap_ratio <= 2.0) {
      } else {
        tap_ratio = 2.0;
      }

      V[i0] = tap_ratio;
    }

    /*  voltage lower than maximum */
    for (i0 = 0; i0 < BUSNUM; i0++) {
      tap_ratio = ang[i0] + delAng->data[i0];
      ang[i0] = tap_ratio;
    }

    /*  calculate the power mismatch and check convergence */
    /* %%%%%%%% calc.m */
    /*  voltage in rectangular coordinate */
    memcpy(&P[0], &ang[0], BUSNUM * sizeof(real_T));
    b_cos(P, BUSNUM);
    memcpy(&Q[0], &ang[0], BUSNUM * sizeof(real_T));
    b_sin(Q, BUSNUM);
    for (i0 = 0; i0 < BUSNUM; i0++) {
      V_rect[i0].re = V[i0] * (P[i0] + 0.0 * Q[i0]);
      V_rect[i0].im = V[i0] * Q[i0];
    }

    /*  bus current injection */
    if (Y->size[1] == 1) {
      i0 = y->size[0];
      y->size[0] = Y->size[0];
      emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(creal_T));
      loop_ub = Y->size[0];
      for (i0 = 0; i0 < loop_ub; i0++) {
        y->data[i0].re = 0.0;
        y->data[i0].im = 0.0;
        ixstart = Y->size[1];
        for (i1 = 0; i1 < ixstart; i1++) {
          tap_ratio = Y->data[i0 + Y->size[0] * i1].re * V_rect[i1].re - Y->
            data[i0 + Y->size[0] * i1].im * V_rect[i1].im;
          V_rect_im = Y->data[i0 + Y->size[0] * i1].re * V_rect[i1].im + Y->
            data[i0 + Y->size[0] * i1].im * V_rect[i1].re;
          y->data[i0].re += tap_ratio;
          y->data[i0].im += V_rect_im;
        }
      }
    } else {
      PQVptr = (uint32_T)Y->size[0];
      i0 = y->size[0];
      y->size[0] = (int32_T)PQVptr;
      emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(creal_T));
      loop_ub = (int32_T)PQVptr;
      for (i0 = 0; i0 < loop_ub; i0++) {
        y->data[i0].re = 0.0;
        y->data[i0].im = 0.0;
      }

      if (Y->size[0] == 0) {
      } else {
        loop_ub = 0;
        while ((Y->size[0] > 0) && (loop_ub <= 0)) {
          i0 = Y->size[0];
          for (ic = 1; ic <= i0; ic++) {
            y->data[ic - 1].re = 0.0;
            y->data[ic - 1].im = 0.0;
          }

          loop_ub = Y->size[0];
        }

        br = 0;
        loop_ub = 0;
        while ((Y->size[0] > 0) && (loop_ub <= 0)) {
          ar = 0;
          i0 = br + Y->size[1];
          for (ib = br; ib + 1 <= i0; ib++) {
            if ((V_rect[ib].re != 0.0) || (V_rect[ib].im != 0.0)) {
              temp.re = V_rect[ib].re - 0.0 * V_rect[ib].im;
              temp.im = V_rect[ib].im + 0.0 * V_rect[ib].re;
              ia = ar;
              i1 = Y->size[0];
              for (ic = 0; ic + 1 <= i1; ic++) {
                ia++;
                tap_ratio = temp.re * Y->data[ia - 1].re - temp.im * Y->data[ia
                  - 1].im;
                temp_im = temp.re * Y->data[ia - 1].im + temp.im * Y->data[ia -
                  1].re;
                y->data[ic].re += tap_ratio;
                y->data[ic].im += temp_im;
              }
            }

            ar += Y->size[0];
          }

          br += Y->size[1];
          loop_ub = Y->size[0];
        }
      }
    }

    /*  power output  */
    for (i0 = 0; i0 < BUSNUM; i0++) {
      y_re = y->data[i0].re;
      temp_im = -y->data[i0].im;
      tap_ratio = V_rect[i0].re;
      V_rect_im = V_rect[i0].im;
      V_rect[i0].re = V_rect[i0].re * y_re - V_rect[i0].im * temp_im;
      V_rect[i0].im = tap_ratio * temp_im + V_rect_im * y_re;
    }

    for (ixstart = 0; ixstart < BUSNUM; ixstart++) {
      delP[ixstart] = (Pg[ixstart] - Pl[ixstart]) - V_rect[ixstart].re;
      delQ[ixstart] = (Qg[ixstart] - Ql[ixstart]) - V_rect[ixstart].im;
      P[ixstart] = V_rect[ixstart].re;
      Q[ixstart] = V_rect[ixstart].im;
    }

    /*  zero out mismatches on swing bus and generation bus */
    for (ixstart = 0; ixstart < (int32_T)BUSNUM; ixstart++) {
      if (bus_type[ixstart] == 1.0) {
        delP[(int32_T)(1.0 + (real_T)ixstart) - 1] = 0.0;
        delQ[(int32_T)(1.0 + (real_T)ixstart) - 1] = 0.0;
      } else {
        if (bus_type[(int32_T)(1.0 + (real_T)ixstart) - 1] == 2.0) {
          delQ[(int32_T)(1.0 + (real_T)ixstart) - 1] = 0.0;
        }
      }
    }

    /*   total mismatch */
    if (norm(delQ, BUSNUM) + norm(delP, BUSNUM) > 1.0E-11) {
    } else {
      conv_flag = 0;
    }

    /* %%%%%%%%calc.m ends */
  }

  emxFree_real_T(&b_J11);
  emxFree_real_T(&a);
  emxFree_creal_T(&b);
  emxFree_real_T(&delV);
  emxFree_real_T(&delAng);
  emxFree_real_T(&b_temp);
  emxFree_real_T(&red_delQ);
  emxFree_real_T(&J22);
  emxFree_real_T(&J21);
  emxFree_real_T(&J12);
  emxFree_real_T(&J11);
  emxFree_creal_T(&XX_1);
  emxFree_creal_T(&V_2);
  emxFree_real_T(&volt_red);
  emxFree_real_T(&ang_red);
  emxFree_real_T(&PQ_no);
  emxFree_real_T(&PQV_no);
  emxFree_creal_T(&y);
  emxFree_real_T(&bus_int);
  emxFree_creal_T(&Y);

  /* % */
  for (ixstart = 0; ixstart < (int32_T)BUSNUM; ixstart++) {
    if (bus_type[ixstart] == 2.0) {
      Pg[(int32_T)(1.0 + (real_T)ixstart) - 1] = P[(int32_T)(1.0 + (real_T)
          ixstart) - 1] + Pl[(int32_T)(1.0 + (real_T)ixstart) - 1];
      Qg[(int32_T)(1.0 + (real_T)ixstart) - 1] = Q[(int32_T)(1.0 + (real_T)
          ixstart) - 1] + Ql[(int32_T)(1.0 + (real_T)ixstart) - 1];
    } else {
      if (bus_type[(int32_T)(1.0 + (real_T)ixstart) - 1] == 3.0) {
        Pl[(int32_T)(1.0 + (real_T)ixstart) - 1] = Pg[(int32_T)(1.0 + (real_T)
            ixstart) - 1] - P[(int32_T)(1.0 + (real_T)ixstart) - 1];
        Ql[(int32_T)(1.0 + (real_T)ixstart) - 1] = Qg[(int32_T)(1.0 + (real_T)
            ixstart) - 1] - Q[(int32_T)(1.0 + (real_T)ixstart) - 1];
      }
    }
  }

  Pg[(int32_T)SB - 1] = P[(int32_T)SB - 1] + Pl[(int32_T)SB - 1];
  Qg[(int32_T)SB - 1] = Q[(int32_T)SB - 1] + Ql[(int32_T)SB - 1];


  for (i0 = 0; i0 < BUSNUM; i0++) {
    bus_sol[i0] = bus[i0];
    bus_sol[BUSNUM + i0] = V[i0];
    bus_sol[BUSNUM*2 + i0] = ang[i0] * 180.0 / 3.1415926535897931;
    bus_sol[BUSNUM*3 + i0] = Pg[i0];
    bus_sol[BUSNUM*4 + i0] = Qg[i0];
    bus_sol[BUSNUM*5 + i0] = Pl[i0];
    bus_sol[BUSNUM*6 + i0] = Ql[i0];
    bus_sol[BUSNUM*7 + i0] = bus[BUSNUM*7 + i0];
    bus_sol[BUSNUM*8 + i0] = bus[BUSNUM*8 + i0];
    bus_sol[BUSNUM*9 + i0] = bus_type[i0];
  }
  return iter;

}

/* End of code generation (loadflow.c) */
