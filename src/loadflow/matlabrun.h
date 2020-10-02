/*
 * matlabrun.h
 *
 * Code generation for function 'matlabrun'
 *
 * C source code generated on: Thu Dec 11 00:48:24 2014
 *
 */

#ifndef __MATLABRUN_H__
#define __MATLABRUN_H__
/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"

#include "rtwtypes.h"
#include "matlabrun_types.h"


extern uint_T matlabrun(real_T *bus, const real_T *line, const uint_T noofgens,
              const uint_T SB, const uint_T noofslack,
              real_T *bus_sol);

#endif
/* End of code generation (matlabrun.h) */
