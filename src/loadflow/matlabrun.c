/*
 * matlabrun.c
 *
 * Code generation for function 'matlabrun'
 *
 * C source code generated on: Thu Dec 11 00:48:24 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "matlabrun.h"
#include "loadflow.h"

#include <stdio.h>

/* Function Definitions */

uint_T matlabrun(real_T *bus, const real_T *line, const uint_T noofgens,
    const uint_T SB, const uint_T noofslack,
    real_T *bus_sol)
{
  return loadflow(bus,line,noofgens, SB, noofslack, bus_sol);
}

/* End of code generation (matlabrun.c) */
