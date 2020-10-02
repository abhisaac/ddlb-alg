/*
 * loadflow.h
 *
 * Code generation for function 'loadflow'
 *
 * C source code generated on: Fri Dec 12 22:50:42 2014
 *
 */

#ifndef __LOADFLOW_H__
#define __LOADFLOW_H__
/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"

#include "rtwtypes.h"
#include "matlabrun_types.h"

/* Function Declarations */
extern uint_T loadflow(const real_T *bus, const real_T *line, const uint_T noofgens,
              const uint_T SB, const uint_T noofslack,
              real_T *bus_sol);
#endif
/* End of code generation (loadflow.h) */
