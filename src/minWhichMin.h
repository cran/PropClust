#ifndef __minWhichMin_h__

#define __minWhichMin_h__

#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/libextern.h>

void minWhichMin(double * matrix, int * nRows, int * nColumns, double * min, double * whichMin);

#endif
