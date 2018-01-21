#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

#include "minWhichMin.h"
#include "PropClustParallelTrials.h"


#define CDEF(name, n, args)  {#name, (DL_FUNC) &name, n, args}

static R_NativePrimitiveArgType 
  // minWhichMin
  minWhich_t[] = { REALSXP, INTSXP, INTSXP, REALSXP, REALSXP },
  // propclusttrial, propclustaccel, propensityclustering
  propclusttrial_t[] = {SINGLESXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP},
  // propdecompaccel, propensitydecomposition
  propdecompaccel_t[] = {SINGLESXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP},
  // singleclusterupdate
  singleclusterupdate_t[] = {SINGLESXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP},
  // spped1 and speed2
  speed_t[] = {REALSXP, INTSXP};
  
  
  

static const R_CMethodDef R_CMethods[] = {
   CDEF(minWhichMin, 5, minWhich_t),
   {NULL, NULL, 0, NULL} };

static const R_FortranMethodDef R_FortranMethods[] = {
   {"propclusttrial", (DL_FUNC) &F77_NAME(propclusttrial), 10, propclusttrial_t},
   {"propclustaccel", (DL_FUNC) &F77_NAME(propclusttrial), 10, propclusttrial_t},
   {"propensityclustering", (DL_FUNC) &F77_NAME(propclusttrial), 10, propclusttrial_t},
   {"propensitydecomposition", (DL_FUNC) &F77_NAME(propensitydecomposition), 9, propdecompaccel_t},
   {"propdecompaccel", (DL_FUNC) &F77_NAME(propdecompaccel), 9, propdecompaccel_t},
   {"singleclusterupdate", (DL_FUNC) & F77_NAME(singleclusterupdate), 6, singleclusterupdate_t},
   {"speedtest1", (DL_FUNC) & F77_NAME(speedtest1), 2, speed_t},
   {"speedtest2", (DL_FUNC) & F77_NAME(speedtest2), 2, speed_t},
   {NULL, NULL, 0, NULL}
};

void R_init_PropClust(DllInfo *dll)
{
    R_registerRoutines(dll, R_CMethods, NULL, R_FortranMethods, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

