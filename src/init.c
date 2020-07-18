#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void mcmcloopM2a1(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void mcmcloopM2b1(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *);

extern void mcmcloopM2a2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void mcmcloopM2b2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *);

static const R_CMethodDef CEntries[] = {
    {"mcmcloopM2a1", (DL_FUNC) &mcmcloopM2a1, 30},
    {"mcmcloopM2b1", (DL_FUNC) &mcmcloopM2b1, 32},
    {"mcmcloopM2a2", (DL_FUNC) &mcmcloopM2a2, 30},
    {"mcmcloopM2b2", (DL_FUNC) &mcmcloopM2b2, 32},
    {NULL, NULL, 0}
};


void R_init_eCAR(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
