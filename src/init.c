#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> //
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void mcmcloop_qva(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *, void *, void *);

extern void mcmcloopM2a1(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *);


extern void mcmcloopM2a2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *);

extern void mcmcloopM2b1(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *, void *);

extern void mcmcloopM2b2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                         void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"mcmcloop_qva", (DL_FUNC) &mcmcloop_qva, 24},
    {"mcmcloopM2a1", (DL_FUNC) &mcmcloopM2a1, 31},
    {"mcmcloopM2a2", (DL_FUNC) &mcmcloopM2a2, 31},
    {"mcmcloopM2b1", (DL_FUNC) &mcmcloopM2b1, 33},
    {"mcmcloopM2b2", (DL_FUNC) &mcmcloopM2b2, 33},
    {NULL, NULL, 0}
};

void R_init_modernVA(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
