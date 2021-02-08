#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP agglom(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gapfill(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP getEIC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP indexed(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP partID(SEXP, SEXP);
extern SEXP pickpeak(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"agglom",   (DL_FUNC) &agglom,    5},
    {"gapfill",  (DL_FUNC) &gapfill,   7},
    {"getEIC",   (DL_FUNC) &getEIC,    9},
    {"indexed",  (DL_FUNC) &indexed,   5},
    {"partID",   (DL_FUNC) &partID,    2},
    {"pickpeak", (DL_FUNC) &pickpeak, 12},
    {NULL, NULL, 0}
};

void R_init_LipidMS(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}