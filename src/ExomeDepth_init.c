#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP C_hmm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP get_loglike_matrix(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"C_hmm",              (DL_FUNC) &C_hmm,              6},
    {"get_loglike_matrix", (DL_FUNC) &get_loglike_matrix, 5},
    {NULL, NULL, 0}
};

void R_init_ExomeDepth(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

