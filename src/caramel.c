#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h> // for NULL
#include <Rmath.h>

void F77_NAME(pareto)(double *X, int nind, int nobj, int *Ft);
extern SEXP c_pareto(SEXP X){
  SEXP dim = PROTECT(getAttrib(X, R_DimSymbol));
  const int nind = INTEGER(dim)[0];
  const int nobj = INTEGER(dim)[1];
  SEXP Ft;
  PROTECT(Ft = allocVector(INTSXP, nind));
  F77_CALL(pareto)(REAL(X), nind, nobj, INTEGER(Ft));
  UNPROTECT(2);
  return(Ft);
}

void F77_NAME(dominate)(double *matobj, int nind, int nobj, int *f);
extern SEXP c_dominate(SEXP matobj){
  SEXP dim = PROTECT(getAttrib(matobj, R_DimSymbol));
  const int nind = INTEGER(dim)[0];
  const int nobj = INTEGER(dim)[1];
  SEXP f;
  PROTECT(f = allocVector(INTSXP, nind));
  F77_CALL(dominate)(REAL(matobj), nind, nobj, INTEGER(f));
  UNPROTECT(2);
  return(f);
}

void F77_NAME(dominated)(double *Xi, double *X, int nind, int nobj, int *is_dominated);
extern SEXP c_dominated(SEXP Xi, SEXP X){
  SEXP dim = PROTECT(getAttrib(X, R_DimSymbol));
  const int nind = INTEGER(dim)[0];
  const int nobj = INTEGER(dim)[1];
  SEXP is_dominated;
  PROTECT(is_dominated = allocVector(INTSXP, nind));
  F77_CALL(dominated)(REAL(Xi), REAL(X), nind, nobj, INTEGER(is_dominated));
  UNPROTECT(2);
  return(is_dominated);
}

static const R_CallMethodDef CallEntries[] = {
    {"c_pareto", (DL_FUNC) &c_pareto, 1},
    {"c_dominate", (DL_FUNC) &c_dominate, 1},
    {"c_dominated", (DL_FUNC) &c_dominated, 2},
    {NULL, NULL, 0}
};

void R_init_caRamel (DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
