#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP windcross(SEXP inSeries1, SEXP inSeries2, SEXP wMax, SEXP tMax, SEXP wInc, SEXP tInc);

static const R_CallMethodDef CallEntries[] = {
    {"windcross", (DL_FUNC) &windcross, 6},
    {NULL, NULL, 0}
};

void R_init_wcc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
