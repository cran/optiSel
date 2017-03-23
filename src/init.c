#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#define ARMA_64BIT_WORD 1

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP optiSel_rcpp_completeness(SEXP, SEXP, SEXP, SEXP);
extern SEXP optiSel_rcpp_haplofreq(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP optiSel_rcpp_nativecont(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP optiSel_rcpp_segBreedComp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP optiSel_rcpp_segIBD(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP optiSel_rcpp_segIBDandN(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP optiSel_rcpp_segIBDandNVersion2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP optiSel_rcpp_segInbreeding(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP optiSel_rcpp_segN(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"optiSel_rcpp_completeness",       (DL_FUNC) &optiSel_rcpp_completeness,        4},
    {"optiSel_rcpp_haplofreq",          (DL_FUNC) &optiSel_rcpp_haplofreq,          22},
    {"optiSel_rcpp_nativecont",         (DL_FUNC) &optiSel_rcpp_nativecont,          6},
    {"optiSel_rcpp_segBreedComp",       (DL_FUNC) &optiSel_rcpp_segBreedComp,        6},
    {"optiSel_rcpp_segIBD",             (DL_FUNC) &optiSel_rcpp_segIBD,             17},
    {"optiSel_rcpp_segIBDandN",         (DL_FUNC) &optiSel_rcpp_segIBDandN,         16},
    {"optiSel_rcpp_segIBDandNVersion2", (DL_FUNC) &optiSel_rcpp_segIBDandNVersion2, 14},
    {"optiSel_rcpp_segInbreeding",      (DL_FUNC) &optiSel_rcpp_segInbreeding,      17},
    {"optiSel_rcpp_segN",               (DL_FUNC) &optiSel_rcpp_segN,                6},
    {NULL, NULL, 0}
};

void R_init_optiSel(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}
