// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// MtM
inline MatrixXd MtM(const MapMatd& M);
RcppExport SEXP _gprcpp_MtM(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MapMatd& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(MtM(M));
    return rcpp_result_gen;
END_RCPP
}
// AtB
inline MatrixXd AtB(const MapMatd& A, const MapMatd& B);
RcppExport SEXP _gprcpp_AtB(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MapMatd& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const MapMatd& >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(AtB(A, B));
    return rcpp_result_gen;
END_RCPP
}
// A_solve
inline MatrixXd A_solve(const MapMatd& M);
RcppExport SEXP _gprcpp_A_solve(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MapMatd& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(A_solve(M));
    return rcpp_result_gen;
END_RCPP
}
// A_solve_B
inline MatrixXd A_solve_B(const MapMatd& A, const MapMatd& B);
RcppExport SEXP _gprcpp_A_solve_B(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MapMatd& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const MapMatd& >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(A_solve_B(A, B));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gprcpp_MtM", (DL_FUNC) &_gprcpp_MtM, 1},
    {"_gprcpp_AtB", (DL_FUNC) &_gprcpp_AtB, 2},
    {"_gprcpp_A_solve", (DL_FUNC) &_gprcpp_A_solve, 1},
    {"_gprcpp_A_solve_B", (DL_FUNC) &_gprcpp_A_solve_B, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_gprcpp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}