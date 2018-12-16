// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// cG_calcGibbsProbY_Simple
NumericVector cG_calcGibbsProbY_Simple(const IntegerMatrix counts, IntegerVector nGbyTS, IntegerMatrix nTSbyC, IntegerVector nbyTS, IntegerVector nbyG, const IntegerVector y, const int L, const int index, const double gamma, const double beta, const double delta);
RcppExport SEXP _celda_cG_calcGibbsProbY_Simple(SEXP countsSEXP, SEXP nGbyTSSEXP, SEXP nTSbyCSEXP, SEXP nbyTSSEXP, SEXP nbyGSEXP, SEXP ySEXP, SEXP LSEXP, SEXP indexSEXP, SEXP gammaSEXP, SEXP betaSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type counts(countsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nGbyTS(nGbyTSSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type nTSbyC(nTSbyCSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nbyTS(nbyTSSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nbyG(nbyGSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< const int >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int >::type index(indexSEXP);
    Rcpp::traits::input_parameter< const double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(cG_calcGibbsProbY_Simple(counts, nGbyTS, nTSbyC, nbyTS, nbyG, y, L, index, gamma, beta, delta));
    return rcpp_result_gen;
END_RCPP
}
// cG_CalcGibbsProbY
NumericVector cG_CalcGibbsProbY(const int index, const IntegerMatrix& counts, const IntegerMatrix& nTSbyC, const IntegerVector& nbyTS, const IntegerVector& nGbyTS, const IntegerVector& nbyG, const IntegerVector& y, const int L, const int nG, const NumericVector& lg_beta, const NumericVector& lg_gamma, const NumericVector& lg_delta, const double delta);
RcppExport SEXP _celda_cG_CalcGibbsProbY(SEXP indexSEXP, SEXP countsSEXP, SEXP nTSbyCSEXP, SEXP nbyTSSEXP, SEXP nGbyTSSEXP, SEXP nbyGSEXP, SEXP ySEXP, SEXP LSEXP, SEXP nGSEXP, SEXP lg_betaSEXP, SEXP lg_gammaSEXP, SEXP lg_deltaSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type index(indexSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type counts(countsSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type nTSbyC(nTSbyCSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type nbyTS(nbyTSSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type nGbyTS(nGbyTSSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type nbyG(nbyGSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const int >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int >::type nG(nGSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lg_beta(lg_betaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lg_gamma(lg_gammaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lg_delta(lg_deltaSEXP);
    Rcpp::traits::input_parameter< const double >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(cG_CalcGibbsProbY(index, counts, nTSbyC, nbyTS, nGbyTS, nbyG, y, L, nG, lg_beta, lg_gamma, lg_delta, delta));
    return rcpp_result_gen;
END_RCPP
}
// eigenMatMultInt
SEXP eigenMatMultInt(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map< Eigen::MatrixXi> B);
RcppExport SEXP _celda_eigenMatMultInt(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< const Eigen::Map< Eigen::MatrixXi> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(eigenMatMultInt(A, B));
    return rcpp_result_gen;
END_RCPP
}
// fastNormProp
SEXP fastNormProp(NumericMatrix R_counts, double R_alpha);
RcppExport SEXP _celda_fastNormProp(SEXP R_countsSEXP, SEXP R_alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type R_counts(R_countsSEXP);
    Rcpp::traits::input_parameter< double >::type R_alpha(R_alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(fastNormProp(R_counts, R_alpha));
    return rcpp_result_gen;
END_RCPP
}
// fastNormPropLog
SEXP fastNormPropLog(NumericMatrix R_counts, double R_alpha);
RcppExport SEXP _celda_fastNormPropLog(SEXP R_countsSEXP, SEXP R_alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type R_counts(R_countsSEXP);
    Rcpp::traits::input_parameter< double >::type R_alpha(R_alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(fastNormPropLog(R_counts, R_alpha));
    return rcpp_result_gen;
END_RCPP
}
// fastNormPropSqrt
SEXP fastNormPropSqrt(NumericMatrix R_counts, double R_alpha);
RcppExport SEXP _celda_fastNormPropSqrt(SEXP R_countsSEXP, SEXP R_alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type R_counts(R_countsSEXP);
    Rcpp::traits::input_parameter< double >::type R_alpha(R_alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(fastNormPropSqrt(R_counts, R_alpha));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _colSumByGroup(SEXP, SEXP);
RcppExport SEXP _colSumByGroup_numeric(SEXP, SEXP);
RcppExport SEXP _colSumByGroupChange(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP _rowSumByGroup(SEXP, SEXP);
RcppExport SEXP _rowSumByGroup_numeric(SEXP, SEXP);
RcppExport SEXP _rowSumByGroupChange(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_celda_cG_calcGibbsProbY_Simple", (DL_FUNC) &_celda_cG_calcGibbsProbY_Simple, 11},
    {"_celda_cG_CalcGibbsProbY", (DL_FUNC) &_celda_cG_CalcGibbsProbY, 13},
    {"_celda_eigenMatMultInt", (DL_FUNC) &_celda_eigenMatMultInt, 2},
    {"_celda_fastNormProp", (DL_FUNC) &_celda_fastNormProp, 2},
    {"_celda_fastNormPropLog", (DL_FUNC) &_celda_fastNormPropLog, 2},
    {"_celda_fastNormPropSqrt", (DL_FUNC) &_celda_fastNormPropSqrt, 2},
<<<<<<< HEAD
    {"_colSumByGroup",                  (DL_FUNC) &_colSumByGroup,                   2},
    {"_colSumByGroup_numeric",          (DL_FUNC) &_colSumByGroup_numeric,           2},
    {"_colSumByGroupChange",            (DL_FUNC) &_colSumByGroupChange,             4},
    {"_rowSumByGroup",                  (DL_FUNC) &_rowSumByGroup,                   2},
    {"_rowSumByGroup_numeric",          (DL_FUNC) &_rowSumByGroup_numeric,           2},
    {"_rowSumByGroupChange",            (DL_FUNC) &_rowSumByGroupChange,             4},
    {NULL, NULL, 0}
};

RcppExport void R_init_celda(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
