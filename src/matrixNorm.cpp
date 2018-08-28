// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
using namespace Rcpp ;

//' Fast normalization for numeric matrix
//' 
//' @param R_counts An integer matrix
//' @param R_alpha A double value to be added to the matrix as a pseudocount
//' @return A numeric matrix where the columns have been normalized to proportions
// [[Rcpp::export]]
SEXP fastNormProp(NumericMatrix R_counts, double R_alpha) {

  // Get colSums and instantiate new matrix
  NumericVector cs = colSums(R_counts);
  NumericMatrix res = NumericMatrix(R_counts.nrow(), R_counts.ncol());
  
  // Normalize cell counts to proportions after adding pseudocount  
  double alpha_tot = R_counts.nrow() * R_alpha;
  for (int i = 0; i < R_counts.ncol(); ++i) {
    res(_,i) = (R_counts(_,i) + R_alpha) / (cs[i] + alpha_tot);
  }
  
  return Rcpp::wrap(res);
}


// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
using namespace Rcpp ;

//' Fast normalization for numeric matrix
//' 
//' @param R_counts An integer matrix
//' @param R_alpha A double value to be added to the matrix as a pseudocount
//' @return A numeric matrix where the columns have been normalized to proportions
// [[Rcpp::export]]
SEXP fastNormPropLog(NumericMatrix R_counts, double R_alpha) {

  // Get colSums and instantiate new matrix
  NumericVector cs = colSums(R_counts);
  NumericMatrix res = NumericMatrix(R_counts.nrow(), R_counts.ncol());
  
  // Normalize cell counts to proportions after adding pseudocount  
  double alpha_tot = R_counts.nrow() * R_alpha;
  for (int i = 0; i < R_counts.ncol(); ++i) {
    res(_,i) = log((R_counts(_,i) + R_alpha) / (cs[i] + alpha_tot));
  }
  
  return Rcpp::wrap(res);
}


// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
using namespace Rcpp ;

//' Fast normalization for numeric matrix
//' 
//' @param R_counts An integer matrix
//' @param R_alpha A double value to be added to the matrix as a pseudocount
//' @return A numeric matrix where the columns have been normalized to proportions
// [[Rcpp::export]]
SEXP fastNormPropSqrt(NumericMatrix R_counts, double R_alpha) {

  // Get colSums and instantiate new matrix
  NumericVector cs = colSums(R_counts);
  NumericMatrix res = NumericMatrix(R_counts.nrow(), R_counts.ncol());
  
  // Normalize cell counts to proportions after adding pseudocount  
  double alpha_tot = R_counts.nrow() * R_alpha;
  for (int i = 0; i < R_counts.ncol(); ++i) {
    res(_,i) = sqrt((R_counts(_,i) + R_alpha) / (cs[i] + alpha_tot));
  }
  
  return Rcpp::wrap(res);
}

