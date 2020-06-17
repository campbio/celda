// [[Rcpp::depends(Rcpp)]]
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
    if (cs[i] + alpha_tot == 0) {
      stop("Division by 0. Make sure colSums of counts does not contain 0 after rounding counts to integers.");
    }
    res(_,i) = (R_counts(_,i) + R_alpha) / (cs[i] + alpha_tot);
  }
  
  return Rcpp::wrap(res);
}


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
    if (cs[i] + alpha_tot == 0) {
      stop("Division by 0. Make sure colSums of counts does not contain 0 after rounding counts to integers.");
    }
    res(_,i) = log((R_counts(_,i) + R_alpha) / (cs[i] + alpha_tot));
  }
  
  return Rcpp::wrap(res);
}


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
    if (cs[i] + alpha_tot == 0) {
      stop("Division by 0. Make sure colSums of counts does not contain 0 after rounding counts to integers.");
    }
    res(_,i) = sqrt((R_counts(_,i) + R_alpha) / (cs[i] + alpha_tot));
  }
  
  return Rcpp::wrap(res);
}



//' get row and column indices of none zero elements in the matrix
//' 
//' @param R_counts A matrix
//' @return An integer matrix where each row is a row, column indices pair 
// [[Rcpp::export]]
SEXP nonzero(NumericMatrix R_counts) {
  
  IntegerVector row(1);
  IntegerVector col(1);
  NumericVector val(1);
  
  int nR = R_counts.nrow();
  int nC = R_counts.ncol();
  double x;
  
  for (int c = 0; c < nC; c++) {
    for (int r = 0; r < nR; r++) {
      x = R_counts[c * nR + r];
      if (x != 0) {
        row.push_back(r + 1);
        col.push_back(c + 1);
        val.push_back(x);
      }
    }
  }
  
  row.erase(0);
  col.erase(0);
  val.erase(0);
  
  List res;
  res["row"] = row;
  res["col"] = col;
  res["val"] = val;
  
  return(res);
}

