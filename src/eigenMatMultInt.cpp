// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>

//' Fast matrix multiplication for double x int
//' 
//' @param A a double matrix
//' @param B an integer matrix
//' @return An integer matrix representing the product of A and B
// [[Rcpp::export]]
SEXP eigenMatMultInt(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map< Eigen::MatrixXi> B){
  Eigen::MatrixXd C = A.transpose() * B.cast<double>();
  return Rcpp::wrap(C);
}
