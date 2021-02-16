#include <RcppEigen.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp; 

  
// [[Rcpp::export]]
Rcpp::List decontXEM(const Eigen::MappedSparseMatrix<double>& counts,
                     const NumericVector& counts_colsums,
                     const NumericVector& theta,
                     const NumericMatrix& eta,
                     const NumericMatrix& phi,
                     const IntegerVector& z,
                     const bool& estimate_delta,
                     const NumericVector& delta,
                     const double& pseudocount) {

  // Perform error checking
  if(counts.cols() != theta.size()) {
    stop("Length of 'theta' must be equal to the number of columns in 'counts'.");
  }
  if(counts.cols() != z.size()) {
    stop("Length of 'z' must be equal to the number of columns in 'counts'.");
  }  
  if(counts.cols() != counts_colsums.size()) {
    stop("Length of 'counts_colsums' must be equal to the number of columns in 'counts'.");
  }    
  if(counts.rows() != phi.nrow()) {
    stop("The number of rows in 'phi' must be equal to the number of rows in 'counts'.");
  }
  if(counts.rows() != eta.nrow()) {
    stop("The number of rows in 'eta' must be equal to the number of rows in 'counts'.");
  }
  if(phi.ncol() != eta.ncol()) {
    stop("The number of columns in 'eta' must be equal to the number of columns in 'phi'.");
  }
  if(min(z) < 1 || max(z) > eta.ncol()) {
    stop("The entries in 'z' need to be between 1 and the number of columns in eta and phi.");
  }
  if(delta.size() != 2 || sum(delta < 0) > 0) {
    stop("'delta' must be a numeric vector of length 2 with positive integers.");
  }
  
  // Declare variables and functions
  NumericVector new_theta(theta.size());
  NumericVector native_total(theta.size());
  NumericMatrix new_phi(phi.nrow(), phi.ncol());
  NumericMatrix new_eta(eta.nrow(), eta.ncol());
    
  // Obtaining 'fit_dirichlet' function from MCMCprecision package
  Environment pkg = Environment::namespace_env("MCMCprecision");
  Function f = pkg["fit_dirichlet"];

  int i;
  int j;
  int k;
  int nr = phi.nrow();
  double x; 
  double pcontamin;
  double pnative;
  double normp;
  double px;
  for(int j = 0; j < counts.cols(); ++j) {
    for (Eigen::MappedSparseMatrix<double>::InnerIterator i_(counts, j); i_; ++i_) {
      i = i_.index();
      x = i_.value();
      k = z[j] - 1;
      
      // Calculate variational probabilities
      // Removing the log/exp speeds it up and produces the same result since
      // there are only 2 probabilities being multiplied
      
      //pnative = log(phi(i,k) + pseudocount) + log(theta(j) + pseudocount);
      //pcontamin = log(eta(i,k) + pseudocount) + log(1 - theta(j) + pseudocount);
      pnative = (phi[nr * k + i] + pseudocount) * (theta[j] + pseudocount);
      pcontamin = (eta[nr * k + i] + pseudocount) * (1 - theta[j] + pseudocount);
      
      // Normalize probabilities and add to proper components
      //normp = exp(pnative) / (exp(pcontamin) + exp(pnative));
      normp = pnative / (pcontamin + pnative);
      px = normp * x;
      new_phi(i,k) += px;
      native_total(j) += px;
    }  
  }
  
  // Calculate Eta using Weights from Phi
  NumericVector phi_rowsum = rowSums(new_phi);
  for(i = 0; i < new_eta.ncol(); i++) {
    for(j = 0; j < new_eta.nrow(); j++) {
      new_eta(j,i) = phi_rowsum[j] - new_phi(j,i);
    }
  }
  
  // Normalize Phi and Eta
  NumericVector phi_colsum = colSums(new_phi);
  NumericVector eta_colsum = colSums(new_eta);  
  for(i = 0; i < new_phi.ncol(); i++) {
    new_phi(_,i) = new_phi(_,i) / phi_colsum[i];
    new_eta(_,i) = new_eta(_,i) / eta_colsum[i];
  }
  
  // Update Theta
  NumericVector contamination_prop = (counts_colsums - native_total) / counts_colsums;
  NumericVector native_prop = 1 - contamination_prop;
  NumericMatrix theta_raw = cbind(native_prop, contamination_prop);
  
  NumericVector new_delta = delta;
  if(estimate_delta == TRUE) {
    Rcpp::List result = f(Named("x", theta_raw)); 
    new_delta = result["alpha"];
  }
  
  // Estimate new theta
  new_theta = (native_total + new_delta[0]) / (counts_colsums + sum(new_delta));
  
  return Rcpp::List::create(Rcpp::Named("phi") = new_phi,
                            Rcpp::Named("eta") = new_eta,
                            Rcpp::Named("theta") = new_theta,
                            Rcpp::Named("delta") = new_delta,
                            Rcpp::Named("contamination") = contamination_prop);
}





// [[Rcpp::export]]
double decontXLogLik(const Eigen::MappedSparseMatrix<double>& counts,
                     const NumericVector& theta,
                     const NumericMatrix& eta,
                     const NumericMatrix& phi,
                     const IntegerVector& z,
                     const double& pseudocount) {
  
  // Perform error checking
  if(counts.cols() != theta.size()) {
    stop("Length of 'theta' must be equal to the number of columns in 'counts'.");
  }
  if(counts.cols() != z.size()) {
    stop("Length of 'z' must be equal to the number of columns in 'counts'.");
  }  
  if(counts.rows() != phi.nrow()) {
    stop("The number of rows in 'phi' must be equal to the number of rows in 'counts'.");
  }
  if(counts.rows() != eta.nrow()) {
    stop("The number of rows in 'eta' must be equal to the number of rows in 'counts'.");
  }
  if(phi.ncol() != eta.ncol()) {
    stop("The number of columns in 'eta' must be equal to the number of columns in 'phi'.");
  }
  if(min(z) < 1 || max(z) > eta.ncol()) {
    stop("The entries in 'z' need to be between 1 and the number of columns in eta and phi.");
  }
  
  // Declare variables and functions
  double loglik = 0;    
  
  int i;
  int k; 
  double x;
  int nr = phi.nrow();
  
  // Original R code:
  // ll <- sum(Matrix::t(counts) * log(theta * t(phi)[z, ] +
  //       (1 - theta) * t(eta)[z, ] + 1e-20))
  
  for(int j = 0; j < counts.cols(); ++j) {
    for (Eigen::MappedSparseMatrix<double>::InnerIterator i_(counts, j); i_; ++i_) {
      i = i_.index();
      x = i_.value();
      k = z[j] - 1;
      
      loglik += x * log((phi[nr * k + i] * theta[j]) + (eta[nr * k + i] * (1 - theta[j])) + pseudocount);
    }  
  }
  
  return loglik;
}  



// [[Rcpp::export]]
Rcpp::List decontXInitialize(const Eigen::MappedSparseMatrix<double>& counts,
                             const NumericVector& theta,
                             const IntegerVector& z,
                             const double& pseudocount) {
  
  // Perform error checking
  if(counts.cols() != theta.size()) {
    stop("Length of 'theta' must be equal to the number of columns in 'counts'.");
  }
  if(counts.cols() != z.size()) {
    stop("Length of 'z' must be equal to the number of columns in 'counts'.");
  }  
  
  // Declare variables and functions
  NumericMatrix new_phi(counts.rows(), max(z));
  NumericMatrix new_eta(counts.rows(), max(z));  
  std::fill(new_phi.begin(), new_phi.end(), pseudocount);
  std::fill(new_eta.begin(), new_eta.end(), pseudocount);
  
  int k;
  int i;
  double x; 
  for(int j = 0; j < counts.cols(); ++j) {
    for (Eigen::MappedSparseMatrix<double>::InnerIterator i_(counts, j); i_; ++i_) {
      i = i_.index();
      x = i_.value();
      k = z[j] - 1;
      
      new_phi(i,k) += x * theta(j);
    }  
  }
  
  // Calculate Eta using Weights from Phi
  NumericVector phi_rowsum = rowSums(new_phi);
  int j;
  for(i = 0; i < new_eta.ncol(); i++) {
    for(j = 0; j < new_eta.nrow(); j++) {
      new_eta(j,i) = phi_rowsum[j] - new_phi(j,i);
    }
  }
  
  // Normalize Phi and Eta
  NumericVector phi_colsum = colSums(new_phi);
  NumericVector eta_colsum = colSums(new_eta);  
  for(i = 0; i < new_phi.ncol(); i++) {
    new_phi(_,i) = new_phi(_,i) / phi_colsum[i];
    new_eta(_,i) = new_eta(_,i) / eta_colsum[i];
  }
  
  return Rcpp::List::create(Rcpp::Named("phi") = new_phi,
                            Rcpp::Named("eta") = new_eta);
  
}  





// [[Rcpp::export]]
Eigen::SparseMatrix<double> calculateNativeMatrix(const Eigen::MappedSparseMatrix<double>& counts,
                             const NumericVector& theta,
                             const NumericMatrix& eta,
                             const NumericMatrix& phi,
                             const IntegerVector& z,
                             const double& pseudocount) {
  
  // Perform error checking
  if(counts.cols() != theta.size()) {
    stop("Length of 'theta' must be equal to the number of columns in 'counts'.");
  }
  if(counts.cols() != z.size()) {
    stop("Length of 'z' must be equal to the number of columns in 'counts'.");
  }  
  if(counts.rows() != phi.nrow()) {
    stop("The number of rows in 'phi' must be equal to the number of rows in 'counts'.");
  }
  if(counts.rows() != eta.nrow()) {
    stop("The number of rows in 'eta' must be equal to the number of rows in 'counts'.");
  }
  if(phi.ncol() != eta.ncol()) {
    stop("The number of columns in 'eta' must be equal to the number of columns in 'phi'.");
  }
  if(min(z) < 1 || max(z) > eta.ncol()) {
    stop("The entries in 'z' need to be between 1 and the number of columns in eta and phi.");
  }

  Eigen::SparseMatrix<double> native_matrix = counts;
  
  int i;
  int k; 
  double x; 
  double pcontamin;
  double pnative;
  double normp;
  for(int j = 0; j < counts.cols(); ++j) {
    for (Eigen::MappedSparseMatrix<double>::InnerIterator i_(counts, j); i_; ++i_) {
      i = i_.index();
      x = i_.value();
      k = z[j] - 1;
      
      // Calculate variational probabilities 
      pnative = log(phi(i,k) + pseudocount) + log(theta(j) + pseudocount);
      pcontamin = log(eta(i,k) + pseudocount) + log(1 - theta(j) + pseudocount);
      
      // Normalize probabilities and add to proper components
      normp = exp(pnative) / (exp(pcontamin) + exp(pnative));
      native_matrix.coeffRef(i, j) *= normp;
    }  
  }

  return native_matrix;
}



// [[Rcpp::export]]
Rcpp::List decontXEM_fixEta(const Eigen::MappedSparseMatrix<double>& counts,
                     const NumericVector& counts_colsums,
                     const NumericVector& theta,
                     const NumericVector& eta,
                     const NumericMatrix& phi,
                     const IntegerVector& z,
                     const bool& estimate_delta,
                     const NumericVector& delta,
                     const double& pseudocount) {

  // Perform error checking
  if(counts.cols() != theta.size()) {
    stop("Length of 'theta' must be equal to the number of columns in 'counts'.");
  }
  if(counts.cols() != z.size()) {
    stop("Length of 'z' must be equal to the number of columns in 'counts'.");
  }  
  if(counts.cols() != counts_colsums.size()) {
    stop("Length of 'counts_colsums' must be equal to the number of columns in 'counts'.");
  }    
  if(counts.rows() != phi.nrow()) {
    stop("The number of rows in 'phi' must be equal to the number of rows in 'counts'.");
  }
  if(counts.rows() != eta.size()) {
    stop("The size of vector 'eta' must be equal to the number of rows in 'counts'.");
  }
  if(min(z) < 1 || max(z) > phi.ncol()) {
    stop("The entries in 'z' need to be between 1 and the number of columns in eta and phi.");
  }
  if(delta.size() != 2 || sum(delta < 0) > 0) {
    stop("'delta' must be a numeric vector of length 2 with positive integers.");
  }
  
  // Declare variables and functions
  NumericVector new_theta(theta.size());
  NumericVector native_total(theta.size());
  NumericMatrix new_phi(phi.nrow(), phi.ncol());
  
    
  // Obtaining 'fit_dirichlet' function from MCMCprecision package
  Environment pkg = Environment::namespace_env("MCMCprecision");
  Function f = pkg["fit_dirichlet"];

  int i;
  int j;
  int k;
  int nr = phi.nrow();
  double x; 
  double pcontamin;
  double pnative;
  double normp;
  double px;
  for(int j = 0; j < counts.cols(); ++j) {
    for (Eigen::MappedSparseMatrix<double>::InnerIterator i_(counts, j); i_; ++i_) {
      i = i_.index();
      x = i_.value();
      k = z[j] - 1;
      
      // Calculate variational probabilities
      // Removing the log/exp speeds it up and produces the same result since
      // there are only 2 probabilities being multiplied
      
      //pnative = log(phi(i,k) + pseudocount) + log(theta(j) + pseudocount);
      //pcontamin = log(eta(i,k) + pseudocount) + log(1 - theta(j) + pseudocount);
      pnative = (phi[nr * k + i] + pseudocount) * (theta[j] + pseudocount);
      pcontamin = (eta[i] + pseudocount) * (1 - theta[j] + pseudocount);
      
      // Normalize probabilities and add to proper components
      //normp = exp(pnative) / (exp(pcontamin) + exp(pnative));
      normp = pnative / (pcontamin + pnative);
      px = normp * x;
      new_phi(i,k) += px;
      native_total(j) += px;
    }  
  }
  
  
  // Normalize Phi
  NumericVector phi_colsum = colSums(new_phi);
  for(i = 0; i < new_phi.ncol(); i++) {
    new_phi(_,i) = new_phi(_,i) / phi_colsum[i];
  }
  
  // Update Theta
  NumericVector contamination_prop = (counts_colsums - native_total) / counts_colsums;
  NumericVector native_prop = 1 - contamination_prop;
  NumericMatrix theta_raw = cbind(native_prop, contamination_prop);
  
  NumericVector new_delta = delta;
  if(estimate_delta == TRUE) {
    Rcpp::List result = f(Named("x", theta_raw)); 
    new_delta = result["alpha"];
  }
  
  // Estimate new theta
  new_theta = (native_total + new_delta[0]) / (counts_colsums + sum(new_delta));
  
  return Rcpp::List::create(Rcpp::Named("phi") = new_phi,
                            Rcpp::Named("eta") = eta,
                            Rcpp::Named("theta") = new_theta,
                            Rcpp::Named("delta") = new_delta,
                            Rcpp::Named("contamination") = contamination_prop);
}

