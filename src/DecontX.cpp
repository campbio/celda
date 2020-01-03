// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List decontXEM(const arma::sp_mat& counts,
                     const NumericVector& counts_colsums,
                     const NumericVector& theta,
                     const NumericMatrix& eta,
                     const NumericMatrix& phi,
                     const IntegerVector& z,
                     const double& pseudocount) {

  // Perform error checking
  if(counts.n_cols != theta.size()) {
    stop("Length of 'theta' must be equal to the number of columns in 'counts'.");
  }
  if(counts.n_cols != z.size()) {
    stop("Length of 'z' must be equal to the number of columns in 'counts'.");
  }  
  if(counts.n_cols != counts_colsums.size()) {
    stop("Length of 'counts_colsums' must be equal to the number of columns in 'counts'.");
  }    
  if(counts.n_rows != phi.nrow()) {
    stop("The number of rows in 'phi' must be equal to the number of rows in 'counts'.");
  }
  if(counts.n_rows != eta.nrow()) {
    stop("The number of rows in 'eta' must be equal to the number of rows in 'counts'.");
  }
  if(phi.ncol() != eta.ncol()) {
    stop("The number of columns in 'eta' must be equal to the number of columns in 'phi'.");
  }
  if(min(z) < 1 || max(z) > eta.ncol()) {
    stop("The entries in 'z' need to be between 1 and the number of columns in eta and phi.");
  }

  // Declare variables and functions
  NumericVector new_theta(theta.size());
  NumericVector native_total(theta.size());
  NumericMatrix new_phi(phi.nrow(), phi.ncol());
  NumericMatrix new_eta(eta.nrow(), eta.ncol());
    
  std::fill(new_phi.begin(), new_phi.end(), pseudocount);
  std::fill(new_eta.begin(), new_eta.end(), pseudocount);
    
  // Obtaining 'fit_dirichlet' function from MCMCprecision package
  Environment pkg = Environment::namespace_env("MCMCprecision");
  Function f = pkg["fit_dirichlet"];

  int i;
  int j;
  int k; 
  double x; 
  double pcontamin;
  double pnative;
  double normp;
  double px;
  for (arma::sp_mat::const_iterator it = counts.begin(); it != counts.end(); ++it) {
    i = it.row();
    j = it.col();
    x = *it;
    k = z[j] - 1;
    
    // Calculate variational probabilities 
    pnative = log(phi(i,k) + pseudocount) + log(theta(j) + pseudocount);
    pcontamin = log(eta(i,k) + pseudocount) + log(1 - theta(j) + pseudocount);
    
    // Normalize probabilities and add to proper components
    normp = exp(pnative) / (exp(pcontamin) + exp(pnative));
    px = normp * x;
    new_phi(i,k) += px;
    native_total(j) += px;
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
  
  Rcpp::List result = f(Named("x", theta_raw)); 
  NumericVector delta = result["alpha"];
  new_theta = (native_total + delta[0]) / (counts_colsums + result["sum"]);
  
  return Rcpp::List::create(Rcpp::Named("phi") = new_phi,
                            Rcpp::Named("eta") = new_eta,
                            Rcpp::Named("theta") = new_theta,
                            Rcpp::Named("delta") = delta,
                            Rcpp::Named("contamination") = contamination_prop);
}




// [[Rcpp::export]]
double decontXLogLik(const arma::sp_mat& counts,
                     const NumericVector& theta,
                     const NumericMatrix& eta,
                     const NumericMatrix& phi,
                     const IntegerVector& z,
                     const double& pseudocount) {

  // Perform error checking
  if(counts.n_cols != theta.size()) {
    stop("Length of 'theta' must be equal to the number of columns in 'counts'.");
  }
  if(counts.n_cols != z.size()) {
    stop("Length of 'z' must be equal to the number of columns in 'counts'.");
  }  
  if(counts.n_rows != phi.nrow()) {
    stop("The number of rows in 'phi' must be equal to the number of rows in 'counts'.");
  }
  if(counts.n_rows != eta.nrow()) {
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
  int j;
  int k; 
  double x; 
  
  // Original R code:
  // ll <- sum(Matrix::t(counts) * log(theta * t(phi)[z, ] +
  //       (1 - theta) * t(eta)[z, ] + 1e-20))
  
  for (arma::sp_mat::const_iterator it = counts.begin(); it != counts.end(); ++it) {
    i = it.row();
    j = it.col();
    x = *it;
    k = z[j] - 1;
    
    loglik += x * log((phi(i,k) * theta(j)) + (eta(i,k) * (1 - theta(j))) + pseudocount);    
  }
    
  return loglik;
}  





// [[Rcpp::export]]
Rcpp::List decontXInitialize(const arma::sp_mat& counts,
                     const NumericVector& theta,
                     const IntegerVector& z,
                     const double& pseudocount) {

  // Perform error checking
  if(counts.n_cols != theta.size()) {
    stop("Length of 'theta' must be equal to the number of columns in 'counts'.");
  }
  if(counts.n_cols != z.size()) {
    stop("Length of 'z' must be equal to the number of columns in 'counts'.");
  }  

  // Declare variables and functions
  NumericMatrix new_phi(counts.n_rows, max(z));
  NumericMatrix new_eta(counts.n_rows, max(z));  
  std::fill(new_phi.begin(), new_phi.end(), pseudocount);
  std::fill(new_eta.begin(), new_eta.end(), pseudocount);

  int i;
  int j;
  int k; 
  double x; 
    
  for (arma::sp_mat::const_iterator it = counts.begin(); it != counts.end(); ++it) {
    i = it.row();
    j = it.col();
    x = *it;
    k = z[j] - 1;
    
    new_phi(i,k) += x * theta(j);
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
    
  return Rcpp::List::create(Rcpp::Named("phi") = new_phi,
                            Rcpp::Named("eta") = new_eta);

}  

