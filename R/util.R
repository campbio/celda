
#' Calculate the marginal likelihood from a single celda chain.
#' Marginal likelihood is estimated as the harmonic mean of the 
#' (non-log) likelihood over all iterations of Gibbs sampling.
#' 
#' @param completeLogLik The complete Gibbs sampling history of log-likelihoods for a single celda chain
#' @return The estimated marginal likelihood as an mpfr number
#' @export
calculate_marginal_likelihood = function(completeLogLik) {
  mpfr_log_lik = Rmpfr::mpfr(completeLogLik, 512)
  complete_likelihood = exp(mpfr_log_lik)
  marginal_likelihood = (Rmpfr::mean((1/complete_likelihood)))^-1
  return(marginal_likelihood)
}


#' Calculate the perplexity from a single celda chain.
#' Perplexity is defined as the inverse of the geometric mean of the 
#' 
#' 
#' @param completeLogLik The complete Gibbs sampling history of log-likelihoods for a single celda chain
#' @return The perplexity for the provided chain as an mpfr number
#' @export
calculate_perplexity = function(completeLogLik) {
  mpfr_log_lik = Rmpfr::mpfr(completeLogLik, 512)
  perplexity = exp(Rmpfr::mean(mpfr_log_lik))^-1
  return(perplexity)
}
