
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


# Convenience function to calculate performance metrics by specifying a method. 
calculate_performance_metric = function(log.likelihoods, method="perplexity") {
    if (method == "perplexity") {
    metric = calculate_perplexity(log.likelihoods)
  } else if (method == "harmonic") {
    metric = calculate_marginal_likelihood(log.likelihoods)
  } else if (method == "loglik") {
     metric = max(log.likelihoods)
  } else stop("Invalid method specified")
  return(metric)
}


# Actually render the plot described in visualize_model_performance.
render_model_performance_plot = function(cluster.scores, cluster.label, metric.type,
                                         title="Model Performance (All Chains") {
  plot = ggplot2::ggplot(cluster.scores, ggplot2::aes(x=factor(size), y=metric)) + 
           ggplot2::geom_boxplot(outlier.color=NA, fill=NA) + 
           ggplot2::geom_point(position=ggplot2::position_jitter(width=0.1)) +
           ggplot2::xlab(cluster.label) + ggplot2::ylab(metric.type) + 
           ggplot2::ggtitle(title) + ggplot2::theme_bw()
  return(plot)
}
