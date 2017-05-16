
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


#' Visualize various performance metrics as a function of K / L to aid parameter choice.
#' 
#' @param celda.list A celda_list object as returned from *celda()*
#' @param metric Which performance metric to visualize. One of ("perplexity", "harmonic", "loglik"). "perplexity" calculates the inverse of the geometric mean of the log likelihoods from each iteration of Gibbs sampling. "harmonic" calculates the marginal likelihood has the harmonic mean of the likelihoods. "loglik" plots the highest log-likelihood during Gibbs iteration.
#' @return A ggplot object containing the requested plot(s), or a list of ggplots if the provided celda_list contains celda_CG models.
#' @export
visualize_performance = function(celda.list, method="perplexity") {
  # TODO use celda_list getter
  log.likelihoods = lapply(celda.list$res.list,
                           function(mod) { return(mod$completeLogLik) })
    
  if (method == "perplexity") {
    metric = lapply(log.likelihoods, calculate_perplexity)
    metric = new("mpfr", unlist(metric))
  } else if (method == "harmonic") {
    metric = lapply(log.likelihoods, calculate_marginal_likelihood)
    metric = new("mpfr", unlist(metric))
  } else if (method == "loglik") {
    # TODO use celda_list getter
    metric = lapply(log.likelihoods, max)
  } else stop("Invalid method specified")
  
  
  #TODO ggplot table building below is vulnerable to error 
  #     if the user modifies the celda_list at all..
  if (celda.list$content.type == "celda_C") {
    Ks = lapply(celda.list$res.list, function(mod) { mod$K })
    plot.df = data.frame(K=as.factor(unlist(Ks)), metric=as.numeric(metric))
    ggplot2::ggplot(plot.df, aes(x=K, y=metric)) + geom_point() +
      xlab("K") + ylab(method)
  }
  
}





