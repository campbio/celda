#' Calculate the marginal likelihood from a single celda model
#' 
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


#' Calculate the perplexity from a single celda chain
#' 
#' Perplexity is defined as the inverse of the geometric mean of the log-likelihoods over all 
#' iterations of Gibbs sampling.
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
           ggplot2::geom_point(position=ggplot2::position_jitter(width=0.1, height=0)) +
           ggplot2::xlab(cluster.label) + ggplot2::ylab(metric.type) + 
           ggplot2::ggtitle(title) + ggplot2::theme_bw()
  return(plot)
}


#' Visualize log likelihood per iteration of Gibbs sampling
#' 
#' @param celda.obj A celda object
#' @return A ggplot object of the likelihood per iteration of Gibbs sampling
#' @export
render_iteration_likelihood_plot = function(celda.obj) {
  df = data.frame(celda.obj$completeLogLik)
  df["iter"] = factor(as.numeric(rownames(df)))
  plot = ggplot2::ggplot(df, ggplot2::aes(x=iter, y=df)) + ggplot2::geom_point() +
    ggplot2::theme_bw() + ggplot2::ylab("Log likelihood") + ggplot2::xlab("Iterations") +
    ggplot2::ggtitle("Log Likelihood Per Iteration")
  return(plot)
}


#' Visualize the performance of celda_CG models
#' 
#' Plot the performance of a list of celda_CG models returned from running celda function.
#' For each number of gene clusters L (cell clusters K), plot the performance of each number of cell
#' clusters K (gene clusters L).
#' @param celda_CG.list A list of celda_CG objects
#' @return A grob object (gridExtra) of the model performance
#' @export
render_celda_CG_list_performance = function(celda_CG.list, by = "L", method = "perplexity") {
  # validate input parameters
  if (class(small.sim.res) != "celda_list") {
    stop("celda_CG.list argument must be of class 'celda_list'")
  } else if (celda_CG.list$content.type != "celda_CG") {
    stop("celda_CG.list must be a 'celda_list' of 'celda_CG' objects")
  } else if (!(by %in% c("K","L"))) {
    stop("'by' has to be either 'L' or 'K'")
  } else if (!(method %in% c("perplexity","harmonic","loglik"))) {
    stop("Invalid method, 'method' has to be either 'perplexity', 'harmonic', or 'loglik'")
  } 
  # if every thing went well
  else {
    # plot title
    title=
    y.lab = method
    
    cluster.sizes = unlist(lapply(celda_CG.list$res.list, function(mod) { getK(mod) }))
    log.likelihoods = lapply(celda_CG.list$res.list,
                             function(mod) { completeLogLikelihood(mod) })
    performance.metric = lapply(log.likelihoods, 
                                calculate_performance_metric,
                                method)
    
    # These methods return Rmpfr numbers that are extremely small and can't be 
    # plotted, so log 'em first
    if (method %in% c("perplexity", "harmonic")) {
      performance.metric = lapply(performance.metric, log)
      performance.metric = new("mpfr", unlist(performance.metric))
      performance.metric = as.numeric(performance.metric)
      y.lab = paste0("Log(",method,")")
    } else {
      performance.metric = as.numeric(performance.metric)
    }
    
    plot.df = data.frame(K=cluster.sizes, L=celda_CG.list$run.params$L,
                         metric=performance.metric)
    
    L.list = sort(unique(celda_CG.list$run.params$L))
    K.list = sort(unique(celda_CG.list$run.params$K))
    #chains = sort(unique(celda_CG.list$run.params$chain))
    
    plots = list()
    
    if (by=="L") {
      nc = round(length(L.list)^.5)
      x.lab = "K"
      
      for (i in L.list) {
        
        plots = c(plots, list(ggplot2::ggplot(subset(plot.df, L==i), 
                                              ggplot2::aes(x=K, y=metric, group=K)) + 
                                ggplot2::geom_boxplot(outlier.color=NA, fill=NA) + 
                                ggplot2::geom_point(position=ggplot2::position_jitter(width=0.1, height=0)) +
                                ggplot2::xlab("K") + ggplot2::ylab(method) + 
                                ggplot2::ggtitle(paste0("L = ", i)) + 
                                ggplot2::theme_bw() + 
                                ggplot2::theme(axis.title.x=element_blank(), 
                                               axis.title.y=element_blank())))
        
      }
      
    } else if (by=="K") {
      nc = round(length(K.list)^.5)
      x.lab = "L"
      
      for (i in K.list) {
        
        plots = c(plots, list(ggplot2::ggplot(subset(plot.df, K==i), 
                                              ggplot2::aes(x=L, y=metric, group=L)) + 
                                ggplot2::geom_boxplot(outlier.color=NA, fill=NA) + 
                                ggplot2::geom_point(position=ggplot2::position_jitter(width=0.1, height=0)) + 
                                ggplot2::ggtitle(paste0("K = ", i)) + 
                                ggplot2::theme_bw() + 
                                ggplot2::theme(axis.title.x=element_blank(), 
                                               axis.title.y=element_blank())))
        
      }
    }
    grid.arrange(grobs=plots, ncol=nc, left = textGrob(y.lab, rot = 90), 
                 top = textGrob("Model Performance (All Chains)"),
                 bottom=textGrob(x.lab))
  }
}

