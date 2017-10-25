#' Calculate the perplexity from a single celda chain
#' 
#' Perplexity is defined as the inverse of the geometric mean of the log-likelihoods over all 
#' iterations of Gibbs sampling.
#' @param completeLogLik The complete Gibbs sampling history of log-likelihoods for a single celda chain
#' @param log Set log to TRUE to visualize the log(perplexity) of Celda_CG objects.
#' @return The perplexity for the provided chain as an mpfr number
#' @export
calculate_perplexity = function(completeLogLik, log = FALSE) {
  if (log) {
    return(-mean(completeLogLik))
  }
  mpfr_log_lik = Rmpfr::mpfr(completeLogLik, 512)
  perplexity = exp(Rmpfr::mean(mpfr_log_lik))^-1
  return(perplexity)
}


# Convenience function to calculate performance metrics by specifying a method. 
calculate_performance_metric = function(log.likelihoods, method="perplexity", log = FALSE) {
  if (method == "perplexity") {
    metric = calculate_perplexity(log.likelihoods, log)
  } else if (method == "loglik") {
     metric = max(log.likelihoods)
  } else stop("Invalid method specified")
  return(metric)
}


# Actually render the plot described in visualizeModelPerformance.
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


#' Plots the log likehood over iteration for each chain 
#' 
#' @param celda.res An object returned from the 'celda' function.
#' @param group.by A vector giving the variables to group the chains by.
#' @param exclude.iter Integers giving the indicies of iterations to excluce in the plotting. Default 1:10.
#' @param line.size Size of the lines for each chain.
#' @param legend.position Position of the color legend for the chains. Default "none".
#' @param scales Are scales shared across all facets (the default, ‘"fixed"’), or do they vary across rows (‘"free_x"’), columns (‘"free_y"’), or both rows and columns (‘"free"’). See ?ggplot2::facet_grid for more information.
#' @export
celdaLogLik = function(celda.res, group.by, exclude.iter = 1:10, line.size=0.5, legend.position = c("none", "left", "right", "bottom", "top"), scales = c("free_y", "free_x", "free", "fixed")) {

  legend.position = match.arg(legend.position)
  scales = match.arg(scales) 
  
  ## Set up groupings
  rp = celda.res$run.params
  chain = rp$chain
    
  if(!all(group.by %in% colnames(rp))) {
    stop("Items in 'group.by' must match column names in 'run.params' object")
  }
  rp = rp[,group.by,drop=FALSE]

  group = sapply(1:ncol(rp), function(i) paste0(colnames(rp)[i], "-", rp[,i]))
  group = apply(group, 1, paste, collapse=";")

  # Extract log likelihoods from each model
  ll = sapply(1:length(celda.res$res.list), function(i) { celda.res$res.list[[i]]$completeLogLik } )  
  rownames(ll) = 0:(nrow(ll)-1)
  if(!is.null(exclude.iter)) {
    ll = ll[-exclude.iter,,drop=FALSE]
  }
  
  ## Combine into data frame and plot
  df = data.frame(Group=group, t(ll), Chain=as.factor(chain), check.names=FALSE)
  df.m = reshape2::melt(df, id=c("Group", "Chain"), value.name="Log_Likelihood", variable.name="Iteration")
  df.m$Iteration = as.numeric(as.character(df.m$Iteration))
  gg = ggplot2::ggplot(df.m, ggplot2::aes(x=Iteration, y=Log_Likelihood, color=Chain)) +
       ggplot2::geom_line(size=line.size) +
       ggplot2::facet_wrap(~ Group, scales=scales) +     
       ggplot2::theme_bw() + ggplot2::theme(panel.spacing = ggplot2::unit(0,"lines"), legend.position=legend.position) + 
       ggplot2::ylab("Log likelihood") + ggplot2::xlab("Iteration") 
  
  return(gg)
}  



