#' Calculate the perplexity from a single celda run
#' 
#' Perplexity can be seen as a measure of how well a provided set of 
#' cluster assignments fit the count data being modeled.
#' 
#' Perplexity is defined in LDA as the exp of the log likelihood of the model
#' divided by the total amount of word tokens. The corresponding perplexity
#' for the celda models are derived in their respective model description 
#' documents. These documents will not be made publicly available until 
#' after publication of celda.
#' 
#' @param celda.mod A single celda run (usually from the _res.list_ property of a celda_list)
#' @param counts The count matrix modeled in the celdaRun parameter
#' @param resample The number of resamplings of the counts matrix to calculate perplexity for
#' @param precision The amount of bits of precision to pass to Rmpfr
#' @return The perplexity for the provided chain as an mpfr number
#' @export
calculatePerplexity = function(celda.mod, counts, resample=1, precision=128) {
  UseMethod("calculatePerplexity", celda.mod)
}


# Convenience function to calculate performance metrics by specifying a method. 
calculatePerformanceMetric = function(celdaRun, counts, log.likelihoods, 
                                      method="perplexity", log = FALSE,
                                      resample=1) {
  if (method == "perplexity") {
    metric = calculatePerplexity(celdaRun, counts, resample)
  } else if (method == "loglik") {
     metric = max(log.likelihoods)
  } else stop("Invalid method specified")
  return(metric)
}


# Actually render the plot described in visualizeModelPerformance.
renderModelPerformancePlot = function(cluster.scores, cluster.label, metric.type,
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
renderIterationLikelihoodPlot = function(celda.obj) {
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


# Resample a counts matrix for evaluating perplexity
#
# Normalizes each column (cell) of a count matrix by the column sum to 
# create a distribution of observing a given number of counts for a given gene in that cell,
# then samples across all cells.
#
# This is primarily used to evaluate the stability of the perplexity for a given K/L combination.
# 
# @param celda.mod A single celda run (usually from the _res.list_ property of a celda_list)
# @return The perplexity for the provided chain as an mpfr number
resampleCountMatrix = function(count.matrix) {
  colsums  = colSums(count.matrix)
  prob     = t(t(count.matrix) / colsums)
  resample = sapply(1:ncol(count.matrix), function(idx){
                      rmultinom(n=1, size=colsums[idx], prob=prob[, idx])
                   })
  return(resample)
}


