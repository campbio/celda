#' Identify features with the highest influence on clustering
#' 
#' @param fm factorized matrix.
#' @param n Maximum number of items returned for each entry. 
#' @param margin 1 for rows, 2 for columns.
#' @param threshold only include entries in the matrix that is greader than the threshold.
#' @param decreasing Logical; specifying if rank should be decreasing. Default to be TRUE. 
#' @export
topRank = function(fm, n=25, margin=2, threshold=0, decreasing=TRUE) {
  if(is.null(threshold) | is.na(threshold)) {
    threshold = min(fm) - 1 
  }
  
  ## Function to sort values in a vector and return 'n' top results
  ## If there are not 'n' top results above 'thresh', then the
  ## number of entries in 'v' that are above 'thresh' will be returned
  topFunction = function(v, n, thresh) {
    v.above.thresh = sum(v > thresh)
    n.to.select = min(v.above.thresh, n)
    
    h = NA
    if(n.to.select > 0) {
      h = head(order(v, decreasing=decreasing), n.to.select)
    }  
    return(h)
  }

  ## Parse top ranked indices from matrix
  top.ix = base::apply(fm, margin, topFunction, thresh=threshold, n=n)
  
  ## Convert to list if apply converted to a matrix because all
  ## elements had the same length
  if(is.matrix(top.ix)) {
    top.ix = lapply(1:ncol(top.ix), function(i) top.ix[,i]) 
    names(top.ix) = dimnames(fm)[[margin]]
  }  
  
  ## Parse names from returned margin
  opposite.margin = ifelse(margin - 1 > 0, margin-1, length(dim(fm)))
  top.names = NULL
  names.to.parse = dimnames(fm)[[opposite.margin]]
  if(!is.null(names.to.parse) & all(!is.na(top.ix))) {
    top.names = lapply(1:length(top.ix),
    			function(i) {
    			  ifelse(is.na(top.ix[[i]]), NA, names.to.parse[top.ix[[i]]])
    			})
    names(top.names) = names(top.ix)
  }
  
  return(list(index=top.ix, names=top.names))
}

#' @title Create a bar chart based on Gini coefficients of transcriptional states
#' @description A bar plot shows the ordered Gini coefficients of transcriptional states. The transcriptional states with high Gini coefficients (i.e. close to 1) 
#'    indicates that only a few cell clusters highly express the genes in these transcriptional states. It is easier for users to define cell types based on these
#'    genes in these significant states.
#'    
#' @param counts A numeric count matrix.
#' @param celda.mod An object of class celda_CG.
#' @param cell_clusters A list of cell clusters. Calculate Gini coefficients using cells from specific cell clusters. If NULL, uses all cells (NULL by default).
#' @param label_size Numeric. The size of labels of transcrptional states on top of bars (4 by default).
#' @param bar_col The color of bars (#FF8080FF by default).
#' @export
GiniPlot <- function(counts, celda.mod, cell_clusters = NULL, label_size = 4, bar_col = NULL) {
  compareCountMatrix(counts, celda.mod)
  
  factorize.matrix <-
    factorizeMatrix(counts = counts, celda.mod = celda.mod)
  if(is.null(bar_col)){
    bar_col <- "#FF8080FF"
  }
  if (is.null(cell_clusters)) {
    cell_clusters <- sort(unique(celda.mod$z))
  }
  gini <-
    apply(factorize.matrix$proportions$population.states[, cell_clusters],
          1,
          ineq::Gini)
  sorted_gini <- sort(gini, decreasing = TRUE)
  sorted_gini_states <- as.numeric(substring(names(sorted_gini), 2))
  df <-
    data.frame(
      Gini_Coefficient = sorted_gini,
      Transcriptional_States = factor(sorted_gini_states, levels = sorted_gini_states)
    )
  ggplot2::ggplot(df, ggplot2::aes(Transcriptional_States, Gini_Coefficient)) +
    ggplot2::geom_col(fill = bar_col) + ggplot2::geom_text(ggplot2::aes(label = Transcriptional_States), position = ggplot2::position_dodge(0.9),
                                                           vjust = 0, size = label_size) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black")
    )
}
