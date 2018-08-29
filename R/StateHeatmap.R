#' @title Transcriptional state heatmap
#' @description Draws a heatmap focusing on a transcriptional state. Both cells and genes are sorted by 
#'    their proportions of counts in a given transcriptional state. Allows for nice visualization of 
#'    co-expression of those genes grouped into transcriptional states by Celda.    
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class "celda_G" or "celda_CG". 
#' @param feature.module Integer. The feature module to display.
#' @param top.cells Integer. Number of cells with the highest and lowest probabilities for this module to plot. For example, if `top.cells` = 50, the 50 cells with the lowest probability and the 50 cells with the highest probability for that feature module will be plotted. If NULL, all cells will be plotted. Default NULL. 
#' @param top.features Integer. Plot `top.features` with the highest probability in the feature module. If NULL, plot all features in the module. Default NULL.
#' @param normalize Logical. Whether to normalize the columns of `counts`. Default TRUE. 
#' @param scale.row Character. Which function to use to scale each individual row. Set to NULL to disable. Occurs after normalization and log transformation. 'scale' will Z-score transform each row. Default 'scale'.
#' @param show_featurenames Logical. Specifies if feature names should be shown. Default TRUE. 
#' @export 
moduleHeatmap <- function(counts, celda.mod, feature.module = 1, top.cells = NULL, top.features = NULL, normalize = TRUE, scale.row = scale, show_featurenames = TRUE){
  #input checks
  if (is.null(counts)) {
    stop("'counts' should be a numeric count matrix")
  }
  if (is.null(celda.mod) || is.null(celda.mod$y)){
    stop("'celda.mod' should be an object of class celda_G or celda_CG")
  }
  compareCountMatrix(counts, celda.mod)
  
  #factorize counts matrix
  factorize.matrix <-
    factorizeMatrix(celda.mod = celda.mod, counts = counts)
  
  #take topRank
  if(!is.null(top.features) && (is.numeric(top.features))|is.integer(top.features)){
    top.rank <- topRank(matrix = factorize.matrix$proportions$gene.states, 
                        n = top.features)
  }else{
    top.rank <- topRank(matrix = factorize.matrix$proportions$gene.states, 
                        n = nrow(factorize.matrix$proportions$gene.states))  
  }
  
  #filter topRank using feature.module into feature.indices
  feature.indices <- c()
  for(module in feature.module){
    feature.indices <- c(feature.indices, top.rank$index[[module]])
  }
  
  
  #Determine cell order from factorize.matrix$proportions$cell.states
  cell.states <- factorize.matrix$proportions$cell.states 
  cell.states <- cell.states[feature.module, , drop = FALSE]
  
  if(is.null(top.cells)){
    top.cells <- ncol(cell.states)
  }
  cell.indices <- c()
  for(modules in 1:nrow(cell.states)){
    single.module <- cell.states[modules, ]
    single.module.ordered <- order(single.module, decreasing = TRUE)
    cell.indices <- c(cell.indices, head(single.module.ordered, n = top.cells), tail(single.module.ordered, n = top.cells))
  }
  cell.indices <- cell.indices[!duplicated(cell.indices)] 
  
  
  #normalize counts matrix
  if(normalize){
    norm.counts <- normalizeCounts(counts, normalize="proportion", transformation.fun=sqrt)
  } else{
    norm.counts <- counts
  }
  
  #filter counts based on feature.indices 
  filtered_norm.counts <- norm.counts[feature.indices, cell.indices]
  
  filtered_norm.counts <- filtered_norm.counts[rowSums(filtered_norm.counts>0)>0, ]
  
  if(!is.null(celda.mod$z)){
    cell <- distinct_colors(length(unique(celda.mod$z)))[sort(unique(celda.mod$z[cell.indices]))]
    names(cell) <- sort(unique(celda.mod$z[cell.indices]))
    anno_cell_colors <- list(cell = cell)
  }else{
    anno_cell_colors <- NULL
  }
  plotHeatmap(
    filtered_norm.counts,
    z = celda.mod$z[cell.indices],
    y = celda.mod$y[feature.indices],
    scale.row = scale.row,
    color_scheme = "divergent",
    show_featurenames = show_featurenames,
    cluster_gene = FALSE,
    cluster_cell = FALSE,
    annotation_color = anno_cell_colors
  )
}
