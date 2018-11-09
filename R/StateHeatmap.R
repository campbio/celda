#' @title Gene module heatmap
#' @description Draws a heatmap focusing on a gene module. Both cells and genes are sorted by 
#'    their proportions of counts in a given gene module. Allows for nice visualization of 
#'    co-expression of those genes grouped into gene modules by Celda.    
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class "celda_G" or "celda_CG". 
#' @param feature.module Integer. The feature module to display.
#' @param top.cells Integer. Number of cells with the highest and lowest probabilities for this module to plot. For example, if `top.cells` = 50, the 50 cells with the lowest probability and the 50 cells with the highest probability for that feature module will be plotted. If NULL, all cells will be plotted. Default NULL. 
#' @param top.features Integer. Plot `top.features` with the highest probability in the feature module. If NULL, plot all features in the module. Default NULL.
#' @param normalize Logical. Whether to normalize the columns of `counts`. Default TRUE. 
#' @param scale.row Character. Which function to use to scale each individual row. Set to NULL to disable. Occurs after normalization and log transformation. 'scale' will Z-score transform each row. Default 'scale'.
#' @param show_featurenames Logical. Specifies if feature names should be shown. Default TRUE. 
#' @return A list containing row and column dendrogram information, as well as a gtable for grob plotting
#' @examples
#' celda.sim = simulateCells("celda_CG")
#' celda.mod = celda_CG(celda.sim$counts, K=celda.sim$K, L=celda.sim$L,
#'                      nchains=1, max.iter=1)
#' moduleHeatmap(celda.sim$counts, celda.mod)
#' @export 
moduleHeatmap <- function(counts, celda.mod, feature.module = 1, top.cells = NULL, top.features = NULL, normalize = TRUE, scale.row = scale, show_featurenames = TRUE){
  #input checks
  if (is.null(counts) || !is.matrix(counts) & !is.data.frame(counts)){
    stop("'counts' should be a numeric count matrix")
  }
  if (is.null(celda.mod) || class(celda.mod) != "celda_G" & class(celda.mod) != "celda_CG"){
    stop("'celda.mod' should be an object of class celda_G or celda_CG")
  }
  compareCountMatrix(counts, celda.mod)
  
  #factorize counts matrix
  factorize.matrix <-
    factorizeMatrix(celda.mod = celda.mod, counts = counts)
  
  #take topRank
  if(!is.null(top.features) && (is.numeric(top.features))|is.integer(top.features)){
    top.rank <- topRank(matrix = factorize.matrix$proportions$module, 
                        n = top.features)
  }else{
    top.rank <- topRank(matrix = factorize.matrix$proportions$module, 
                        n = nrow(factorize.matrix$proportions$module))  
  }
  
  #filter topRank using feature.module into feature.indices
  feature.indices <- c()
  for(module in feature.module){
    feature.indices <- c(feature.indices, top.rank$index[[module]])
  }
  
  
  #Determine cell order from factorize.matrix$proportions$cell
  cell.states <- factorize.matrix$proportions$cell
  cell.states <- cell.states[feature.module, ,drop = FALSE]
  
  single.module <- cell.states[1, ]
  single.module.ordered <- order(single.module, decreasing = TRUE)
  
  if(!is.null(top.cells)){
    if(top.cells * 2 < ncol(cell.states)){
      cell.indices <- c(utils::head(single.module.ordered, n = top.cells), 
                        utils::tail(single.module.ordered, n = top.cells))
    }else{
      cell.indices <- single.module.ordered 
    }
  }else{
    cell.indices <- single.module.ordered
  }
  
  cell.indices <- rev(cell.indices)
  if(normalize){
    norm.counts <- normalizeCounts(counts, normalize="proportion", transformation.fun=sqrt)
  } else{
    norm.counts <- counts
  }
  
  #filter counts based on feature.indices 
  filtered_norm.counts <- norm.counts[feature.indices, cell.indices]
  
  filtered_norm.counts <- filtered_norm.counts[rowSums(filtered_norm.counts>0)>0,]
  
  gene_ix = match(rownames(filtered_norm.counts), celda.mod$names$row)
  cell_ix = match(colnames(filtered_norm.counts), celda.mod$names$column)
  if(!is.null(celda.mod$z)){
    cell <- distinct_colors(length(unique(celda.mod$z)))[sort(unique(celda.mod$z[cell_ix]))]
    names(cell) <- sort(unique(celda.mod$z[cell_ix]))
    anno_cell_colors <- list(cell = cell)
  }else{
    anno_cell_colors <- NULL
  }
  plotHeatmap(
    filtered_norm.counts,
    z = celda.mod$z[cell.indices],
    y = celda.mod$y[gene_ix],
    scale.row = scale.row,
    color.scheme = "divergent",
    show.names.feature = show_featurenames,
    cluster.feature = FALSE,
    cluster.cell = FALSE,
    annotation.color = anno_cell_colors
  )
}

