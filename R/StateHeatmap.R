#' @title Heatmap for feature modules
#' @description Renders a heatmap for selected feature modules. Cells are ordered from those with the lowest probability of the module on the left to the highest probability on the right. If more than one module is used, then cells will be ordered by the probabilities of the first module only. Features are ordered from those with the highest probability in the module on the top to the lowest probability on the bottom.
#'
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class `celda_G` or `celda_CG`. 
#' @param feature.module Integer Vector. The feature module(s) to display. Multiple modules can be included in a vector. 
#' @param top.cells Integer. Number of cells with the highest and lowest probabilities for this module to include in the heatmap. For example, if `top.cells` = 50, the 50 cells with the lowest probability and the 50 cells with the highest probability for that feature module will be included. If NULL, all cells will be plotted. Default 100. 
#' @param top.features Integer. Plot `top.features` with the highest probability in the feature module. If NULL, plot all features in the module. Default NULL.
#' @param normalize Logical. Whether to normalize the columns of `counts`. Default TRUE. 
#' @param scale.row Character. Which function to use to scale each individual row. Set to NULL to disable. Occurs after normalization and log transformation. For example, `scale` will Z-score transform each row. Default `scale`.
#' @param show_featurenames Logical. Wheter feature names should be displayed. Default TRUE. 
#' @return A list containing row and column dendrograms as well as a gtable for grob plotting
#' @examples
#' moduleHeatmap(celda.CG.sim$counts, celda.CG.mod)
#' @export 
moduleHeatmap <- function(counts, celda.mod, feature.module = 1, top.cells = 100, top.features = NULL, normalize = TRUE, scale.row = scale, show_featurenames = TRUE){
  #input checks
  if (is.null(counts) || !is.matrix(counts) & !is.data.frame(counts)){
    stop("'counts' should be a numeric count matrix")
  }
  if (is.null(celda.mod) || !methods::is(celda.mod, "celda_G") & !methods::is(celda.mod, "celda_CG")){
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
  feature.indices <- lapply(feature.module, function(module) {
    top.rank$index[[module]]
  })
  feature.indices <- unlist(feature.indices)

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
  filtered_norm.counts <- norm.counts[feature.indices, cell.indices, drop = FALSE]
  
  filtered_norm.counts <- filtered_norm.counts[rowSums(filtered_norm.counts>0)>0,,drop = FALSE]
  
  gene_ix = match(rownames(filtered_norm.counts), celda.mod@names$row)
  cell_ix = match(colnames(filtered_norm.counts), celda.mod@names$column)
  z.to.plot = c()
  if(methods::.hasSlot(celda.mod, "z")){
    cell <- distinct_colors(length(unique(celda.mod@clusters$z)))[sort(unique(celda.mod@clusters$z[cell_ix]))]
    names(cell) <- sort(unique(celda.mod@clusters$z[cell_ix]))
    anno_cell_colors <- list(cell = cell)
    z.to.plot = celda.mod@clusters$z[cell.indices]
  }else{
    anno_cell_colors <- NULL
  }
  plotHeatmap(
    filtered_norm.counts,
    z = z.to.plot,
    y = celda.mod@clusters$y[gene_ix],
    scale.row = scale.row,
    color.scheme = "divergent",
    show.names.feature = show_featurenames,
    cluster.feature = FALSE,
    cluster.cell = FALSE,
    annotation.color = anno_cell_colors
  )
}

