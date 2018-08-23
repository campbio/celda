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
  if (is.null(counts)) {
    stop("'counts' should be a numeric count matrix")
  }
  if (is.null(celda.mod) || is.null(celda.mod$y)){
    stop("'celda.mod' should be an object of class celda_G or celda_CG")
  }
  compareCountMatrix(counts, celda.mod)
  
  factorize.matrix <-
    factorizeMatrix(celda.mod = celda.mod, counts = counts)
  genes <- celda.mod$names$row[celda.mod$y %in% feature.module]
  ascending_ordered_matrix <-
    factorize.matrix$proportions$gene.states[, feature.module, drop = FALSE][which(rownames(factorize.matrix$proportions$gene.states[, feature.module, drop = FALSE]) %in% genes),]
  if(class(ascending_ordered_matrix) == "numeric"){
    ascending_ordered_genes <- names(sort(ascending_ordered_matrix))
  }else{
    ascending_ordered_genes <- names(sort(rowMeans(ascending_ordered_matrix)))
  }
  if (class(top.features) == "character") {
    if (setequal(intersect(top.features, ascending_ordered_genes), top.features)) {
      filtered_genes <- names(sort(factorize.matrix$proportions$gene.states[, feature.module][which(rownames(factorize.matrix$proportions$gene.states[, feature.module, drop = FALSE]) %in% top.features)]))
    } else {
      miss_genes <- setdiff(top.features, intersect(top.features, ascending_ordered_genes))
      miss_genes[-length(miss_genes)] <- paste0(miss_genes[-length(miss_genes)], ', ')
      stop("'top.features': gene(s) (", miss_genes, ") is/are not in L", feature.module)
    }
  }else{
    if(is.null(top.features) || top.features >= ceiling(length(ascending_ordered_genes) / 2)){
      filtered_genes <- ascending_ordered_genes
      message(paste0("Use all genes in L", feature.module, " "))
    } else{
      filtered_genes <-
        c(
          head(ascending_ordered_genes, top.features),
          tail(ascending_ordered_genes, top.features)
        )
    }
  }
  ascending_ordered_cells <-
    names(sort(colMeans(factorize.matrix$proportions$cell.states[feature.module,,drop = FALSE])))
  if (is.null(top.cells)) {
    filtered_cells <-
      ascending_ordered_cells
    message("Use all cells")
  } else if (is.numeric(top.cells)){
    if (top.cells >= ceiling(length(ascending_ordered_cells) / 2)) {
      filtered_cells <- ascending_ordered_cells
      message("Use all cells")
    } else{
      filtered_cells <-
        c(
          head(ascending_ordered_cells, top.cells),
          tail(ascending_ordered_cells, top.cells)
        )
    }
  }else{
    filtered_cells <-
      names(sort(factorize.matrix$proportions$cell.states[feature.module, top.cells]))
  }
  if(normalize){
    norm.counts <- normalizeCounts(counts, normalize="proportion", transformation.fun=sqrt)
  } else{
    norm.counts <- counts
  }
  filtered_norm.counts <- norm.counts[rev(filtered_genes), filtered_cells]
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
  renderCeldaHeatmap(
    filtered_norm.counts,
    z = celda.mod$z[cell_ix],
    y = celda.mod$y[gene_ix],
    scale.row = scale.row,
    color_scheme = "divergent",
    show_featurenames = show_featurenames,
    cluster_gene = FALSE,
    cluster_cell = FALSE,
    annotation_color = anno_cell_colors
  )
}
