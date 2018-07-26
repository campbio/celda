#' @title Transcriptional state heatmap
#' @description Draws a heatmap focusing on a transcriptional state. Both cells and genes are sorted by 
#'    their proportions of counts in a given transcriptional state. Allows for nice visualization of 
#'    co-expression of those genes grouped into transcriptional states by Celda.
#'    
#' @param counts A numeric count matrix.
#' @param celda.mod An object of class celda_G or celda_CG.
#' @param state.use Numeric. A transcriptional state to plot.
#' @param cells.use Numeric or character. If a number, plot only this number of cells with respectively 
#'    highest and lowest proportions of counts in the transcriptional state. If a list of cell names, 
#'    plot only these cells orderd by their proportions of counts in the transcriptional state from 
#'    left to right. If NULL, plot all cells (NULL by default).
#' @param genes.use Numeric or character. If a number, plot only this number of genes with respectively 
#'    highest and lowest proportions of counts in the transcriptional state. If a list of gene names, 
#'    plot only these genes orderd by their proportions of counts in the transcriptional state from 
#'    the bottom up. If NULL, plot all genes in the transcriptional state (NULL by default).
#' @param normalize Logical. Should the columns be normmalized? Set to FALSE to disable. 'normalizeCounts' 
#'    normalizes to counts per million (CPM) (default by TRUE).
#' @param scale_row Function; A function to scale each individual row. Set to NULL to disable. Occurs 
#'    after normalization and log transformation. Defualt is 'scale' and thus will Z-score transform each row.
#' @param show_genenames Logical; specifying if gene names should be shown. Default to be TRUE.
#' @export 
stateHeatmap <- function(counts, celda.mod, state.use = 1, cells.use = NULL, genes.use = NULL, normalize = TRUE, scale_row = scale, show_genenames = TRUE){
  if (is.null(counts)) {
    stop("'counts' should be a numeric count matrix")
  }
  if (is.null(celda.mod) || is.null(celda.mod$y)){
    stop("'celda.mod' should be an object of class celda_G or celda_CG")
  }
  compareCountMatrix(counts, celda.mod)
  
  factorize.matrix <-
    factorizeMatrix(celda.mod = celda.mod, counts = counts)
  genes <- celda.mod$names$row[celda.mod$y %in% state.use]
  ascending_ordered_matrix <-
    factorize.matrix$proportions$gene.states[, state.use, drop = FALSE][which(rownames(factorize.matrix$proportions$gene.states[, state.use, drop = FALSE]) %in% genes),]
  if(class(ascending_ordered_matrix) == "numeric"){
    ascending_ordered_genes <- names(sort(ascending_ordered_matrix))
  }else{
    ascending_ordered_genes <- names(sort(rowMeans(ascending_ordered_matrix)))
  }
  if (class(genes.use) == "character") {
    if (setequal(intersect(genes.use, ascending_ordered_genes), genes.use)) {
      filtered_genes <- names(sort(factorize.matrix$proportions$gene.states[, state.use][which(rownames(factorize.matrix$proportions$gene.states[, state.use, drop = FALSE]) %in% genes.use)]))
    } else {
      miss_genes <- setdiff(genes.use, intersect(genes.use, ascending_ordered_genes))
      miss_genes[-length(miss_genes)] <- paste0(miss_genes[-length(miss_genes)], ', ')
      stop("'genes.use': gene(s) (", miss_genes, ") is/are not in L", state.use)
    }
  }else{
    if(is.null(genes.use) || genes.use >= ceiling(length(ascending_ordered_genes) / 2)){
      filtered_genes <- ascending_ordered_genes
      message(paste0("Use all genes in L", state.use, " "))
    } else{
      filtered_genes <-
        c(
          head(ascending_ordered_genes, genes.use),
          tail(ascending_ordered_genes, genes.use)
        )
    }
  }
  ascending_ordered_cells <-
    names(sort(colMeans(factorize.matrix$proportions$cell.states[state.use,,drop = FALSE])))
  if (is.null(cells.use)) {
    filtered_cells <-
      ascending_ordered_cells
    message("Use all cells")
  } else if (is.numeric(cells.use)){
    if (cells.use >= ceiling(length(ascending_ordered_cells) / 2)) {
      filtered_cells <- ascending_ordered_cells
      message("Use all cells")
    } else{
      filtered_cells <-
        c(
          head(ascending_ordered_cells, cells.use),
          tail(ascending_ordered_cells, cells.use)
        )
    }
  }else{
    filtered_cells <-
      names(sort(factorize.matrix$proportions$cell.states[state.use, cells.use]))
  }
  if(normalize){
    norm.counts <- normalizeCounts(counts)
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
    normalize = NULL,
    scale_row = scale_row,
    color_scheme = "divergent",
    show_genenames = show_genenames,
    cluster_gene = FALSE,
    cluster_cell = FALSE,
    annotation_color = anno_cell_colors
  )
}
