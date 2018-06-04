#' Render a heatmap based on a matrix of counts where rows are genes and columns are cells.
#' 
#' @param counts A count matrix where rows are genes and columns are cells. 
#' @param z A numeric vector of cluster assignments for cell. 
#' @param y A numeric vector of cluster assignments for gene.
#' @param scale_log Function; Applys a scale function such as log, log2, log10. Set to NULL to disable. Occurs after normalization. Default NULL.
#' @param pseudocount_log Numeric; A pseudocount to add to data before log transforming. Default  0. 
#' @param pseudocount_normalize Numeric; A pseudocount to add to data before normalization. Default  1. 
#' @param gene.ix Index of genes to pull out of the matrix after normalization. If NULL, no subsettig will be performed. Default NULL.
#' @param cell.ix Index of genes to pull out of the matrix after normalization. If NULL, no subsettig will be performed. Default NULL.
#' @param scale_row Function; A function to scale each individual row. Set to NULL to disable. Occurs after normalization and log transformation. Defualt is 'scale' and thus will Z-score transform each row. 
#' @param trim A two element vector to specify the lower and upper cutoff for the data. Occurs after normalization, log transformation, and row scaling. Set to NULL to disable. Default c(-2,2).
#' @param normalize A function to normalize the columns. Set to NULL to disable. Default is 'normalizeCounts', which normalizes to counts per million (CPM). 
#' @param cluster_gene Logical; determining if rows should be clustered.
#' @param cluster_cell Logical; determining if columns should be clustered.
#' @param annotation_cell A data frame for the cell annotations (columns).
#' @param annotation_gene A data frame for the gene annotations (rows).
#' @param annotation_color A list containing color scheme for cell and/or gene annotation. See '?pheatmap' for more details.
#' @param color_scheme One of "divergent" or "sequential". A "divergent" scheme is best for highlighting either end of the data with a break in the middle (denoted by 'color_scheme_center') such as gene expression data that has been normalized and centered. 
#' A "sequential" scheme is best for data that are ordered low to high such as raw counts or probabilities.
#' @param color_scheme_symmetric Boolean; When the color_scheme is "divergent" and the data contains both positive and negative numbers, TRUE indicates that the color scheme should be symmetric about the center. For example, if the data ranges goes from -1.5 to 2, then setting this to TRUE will force the colors to range from -2 to 2.
#' @param color_scheme_center Numeric. Indicates the center of a "divergent" color_scheme. Default 0.
#' @param col color for the heatmap. 
#' @param breaks a sequence of numbers that covers the range of values in mat and is one 
#' element longer than color vector. Used for mapping values to colors. Useful, if needed 
#' to map certain values to certain colors, to certain values. If value is NA then the 
#' breaks are calculated automatically.
#' @param legend logical; determining if legend should be drawn or not. Default to be TRUE.
#' @param annotation_legend Logical; showing if the legend for annotation tracks should be drawn.
#' @param annotation_names_gene Logical; showing if the names for gene annotation tracks should be drawn. Default to be TRUE.
#' @param annotation_names_cell Logical; showing if the names for cell annotation tracks should be drawn. Default to be TRUE. 
#' @param show_genenames Logical; specifying if gene names should be shown. Default to be FALSE. 
#' @param show_cellnames Logical; specifying if cell names should be shown. Default to be FALSE. 
#' @param hclust_method Character; Specifies the method to use for the 'hclust' function. Default is "ward.D2". See ?hclust for possible values. 
#' @param treeheight_gene Numeric; Width of the gene dedrogram. Set to 0 to disable plotting of this dendrogram. 
#' @param treeheight_cell Numeric; Hieght of the cell dedrogram. Set to 0 to disable plotting of this dendrogram. 
#' @param ... Other arguments to be passed to underlying pheatmap function
#' @import gtable
#' @import grid
#' @import scales
#' @import RColorBrewer
#' @import grDevices
#' @import graphics
#' @export 
renderCeldaHeatmap <- function(counts, z = NULL, y = NULL, 
                                 scale_log = NULL,
                                 scale_row = scale,
                                 normalize = normalizeCounts,
                                 trim=c(-2,2), 
                                 pseudocount_normalize=0,
                                 pseudocount_log=0,
                                 gene.ix = NULL,
                                 cell.ix = NULL,
                                 cluster_gene = TRUE, cluster_cell = TRUE,
                                 color_scheme = c("divergent", "sequential"),
                                 color_scheme_symmetric = TRUE,
                                 color_scheme_center = 0,
                                 col= NULL,
                                 annotation_cell = NULL, annotation_gene = NULL, 
                                 annotation_color = NULL,
                                 breaks = NULL, 
                                 legend = TRUE,
                                 annotation_legend = TRUE,
                                 annotation_names_gene = TRUE, 
                                 annotation_names_cell = TRUE,
                                 show_genenames = FALSE, 
                                 show_cellnames = FALSE,
                                 hclust_method = "ward.D2",
                                 treeheight_gene = ifelse(cluster_gene, 50, 0), 
								 treeheight_cell = ifelse(cluster_cell, 50, 0),
                                 ...) {
  
  # Check for same lengths for z and y group variables
  if (!is.null(z) & length(z) != ncol(counts)) stop("Length of z must match number of columns in counts matrix")
  if (!is.null(y) & length(y) != nrow(counts)) stop("Length of y must match number of rows in counts matrix")
  color_scheme = match.arg(color_scheme)

  ## Normalize, transform, row scale, and then trim data
  if(!is.null(normalize)) {  
    counts <- do.call(normalize, list(counts + pseudocount_normalize))
  }
  if(!is.null(scale_log)){
    counts <- do.call(scale_log, list(counts + pseudocount_log))
  }
  if(!is.null(scale_row)) {
    counts <- t(base::apply(counts, 1, scale_row))
  }  
  if(!is.null(trim)){
    if(length(trim) != 2) {
      stop("'trim' should be a 2 element vector specifying the lower and upper boundaries")
    }
    trim<-sort(trim)
    counts[counts < trim[1]] <- trim[1]
    counts[counts > trim[2]] <- trim[2]
  }
  
  ## Create cell annotation  
  if(!is.null(annotation_cell) & !is.null(z)){
  
    if(is.null(rownames(annotation_cell))) {
      rownames(annotation_cell) = colnames(counts)
    } else {
      if(any(rownames(annotation_cell) != colnames(counts))) {
        stop("Row names of 'annotation_cell' are different than the column names of 'counts'")
      }
    }
    annotation_cell = data.frame(cell = as.factor(z), annotation_cell)
    
  } else if(is.null(annotation_cell) & !is.null(z)) { 
	annotation_cell <- data.frame(cell = as.factor(z)) 
    rownames(annotation_cell) <- colnames(counts)  
  } else {
    annotation_cell = NA
  }
  
  
  # Set gene annotation
  if(!is.null(annotation_gene) & !is.null(y)){
  
    if(is.null(rownames(annotation_gene))) {
      rownames(annotation_cell) = rownames(counts)
    } else {
      if(any(rownames(annotation_cell) != rownames(counts))) {
        stop("Row names of 'annotation_gene' are different than the row names of 'counts'")
      }
    }
    annotation_gene = data.frame(gene = as.factor(y), annotation_gene)   
    
  } else if(is.null(annotation_gene) & !is.null(y)) { 
	annotation_gene <- data.frame(gene = as.factor(y)) 
    rownames(annotation_gene) <- rownames(counts)  
  } else {
    annotation_gene = NA
  }
  
  ## Set annotation colors
  if(!is.null(z)) {
    K = sort(unique(z))
    K.col = distinct_colors(length(K))
    names(K.col) = K

    if(!is.null(annotation_color)) {
      if(!("cell" %in% names(annotation_color))) {
        annotation_color = c(list(cell = K.col), annotation_color)
      }
    } else {
      annotation_color = list(cell = K.col)
    }
  }
     
  if(!is.null(y)) {
    L = sort(unique(y))
    L.col = distinct_colors(length(L))
    names(L.col) = L

    if(!is.null(annotation_color)) {
      if(!("gene" %in% names(annotation_color))) {
        annotation_color = c(list(gene = L.col), annotation_color)
      }
    } else {
      annotation_color = list(gene = L.col)
    }
  }


  ## Select subsets of genes/cells
  if(!is.null(gene.ix)) {
    counts = counts[gene.ix,,drop=FALSE]
    if(length(annotation_gene) > 1 || (length(annotation_gene) == 1 & !is.na(annotation_gene))) {
      annotation_gene = annotation_gene[gene.ix,,drop=FALSE]
    }
    if(!is.null(y)) {
      y = y[gene.ix]
    }
  }
  if(!is.null(cell.ix)) {
    counts = counts[,cell.ix,drop=FALSE]
    if(length(annotation_cell) > 1 || (length(annotation_cell) == 1 & !is.na(annotation_cell))) {
      annotation_cell = annotation_cell[cell.ix,,drop=FALSE]
    }
    if(!is.null(z)) {  
      z = z[cell.ix]
    }  
  }

  ## Set color scheme and breaks
  ubound.range <- max(counts)
  lbound.range <- min(counts)

  if(color_scheme == "divergent"){  
    if(color_scheme_symmetric == TRUE){
      ubound.range <- max(abs(ubound.range), abs(lbound.range))
      lbound.range <- -ubound.range
    }
    if(is.null(col)){
      col <- colorRampPalette(c("#1E90FF","#FFFFFF","#CD2626"),space = "Lab")(100)
    }
    col.len <- length(col)
    if(is.null(breaks)){
      breaks <- c(seq(lbound.range, color_scheme_center, length.out = round(col.len/2) + 1  ),
                  seq(color_scheme_center+1e-6, ubound.range, length.out = col.len-round(col.len/2) ))
    }
  } else {  # Sequential color scheme
    if(is.null(col)){
      col <- colorRampPalette(c("#FFFFFF", brewer.pal(n = 9, name = "Blues")))(100)
      col.len = length(col)
    }
    if(is.null(breaks)){
        breaks <- seq(lbound.range, ubound.range, length.out = col.len)
    }
  }

  semi_pheatmap(mat = counts,
  	color = col,
    breaks = breaks, 
    cluster_cols = cluster_cell,
    cluster_rows = cluster_gene,
    annotation_row = annotation_gene,
    annotation_col = annotation_cell,
    annotation_colors = annotation_color,
    legend = legend,
    annotation_legend = annotation_legend, 
    annotation_names_row = annotation_names_gene, 
    annotation_names_col = annotation_names_cell,
    show_rownames = show_genenames,
    show_colnames = show_cellnames,
    clustering_method =  hclust_method,
    treeheight_row = treeheight_cell,
    treeheight_col = treeheight_cell,
    row_label = y,
    col_label = z,
    ...)
}



#' Renders a heatmap based on a population matrix from the factorized counts matrix.
#' 
#' @param counts A count matrix where rows are genes and columns are cells, used to generate celda model. 
#' @param celda.mod celda model of class "celda_CG". 
#' @param main The title of the plot; default = NA. 
#' @import gtable
#' @import grid
#' @import scales
#' @import RColorBrewer
#' @import grDevices
#' @import graphics
#' @export 
absoluteProbabilityHeatmap <- function(counts, celda.mod, main = NA){
  factorized <- factorizeMatrix(celda.mod = celda.mod, counts = counts)
  pop <- factorized$proportions$population.states
  z <- 1:ncol(pop)
  y <- 1:nrow(pop)
  
  K = sort(unique(z))
  K.col = distinct_colors(length(K))
  names(K.col) = K
  
  annotation_color = list(cell = K.col)
  
  L = sort(unique(y))
  L.col = distinct_colors(length(L))
  names(L.col) = L
  
  percentile.9 <- round(quantile(pop,.9), digits = 2) * 100
  col1 <- colorRampPalette(c("#FFFFFF", brewer.pal(n = 9, name = "Blues")))(percentile.9)
  col2 <- colorRampPalette(c("#08306B", c("#006D2C","Yellowgreen","Yellow","Orange","Red")))(100-percentile.9)
  col <- c(col1,col2)
  
  breaks <-  seq(0, 1, length.out = length(col)) 
  
  semi_pheatmap(pop, row_label = NULL, col_label = NULL, col = col, breaks = breaks, cluster_cols = FALSE, cluster_rows = FALSE, main = main)
}



#' Renders a heatmap based on a population matrix from the factorized counts matrix. The relative probability of each transcriptional state in each cell subpopulation is visualized.
#' 
#' @param counts A count matrix where rows are genes and columns are cells, used to generate celda model. 
#' @param celda.mod celda model of class "celda_CG". 
#' @param main The title of the plot; default = NA. 
#' @import gtable
#' @import grid
#' @import scales
#' @import RColorBrewer
#' @import grDevices
#' @import graphics
#' @export 
relativeProbabilityHeatmap <- function(counts, celda.mod, main = NA){
  factorized <- factorizeMatrix(celda.mod = celda.mod, counts = counts)
  pop <- factorized$proportions$population.states
  z <- 1:ncol(pop)
  y <- 1:nrow(pop)
  
  K = sort(unique(z))
  K.col = distinct_colors(length(K))
  names(K.col) = K
  
  annotation_color = list(cell = K.col)
  
  L = sort(unique(y))
  L.col = distinct_colors(length(L))
  names(L.col) = L
  
  pop <- sweep(pop, 1, rowSums(pop), "/")
  
  percentile.9 <- round(quantile(pop,.9), digits = 2) * 100
  col1 <- colorRampPalette(c("#FFFFFF", brewer.pal(n = 9, name = "Blues")))(percentile.9)
  col2 <- colorRampPalette(c("#08306B", c("#006D2C","Yellowgreen","Yellow","Orange","Red")))(100-percentile.9)
  col <- c(col1,col2)
  
  breaks <-  seq(0, 1, length.out = length(col)) 
  
  semi_pheatmap(pop, row_label = NULL, col_label = NULL, col = col, breaks = breaks, cluster_cols = FALSE, cluster_rows = FALSE, main = main)
}
