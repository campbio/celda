#' Render a heatmap based on a matrix of counts where rows are genes and columns are cells.
#' 
#' @param counts A count matrix where rows are genes and columns are cells. 
#' @param z A numeric vector of cluster assignments for cell. 
#' @param y A numeric vector of cluster assignments for gene.
#' @param scale_log Function; Applys a scale function such as log, log2, log10. Set to NULL to disable. Occurs after normalization. Default NULL.
#' @param pseudocount_log Numeric; A pseudocount to add to data before log transforming. Default  0. 
#' @param pseudocount_normalize Numeric; A pseudocount to add to data before normalization. Default  1. 
#' @param scale_row Function; A function to scale each individual row. Set to NULL to disable. Occurs after normalization and log transformation. Defualt is 'scale' and thus will Z-score transform each row. 
#' @param trim A two element vector to specify the lower and upper cutoff for the data. Occurs after normalization, log transformation, and row scaling. Set to NULL to disable. Default c(-2,2).
#' @param normalize A function to normalize the columns. Set to NULL to disable. Default is 'normalizeCounts', which normalizes to counts per million (CPM). 
#' @param cluster_row Logical; determining if rows should be clustered.
#' @param cluster_column Logical; determining if columns should be clustered.
#' @param annotation_cell a dataframe for the cell annotations (columns).
#' @param annotation_gene a dataframe for the gene annotations (rows).
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
#' @import gtable
#' @import grid
#' @import scales
#' @import RColorBrewer
#' @import grDevices
#' @import graphics
#' @export 

render_celda_heatmap <- function(counts, z = NULL, y = NULL, 
                                 scale_log = NULL,
                                 scale_row = scale,
                                 normalize = normalizeCounts,
                                 trim=c(-2,2), 
                                 pseudocount_normalize=0,
                                 pseudocount_log=0,
                                 cluster_row = TRUE, cluster_column = TRUE,
                                 annotation_cell = NULL, annotation_gene = NULL, 
                                 col= NULL,
                                 breaks = NULL, 
                                 legend = TRUE,
                                 annotation_legend = TRUE,
                                 annotation_names_gene = TRUE, 
                                 annotation_names_cell = TRUE,
                                 show_genenames = FALSE, 
                                 show_cellnames = FALSE,
                                 hclust_method = "ward.D2",
                                 ...) {
  
  ## Check length of z/y variables
  if (!is.null(z) & z != ncol(counts)) stop("Length of z must match number of columns in counts matrix")
  if (!is.null(y) & y != nrow(counts)) stop("Length of y must match number of rows in counts matrix")
  
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

  #if null y or z 
  if(is.null(y)){
    y <- rep(1, nrow(counts))
  }
  if(is.null(z)){
    z <- rep(1, ncol(counts))
  }
  
  
  K <- length(unique(z))
  L <- length(unique(y))
  
  
  ## Set row and column name to counts matrix 
  if(is.null(rownames(counts))){
    rownames(counts) <- 1:nrow(counts)
  } 
  if(is.null(colnames(counts))) {
    colnames(counts) <- 1:ncol(counts)
  }  
  
  #  Set cell annotation  
  if(!is.null(annotation_cell)){
    annotation_cell$cell_label <- as.factor(z)   # has to be the original cell label order
    rownames(annotation_cell) <- colnames(counts)
  }else{  # annotation_cell is null
	annotation_cell <- data.frame(cell_label = as.factor(z)) 
    rownames(annotation_cell) <- colnames(counts)  # rowname should correspond to counts matrix's col(cell) name
  }
  #  Set gene annotation
  if(!is.null(annotation_gene)){
    annotation_gene$gene_label <- as.factor(y)  # has to be the original gene label order 
    rownames(annotation_gene) <- rownames(counts)
  } else{  #annotation_gene is  null 
	annotation_gene <- data.frame(gene_label = as.factor(y))
    rownames(annotation_gene) <- rownames(counts)  # rowname should correspond to counts matrix's row(gene) name
  }
    
  ## set breaks & color
  ubound.range <- max(counts)
  lbound.range <- min(counts)

  total.range <- max(abs(c(ubound.range, lbound.range)))

  if(lbound.range < 0 & 0 < ubound.range){  # both sides of zero for the counts values
    if(is.null(col)){
      col <- colorRampPalette(c("#1E90FF","#FFFFFF","#CD2626"),space = "Lab")(100)
    }
    col.len <- length(col)
    if(is.null(breaks)){
      breaks <- c(seq(-total.range, 0,  length.out = round(col.len/2) + 1  ),
                  seq(0+1e-6, total.range, length.out = col.len-round(col.len/2) ))
    }
  } else {  # only one side for the counts values (either positive or negative, for probabilities)
    if(is.null(col)){
      col <- colorRampPalette(c("#FFFFFF", brewer.pal(n = 9, name = "Reds")))(100)
    }
  }

  semi_pheatmap(mat = counts,
  	color = col,
    breaks = breaks, 
    cluster_cols = cluster_column,
    cluster_rows = cluster_row,
    annotation_row = annotation_gene,
    annotation_col = annotation_cell,
    legend = legend,
    annotation_legend = annotation_legend, 
    annotation_names_row = annotation_names_gene, 
    annotation_names_col = annotation_names_cell,
    show_rownames = show_genenames,
    show_colnames = show_cellnames,
    clustering_method =  hclust_method,
    row_label = y,
    col_label = z,
    ...)
}
