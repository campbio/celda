#' Renders a heatmap based on a matrix of counts where rows are features and columns are cells.
#' 
#' @param counts Numeric matrix. Normalized counts matrix where rows represent features and columns represent cells. . 
#' @param z Numeric vector. Denotes cell population labels.  
#' @param y Numeric vector. Denotes feature module labels. 
#' @param feature.ix Integer vector. Select features for display in heatmap. If NULL, no subsetting will be performed. Default NULL.
#' @param cell.ix Integer vector. Select cells for display in heatmap. If NULL, no subsetting will be performed. Default NULL.
#' @param scale.row Function; A function to scale each individual row. Set to NULL to disable. Occurs after normalization and log transformation. Defualt is 'scale' and thus will Z-score transform each row. 
#' @param trim Numeric vector. Vector of length two that specifies the lower and upper bounds for the data. This threshold is applied after row scaling. Set to NULL to disable. Default c(-2,2). 
#' @param cluster.feature Logical. Determines whether rows should be clustered. Default TRUE. 
#' @param cluster.cell Logical. Determines whether columns should be clustered. Default TRUE. 
#' @param annotation.cell Data frame. Additional annotations for each cell will be shown in the column color bars. The format of the data frame should be one row for each cell and one column for each annotation. Numeric variables will be displayed as continuous color bars and factors will be displayed as discrete color bars. Default NULL. 
#' @param annotation.feature A data frame for the feature annotations (rows).
#' @param annotation.color List. Contains color scheme for all annotations. See `?pheatmap` for more details. 
#' @param color.scheme "Character. One of ""divergent"" or ""sequential"". A ""divergent"" scheme is best for highlighting relative data (denoted by 'color.scheme.center') such as gene expression data that has been normalized and centered. A ""sequential"" scheme is best for highlighting data that are ordered low to high such as raw counts or probabilities. Default "divergent".
#' @param color.scheme.symmetric Logical. When the color.scheme is "divergent" and the data contains both positive and negative numbers, TRUE indicates that the color scheme should be symmetric from [-max(abs(data)),max(abs(data))]. For example, if the data ranges goes from -1.5 to 2, then setting this to TRUE will force the color scheme to range from -2 to 2. Default TRUE.
#' @param color.scheme.center Numeric. Indicates the center of a "divergent" color.scheme. Default 0.
#' @param col Color for the heatmap. 
#' @param breaks Numeric vector. A sequence of numbers that covers the range of values in the normalized `counts`. Values in the normalized `matrix` are assigned to each bin in `breaks`. Each break is assigned to a unique color from `col`. If NULL, then breaks are calculated automatically. Default NULL. 
#' @param legend Logical. Determines whether legend should be drawn. Default TRUE.
#' @param annotation.legend Logical. Whether legend for all annotations should be drawn. Default TRUE. 
#' @param annotation.names.feature Logical. Whether the names for features should be shown. Default TRUE.
#' @param annotation.names.cell Logical. Whether the names for cells should be shown. Default TRUE. 
#' @param show.names.feature Logical. Specifies if feature names should be shown. Default TRUE.  
#' @param show.names.cell Logical. Specifies if cell names should be shown. Default FALSE. 
#' @param hclust.method Character. Specifies the method to use for the 'hclust' function. See `?hclust` for possible values. Default "ward.D2".  
#' @param treeheight.feature Numeric. Width of the feature dendrogram. Set to 0 to disable plotting of this dendrogram. Default: if cluster.feature == TRUE, then treeheight.feature = 50, else treeheight.feature = 0.  
#' @param treeheight.cell Numeric. Height of the cell dendrogram. Set to 0 to disable plotting of this dendrogram. Default: if cluster.cell == TRUE, then treeheight.cell = 50, else treeheight.cell = 0.  
#' @param silent Logical. Whether to plot the heatmap.
#' @param ... Other arguments to be passed to underlying pheatmap function.
#' @examples 
#' plotHeatmap(celda.CG.sim$counts, z=clusters(celda.CG.mod)$z, y=clusters(celda.CG.mod)$y)
#' @return list A list containing dendrogram information and the heatmap grob
#' @import gtable
#' @import grid
#' @import scales
#' @import RColorBrewer
#' @import grDevices
#' @import graphics
#' @export 
plotHeatmap <- function(counts, 
                        z = NULL, 
                        y = NULL, 
                        scale.row = scale,
                        trim=c(-2,2), 
                        feature.ix = NULL,
                        cell.ix = NULL,
                        cluster.feature = TRUE, 
                        cluster.cell = TRUE,
                        color.scheme = c("divergent", "sequential"),
                        color.scheme.symmetric = TRUE,
                        color.scheme.center = 0,
                        col= NULL,
                        annotation.cell = NULL, 
                        annotation.feature = NULL, 
                        annotation.color = NULL,
                        breaks = NULL, 
                        legend = TRUE,
                        annotation.legend = TRUE,
                        annotation.names.feature = TRUE, 
                        annotation.names.cell = TRUE,
                        show.names.feature = FALSE, 
                        show.names.cell = FALSE,
                        hclust.method = "ward.D2",
                        treeheight.feature = ifelse(cluster.feature, 50, 0), 
      								  treeheight.cell = ifelse(cluster.cell, 50, 0),
      								  silent = FALSE,
                        ...) {
  
  
  # Check for same lengths for z and y group variables
  if (!is.null(z) & length(z) != ncol(counts)) stop("Length of z must match number of columns in counts matrix")
  if (!is.null(y) & length(y) != nrow(counts)) stop("Length of y must match number of rows in counts matrix")
  color.scheme = match.arg(color.scheme)

  if(!is.null(scale.row)) {
    if(is.function(scale.row)) {
      cn = colnames(counts)
      counts <- t(base::apply(counts, 1, scale.row))
      colnames(counts) = cn
    } else {
      stop("'scale.row' needs to be of class 'function'")
    }  
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
  if(!is.null(annotation.cell) & !is.null(z)){
  
    if(is.null(rownames(annotation.cell))) {
      rownames(annotation.cell) = colnames(counts)
    } else {
      if(any(rownames(annotation.cell) != colnames(counts))) {
        stop("Row names of 'annotation.cell' are different than the column names of 'counts'")
      }
    }
    annotation.cell = data.frame(cell = as.factor(z), annotation.cell)
    
  } else if(is.null(annotation.cell) & !is.null(z)) { 
	annotation.cell <- data.frame(cell = as.factor(z)) 
    rownames(annotation.cell) <- colnames(counts)  
  } else {
    annotation.cell = NA
  }
  
  
  # Set feature annotation
  if(!is.null(annotation.feature) & !is.null(y)){
  
    if(is.null(rownames(annotation.feature))) {
      rownames(annotation.feature) = rownames(counts)
    } else {
      if(any(rownames(annotation.feature) != rownames(counts))) {
        stop("Row names of 'annotation.feature' are different than the row names of 'counts'")
      }
    }
    annotation.feature = data.frame(module = as.factor(y), annotation.feature)   
    
  } else if(is.null(annotation.feature) & !is.null(y)) { 
	annotation.feature <- data.frame(module = as.factor(y)) 
    rownames(annotation.feature) <- rownames(counts)  
  } else {
    annotation.feature = NA
  }
  
  ## Set annotation colors
  if(!is.null(z)) {
    K = sort(unique(z))
    K.col = distinct_colors(length(K))
    names(K.col) = K

    if(!is.null(annotation.color)) {
      if(!("cell" %in% names(annotation.color))) {
        annotation.color = c(list(cell = K.col), annotation.color)
      }
    } else {
      annotation.color = list(cell = K.col)
    }
  }
     
  if(!is.null(y)) {
    L = sort(unique(y))
    L.col = distinct_colors(length(L))
    names(L.col) = L

    if(!is.null(annotation.color)) {
      if(!("module" %in% names(annotation.color))) {
        annotation.color = c(list(module = L.col), annotation.color)
      }
    } else {
      annotation.color = list(module = L.col)
    }
  }


  ## Select subsets of features/cells
  if(!is.null(feature.ix)) {
    counts = counts[feature.ix,,drop=FALSE]
    if(length(annotation.feature) > 1 || (length(annotation.feature) == 1 & !is.na(annotation.feature))) {
      annotation.feature = annotation.feature[feature.ix,,drop=FALSE]
    }
    if(!is.null(y)) {
      y = y[feature.ix]
    }
  }
  if(!is.null(cell.ix)) {
    counts = counts[,cell.ix,drop=FALSE]
    if(length(annotation.cell) > 1 || (length(annotation.cell) == 1 & !is.na(annotation.cell))) {
      annotation.cell = annotation.cell[cell.ix,,drop=FALSE]
    }
    if(!is.null(z)) {  
      z = z[cell.ix]
    }  
  }

  ## Set color scheme and breaks
  ubound.range <- max(counts)
  lbound.range <- min(counts)

  if(color.scheme == "divergent"){  
    if(color.scheme.symmetric == TRUE){
      ubound.range <- max(abs(ubound.range), abs(lbound.range))
      lbound.range <- -ubound.range
    }
    if(is.null(col)){
      col <- colorRampPalette(c("#1E90FF","#FFFFFF","#CD2626"),space = "Lab")(100)
    }
    col.len <- length(col)
    if(is.null(breaks)){
      breaks <- c(seq(lbound.range, color.scheme.center, length.out = round(col.len/2) + 1  ),
                  seq(color.scheme.center+1e-6, ubound.range, length.out = col.len-round(col.len/2) ))
    }
  } else {  # Sequential color scheme
    if(is.null(col)){
      col <- colorRampPalette(c("#FFFFFF", brewer.pal(n = 9, name = "Blues")))(100)
    }
    col.len = length(col)
    if(is.null(breaks)){
        breaks <- seq(lbound.range, ubound.range, length.out = col.len)
    }
  }

  sp = semi_pheatmap(mat = counts,
  	color = col,
    breaks = breaks, 
    cluster_cols = cluster.cell,
    cluster_rows = cluster.feature,
    annotation_row = annotation.feature,
    annotation_col = annotation.cell,
    annotation_colors = annotation.color,
    legend = legend,
    annotation_legend = annotation.legend, 
    annotation_names_row = annotation.names.feature, 
    annotation_names_col = annotation.names.cell,
    show_rownames = show.names.feature,
    show_colnames = show.names.cell,
    clustering_method =  hclust.method,
    treeheight_row = treeheight.feature,
    treeheight_col = treeheight.cell,
    row_label = y,
    col_label = z,
    silent = TRUE,
    ...)
  
  if(!isTRUE(silent)) {
    grid::grid.newpage() 
    grid::grid.draw(sp$gtable)  
  }  
  
  invisible(sp)  
}


