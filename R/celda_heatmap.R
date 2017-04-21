order_index <- function(mat, class.label, col=F) {   # order index : gather in group then order group
  if(col==T){
    mat <- t(mat)
  }
  
  rowmean.Mat  <- rowMeans(mat)
  df  <- aggregate(x = rowmean.Mat, by=list(class.label) , mean)   # a dataframe with label and group mean 
  row.Gmean.mat <- df[,2][match(class.label, df[,1])]
  
  # order matrix 
  mat <- mat[order( row.Gmean.mat, class.label, decreasing = TRUE), ]
  class.label <- class.label[order(row.Gmean.mat, class.label, decreasing = TRUE)]
  
  if(col==T){
    mat <- t(mat)
  }
  
  return(list(mat=mat, class.label=class.label))
}


##ToDo:  Need to (1)  cluster by group(row/col label)
##       (2) within each group(row/col label) --> hierarchical clustering 
#' plot the heatmap of the counts data
#' @param counts the counts matrix 
#' @param K The number of clusters being considered  (Question1)or: Total number of cell populations??
#' @param z A numeric vector of cluster assignments for cell 
#' @param L Total number of transcriptional states
#' @param y A numeric vector of cluster assignments for gene
#' @param scale_type specify the transformation type of the matrix for (semi-)heatmap, can be "log","row"(z-acore by row),"col"(z-score by column), etc. #To be completed 
#' @example TODO
#' @export 
celda_heatmap <- function(counts, K, z, L, y, scale_type="row") {
  
  # matrix transformation for heatmap
  if(scale_type=="log"){
    counts <- log(counts+1)
  }else if(scale_type=="row"){
    counts <- t(apply(counts, 1, scale, center=T))
  } # else use the raw matrix 
  
  # order gene (row) 
  order.gene <- order_index(mat = counts, class.label = y, col = FALSE)
  # order cell (col)
  order.gene_cell <- order_index(mat = order.gene$mat, class.label = z, col = TRUE)
  
  counts <- order.gene_cell$mat
  y <- order.gene$class.label
  z <- order.gene_cell$class.label
  
  
  ## Set row and column name to counts matrix 
  if(is.null(rownames(counts))){
    rownames(counts) <- 1:nrow(counts)
  } 
  if(is.null(colnames(counts))) {
    colnames(counts) <- 1:ncol(counts)
  }
  
  
  
  #  Set cell annotation    
  annotation_cell <- data.frame(cell_lable = as.factor(order.gene_cell$class.label)   ) 
  rownames(annotation_cell) <- colnames(counts)  # rowname should correspond to counts matrix's col(cell) name
  #  Set gene annotation
  annotation_gene <- data.frame(gene_label = as.factor(order.gene$class.label)  )
  rownames(annotation_gene) <- rownames(counts)  # rowname should correspond to counts matrix's row(gene) name
  
  ## Set color 
  #col.pal <- colorRampPalette(RColorBrewer::brewer.pal(n = 2, name = col))(100)  # ToDo: need to be more flexible or fixed to a better color list
  #col.pal <- gplots::bluered(100)
  
  celda::semi_pheatmap(mat = counts, 
                       #color = col.pal,
                       cutree_rows = L,
                       cutree_cols = K,
                       annotation_row = annotation_gene,
                       annotation_col = annotation_cell,
                       row_label = y,
                       col_label = z
                       )

}