##ToDo:  Need to (1)  cluster by group(row/col label)
##       (2) within each group(row/col label) --> hierarchical clustering 
#' plot the heatmap of the counts data
#' @param counts the counts matrix 
#' @param K The number of clusters being considered  (Question1)or: Total number of cell populations??
#' @param z A numeric vector of cluster assignments for cell 
#' @param L Total number of transcriptional states
#' @param y A numeric vector of cluster assignments for gene
#' @param scale_type specify the transformation type of the matrix for (semi-)heatmap, can be "log","row"(z-acore by row),"col"(z-score by column), etc. #To be completed 
#' @param cluster_gene ##  has to set to be FALSE in pheatmap --> need to improve this ourself 
#' @param cluster_cell ##  has to set to be FALSE in pheatmap --> need to improve this ourself 
#' @example TODO
#' @export 
celda_heatmap <- function(counts, K, z, L, y, scale_type="log") {
  
  # matrix transformation for heatmap
  if(scale_type=="log"){
    counts <- log(counts+1)
  }else if(scale_type=="row"){
    counts <- t(apply(counts, 1, scale, center=T))
  }else if(scale_type=="column"){
    counts <- apply(counts, 2, scale, center=T)
  }
  
  ## Set row and column name to counts matrix 
  if(is.null(rownames(counts))){
    rownames(counts) <- 1:nrow(counts)
  } 
  if(is.null(colnames(counts))) {
    colnames(counts) <- 1:ncol(counts)
  }
  

  
  ## order(cluster) count matrix   
  counts.orderCell_gene <- counts[order(y),order(z) ] 
  
  
  #  Set cell annotation    
  annotation_cell <- data.frame(cell_lable = as.factor(z[order(z)])  ) 
  rownames(annotation_cell) <- colnames(counts.orderCell_gene)  # rowname should correspond to counts matrix's col(cell) name
  #  Set gene annotation
  annotation_gene <- data.frame(gene_label = as.factor(y[order(y)])  )
  rownames(annotation_gene) <- rownames(counts.orderCell_gene)  # rowname should correspond to counts matrix's row(gene) name
  
  
  ## Set color 
  #col.pal <- colorRampPalette(RColorBrewer::brewer.pal(n = 2, name = col))(100)  # ToDo: need to be more flexible or fixed to a better color list
  #col.pal <- gplots::bluered(100)
  
  celda::semi_pheatmap(counts.orderCell_gene, 
                       #color = col.pal,
                       cutree_rows = L,
                       cutree_cols = K,
                       annotation_row = annotation_gene,
                       annotation_col = annotation_cell,
                       row_label = y[order(y)],
                       col_label = z[order(z)]
                       )
}