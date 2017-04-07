##ToDo:  Need to (1)  cluster by group(row/col label)
##       (2) within each group(row/col label) --> hierarchical clustering 
#' plot the heatmap of the counts data
#' @param counts the counts matrix 
#' @param K The number of clusters being considered  (Question1)or: Total number of cell populations??
#' @param z A numeric vector of cluster assignments for cell 
#' @param L Total number of transcriptional states
#' @param y A numeric vector of cluster assignments for gene
#' @param col vector of colors used in heatmap
#' @param cluster_gene ##  has to set to be FALSE in pheatmap --> need to improve this ourself 
#' @param cluster_cell ##  has to set to be FALSE in pheatmap --> need to improve this ourself 
#' @example TODO
#' @export 
celda_heatmap <- function(counts, K, z, L, y,  col="YlOrBr") {
  ## Set row and column name to counts matrix 
  if(is.null(rownames(counts))){
    rownames(counts) <- 1:nrow(counts)
  } 
  if(is.null(colnames(counts))) {
    colnames(counts) <- 1:ncol(counts)
  }
  
  
  ## order(cluster) count matrix   
  counts.orderCell_gene <- counts[order(y),order(z) ] 
  
  
  lLenth.cell = table(z)   # cell : label length 
  lLenth.gene = table(y)   # gene : label length 
  
  pos.lcell = cumsum(lLenth.cell)  # cell : position of the label gap
  pos.lgene = cumsum(lLenth.gene)  # gene : position of the label gap
  
  
  
  #  Set cell annotation    
  annotation_cell <- data.frame(cell_lable = as.factor(z[order(z)])  ) 
  rownames(annotation_cell) <- colnames(counts.orderCell_gene)  # rowname should correspond to counts matrix's col(cell) name
  #  Set gene annotation
  annotation_gene <- data.frame(gene_label = as.factor(y[order(y)])  )
  rownames(annotation_gene) <- rownames(counts.orderCell_gene)  # rowname should correspond to counts matrix's row(gene) name
  
  
  ## Set color 
  #col.pal <- colorRampPalette(RColorBrewer::brewer.pal(n = 2, name = col))(100)  # ToDo: need to be more flexible or fixed to a better color list
  col.pal <- gplots::bluered(100)
  
  pheatmap::pheatmap(counts.orderCell_gene, 
                     color = col.pal,
                     gaps_row = pos.lgene[-L],
                     gaps_col = pos.lcell[-K],
                     annotation_row = annotation_gene,
                     annotation_col = annotation_cell,
                     cluster_rows = FALSE,   # has to set to be FALSE
                     cluster_cols = FALSE,   # has to set to be FALSE 
                     fontsize = 6.5,
                     fontsize_col = 5
  )
}