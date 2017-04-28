order_index <- function(mat, class.label, col=F) {   # order index : gather in group then order group
  if(col==T){
    mat <- t(mat)
  }
  
  rowmean.Mat  <- rowMeans(mat)
  df  <- aggregate(x = rowmean.Mat, by=list(class.label) , mean)   # a dataframe with label and group mean 
  row.Gmean.mat <- df[,2][match(class.label, df[,1])]
  
  # order matrix 
  ordlab <-order( row.Gmean.mat, class.label, decreasing = TRUE)
  mat <- mat[ordlab, ]
  class.label <- class.label[ordlab]
  
  if(col==T){
    mat <- t(mat)
  }
  
  return(list(mat=mat, class.label=class.label))
}



order_index.median <- function(mat, class.label, col=F) {   # order index : gather in group then order group
  if(col==T){
    mat <- t(mat)
  }
  
  unique.label <- unique(class.label)
  group.sd.median <- unlist(lapply(unique.label, function(l) length(which(class.label %in% l))*mean(mat[which(class.label %in% l), ])/sd(mat[which(class.label %in% l), ])))
  df <- data.frame(Group.1=unique.label, x=group.sd.median)  # a dataframe with label and group mean 
  row.Gmean.mat <- df[,2][match(class.label, df[,1])]
  
  # order matrix
  ordlab <-order( row.Gmean.mat, class.label, decreasing = TRUE)
  mat <- mat[ordlab, ]
  class.label <- class.label[ordlab]
  
  if(col==T){
    mat <- t(mat)
  }
  
  return(list(mat=mat, class.label=class.label))
}

cpm<- function(x){
  t(t(x)/colSums(x)*1000000)
}

##ToDo:  Need to (1)  cluster by group(row/col label)
##       (2) within each group(row/col label) --> hierarchical clustering 
#' plot the heatmap of the counts data
#' @param counts the counts matrix 
#' @param K The number of clusters being considered  (Question1)or: Total number of cell populations??
#' @param z A numeric vector of cluster assignments for cell 
#' @param L Total number of transcriptional states
#' @param y A numeric vector of cluster assignments for gene
#' @param scale.log specify the transformation type of the matrix for (semi-)heatmap, can be "log","row"(z-acore by row),"col"(z-score by column), etc. #To be completed
#' @param scale.row specify the transformation type of the matrix for (semi-)heatmap, can be "log","row"(z-acore by row),"col"(z-score by column), etc. #To be completed
#' @param z.trim two element vector to specify the lower and upper cutoff of the z-score normalization result by default it is set to NULL so no trimming will be done.
#' @example TODO
#' @export 
celda_heatmap <- function(counts, K, z, L, y, scale.log=FALSE, scale.row=FALSE, normalize = c("none","cpm"),
                           z.trim=NULL) {
  require(gtable)
  require(grid)
  require(scales)
  require(stats)
  require(RColorBrewer)
  require(grDevices)
  require(graphics)
  
  if(length(z.trim!=2)) {
    stop("z.trim should be a 2 element vector specifying the lower and upper cutoffs")
  }
  
  if(normalize =="cpm"){
    counts <- cpm(counts)
  }
  
  # matrix transformation for heatmap
  if(scale.log){
    counts <- log(counts+1)
  }
  
  if(scale.row){
    counts <- t(apply(counts, 1, scale, center=T))
  }
  
  if(!is.null(z.trim)){
    z.trim<-sort(z.trim)
    counts[counts < z.trim[1]] <- z.trim[1]
    counts[counts > z.trim[2]] <- z.trim[2]
  }
  
  # order gene (row) 
  order.gene <- order_index.median(mat = counts, class.label = y, col = FALSE)
  # order cell (col)
  order.gene_cell <- order_index.median(mat = order.gene$mat, class.label = z, col = TRUE)
  
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
                       col_label = z,
                       clustering_method =  "ward.D"   # should also add this parameter into celda_pheatmap 
                       )

}