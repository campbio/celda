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
  group.sd.median <- unlist(lapply(unique.label, function(l) mean(mat[which(class.label %in% l), ])/sd(mat[which(class.label %in% l), ])))
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
  t(t(x)/colSums(x)*1e6)
}

robust_scale <- function(x){
  madx<-ifelse(mad(x)==0,1,mad(x))
  median(x)/mad(x)
}
##ToDo:  Need to (1)  cluster by group(row/col label)
##       (2) within each group(row/col label) --> hierarchical clustering 
#' plot the heatmap of the counts data
#' @param counts the counts matrix 
#' @param z A numeric vector of cluster assignments for cell 
#' @param y A numeric vector of cluster assignments for gene
#' @param scale.log specify the transformation type of the matrix for (semi-)heatmap, can be "log","row"(z-acore by row),"col"(z-score by column), etc. #To be completed
#' @param scale.row specify the transformation type of the matrix for (semi-)heatmap, can be "log","row"(z-acore by row),"col"(z-score by column), etc. #To be completed
#' @param z.trim two element vector to specify the lower and upper cutoff of the z-score normalization result by default it is set to NULL so no trimming will be done.
#' @param scale_fun specify the function for scaling 
#' @param cluster.row boolean values determining if rows should be clustered
#' @param cluster.column boolean values determining if columns should be clustered
#' @example TODO
#' @export 
render_celda_heatmap <- function(counts, z=NULL, y=NULL, 
                                 scale.log=FALSE, scale.row=TRUE,
                                 scale_function=scale, normalize = "none",
                                 z.trim=c(-2,2), 
                                 cluster.row=TRUE, cluster.column = TRUE) {
  require(gtable)
  require(grid)
  require(scales)
  require(stats)
  require(RColorBrewer)
  require(grDevices)
  require(graphics)
  
  
  if(normalize =="cpm"){
    counts <- cpm(counts)
  }
  
  # matrix transformation for heatmap
  if(scale.log){
    counts <- log(counts+1)
  }
  
  if(scale.row){
    counts <- t(apply(counts, 1, scale_function))
    
    if(!is.null(z.trim)){
      if(length(z.trim)!=2) {
        stop("z.trim should be a 2 element vector specifying the lower and upper cutoffs")
      }
      z.trim<-sort(z.trim)
      counts[counts < z.trim[1]] <- z.trim[1]
      counts[counts > z.trim[2]] <- z.trim[2]
    }
  }
  
  
  #if null y or z 
  if(is.null(y)){
    y <- rep(1, nrow(counts))
  }
  if(is.null(z)){
    z <- rep(1, ncol(counts))
  }
  
  # order gene (row) 
  order.gene <- order_index(mat = counts, class.label = y, col = FALSE)
  # order cell (col)
  order.gene_cell <- order_index(mat = order.gene$mat, class.label = z, col = TRUE)
  
  counts <- order.gene_cell$mat
  y <- order.gene$class.label
  z <- order.gene_cell$class.label
  
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
  annotation_cell <- data.frame(cell_lable = as.factor(order.gene_cell$class.label)   ) 
  rownames(annotation_cell) <- colnames(counts)  # rowname should correspond to counts matrix's col(cell) name
  #  Set gene annotation
  annotation_gene <- data.frame(gene_label = as.factor(order.gene$class.label)  )
  rownames(annotation_gene) <- rownames(counts)  # rowname should correspond to counts matrix's row(gene) name
  
  ## Set color 
  #col.pal <- colorRampPalette(RColorBrewer::brewer.pal(n = 3, name = "RdYlBu"))(100)  # ToDo: need to be more flexible or fixed to a better color list
  #col.pal <- gplots::bluered(200)
  
  if(cluster.row & cluster.column){
    celda::semi_pheatmap(mat = counts, 
                         #color = colorRampPalette(c( "blue", "red"))(length(-12:12)),breaks=c(seq(0, 8.871147e-10, length.out = 11), 5.323741e-08 ),
                         #color = col.pal, 
                         cutree_rows = L,
                         cutree_cols = K,
                         annotation_row = annotation_gene,
                         annotation_col = annotation_cell,
                         row_label = y,
                         col_label = z,
                         scale = "none" , 
                         clustering_method =  "ward.D"   # should also add this parameter into celda_pheatmap 
    )
  }
  
  if(cluster.row & (!cluster.column)){
    celda::semi_pheatmap(mat = counts, 
                         cutree_rows = L,
                         cluster_cols = FALSE,
                         annotation_row = annotation_gene,
                         annotation_col = annotation_cell,
                         row_label = y,
                         clustering_method =  "ward.D"   # should also add this parameter into celda_pheatmap 
    )
    }
    
    
    if((!cluster.row) & cluster.column){
      celda::semi_pheatmap(mat = counts, 
                           cluster_rows = FALSE,
                           cutree_cols = K,
                           annotation_row = annotation_gene,
                           annotation_col = annotation_cell,
                           col_label = z,
                           clustering_method =  "ward.D"   # should also add this parameter into celda_pheatmap 
      )
      }
    
    if((!cluster.row) & (!cluster.column) ){
      celda::semi_pheatmap(mat = counts, 
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           annotation_row = annotation_gene,
                           annotation_col = annotation_cell,
                           clustering_method =  "ward.D"   # should also add this parameter into celda_pheatmap 
      )
      }
}
