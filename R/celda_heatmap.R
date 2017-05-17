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
  
  return(list(mat=mat, class.label=class.label, ordlab=ordlab))
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
  
  return(list(mat=mat, class.label=class.label, ordlab=ordlab))
}

order_index.label <- function(mat, class.label, col=F) {   # order index : gather in group then order group
  if(col==T){
    mat <- t(mat)
  }
  
  # order matrix
  ordlab <-order( class.label, decreasing = FALSE)
  mat <- mat[ordlab, ]
  class.label <- class.label[ordlab]
  
  if(col==T){
    mat <- t(mat)
  }
  
  return(list(mat=mat, class.label=class.label, ordlab=ordlab))
}


cpm<- function(x){
  t(t(x)/colSums(x)*1e6)
}

robust_scale <- function(x){
  madx<-ifelse(mad(x)==0,1,mad(x))
  median(x)/mad(x)
}


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
#' @param annotation_cell a dataframe for the cell annotations (columns)
#' @param annotation_gene a dataframe for the gene annotations (rows)
#' @param col color for the heatmap
#' @param legend logical to determine if legend should be drawn or not
#' @param annotation_legend boolean value showing if the legend for annotation tracks should be drawn 
#' @param annotation_names_gene boolean value showing if the names for gene annotation tracks should be drawn 
#' @param annotation_names_cell boolean value showing if the names for cell annotation tracks should be drawn 
#' @param show_genenames boolean specifying if gene names are be shown
#' @param show_cellnames boolean specifying if cell names are be shown
#' @example TODO
#' @export 
render_celda_heatmap <- function(counts, z=NULL, y=NULL, 
                                 scale.log=FALSE, scale.row=TRUE,
                                 scale_function=scale, normalize = "none",
                                 z.trim=c(-2,2), 
                                 cluster.row=TRUE, cluster.column = TRUE,
                                 annotation_cell = NULL, annotation_gene = NULL, 
                                 col=colorRampPalette(c("#1E90FF","#FFFFFF","#CD2626"),space = "Lab")(100),
                                 legend = TRUE,
                                 annotation_legend = TRUE,
                                 annotation_names_gene = TRUE, 
                                 annotation_names_cell = TRUE,
                                 show_genenames = FALSE, 
                                 show_cellnames = FALSE) {
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
  if(cluster.row){
    order.gene <- order_index.median(mat = counts, class.label = y, col = FALSE) 
  }else{
    order.gene <- order_index.label(mat = counts, class.label = y, col = FALSE) 
  }
  # order cell (col)
  if(cluster.column){
    order.gene_cell <- order_index.median(mat = order.gene$mat, class.label = z, col = TRUE) 
  }else{
    order.gene_cell <- order_index.label(mat = order.gene$mat, class.label = z, col = TRUE)
  }
  
  counts <- order.gene_cell$mat
  #y <- order.gene$class.label
  #z <- order.gene_cell$class.label
  
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
    annotation_cell <- annotation_cell[order.gene_cell$ordlab, ]
    rownames(annotation_cell) <- colnames(counts)
  }else{  # annotation_cell is null
    annotation_cell <- data.frame(cell_label = as.factor(order.gene_cell$class.label)   ) 
    rownames(annotation_cell) <- colnames(counts)  # rowname should correspond to counts matrix's col(cell) name
  }
  #  Set gene annotation
  if(!is.null(annotation_gene)){
    annotation_gene$gene_label <- as.factor(y)  # has to be the original gene label order 
    annotation_gene <- annotation_gene[order.gene$ordlab,]
    rownames(annotation_gene) <- rownames(counts)
  }else{  #annotation_gene is  null 
    annotation_gene <- data.frame(gene_label = as.factor(order.gene$class.label)  )
    rownames(annotation_gene) <- rownames(counts)  # rowname should correspond to counts matrix's row(gene) name
  }
  
  # update the label 
  y <- order.gene$class.label
  z <- order.gene_cell$class.label
  
  
  ## Set color 
  #col.pal <- colorRampPalette(RColorBrewer::brewer.pal(n = 3, name = "RdYlBu"))(100)  # ToDo: need to be more flexible or fixed to a better color list
  #col.pal <- gplots::bluered(200)
  
  if(cluster.row & cluster.column){
    celda::semi_pheatmap(mat = counts, 
                         #color = colorRampPalette(c( "blue", "red"))(length(-12:12)),breaks=c(seq(0, 8.871147e-10, length.out = 11), 5.323741e-08 ),
                         color = col, 
                         cutree_rows = L,
                         cutree_cols = K,
                         annotation_row = annotation_gene,
                         annotation_col = annotation_cell,
                         row_label = y,
                         col_label = z,
                         legend = legend,
                         annotation_legend = annotation_legend, 
                         annotation_names_row = annotation_names_gene, 
                         annotation_names_col = annotation_names_cell,
                         show_rownames = show_genenames,
                         show_colnames = show_cellnames,
                         clustering_method =  "ward.D"   # should also add this parameter into celda_pheatmap 
    )
  }
  
  if(cluster.row & (!cluster.column)){
    celda::semi_pheatmap(mat = counts, 
                         color = col,
                         cutree_rows = L,
                         cluster_cols = FALSE,
                         annotation_row = annotation_gene,
                         annotation_col = annotation_cell,
                         row_label = y,
                         legend = legend,
                         annotation_legend = annotation_legend, 
                         annotation_names_row = annotation_names_gene, 
                         annotation_names_col = annotation_names_cell,
                         show_rownames = show_genenames,
                         show_colnames = show_cellnames,
                         clustering_method =  "ward.D"   # should also add this parameter into celda_pheatmap 
    )
    }
    
    
    if((!cluster.row) & cluster.column){
      celda::semi_pheatmap(mat = counts, 
                           color = col,
                           cluster_rows = FALSE,
                           cutree_cols = K,
                           annotation_row = annotation_gene,
                           annotation_col = annotation_cell,
                           col_label = z,
                           legend = legend,
                           annotation_legend = annotation_legend, 
                           annotation_names_row = annotation_names_gene, 
                           annotation_names_col = annotation_names_cell,
                           show_rownames = show_genenames,
                           show_colnames = show_cellnames,
                           clustering_method =  "ward.D"   # should also add this parameter into celda_pheatmap 
      )
      }
    
    if((!cluster.row) & (!cluster.column) ){
      celda::semi_pheatmap(mat = counts,
                           color = col,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           annotation_row = annotation_gene,
                           annotation_col = annotation_cell,
                           legend = legend,
                           annotation_legend = annotation_legend, 
                           annotation_names_row = annotation_names_gene, 
                           annotation_names_col = annotation_names_cell,
                           show_rownames = show_genenames,
                           show_colnames = show_cellnames,
                           clustering_method =  "ward.D"   # should also add this parameter into celda_pheatmap 
      )
      }
}
