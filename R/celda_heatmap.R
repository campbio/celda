##ToDo:  Need to (1)scale the row height accordingly; (2) pick more contradictory color;
##       (3)Aanotation need to change and; (4)  what else need to do?
 #' plot the heatmap of the counts data
 #' @param counts the counts matrix 
 #' @param K The number of clusters being considered  (Question1)or: Total number of cell populations??
 #' @param z A numeric vector of cluster assignments
 #' @param L Total number of transcriptional states
 #' @param col vector of colors used in heatmap
 #' @param cluster_gene boolean values determining if genes should be clustered
 #' @param cluster_cell boolean values determining if cells should be clustered
 #' @param annotation_gene data frame that specifies the annotations for genes 
 #' @param annotation_cell data frame that specifies the annotations for cells 
 #' @example TODO
 #' @export 
  celda_heatmap <- function(counts, K, z, L, col="YlOrBr", cluster_gene = TRUE, cluster_cell = FALSE, 
                            annotation_gene, annotation_cell) {
    ## Set row name to counts matrix 
    if(is.null(rownames(counts))){
      rownames(counts) <- 1:nrow(counts)
    } 
    else if(is.null(colnames(counts))) {
      colnames(counts) <- 1:ncol(counts)
    }
    ##-- Set cell annotation    # need to do 
    #annotaion_cell <- data.frame(pseudoanno = sample(c("Tcell","BCell"), nrow(counts), replace = T))   # ToDo: need to change 
    ##-- Set gene annotation
    #annotation_gene <- data.frame(pseudo=sample(1:6,nrow(counts),replace=T))
    
    ## Set color 
    col.pal <- colorRampPalette(RColorBrewer::brewer.pal(n = 9, name =col))(100)  # ToDo: need to be more flexible or fixed to a better color list
    pheatmap::pheatmap(counts, 
                       color = col.pal,
                       cluster_rows = cluster_gene,
                       cluster_cols = cluster_cell,
                       annotation_row = annotation_gene,
                       annotation_col = annotation_cell,
                       cutree_rows = L,   # Question1: not sure about this
                       cutree_cols = L,   # Question2: not sure about this either 
                       fontsize = 6.5,
                       fontsize_col = 5
                       )
  }