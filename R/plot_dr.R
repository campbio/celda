#' Create a scatterplot given two dimensions from a data dimensionality reduction tool (e.g tSNE)
#' 
#' @param dim1 Numeric vector; first dimension from data dimensionality reduction output.
#' @param dim2 Numeric vector; second dimension from data dimensionality reduction output.
#' @param matrix Matrix, will contain cells/samples as columns and variable of interest as rows.
#' @param size Numeic vector; size of point on plot. 
#' @param xlab Character vector, used as label for x axis. 
#' @param ylab Character vector, used as label for y axis. 
#' @param color_low Character vector of R colors available from the colors() function. The color will be used to signify the lowest values on the scale. 
#' @param color_mid Character vector of R colors available from the colors() function. The color will be used to signify the midpoint on the scale. 
#' @param color_high Character vector of R colors available from the colors() function. The color will be used to signify the highest values on the scale. 
#' @param var_label Character vector, used as label for the scale.
#' @export
plotDrGrid <- function(dim1, dim2, matrix, size, xlab, ylab, color_low, color_mid, color_high, var_label){
  df <- data.frame(dim1,dim2,t(as.data.frame(matrix)))
  m <- reshape2::melt(df, id.vars = c("dim1","dim2"))
  colnames(m) <- c(xlab,ylab,"facet",var_label)
  ggplot2::ggplot(m, ggplot2::aes_string(x=xlab, y=ylab)) + ggplot2::geom_point(stat = "identity", size = size, ggplot2::aes_string(color = var_label)) + 
    ggplot2::facet_wrap(~facet) + ggplot2::theme_bw() + ggplot2::scale_colour_gradient2(low = color_low, high = color_high, mid = color_mid, midpoint = (max(m[,4])-min(m[,4]))/2 ,name = gsub("_"," ",var_label)) + 
    ggplot2::theme(strip.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.spacing = unit(0,"lines"),
                   panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))
}

#' Create a scatterplot for each row of a normalized gene expression matrix where x and y axis are from a data dimensionality reduction tool.  
#' 
#' @param dim1 Numeric vector; first dimension from data dimensionality reduction output.
#' @param dim2 Numeric vector; second dimension from data dimensionality reduction output.
#' @param matrix Counts matrix, will have cell name for column name and gene name for row name.
#' @param trim A two element vector to specify the lower and upper cutoff for the data. Occurs after row scaling. Set to NULL to disable. Default c(-2,2).A two element vector to specify the lower and upper cutoff for the data. Occurs after normalization, log transformation, and row scaling. Set to NULL to disable. Default c(-2,2).
#' @param rescale Boolean. If TRUE, will rescale counts matrix on a 0~1 scale. Default TRUE.
#' @param size Numeic vector; size of point on plot. Default = 1.
#' @param xlab Character vector, used as label for x axis. Default "Dimension_1".
#' @param ylab Character vector, used as label for y axis. Default "Dimension_2".
#' @param color_low Character vector of R colors available from the colors() function. The color will be used to signify the lowest values on the scale. Default: 'grey'
#' @param color_mid Character vector of R colors available from the colors() function. The color will be used to signify the midpoint on the scale. 
#' @param color_high Character vector of R colors available from the colors() function. The color will be used to signify the highest values on the scale. Default: 'blue'
#' @export 
plotDrGene <- function(dim1, dim2, matrix, trim = c(-2,2), rescale = TRUE, size = 1, xlab = "Dimension_1", ylab = "Dimension_2", color_low = "grey", color_mid = NULL, color_high = "blue"){
  if(rescale == TRUE){
    matrix <- t(scale(t(matrix)))
    if(!is.null(trim)){
      if(length(trim) != 2) {
        stop("'trim' should be a 2 element vector specifying the lower and upper boundaries")
      }
      trim<-sort(trim)
      matrix[matrix < trim[1]] <- trim[1]
      matrix[matrix > trim[2]] <- trim[2]
    }
  }
  var_label = "Expression"
  plotDrGrid(dim1,dim2,matrix,size,xlab,ylab,color_low,color_mid,color_high, var_label)
}

#' Create a scatterplot based off of a matrix containing the celda state probabilities per cell.
#' 
#' @param dim1 Numeric vector; first dimension from data dimensionality reduction output.
#' @param dim2 Numeric vector; second dimension from data dimensionality reduction output.
#' @param matrix Cell state probability matrix, will have cell name for column name and state probability for row name.
#' @param rescale Boolean. If TRUE z-score normalize the matrix. Default TRUE.
#' @param size Numeic vector; size of point on plot. Default = 1.
#' @param xlab Character vector, used as label for x axis. Default "Dimension_1".
#' @param ylab Character vector, used as label for y axis. Default "Dimension_2".
#' @param color_low Character vector of R colors available from the colors() function. The color will be used to signify the lowest values on the scale. Default: 'grey'
#' @param color_mid Character vector of R colors available from the colors() function. The color will be used to signify the midpoint on the scale. 
#' @param color_high Character vector of R colors available from the colors() function. The color will be used to signify the highest values on the scale. Default: 'blue'
#' @export 
plotDrState <- function(dim1, dim2, matrix, rescale = TRUE, size = 1, xlab = "Dimension_1", ylab = "Dimension_2", color_low = "grey", color_mid = NULL, color_high = "blue"){
  if(rescale == TRUE){
    for(x in 1:nrow(matrix)){ 
      matrix[x,] <- matrix[x,]-min(matrix[x,])
      matrix[x,] <- matrix[x,]/max(matrix[x,])
      var_label = "Scaled_Probability"
    }
  }else{
    var_label = "Probability"
  }
  plotDrGrid(dim1,dim2,matrix,size,xlab,ylab,color_low,color_mid,color_high, var_label)
}

#' Create a scatterplot based on celda cluster labels.
#' 
#' @param dim1 Numeric vector; first dimension from data dimensionality reduction output.
#' @param dim2 Numeric vector; second dimension from data dimensionality reduction output.
#' @param cluster Vector; Contains cluster labels (e.g. from celda_C or celda_CG).
#' @param size Numeic vector; size of point on plot.Default 1.
#' @param xlab Character vector, used as label for rows. Default "Dimension_1".
#' @param ylab Character vector, used as label for columns. Default "Dimension_2".
#' @param specific_clusters Numeic vector; Contains specific cluster labels.
#' @export 
plotDrCluster <- function(dim1, dim2, cluster, size = 1, xlab = "Dimension_1", ylab = "Dimension_2", specific_clusters = NULL){
  df <- data.frame(dim1, dim2, cluster)
  colnames(df) <- c(xlab,ylab,"Cluster")
  if(!is.null(specific_clusters)){
    df[3][!(df[[3]] %in% specific_clusters),] <- 0
    df <- df[order(df[[3]]),]
    df[3] <- as.factor(df[[3]])
    cluster_colors <- c('grey',distinct_colors(nlevels(as.factor(cluster)))[sort(specific_clusters)])
  } else{
    cluster_colors <- distinct_colors(nlevels(as.factor(cluster)))
    df[3] <- as.factor(df[[3]])
  }
  ggplot2::ggplot(df, ggplot2::aes_string(x = xlab, y = ylab)) +
    ggplot2::geom_point(stat = "identity", size = size, ggplot2::aes(color = Cluster)) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(color = "black")) +
    ggplot2::scale_color_manual(values = cluster_colors) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 1)))
}

#' Runs tSNE via Rtsne based on the CELDA model and specified cell states.
#' 
#' @param counts Counts matrix, should have cell name for column name and gene name for row name.
#' @param celda.mod Celda model to use for tsne. class "celda_C","celda_G" or "celda_CG".
#' @param states Numeric vector; determines which cell populations to use for tsne. If none are defined, all states will be used.
#' @param perplexity Numeric vector; determines perplexity for tsne. Default 20.
#' @param max.iter Numeric vector; determines iterations for tsne. Default 1000.
#' @param distance Character vector; determines which distance metric to use for tsne. Options: cosine, hellinger, spearman.
#' @export
celdaTsne = function(counts, celda.mod, states=NULL, perplexity=20, max.iter=2500, distance="hellinger") {
  if (!isTRUE(class(celda.mod) %in% c("celda_CG","celda_C","celda_G"))) {
    stop("celda.mod argument is not of class celda_C, celda_G or celda_CG")
  } 
  
  if (class(celda.mod) == "celda_CG") {
    fm = factorizeMatrix(counts=counts, celda.mod=celda.mod, type="counts")
    
    states.to.use = 1:nrow(fm$counts$cell.states)
    if (!is.null(states)) {
      if (!all(states %in% states.to.use)) {
        stop("'states' must be a vector of numbers between 1 and ", states.to.use, ".")
      }
      states.to.use = states 
    } 
    new.counts = fm$counts$cell.states[states.to.use,]
    norm = normalizeCounts(new.counts, scale.factor=1)
  } else {
    norm = normalizeCounts(counts = counts, scale.factor = 1)
  }
  
  distance = match.arg(distance, choices = c("hellinger","cosine","spearman"))
  if (distance == "cosine") {
    d = cosineDist(norm)  
  } else if(distance == "hellinger") {
    d = hellingerDist(norm)  
  } else if(distance == "spearman") {
    d = spearmanDist(norm)
  } else {
    stop("distances must be either 'cosine' or 'hellinger' or 'spearman")
  }
  
  do.pca = class(celda.mod) == "celda_C"
  res = Rtsne::Rtsne(d, pca=do.pca, max_iter=max.iter, perplexity = perplexity, 
                     check_duplicates = FALSE, is_distance = TRUE)$Y
  colnames(res) = c("tsne_1", "tsne_2")
  return(res)
}

