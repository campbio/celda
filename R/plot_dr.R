#' Create a scatterplot given two dimensions from a data dimensionality reduction tool (e.g tSNE)
#' 
#' @param dim1 Numeric vector; first dimension from data dimensionality reduction output.
#' @param dim2 Numeric vector; second dimension from data dimensionality reduction output.
#' @param matrix Matrix, will contain cells/samples as columns and variable of interest as rows.
#' @param xlab Character vector, used as label for x axis. 
#' @param ylab Character vector, used as label for y axis. 
#' @param color_low Character vector of R colors available from the colors() function. The color will be used to signify the lowest values on the scale. 
#' @param color_mid Character vector of R colors available from the colors() function. The color will be used to signify the midpoint on the scale. 
#' @param color_high Character vector of R colors available from the colors() function. The color will be used to signify the highest values on the scale. 
#' @param var_label Character vector, used as label for the scale.
#' @export
plot_dr_grid <- function(dim1, dim2, matrix, xlab, ylab, color_low, color_mid, color_high, scale_label){
  df <- data.frame(dim1,dim2,t(as.data.frame(matrix)))
  m <- reshape2::melt(df, id.vars = c("dim1","dim2"))
  colnames(m) <- c(xlab,ylab,"facet",scale_label)
  ggplot2::ggplot(m, aes_string(x=xlab, y=ylab)) + ggplot2::geom_point(stat = "identity",ggplot2::aes_string(color = scale_label)) + 
    ggplot2::facet_wrap(~facet) + ggplot2::theme_bw() + ggplot2::scale_colour_gradient2(low = color_low, mid = color_mid, high = color_high, midpoint = (max(m[,4])-min(m[,4]))/2) + 
    ggplot2::theme(strip.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0,"lines"),
                   panel.background = element_blank(), axis.line = ggplot2::element_line(colour = "black"))
}

#' Create a scatterplot for each row of a normalized gene expression matrix where x and y axis are from a data dimensionality reduction tool.  
#' 
#' @param dim1 Numeric vector; first dimension from data dimensionality reduction output.
#' @param dim2 Numeric vector; second dimension from data dimensionality reduction output.
#' @param matrix Counts matrix, will have cell name for column name and gene name for row name.
#' @param trim A two element vector to specify the lower and upper cutoff for the data. Occurs after row scaling. Set to NULL to disable. Default c(-2,2).A two element vector to specify the lower and upper cutoff for the data. Occurs after normalization, log transformation, and row scaling. Set to NULL to disable. Default c(-2,2).
#' @param rescale Boolean. If TRUE, will rescale counts matrix on a 0~1 scale. Default FALSE.
#' @param xlab Character vector, used as label for x axis. Default "Dimension_1".
#' @param ylab Character vector, used as label for y axis. Default "Dimension_2".
#' @param color_low Character vector of R colors available from the colors() function. The color will be used to signify the lowest values on the scale. Default: 'blue'
#' @param color_mid Character vector of R colors available from the colors() function. The color will be used to signify the midpoint on the scale. Default: 'yellow'
#' @param color_high Character vector of R colors available from the colors() function. The color will be used to signify the highest values on the scale. Default: 'red'
#' @export 
plot_dr_gene <- function(dim1, dim2, matrix, trim = c(-2,2), rescale = FALSE, xlab = "Dimension_1", ylab = "Dimension_2", color_low = "blue", color_mid = "yellow", color_high = "red"){
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
  scale_label = "Expression"
  plot_dr_grid(dim1,dim2,matrix,xlab,ylab,color_low,color_mid,color_high, scale_label)
}

#' Create a scatterplot based off of a matrix containing the celda state probabilities per cell.
#' 
#' @param dim1 Numeric vector; first dimension from data dimensionality reduction output.
#' @param dim2 Numeric vector; second dimension from data dimensionality reduction output.
#' @param matrix Cell state probability matrix, will have cell name for column name and state probability for row name.
#' @param rescale Boolean. If TRUE z-score normalize the matrix. Default FALSE.
#' @param xlab Character vector, used as label for x axis. Default "Dimension_1".
#' @param ylab Character vector, used as label for y axis. Default "Dimension_2".
#' @param color_low Character vector of R colors available from the colors() function. The color will be used to signify the lowest values on the scale. Default: 'blue'
#' @param color_mid Character vector of R colors available from the colors() function. The color will be used to signify the midpoint on the scale. Default: 'yellow'
#' @param color_high Character vector of R colors available from the colors() function. The color will be used to signify the highest values on the scale. Default: 'red'
#' @export 
plot_dr_state <- function(dim1, dim2, matrix, rescale = FALSE, xlab = "Dimension_1", ylab = "Dimension_2", color_low = "blue", color_mid = "yellow", color_high = "red"){
  if(rescale == TRUE){
    for(x in 1:nrow(matrix)){ 
      matrix[x,] <- matrix[x,]-min(matrix[x,])
      matrix[x,] <- matrix[x,]/max(matrix[x,])
    }
  }
  scale_label = "Normalized_Prob"
  plot_dr_grid(dim1,dim2,matrix,xlab,ylab,color_low,color_mid,color_high, scale_label)
}

#' Create a scatterplot based on celda cluster labels.
#' 
#' @param dim1 Numeric vector; first dimension from data dimensionality reduction output.
#' @param dim2 Numeric vector; second dimension from data dimensionality reduction output.
#' @param cluster Vector; Contains cluster labels (e.g. from celda_C or celda_CG).
#' @param xlab Character vector, used as label for rows. Default "Dimension_1".
#' @param ylab Character vector, used as label for columns. Default "Dimension_2".
#' @param color_low Character vector of R colors available from the colors() function. The color will be used to signify the lowest values on the scale. Default: 'blue'
#' @param color_mid Character vector of R colors available from the colors() function. The color will be used to signify the midpoint on the scale. Default: 'yellow'
#' @param color_high Character vector of R colors available from the colors() function. The color will be used to signify the highest values on the scale. Default: 'red'
#' @export 
plot_dr_cluster <- function(dim1, dim2, cluster, xlab = "Dimension_1", ylab = "Dimension_2"){
  df <- data.frame(dim1,dim2,cluster)
  colnames(df) <- c(xlab,ylab,"Cluster") 
  ggplot2::ggplot(df, aes_string(x = xlab,y = ylab)) + ggplot2::geom_point(stat = "identity", ggplot2::aes(color = Cluster))
}
