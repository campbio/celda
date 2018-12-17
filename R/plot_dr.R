#' @title Mapping the dimensionality reduction plot 
#' @description Creates a scatterplot given two dimensions from a data dimensionality reduction tool (e.g tSNE) output.
#' 
#' @param dim1 Numeric vector. First dimension from data dimensionality reduction output.
#' @param dim2 Numeric vector. Second dimension from data dimensionality reduction output.
#' @param matrix Numeric matrix. Each row of the matrix will be plotted as a separate facet. 
#' @param size Numeric. Sets size of point on plot. Default 1. 
#' @param xlab Character vector. Label for the x-axis. Default 'Dimension_1'. 
#' @param ylab Character vector. Label for the y-axis. Default 'Dimension_2'. 
#' @param color_low Character. A color available from `colors()`. The color will be used to signify the lowest values on the scale. Default 'grey'. 
#' @param color_mid Character. A color available from `colors()`. The color will be used to signify the midpoint on the scale. 
#' @param color_high Character. A color available from `colors()`. The color will be used to signify the highest values on the scale. Default 'blue'.
#' @param var_label Character vector. Title for the color legend. 
#' @return The plot as a ggplot object
#' @examples
#' \donttest{
#' celda.tsne <- celdaTsne(counts = celda.CG.sim$counts, celda.mod = celda.CG.mod)
#' plotDimReduceGrid(celda.tsne[,1], celda.tsne[,2], matrix = celda.CG.sim$counts, 
#'                   xlab = "Dimension1", ylab = "Dimension 2", var_label = "tsne", 
#'                   size = 1, color_low = "grey", color_mid = NULL, color_high = "blue")}
#' @export
plotDimReduceGrid = function(dim1, dim2, matrix, size, xlab, ylab, color_low, color_mid, color_high, var_label){
  df = data.frame(dim1,dim2,t(as.data.frame(matrix)))
  na.ix = is.na(dim1) | is.na(dim2)
  df = df[!na.ix,]
  
  m = reshape2::melt(df, id.vars = c("dim1","dim2"))
  colnames(m) = c(xlab,ylab,"facet",var_label)
  ggplot2::ggplot(m, ggplot2::aes_string(x=xlab, y=ylab)) + ggplot2::geom_point(stat = "identity", size = size, ggplot2::aes_string(color = var_label)) + 
    ggplot2::facet_wrap(~facet) + ggplot2::theme_bw() + ggplot2::scale_colour_gradient2(low = color_low, high = color_high, mid = color_mid, midpoint = (max(m[,4])+min(m[,4]))/2 ,name = gsub("_"," ",var_label)) + 
    ggplot2::theme(strip.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.spacing = unit(0,"lines"),
                   panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))
}

#' @title Plotting feature expression on a dimensionality reduction plot
#' @description Create a scatterplot for each row of a normalized gene expression matrix where x and y axis are from a data dimensionality reduction tool. The cells are colored by expression of the specified feature.
#' 
#' @param dim1 Numeric vector. First dimension from data dimensionality reduction output.
#' @param dim2 Numeric vector. Second dimension from data dimensionality reduction output.
#' @param counts Integer matrix. Rows represent features and columns represent cells. 
#' @param features Character vector. Uses these genes for plotting.
#' @param normalize Logical. Whether to normalize the columns of `counts`. Default TRUE.
#' @param exact.match Logical. Whether an exact match or a partial match using `grep()` is used to look up the feature in the rownames of the counts matrix. Default TRUE. 
#' @param trim Numeric vector. Vector of length two that specifies the lower and upper bounds for the data. This threshold is applied after row scaling. Set to NULL to disable. Default c(-2,2). 
#' @param size Numeric. Sets size of point on plot. Default 1.
#' @param xlab Character vector. Label for the x-axis. Default "Dimension_1".
#' @param ylab Character vector. Label for the y-axis. Default "Dimension_2".
#' @param color_low Character. A color available from `colors()`. The color will be used to signify the lowest values on the scale. Default 'grey'.
#' @param color_mid Character. A color available from `colors()`. The color will be used to signify the midpoint on the scale. 
#' @param color_high Character. A color available from `colors()`. The color will be used to signify the highest values on the scale. Default 'blue'.
#' @return The plot as a ggplot object
#' @examples
#' \donttest{
#' celda.tsne <- celdaTsne(counts = celda.CG.sim$counts,
#'                         celda.mod = celda.CG.mod)
#' plotDimReduceFeature(dim1 = celda.tsne[,1], dim2 = celda.tsne[,2],
#'                      counts = celda.CG.sim$counts,
#'                      features = c("Gene_99"), exact.match = TRUE)
#'}
#' @export 
plotDimReduceFeature = function(dim1, dim2, counts, features, normalize = TRUE, exact.match = TRUE, trim = c(-2,2), size = 1, xlab = "Dimension_1", ylab = "Dimension_2", color_low = "grey", color_mid = NULL, color_high = "blue"){
  if(isTRUE(normalize)){
    counts = normalizeCounts(counts, transformation.fun = sqrt, scale.fun = base::scale) 
  }
  if(is.null(features)){
    stop("at least one feature is required to create a plot")
  }
  if(!is.null(trim)){
    if(length(trim) != 2) {
      stop("'trim' should be a 2 element vector specifying the lower and upper boundaries")
    }
    trim = sort(trim)
    counts[counts < trim[1]] = trim[1]
    counts[counts > trim[2]] = trim[2]
  }  
  var_label = "Expression"
  
  if(!isTRUE(exact.match)){
    features.indices = c()  
    for(gene in features){
      features.indices = c(features.indices, grep(gene, rownames(counts)))
    }
    counts = counts[features.indices, , drop = FALSE]
  } else {
    features.not.found = setdiff(features, intersect(features, rownames(counts)))
    if (length(features.not.found) > 0) {
      if (length(features.not.found) == length(features)) {
        stop("None of the provided features had matching rownames in the provided counts matrix.")
      }
      warning(paste0("The following features were not present in the provided count matrix: ",
                     paste0(features.not.found, ",")))
    }
    features.found = setdiff(features, features.not.found)
    counts = counts[features.found, , drop = FALSE]
  }
  plotDimReduceGrid(dim1, dim2, counts, size, xlab, ylab, color_low, color_mid, color_high, var_label)
}

#' @title Plotting the Celda module probability on a dimensionality reduction plot
#' @description Create a scatterplot for each row of a normalized gene expression matrix where x and y axis are from a data dimensionality reduction tool. The cells are colored by the module probability(s).
#' 
#' @param dim1 Numeric vector. First dimension from data dimensionality reduction output.
#' @param dim2 Numeric vector. Second dimension from data dimensionality reduction output.
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`. 
#' @param celda.mod Celda object of class "celda_G" or "celda_CG".
#' @param modules Character vector. Module(s) from celda model to be plotted.
#' @param rescale Logical. Whether rows of the matrix should be rescaled to [0,1]. Default TRUE.
#' @param size Numeric. Sets size of point on plot. Default 1.
#' @param xlab Character vector. Label for the x-axis. Default "Dimension_1".
#' @param ylab Character vector. Label for the y-axis. Default "Dimension_2".
#' @param color_low Character. A color available from `colors()`. The color will be used to signify the lowest values on the scale. Default 'grey'.
#' @param color_mid Character. A color available from `colors()`. The color will be used to signify the midpoint on the scale. 
#' @param color_high Character. A color available from `colors()`. The color will be used to signify the highest values on the scale. Default 'blue'.
#' @return The plot as a ggplot object
#' @examples
#' \donttest{
#' celda.tsne <- celdaTsne(counts = celda.CG.sim$counts, 
#'                         celda.mod = celda.CG.mod)
#' plotDimReduceModule(dim1 = celda.tsne[,1], dim2 = celda.tsne[,2], 
#'                     counts = celda.CG.sim$counts, celda.mod = celda.CG.mod,
#'                     modules = c("L1","L2"))
#'}
#' @export 
plotDimReduceModule = function(dim1, dim2, counts, celda.mod, modules = NULL, rescale = TRUE, size = 1, xlab = "Dimension_1", ylab = "Dimension_2", color_low = "grey", color_mid = NULL, color_high = "blue"){
  
  factorized = factorizeMatrix(celda.mod = celda.mod, counts = counts)
  matrix = factorized$proportions$cell
  if(rescale == TRUE){
    for(x in 1:nrow(matrix)){ 
      matrix[x,] = matrix[x,]-min(matrix[x,])
      matrix[x,] = matrix[x,]/max(matrix[x,])
      var_label = "Scaled_Probability"
    }
  }else{
    var_label = "Probability"
  }
  
  if(!is.null(modules)){
    if(length(rownames(matrix)[rownames(matrix) %in% modules]) < 1){
      stop("All modules selected do not exist in the model.")
    }
    
    matrix = matrix[which(rownames(matrix) %in% modules), , drop=FALSE]
    matrix = matrix[match(rownames(matrix), modules), , drop=FALSE]
  }
  plotDimReduceGrid(dim1,dim2,matrix,size,xlab,ylab,color_low,color_mid,color_high, var_label)
}

#' @title Plotting the cell labels on a dimensionality reduction plot
#' @description Create a scatterplot for each row of a normalized gene expression matrix where x and y axis are from a data dimensionality reduction tool. The cells are colored by its given `cluster` label.
#' 
#' @param dim1 Numeric vector. First dimension from data dimensionality reduction output.
#' @param dim2 Numeric vector. Second dimension from data dimensionality reduction output.
#' @param cluster Integer vector. Contains cluster labels for each cell. 
#' @param size Numeric. Sets size of point on plot. Default 1.
#' @param xlab Character vector. Label for the x-axis. Default "Dimension_1".
#' @param ylab Character vector. Label for the y-axis. Default "Dimension_2".
#' @param specific_clusters Numeric vector. Only color cells in the specified clusters. All other cells will be grey. If NULL, all clusters will be colored. Default NULL. 
#' @return The plot as a ggplot object
#' @examples
#' \donttest{
#' celda.tsne <- celdaTsne(counts = celda.CG.sim$counts, celda.mod = celda.CG.mod)
#' plotDimReduceCluster(dim1 = celda.tsne[,1], dim2 = celda.tsne[,2],
#'                      cluster = as.factor(z(celda.CG.mod)),
#'                      specific_clusters = c(1,2,3))
#' }
#' @export 
plotDimReduceCluster = function(dim1, dim2, cluster, size = 1, xlab = "Dimension_1", ylab = "Dimension_2", specific_clusters = NULL){
  df = data.frame(dim1, dim2, cluster)
  colnames(df) = c(xlab, ylab, "Cluster")
  na.ix = is.na(dim1) | is.na(dim2)
  df = df[!na.ix,]
  df[3] = as.factor(df[[3]])
  cluster_colors = distinct_colors(nlevels(as.factor(cluster)))
  if(!is.null(specific_clusters)){
    cluster_colors[!levels(df[[3]]) %in% specific_clusters] = "gray92"
  }
  ggplot2::ggplot(df, ggplot2::aes_string(x = xlab, y = ylab)) +
    ggplot2::geom_point(stat = "identity", size = size, 
                        ggplot2::aes_string(color = "Cluster")) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(color = "black")) +
    ggplot2::scale_color_manual(values = cluster_colors) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 1)))
}





# Run the t-SNE algorithm for dimensionality reduction
# 
# @param norm Normalized count matrix.
# @param perplexity Numeric vector. Determines perplexity for tsne. Default 20.
# @param max.iter Numeric vector. Determines iterations for tsne. Default 1000.
# @param seed Integer. Seed for random number generation. Defaults to 12345.
# @param do.pca Logical. Whether to perform dimensionality reduction with PCA before tSNE.
# @param initial.dims Integer. Number of dimensions from PCA to use as input in tSNE.
calculateTsne = function(norm, perplexity=20, max.iter=2500, seed=12345, do.pca=FALSE, initial.dims = 20) {

  set.seed(seed)
  res = Rtsne::Rtsne(norm, pca=do.pca, max_iter=max.iter, perplexity = perplexity, 
                     check_duplicates = FALSE, is_distance = FALSE, initial_dims=initial.dims)$Y
  return(res)                     
}

