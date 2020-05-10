#' @title Mapping the dimensionality reduction plot
#' @description Creates a scatterplot given two dimensions from a data
#'  dimensionality reduction tool (e.g tSNE) output.
#' @param dim1 Numeric vector. First dimension from data dimensionality
#'  reduction output.
#' @param dim2 Numeric vector. Second dimension from data dimensionality
#'  reduction output.
#' @param matrix Numeric matrix. Each row of the matrix will be plotted as
#'  a separate facet.
#' @param size Numeric. Sets size of point on plot. Default 1.
#' @param xlab Character vector. Label for the x-axis. Default 'Dimension_1'.
#' @param ylab Character vector. Label for the y-axis. Default 'Dimension_2'.
#' @param colorLow Character. A color available from `colors()`.
#'  The color will be used to signify the lowest values on the scale.
#'  Default 'grey'.
#' @param colorMid Character. A color available from `colors()`.
#'  The color will be used to signify the midpoint on the scale.
#' @param colorHigh Character. A color available from `colors()`.
#'  The color will be used to signify the highest values on the scale.
#'  Default 'blue'.
#' @param midpoint Numeric. The value indicating the midpoint of the
#' diverging color scheme. If \code{NULL}, defaults to the mean
#' with 10 percent of values trimmed. Default \code{NULL}.
#' @param varLabel Character vector. Title for the color legend.
#' @param ncol Integer. Passed to \link[ggplot2]{facet_wrap}. Specify the
#'  number of columns for facet wrap.
#' @param headers Character vector. If `NULL`, the corresponding rownames are
#'  used as labels. Otherwise, these headers are used to label the genes.
#' @return The plot as a ggplot object
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' celdaTsne <- celdaTsne(
#'   counts = celdaCGSim$counts,
#'   celdaMod = celdaCGMod
#' )
#' plotDimReduceGrid(celdaTsne[, 1],
#'   celdaTsne[, 2],
#'   matrix = celdaCGSim$counts,
#'   xlab = "Dimension1",
#'   ylab = "Dimension2",
#'   varLabel = "tsne",
#'   size = 1,
#'   colorLow = "grey",
#'   colorMid = NULL,
#'   colorHigh = "blue"
#' )
#' @importFrom reshape2 melt
#' @export
plotDimReduceGrid <- function(dim1,
                              dim2,
                              matrix,
                              size = 1,
                              xlab = "Dimension_1",
                              ylab = "Dimension_2",
                              colorLow = "grey",
                              colorMid = NULL,
                              colorHigh = "blue",
                              midpoint = NULL,
                              varLabel = NULL,
                              ncol = NULL,
                              headers = NULL) {
  df <- data.frame(dim1, dim2, t(as.data.frame(matrix)))
  naIx <- is.na(dim1) | is.na(dim2)
  df <- df[!naIx, ]

  m <- reshape2::melt(df, id.vars = c("dim1", "dim2"))
  colnames(m) <- c(xlab, ylab, "facet", "Expression")

  if (is.null(midpoint)) {
    midpoint <- mean(m[, 4], trim = 0.1)
  }

  varLabel <- gsub("_", " ", varLabel)

  if (isFALSE(is.null(headers))) {
    names(headers) <- levels(m$facet)
    headers <- ggplot2::as_labeller(headers)

    g <- ggplot2::ggplot(
      m,
      ggplot2::aes_string(x = xlab, y = ylab)
    ) +
      ggplot2::geom_point(
        stat = "identity",
        size = size,
        ggplot2::aes_string(color = m$Expression)
      ) +
      ggplot2::theme_bw() +
      ggplot2::scale_colour_gradient2(
        low = colorLow,
        high = colorHigh,
        mid = colorMid,
        midpoint = midpoint,
        name = varLabel
      ) +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.spacing = unit(0, "lines"),
        panel.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(colour = "black")
      )
    if (isFALSE(is.null(ncol))) {
      g <- g + ggplot2::facet_wrap(~facet,
        labeller = headers,
        ncol = ncol
      )
    } else {
      g <- g + ggplot2::facet_wrap(~facet, labeller = headers)
    }
  } else {
    g <- ggplot2::ggplot(
      m,
      ggplot2::aes_string(x = xlab, y = ylab)
    ) +
      ggplot2::geom_point(
        stat = "identity",
        size = size,
        ggplot2::aes_string(color = m$Expression)
      ) +
      ggplot2::facet_wrap(~facet) +
      ggplot2::theme_bw() +
      ggplot2::scale_colour_gradient2(
        low = colorLow,
        high = colorHigh,
        mid = colorMid,
        midpoint = midpoint,
        name = varLabel
      ) +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.spacing = unit(0, "lines"),
        panel.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(colour = "black")
      )
    if (isFALSE(is.null(ncol))) {
      g <- g + ggplot2::facet_wrap(~facet, ncol = ncol)
    } else {
      g <- g + ggplot2::facet_wrap(~facet)
    }
  }
  return(g)
}


#' @title Plotting feature expression on a dimensionality reduction plot
#' @description Create a scatterplot for each row of a normalized gene
#'  expression matrix where x and y axis are from a data dimensionality
#'  reduction tool. The cells are colored by expression of
#'  the specified feature.
#' @param dim1 Numeric vector. First dimension from data
#'  dimensionality reduction output.
#' @param dim2 Numeric vector. Second dimension from data dimensionality
#'  reduction output.
#' @param counts Integer matrix. Rows represent features and columns
#'  represent cells.
#' @param features Character vector. Features in the rownames of counts to plot.
#' @param headers Character vector. If `NULL`, the corresponding rownames are
#'  used as labels. Otherwise, these headers are used to label the features.
#' @param normalize Logical. Whether to normalize the columns of `counts`.
#'  Default \code{FALSE}.
#' @param zscore Logical. Whether to scale each feature to have a mean 0
#' and standard deviation of 1. Default \code{TRUE}.
#' @param exactMatch Logical. Whether an exact match or a partial match using
#'  \code{grep()} is used to look up the feature in the rownames of the counts
#'   matrix. Default TRUE.
#' @param trim Numeric vector. Vector of length two that specifies the lower
#'  and upper bounds for the data. This threshold is applied after row scaling.
#'  Set to NULL to disable. Default \code{c(-1,1)}.
#' @param size Numeric. Sets size of point on plot. Default 1.
#' @param xlab Character vector. Label for the x-axis. Default "Dimension_1".
#' @param ylab Character vector. Label for the y-axis. Default "Dimension_2".
#' @param colorLow Character. A color available from `colors()`. The color
#'  will be used to signify the lowest values on the scale. Default 'blue'.
#' @param colorMid Character. A color available from `colors()`. The color
#'  will be used to signify the midpoint on the scale. Default 'white'.
#' @param colorHigh Character. A color available from `colors()`. The color
#'  will be used to signify the highest values on the scale. Default 'red'.
#' @param midpoint Numeric. The value indicating the midpoint of the
#' diverging color scheme. If \code{NULL}, defaults to the mean
#' with 10 percent of values trimmed. Default \code{NULL}.
#' @param ncol Integer. Passed to \link[ggplot2]{facet_wrap}. Specify the
#'  number of columns for facet wrap.
#' @return The plot as a ggplot object
#' @examples
#' \donttest{
#' data(celdaCGSim, celdaCGMod)
#' celdaTsne <- celdaTsne(
#'   counts = celdaCGSim$counts,
#'   celdaMod = celdaCGMod
#' )
#' plotDimReduceFeature(
#'   dim1 = celdaTsne[, 1],
#'   dim2 = celdaTsne[, 2],
#'   counts = celdaCGSim$counts,
#'   normalize = TRUE,
#'   features = c("Gene_99"),
#'   exactMatch = TRUE
#' )
#' }
#' @export
plotDimReduceFeature <- function(dim1,
                                 dim2,
                                 counts,
                                 features,
                                 headers = NULL,
                                 normalize = FALSE,
                                 zscore = TRUE,
                                 exactMatch = TRUE,
                                 trim = c(-1, 1),
                                 size = 1,
                                 xlab = "Dimension_1",
                                 ylab = "Dimension_2",
                                 colorLow = "blue4",
                                 colorMid = "white",
                                 colorHigh = "firebrick1",
                                 midpoint = NULL,
                                 ncol = NULL) {

  # Perform checks
  if (is.null(features)) {
    stop("at least one feature is required to create a plot")
  }

  if (isFALSE(is.null(headers))) {
    if (length(headers) != length(features)) {
      stop(
        "Headers ",
        headers,
        " should be the same length as features ",
        features
      )
    }

    if (isFALSE(exactMatch)) {
      warning("exactMatch is FALSE. headers will not be used!")
      headers <- NULL
    }
  }

  ## Normalize data if needed
  if (isTRUE(normalize)) {
    counts <- normalizeCounts(counts, transformationFun = sqrt)
  }

  # After normalization, features can be selected
  featuresIx <- retrieveFeatureIndex(features,
    counts,
    by = "rownames",
    exactMatch = exactMatch
  )
  counts <- as.matrix(counts[featuresIx, , drop = FALSE])

  # Scale/zscore data if needed
  varLabel <- "Expression"
  if (isTRUE(zscore)) {
    counts <- t(scale(t(counts)))
    varLabel <- "Scaled\nExpression"
  }

  if (!is.null(trim)) {
    if (length(trim) != 2) {
      stop(
        "'trim' should be a 2 element vector",
        "specifying the lower and upper boundaries"
      )
    }
    trim <- sort(trim)
    counts[counts < trim[1]] <- trim[1]
    counts[counts > trim[2]] <- trim[2]
  }


  plotDimReduceGrid(
    dim1 = dim1,
    dim2 = dim2,
    matrix = counts,
    size = size,
    xlab = xlab,
    ylab = ylab,
    colorLow = colorLow,
    colorMid = colorMid,
    colorHigh = colorHigh,
    varLabel = varLabel,
    midpoint = midpoint,
    ncol = ncol,
    headers = headers
  )
}


#' @title Plotting the Celda module probability on a
#'  dimensionality reduction plot
#' @description Create a scatterplot for each row of a normalized
#'  gene expression matrix where x and y axis are from a data
#'  dimensionality reduction tool.
#'  The cells are colored by the module probability(s).
#' @param dim1 Numeric vector.
#'  First dimension from data dimensionality reduction output.
#' @param dim2 Numeric vector.
#'  Second dimension from data dimensionality reduction output.
#' @param counts Integer matrix.
#'  Rows represent features and columns represent cells.
#'  This matrix should be the same as the one used to generate `celdaMod`.
#' @param celdaMod Celda object of class "celda_G" or "celda_CG".
#' @param modules Character vector. Module(s) from celda model to be plotted.
#' e.g. c("1", "2").
#' @param rescale Logical.
#'  Whether rows of the matrix should be rescaled to [0, 1]. Default TRUE.
#' @param size Numeric. Sets size of point on plot. Default 1.
#' @param xlab Character vector. Label for the x-axis. Default "Dimension_1".
#' @param ylab Character vector. Label for the y-axis. Default "Dimension_2".
#' @param colorLow Character. A color available from `colors()`.
#'  The color will be used to signify the lowest values on the scale.
#'  Default 'grey'.
#' @param colorMid Character. A color available from `colors()`.
#'  The color will be used to signify the midpoint on the scale.
#' @param colorHigh Character. A color available from `colors()`.
#'  The color will be used to signify the highest values on the scale.
#'  Default 'blue'.
#' @param ncol Integer. Passed to \link[ggplot2]{facet_wrap}. Specify the
#'  number of columns for facet wrap.
#' @return The plot as a ggplot object
#' @examples
#' \donttest{
#' data(celdaCGSim, celdaCGMod)
#' celdaTsne <- celdaTsne(
#'   counts = celdaCGSim$counts,
#'   celdaMod = celdaCGMod
#' )
#' plotDimReduceModule(
#'   dim1 = celdaTsne[, 1], dim2 = celdaTsne[, 2],
#'   counts = celdaCGSim$counts, celdaMod = celdaCGMod,
#'   modules = c("1", "2")
#' )
#' }
#' @export
plotDimReduceModule <-
  function(dim1,
           dim2,
           counts,
           celdaMod,
           modules = NULL,
           rescale = TRUE,
           size = 1,
           xlab = "Dimension_1",
           ylab = "Dimension_2",
           colorLow = "grey",
           colorMid = NULL,
           colorHigh = "blue",
           ncol = NULL) {
    factorized <- factorizeMatrix(
      celdaMod = celdaMod,
      counts = counts
    )
    matrix <- factorized$proportions$cell

    if (rescale == TRUE) {
      for (x in seq(nrow(matrix))) {
        matrix[x, ] <- matrix[x, ] - min(matrix[x, ])
        matrix[x, ] <- matrix[x, ] / max(matrix[x, ])
        varLabel <- "Scaled_Probability"
      }
    } else {
      varLabel <- "Probability"
    }

    rownames(matrix) <- gsub("L", "", rownames(matrix))
    if (!is.null(modules)) {
      if (length(rownames(matrix)[rownames(matrix) %in% modules]) < 1) {
        stop("All modules selected do not exist in the model.")
      }
      matrix <- matrix[which(rownames(matrix) %in% modules), ,
        drop = FALSE
      ]
      matrix <- matrix[match(rownames(matrix), modules), ,
        drop = FALSE
      ]
    }

    rownames(matrix) <- paste0("L", rownames(matrix))
    plotDimReduceGrid(
      dim1 = dim1,
      dim2 = dim2,
      matrix = matrix,
      size = size,
      xlab = xlab,
      ylab = ylab,
      colorLow = colorLow,
      colorMid = colorMid,
      colorHigh = colorHigh,
      varLabel = varLabel,
      ncol = ncol
    )
  }


# Labeling code adapted from Seurat (https://github.com/satijalab/seurat)
#' @title Plotting the cell labels on a dimensionality reduction plot
#' @description Create a scatterplot for each row of a normalized
#'  gene expression matrix where x and y axis are from a
#'  data dimensionality reduction tool.
#'  The cells are colored by its given `cluster` label.
#' @param dim1 Numeric vector. First dimension from data
#'  dimensionality reduction output.
#' @param dim2 Numeric vector. Second dimension from data
#'  dimensionality reduction output.
#' @param cluster Integer vector. Contains cluster labels for each cell.
#' @param size Numeric. Sets size of point on plot. Default 1.
#' @param xlab Character vector. Label for the x-axis. Default "Dimension_1".
#' @param ylab Character vector. Label for the y-axis. Default "Dimension_2".
#' @param specificClusters Numeric vector.
#'  Only color cells in the specified clusters.
#'  All other cells will be grey.
#'  If NULL, all clusters will be colored. Default NULL.
#' @param labelClusters Logical. Whether the cluster labels are plotted.
#'  Default FALSE.
#' @param groupBy Character vector. Contains sample labels for each cell.
#'  If NULL, all samples will be plotted together. Default NULL.
#' @param labelSize Numeric. Sets size of label if labelClusters is TRUE.
#'  Default 3.5.
#' @return The plot as a ggplot object
#' @importFrom ggrepel geom_text_repel
#' @examples
#' \donttest{
#' data(celdaCGSim, celdaCGMod)
#' celdaTsne <- celdaTsne(
#'   counts = celdaCGSim$counts,
#'   celdaMod = celdaCGMod
#' )
#' plotDimReduceCluster(
#'   dim1 = celdaTsne[, 1],
#'   dim2 = celdaTsne[, 2],
#'   cluster = as.factor(celdaCGMod@clusters$z),
#'   specificClusters = c(1, 2, 3)
#' )
#' }
#' @export
plotDimReduceCluster <- function(dim1,
                                 dim2,
                                 cluster,
                                 size = 1,
                                 xlab = "Dimension_1",
                                 ylab = "Dimension_2",
                                 specificClusters = NULL,
                                 labelClusters = FALSE,
                                 groupBy = NULL,
                                 labelSize = 3.5) {
  if (!is.null(groupBy)) {
    df <- data.frame(dim1, dim2, cluster, groupBy)
    colnames(df) <- c(xlab, ylab, "Cluster", "Sample")
  } else {
    df <- data.frame(dim1, dim2, cluster)
    colnames(df) <- c(xlab, ylab, "Cluster")
  }

  naIx <- is.na(dim1) | is.na(dim2)
  df <- df[!naIx, ]
  df[3] <- as.factor(df[[3]])
  clusterColors <- distinctColors(nlevels(as.factor(cluster)))

  if (!is.null(specificClusters)) {
    clusterColors[!levels(df[[3]]) %in% specificClusters] <- "gray92"
  }

  g <- ggplot2::ggplot(df, ggplot2::aes_string(x = xlab, y = ylab)) +
    ggplot2::geom_point(
      stat = "identity",
      size = size,
      ggplot2::aes_string(color = "Cluster")
    ) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black")
    ) +
    ggplot2::scale_color_manual(values = clusterColors) +
    ggplot2::guides(
      color =
        ggplot2::guide_legend(override.aes = list(size = 1))
    )

  if (isTRUE(labelClusters)) {
    # centroidList <- lapply(seq(length(unique(cluster))), function(x) {
    centroidList <- lapply(unique(cluster), function(x) {
      df.sub <- df[df$Cluster == x, ]
      median.1 <- stats::median(df.sub[, xlab])
      median.2 <- stats::median(df.sub[, ylab])
      cbind(median.1, median.2, x)
    })
    centroid <- do.call(rbind, centroidList)
    centroid <- data.frame(
      Dimension_1 = as.numeric(centroid[, 1]),
      Dimension_2 = as.numeric(centroid[, 2]),
      Cluster = centroid[, 3]
    )

    colnames(centroid)[seq(2)] <- c(xlab, ylab)
    g <- g + ggplot2::geom_point(
      data = centroid,
      mapping = ggplot2::aes_string(
        x = xlab,
        y = ylab
      ),
      size = 0,
      alpha = 0
    ) +
      ggrepel::geom_text_repel(
        data = centroid,
        mapping = ggplot2::aes(label = centroid$Cluster),
        size = labelSize
      )
  }
  if (!is.null(x = groupBy)) {
    g <- g + ggplot2::facet_wrap(
      facets = ggplot2::vars(!!ggplot2::sym(x = "Sample"))) +
      ggplot2::theme(strip.background = ggplot2::element_blank())
  }
  return(g)
}


# Run the t-SNE algorithm for dimensionality reduction
# @param norm Normalized count matrix.
# @param perplexity Numeric vector. Determines perplexity for tsne. Default 20.
# @param maxIter Numeric vector. Determines iterations for tsne. Default 1000.
# @param doPca Logical. Whether to perform
# dimensionality reduction with PCA before tSNE.
# @param initialDims Integer. Number of dimensions from PCA to use as
# input in tSNE. Default 50.
#' @importFrom Rtsne Rtsne
.calculateTsne <- function(norm,
    perplexity,
    maxIter,
    doPca,
    initialDims) {

    res <- Rtsne::Rtsne(
        norm,
        pca = doPca,
        max_iter = maxIter,
        perplexity = perplexity,
        check_duplicates = FALSE,
        is_distance = FALSE,
        initial_dims = initialDims)$Y

    return(res)
}


# Run the UMAP algorithm for dimensionality reduction
# @param norm Normalized count matrix.
# @param nNeighbors The size of local neighborhood used for
#   manifold approximation. Larger values result in more global
#   views of the manifold, while smaller values result in more
#   local data being preserved. Default 30.
#    See `?uwot::umap` for more information.
# @param minDist The effective minimum distance between embedded points.
#    Smaller values will result in a more clustered/clumped
#    embedding where nearby points on the manifold are drawn
#    closer together, while larger values will result on a more
#    even dispersal of points. Default 0.2.
#    See `?uwot::umap` for more information.
# @param spread The effective scale of embedded points. In combination with
#    ‘min_dist’, this determines how clustered/clumped the
#    embedded points are. Default 1.
#    See `?uwot::umap` for more information.
# @param pca Logical. Whether to perform
# dimensionality reduction with PCA before UMAP.
# @param initialDims Integer. Number of dimensions from PCA to use as
# input in UMAP. Default 50.
# @param cores Number of threads to use. Default 1.
# @param ... Other parameters to pass to `uwot::umap`.
#' @import uwot
.calculateUmap <- function(norm,
                           nNeighbors = 30,
                           minDist = 0.75,
                           spread = 1,
                           pca = FALSE,
                           initialDims = 50,
                           cores = 1,
                           ...) {
  if (isTRUE(pca)) {
    doPCA <- initialDims
  } else {
    doPCA <- NULL
  }

  res <- uwot::umap(norm,
    n_neighbors = nNeighbors,
    min_dist = minDist, spread = spread,
    n_threads = cores, n_sgd_threads = 1, pca = doPCA, ...
  )
  return(res)
}



#' @title Feature Expression Violin Plot
#' @description Outputs a violin plot for feature expression data.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells.
#' @param celdaMod Celda object of class "celda_G" or "celda_CG".
#' @param features Character vector. Uses these genes for plotting.
#' @param exactMatch Logical. Whether an exact match or a partial match using
#'  \code{grep()} is used to look up the feature in the rownames of the counts
#'   matrix. Default \code{TRUE}.
#' @param plotDots Boolean. If \code{TRUE}, the
#'  expression of features will be plotted as points in addition to the violin
#'  curve. Default \code{TRUE}.
#' @param dotSize Numeric. Size of points if \code{plotDots = TRUE}.
#' Default \code{0.1}.
#' @return Violin plot for each feature, grouped by celda cluster
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' plotCeldaViolin(
#'   counts = celdaCGSim$counts,
#'   celdaMod = celdaCGMod, features = "Gene_1"
#' )
#' @export
plotCeldaViolin <- function(counts,
                            celdaMod,
                            features,
                            exactMatch = TRUE,
                            plotDots = TRUE,
                            dotSize = 0.1) {
  cluster <- celdaMod@clusters$z

  featuresIx <- retrieveFeatureIndex(features,
    counts,
    by = "rownames",
    exactMatch = exactMatch
  )
  dataFeature <- as.matrix(counts[featuresIx, , drop = FALSE])
  dataFeature <- as.data.frame(t(dataFeature))
  df <- cbind(cluster, dataFeature)
  df$cluster <- as.factor(df$cluster)
  m <- reshape2::melt(df, id.vars = c("cluster"))
  colnames(m) <- c("Cluster", "Feature", "Expression")
  colorPal <- distinctColors(length(unique(cluster)))

  p <- ggplot2::ggplot(
    m,
    ggplot2::aes_string(
      x = "Cluster",
      y = "Expression",
      fill = "Cluster"
    )
  ) +
    ggplot2::facet_wrap(~Feature) +
    ggplot2::geom_violin(trim = TRUE, scale = "width") +
    ggplot2::scale_fill_manual(values = colorPal) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing = grid::unit(0, "lines"),
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black")
    )

  if (isTRUE(plotDots)) {
    p <- p + ggplot2::geom_jitter(height = 0, size = dotSize)
  }

  return(p)
}
