#' @title Mapping the dimension reduction plot
#' @description Creates a scatterplot given two dimensions from a data
#'  dimension reduction tool (e.g tSNE) output.
#' @param x Numeric matrix or a \linkS4class{SingleCellExperiment} object
#'  with the matrix located in the assay slot under \code{useAssay}. Each
#'  row of the matrix will be plotted as a separate facet.
#' @param reducedDimName The name of the dimension reduction slot in
#'  \code{reducedDimNames(x)} if \code{x} is a
#'  \linkS4class{SingleCellExperiment} object. Ignored if both \code{dim1} and
#'  \code{dim2} are set.
#' @param dim1 Numeric vector. First dimension from data dimension
#'  reduction output.
#' @param dim2 Numeric vector. Second dimension from data dimension
#'  reduction output.
#' @param useAssay A string specifying which \link{assay}
#'  slot to use if \code{x} is a
#'  \linkS4class{SingleCellExperiment} object. Default "counts".
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default "featureSubset".
#' @param size Numeric. Sets size of point on plot. Default 1.
#' @param xlab Character vector. Label for the x-axis. Default 'Dimension_1'.
#' @param ylab Character vector. Label for the y-axis. Default 'Dimension_2'.
#' @param limits Passed to \link{scale_colour_gradient2}. The range
#'  of color scale.
#' @param colorLow Character. A color available from `colors()`.
#'  The color will be used to signify the lowest values on the scale.
#'  Default "blue4".
#' @param colorMid Character. A color available from `colors()`.
#'  The color will be used to signify the midpoint on the scale. Default
#'  "grey90".
#' @param colorHigh Character. A color available from `colors()`.
#'  The color will be used to signify the highest values on the scale.
#'  Default "firebrick1".
#' @param midpoint Numeric. The value indicating the midpoint of the
#' diverging color scheme. If \code{NULL}, defaults to the mean
#' with 10 percent of values trimmed. Default \code{0}.
#' @param varLabel Character vector. Title for the color legend.
#' @param ncol Integer. Passed to \link[ggplot2]{facet_wrap}. Specify the
#'  number of columns for facet wrap.
#' @param headers Character vector. If `NULL`, the corresponding rownames are
#'  used as labels. Otherwise, these headers are used to label the genes.
#' @param ... Ignored. Placeholder to prevent check warning.
#' @param decreasing logical. Specifies the order of plotting the points.
#'  If \code{FALSE}, the points will be plotted in increasing order where
#'  the points with largest values will be on top. \code{TRUE} otherwise.
#'  If \code{NULL}, no sorting is performed. Points will be plotted in their
#'  current order in \code{x}. Default \code{FALSE}.
#' @return The plot as a ggplot object
#' @export
setGeneric("plotDimReduceGrid", function(x, ...) {
    standardGeneric("plotDimReduceGrid")})


#' @rdname plotDimReduceGrid
#' @examples
#' data(sceCeldaCG)
#' sce <- celdaTsne(sceCeldaCG)
#' plotDimReduceGrid(x = sce,
#'   reducedDimName = "celda_tSNE",
#'   xlab = "Dimension1",
#'   ylab = "Dimension2",
#'   varLabel = "tSNE")
#' @export
setMethod("plotDimReduceGrid",
    signature(x = "SingleCellExperiment"),
    function(x,
        reducedDimName,
        dim1 = NULL,
        dim2 = NULL,
        useAssay = "counts",
        altExpName = "featureSubset",
        size = 1,
        xlab = "Dimension_1",
        ylab = "Dimension_2",
        limits = c(-2, 2),
        colorLow = "blue4",
        colorMid = "grey90",
        colorHigh = "firebrick1",
        midpoint = 0,
        varLabel = NULL,
        ncol = NULL,
        headers = NULL,
        decreasing = FALSE) {

        altExp <- SingleCellExperiment::altExp(x, altExpName)
        matrix <- SummarizedExperiment::assay(x, i = useAssay)

        if (is.null(dim1)) {
            dim1 <- SingleCellExperiment::reducedDim(altExp,
                reducedDimName)[, 1]
        }

        if (is.null(dim2)) {
            dim2 <- SingleCellExperiment::reducedDim(altExp,
                reducedDimName)[, 2]
        }

        g <- .plotDimReduceGrid(dim1 = dim1,
            dim2 = dim2,
            matrix = matrix,
            size = size,
            xlab = xlab,
            ylab = ylab,
            limits = limits,
            colorLow = colorLow,
            colorMid = colorMid,
            colorHigh = colorHigh,
            midpoint = midpoint,
            varLabel = varLabel,
            ncol = ncol,
            headers = headers,
            decreasing = decreasing)
        return(g)
    }
)


#' @rdname plotDimReduceGrid
#' @examples
#' library(SingleCellExperiment)
#' data(sceCeldaCG)
#' sce <- celdaTsne(sceCeldaCG)
#' plotDimReduceGrid(x = counts(sce),
#'   dim1 = reducedDim(altExp(sce), "celda_tSNE")[, 1],
#'   dim2 = reducedDim(altExp(sce), "celda_tSNE")[, 2],
#'   xlab = "Dimension1",
#'   ylab = "Dimension2",
#'   varLabel = "tSNE")
#' @export
setMethod("plotDimReduceGrid",
    signature(x = "matrix"),
    function(x,
        dim1,
        dim2,
        size = 1,
        xlab = "Dimension_1",
        ylab = "Dimension_2",
        limits = c(-2, 2),
        colorLow = "blue4",
        colorMid = "grey90",
        colorHigh = "firebrick1",
        midpoint = NULL,
        varLabel = NULL,
        ncol = NULL,
        headers = NULL,
        decreasing = FALSE) {

        g <- .plotDimReduceGrid(dim1 = dim1,
            dim2 = dim2,
            matrix = x,
            size = size,
            xlab = xlab,
            ylab = ylab,
            limits = limits,
            colorLow = colorLow,
            colorMid = colorMid,
            colorHigh = colorHigh,
            midpoint = midpoint,
            varLabel = varLabel,
            ncol = ncol,
            headers = headers,
            decreasing = decreasing)
        return(g)
    }
)


#' @importFrom reshape2 melt
.plotDimReduceGrid <- function(dim1,
    dim2,
    matrix,
    size,
    xlab,
    ylab,
    limits,
    colorLow,
    colorMid,
    colorHigh,
    midpoint,
    varLabel,
    ncol,
    headers,
    decreasing) {

    df <- data.frame(dim1, dim2, t(as.data.frame(matrix)), check.names = FALSE)
    naIx <- is.na(dim1) | is.na(dim2)
    df <- df[!naIx, ]

    m <- reshape2::melt(df, id.vars = c("dim1", "dim2"))
    colnames(m) <- c(xlab, ylab, "facet", "Expression")

    if (!is.null(decreasing)) {
        m <- m[order(m$facet, m$Expression, decreasing = decreasing), ]
    }

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
                limits = limits,
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
                limits = limits,
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


#' @title Plotting feature expression on a dimension reduction plot
#' @description Create a scatterplot for each row of a normalized gene
#'  expression matrix where x and y axis are from a data dimension
#'  reduction tool. The cells are colored by expression of
#'  the specified feature.
#' @param x Numeric matrix or a \linkS4class{SingleCellExperiment} object
#'  with the matrix located in the assay slot under \code{useAssay}. Rows
#'  represent features and columns represent cells.
#' @param reducedDimName The name of the dimension reduction slot in
#'  \code{reducedDimNames(x)} if \code{x} is a
#'  \linkS4class{SingleCellExperiment} object. Ignored if both \code{dim1} and
#'  \code{dim2} are set.
#' @param dim1 Numeric vector. First dimension from data
#'  dimension reduction output.
#' @param dim2 Numeric vector. Second dimension from data dimension
#'  reduction output.
#' @param useAssay A string specifying which \link{assay}
#'  slot to use if \code{x} is a
#'  \linkS4class{SingleCellExperiment} object. Default "counts".
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default "featureSubset".
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
#' @param limits Passed to \link{scale_colour_gradient2}. The range
#'  of color scale.
#' @param size Numeric. Sets size of point on plot. Default 1.
#' @param xlab Character vector. Label for the x-axis. Default "Dimension_1".
#' @param ylab Character vector. Label for the y-axis. Default "Dimension_2".
#' @param colorLow Character. A color available from `colors()`. The color
#'  will be used to signify the lowest values on the scale.
#' @param colorMid Character. A color available from `colors()`. The color
#'  will be used to signify the midpoint on the scale.
#' @param colorHigh Character. A color available from `colors()`. The color
#'  will be used to signify the highest values on the scale.
#' @param midpoint Numeric. The value indicating the midpoint of the
#' diverging color scheme. If \code{NULL}, defaults to the mean
#' with 10 percent of values trimmed. Default \code{0}.
#' @param ncol Integer. Passed to \link[ggplot2]{facet_wrap}. Specify the
#'  number of columns for facet wrap.
#' @param decreasing logical. Specifies the order of plotting the points.
#'  If \code{FALSE}, the points will be plotted in increasing order where
#'  the points with largest values will be on top. \code{TRUE} otherwise.
#'  If \code{NULL}, no sorting is performed. Points will be plotted in their
#'  current order in \code{x}. Default \code{FALSE}.
#' @param ... Ignored. Placeholder to prevent check warning.
#' @return The plot as a ggplot object
#' @export
setGeneric("plotDimReduceFeature", function(x, ...) {
    standardGeneric("plotDimReduceFeature")})


#' @rdname plotDimReduceFeature
#' @examples
#' data(sceCeldaCG)
#' sce <- celdaTsne(sceCeldaCG)
#' plotDimReduceFeature(x = sce,
#'   reducedDimName = "celda_tSNE",
#'   normalize = TRUE,
#'   features = c("Gene_99"),
#'   exactMatch = TRUE)
#' @export
setMethod("plotDimReduceFeature",
    signature(x = "SingleCellExperiment"),
    function(x,
        reducedDimName,
        dim1 = NULL,
        dim2 = NULL,
        useAssay = "counts",
        altExpName = "featureSubset",
        features,
        headers = NULL,
        normalize = FALSE,
        zscore = TRUE,
        exactMatch = TRUE,
        trim = c(-2, 2),
        limits = c(-2, 2),
        size = 1,
        xlab = "Dimension_1",
        ylab = "Dimension_2",
        colorLow = "blue4",
        colorMid = "grey90",
        colorHigh = "firebrick1",
        midpoint = 0,
        ncol = NULL,
        decreasing = FALSE) {

        altExp <- SingleCellExperiment::altExp(x, altExpName)
        counts <- SummarizedExperiment::assay(x, i = useAssay)

        if (is.null(dim1)) {
            dim1 <- SingleCellExperiment::reducedDim(altExp,
                reducedDimName)[, 1]
        }

        if (is.null(dim2)) {
            dim2 <- SingleCellExperiment::reducedDim(altExp,
                reducedDimName)[, 2]
        }

        g <- .plotDimReduceFeature(dim1 = dim1,
            dim2 = dim2,
            counts = counts,
            features = features,
            headers = headers,
            normalize = normalize,
            zscore = zscore,
            exactMatch = exactMatch,
            trim = trim,
            limits = limits,
            size = size,
            xlab = xlab,
            ylab = ylab,
            colorLow = colorLow,
            colorMid = colorMid,
            colorHigh = colorHigh,
            midpoint = midpoint,
            ncol = ncol,
            decreasing = decreasing)
        return(g)
    }
)


#' @rdname plotDimReduceFeature
#' @examples
#' library(SingleCellExperiment)
#' data(sceCeldaCG)
#' sce <- celdaTsne(sceCeldaCG)
#' plotDimReduceFeature(x = counts(sce),
#'   dim1 = reducedDim(altExp(sce), "celda_tSNE")[, 1],
#'   dim2 = reducedDim(altExp(sce), "celda_tSNE")[, 2],
#'   normalize = TRUE,
#'   features = c("Gene_99"),
#'   exactMatch = TRUE)
#' @export
setMethod("plotDimReduceFeature",
    signature(x = "matrix"),
    function(x,
        dim1,
        dim2,
        features,
        headers = NULL,
        normalize = FALSE,
        zscore = TRUE,
        exactMatch = TRUE,
        trim = c(-2, 2),
        limits = c(-2, 2),
        size = 1,
        xlab = "Dimension_1",
        ylab = "Dimension_2",
        colorLow = "blue4",
        colorMid = "grey90",
        colorHigh = "firebrick1",
        midpoint = 0,
        ncol = NULL,
        decreasing = FALSE) {

        g <- .plotDimReduceFeature(dim1 = dim1,
            dim2 = dim2,
            counts = x,
            features = features,
            headers = headers,
            normalize = normalize,
            zscore = zscore,
            exactMatch = exactMatch,
            trim = trim,
            limits = limits,
            size = size,
            xlab = xlab,
            ylab = ylab,
            colorLow = colorLow,
            colorMid = colorMid,
            colorHigh = colorHigh,
            midpoint = midpoint,
            ncol = ncol,
            decreasing = decreasing)
        return(g)
    }
)


.plotDimReduceFeature <- function(dim1,
                                 dim2,
                                 counts,
                                 features,
                                 headers,
                                 normalize,
                                 zscore,
                                 exactMatch,
                                 trim,
                                 limits,
                                 size,
                                 xlab,
                                 ylab,
                                 colorLow,
                                 colorMid,
                                 colorHigh,
                                 midpoint,
                                 ncol,
                                 decreasing) {

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
    exactMatch = exactMatch)
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

  .plotDimReduceGrid(
    dim1 = dim1,
    dim2 = dim2,
    matrix = counts,
    size = size,
    xlab = xlab,
    ylab = ylab,
    limits = limits,
    colorLow = colorLow,
    colorMid = colorMid,
    colorHigh = colorHigh,
    varLabel = varLabel,
    midpoint = midpoint,
    ncol = ncol,
    headers = headers,
    decreasing = decreasing
  )
}


#' @title Plotting Celda module probability on a
#'  dimension reduction plot
#' @description Create a scatterplot for each row of a normalized
#'  gene expression matrix where x and y axis are from a data
#'  dimension reduction tool.
#'  The cells are colored by the module probability.
#' @param x Numeric matrix or a \linkS4class{SingleCellExperiment} object
#'  with the matrix located in the assay slot under \code{useAssay}. Rows
#'  represent features and columns represent cells.
#' @param reducedDimName The name of the dimension reduction slot in
#'  \code{reducedDimNames(x)} if \code{x} is a
#'  \linkS4class{SingleCellExperiment} object. Ignored if both \code{dim1} and
#'  \code{dim2} are set.
#' @param dim1 Numeric vector.
#'  First dimension from data dimension reduction output.
#' @param dim2 Numeric vector.
#'  Second dimension from data dimension reduction output.
#' @param useAssay A string specifying which \link{assay}
#'  slot to use if \code{x} is a
#'  \linkS4class{SingleCellExperiment} object. Default "counts".
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default "featureSubset".
#' @param celdaMod Celda object of class "celda_G" or "celda_CG". Used only if
#'  \code{x} is a matrix object.
#' @param modules Character vector. Module(s) from celda model to be plotted.
#'  e.g. c("1", "2").
#' @param rescale Logical.
#'  Whether rows of the matrix should be rescaled to [0, 1]. Default TRUE.
#' @param limits Passed to \link{scale_colour_gradient}. The range
#'  of color scale.
#' @param size Numeric. Sets size of point on plot. Default 1.
#' @param xlab Character vector. Label for the x-axis. Default "Dimension_1".
#' @param ylab Character vector. Label for the y-axis. Default "Dimension_2".
#' @param colorLow Character. A color available from `colors()`.
#'  The color will be used to signify the lowest values on the scale.
#' @param colorHigh Character. A color available from `colors()`.
#'  The color will be used to signify the highest values on the scale.
#' @param ncol Integer. Passed to \link[ggplot2]{facet_wrap}. Specify the
#'  number of columns for facet wrap.
#' @param decreasing logical. Specifies the order of plotting the points.
#'  If \code{FALSE}, the points will be plotted in increasing order where
#'  the points with largest values will be on top. \code{TRUE} otherwise.
#'  If \code{NULL}, no sorting is performed. Points will be plotted in their
#'  current order in \code{x}. Default \code{FALSE}.
#' @param ... Ignored. Placeholder to prevent check warning.
#' @return The plot as a ggplot object
#' @export
setGeneric("plotDimReduceModule", function(x, ...) {
    standardGeneric("plotDimReduceModule")})


#' @rdname plotDimReduceModule
#' @examples
#' data(sceCeldaCG)
#' sce <- celdaTsne(sceCeldaCG)
#' plotDimReduceModule(x = sce,
#'   reducedDimName = "celda_tSNE",
#'   modules = c("1", "2"))
#' @export
setMethod("plotDimReduceModule",
    signature(x = "SingleCellExperiment"),
    function(x,
        reducedDimName,
        dim1 = NULL,
        dim2 = NULL,
        useAssay = "counts",
        altExpName = "featureSubset",
        modules = NULL,
        rescale = TRUE,
        limits = c(0, 1),
        size = 1,
        xlab = "Dimension_1",
        ylab = "Dimension_2",
        colorLow = "grey90",
        colorHigh = "firebrick1",
        ncol = NULL,
        decreasing = FALSE) {

        altExp <- SingleCellExperiment::altExp(x, altExpName)

        if (is.null(dim1)) {
            dim1 <- SingleCellExperiment::reducedDim(altExp,
                reducedDimName)[, 1]
        }

        if (is.null(dim2)) {
            dim2 <- SingleCellExperiment::reducedDim(altExp,
                reducedDimName)[, 2]
        }

        counts <- SummarizedExperiment::assay(x, i = useAssay)
        factorized <- factorizeMatrix(x,
            useAssay = useAssay,
            altExpName = altExpName,
            type = "proportion")

        g <- .plotDimReduceModule(dim1 = dim1,
            dim2 = dim2,
            counts = counts,
            factorized = factorized,
            modules = modules,
            rescale = rescale,
            limits = limits,
            size = size,
            xlab = xlab,
            ylab = ylab,
            colorLow = colorLow,
            colorHigh = colorHigh,
            ncol = ncol,
            decreasing = decreasing)
        return(g)
    }
)


#' @rdname plotDimReduceModule
#' @examples
#' library(SingleCellExperiment)
#' data(sceCeldaCG, celdaCGMod)
#' sce <- celdaTsne(sceCeldaCG)
#' plotDimReduceModule(x = counts(sce),
#'   dim1 = reducedDim(altExp(sce), "celda_tSNE")[, 1],
#'   dim2 = reducedDim(altExp(sce), "celda_tSNE")[, 2],
#'   celdaMod = celdaCGMod,
#'   modules = c("1", "2"))
#' @export
setMethod("plotDimReduceModule",
    signature(x = "matrix"),
    function(x,
        dim1,
        dim2,
        celdaMod,
        modules = NULL,
        rescale = TRUE,
        limits = c(0, 1),
        size = 1,
        xlab = "Dimension_1",
        ylab = "Dimension_2",
        colorLow = "blue4",
        colorHigh = "firebrick1",
        ncol = NULL,
        decreasing = FALSE) {

        factorized <- factorizeMatrix(x = x, celdaMod = celdaMod)
        g <- .plotDimReduceModule(dim1 = dim1,
            dim2 = dim2,
            counts = x,
            factorized = factorized,
            modules = modules,
            rescale = rescale,
            limits = limits,
            size = size,
            xlab = xlab,
            ylab = ylab,
            colorLow = colorLow,
            colorHigh = colorHigh,
            ncol = ncol,
            decreasing = decreasing)
        return(g)
    }
)


.plotDimReduceModule <- function(dim1,
    dim2,
    counts,
    factorized,
    modules,
    rescale,
    limits,
    size,
    xlab,
    ylab,
    colorLow,
    colorHigh,
    ncol,
    decreasing) {

    matrix <- factorized$proportions$cell

    if (rescale == TRUE) {
        for (x in seq(nrow(matrix))) {
            matrix[x, ] <- matrix[x, ] - min(matrix[x, ])
            matrix[x, ] <- matrix[x, ] / max(matrix[x, ])
            varLabel <- "Scaled Probability"
        }
    } else {
        varLabel <- "Probability"
    }

    rownames(matrix) <- gsub("L", "", rownames(matrix))
    if (!is.null(modules)) {
        if (length(rownames(matrix)[rownames(matrix) %in% modules]) < 1) {
            stop("All modules selected do not exist in the model.")
        }
        matrix <- matrix[which(rownames(matrix) %in% modules), , drop = FALSE]
        matrix <- matrix[match(rownames(matrix), modules), , drop = FALSE]
    }

    rownames(matrix) <- paste0("L", rownames(matrix))

    df <- data.frame(dim1, dim2, t(as.data.frame(matrix)), check.names = FALSE)
    naIx <- is.na(dim1) | is.na(dim2)
    df <- df[!naIx, ]

    m <- reshape2::melt(df, id.vars = c("dim1", "dim2"))
    colnames(m) <- c(xlab, ylab, "facet", "Expression")

    if (!is.null(decreasing)) {
        m <- m[order(m$facet, m$Expression, decreasing = decreasing), ]
    }

    g <- ggplot2::ggplot(m, ggplot2::aes_string(x = xlab, y = ylab)) +
        ggplot2::geom_point(stat = "identity",
            size = size,
            ggplot2::aes_string(color = m$Expression)) +
        ggplot2::facet_wrap(~facet) +
        ggplot2::theme_bw() +
        ggplot2::scale_colour_gradient(limits = limits,
            low = colorLow,
            high = colorHigh,
            name = varLabel) +
        ggplot2::theme(strip.background = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.spacing = unit(0, "lines"),
            panel.background = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(colour = "black"))
    if (isFALSE(is.null(ncol))) {
        g <- g + ggplot2::facet_wrap(~facet, ncol = ncol)
    } else {
        g <- g + ggplot2::facet_wrap(~facet)
    }

    return(g)
}


# Labeling code adapted from Seurat (https://github.com/satijalab/seurat)
#' @title Plotting the cell labels on a dimension reduction plot
#' @description Create a scatterplot for each row of a normalized
#'  gene expression matrix where x and y axis are from a
#'  data dimension reduction tool.
#'  The cells are colored by "celda_cell_cluster" column in
#'   \code{colData(altExp(x, altExpName))} if \code{x} is a
#'   \linkS4class{SingleCellExperiment} object, or \code{x} if \code{x} is
#'   a integer vector of cell cluster labels.
#' @param x Integer vector of cell cluster labels or a
#'  \linkS4class{SingleCellExperiment} object
#'  containing cluster labels for each cell in \code{"celda_cell_cluster"}
#'  column in \code{colData(x)}.
#' @param reducedDimName The name of the dimension reduction slot in
#'  \code{reducedDimNames(x)} if \code{x} is a
#'  \linkS4class{SingleCellExperiment} object. Ignored if both \code{dim1} and
#'  \code{dim2} are set.
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default "featureSubset".
#' @param dim1 Numeric vector. First dimension from data
#'  dimension reduction output.
#' @param dim2 Numeric vector. Second dimension from data
#'  dimension reduction output.
#' @param size Numeric. Sets size of point on plot. Default 1.
#' @param xlab Character vector. Label for the x-axis. Default "Dimension_1".
#' @param ylab Character vector. Label for the y-axis. Default "Dimension_2".
#' @param specificClusters Numeric vector.
#'  Only color cells in the specified clusters.
#'  All other cells will be grey.
#'  If NULL, all clusters will be colored. Default \code{NULL}.
#' @param labelClusters Logical. Whether the cluster labels are plotted.
#'  Default FALSE.
#' @param groupBy Character vector. Contains sample labels for each cell.
#'  If NULL, all samples will be plotted together. Default NULL.
#' @param labelSize Numeric. Sets size of label if labelClusters is TRUE.
#'  Default 3.5.
#' @param ... Ignored. Placeholder to prevent check warning.
#' @return The plot as a ggplot object
#' @importFrom ggrepel geom_text_repel
#' @export
setGeneric("plotDimReduceCluster", function(x, ...) {
    standardGeneric("plotDimReduceCluster")})


#' @rdname plotDimReduceCluster
#' @examples
#' data(sceCeldaCG)
#' sce <- celdaTsne(sceCeldaCG)
#' plotDimReduceCluster(x = sce,
#'   reducedDimName = "celda_tSNE",
#'   specificClusters = c(1, 2, 3))
#' @export
setMethod("plotDimReduceCluster",
    signature(x = "SingleCellExperiment"),
    function(x,
        reducedDimName,
        altExpName = "featureSubset",
        dim1 = NULL,
        dim2 = NULL,
        size = 1,
        xlab = "Dimension_1",
        ylab = "Dimension_2",
        specificClusters = NULL,
        labelClusters = FALSE,
        groupBy = NULL,
        labelSize = 3.5) {

        altExp <- SingleCellExperiment::altExp(x, altExpName)

        if (!("celda_cell_cluster" %in%
                colnames(SummarizedExperiment::colData(altExp)))) {
            stop("Must have column 'celda_cell_cluster' in",
                " colData(altExp(x, altExpName))!")
        }
        cluster <- SummarizedExperiment::colData(altExp)[["celda_cell_cluster"]]

        if (is.null(dim1)) {
            dim1 <- SingleCellExperiment::reducedDim(altExp,
                reducedDimName)[, 1]
        }

        if (is.null(dim2)) {
            dim2 <- SingleCellExperiment::reducedDim(altExp,
                reducedDimName)[, 2]
        }

        g <- .plotDimReduceCluster(dim1 = dim1,
            dim2 = dim2,
            cluster = cluster,
            size = size,
            xlab = xlab,
            ylab = ylab,
            specificClusters = specificClusters,
            labelClusters = labelClusters,
            groupBy = groupBy,
            labelSize = labelSize)
        return(g)
    }
)


#' @rdname plotDimReduceCluster
#' @examples
#' library(SingleCellExperiment)
#' data(sceCeldaCG, celdaCGMod)
#' sce <- celdaTsne(sceCeldaCG)
#' plotDimReduceCluster(x = celdaClusters(celdaCGMod)$z,
#'   dim1 = reducedDim(altExp(sce), "celda_tSNE")[, 1],
#'   dim2 = reducedDim(altExp(sce), "celda_tSNE")[, 2],
#'   specificClusters = c(1, 2, 3))
#' @export
setMethod("plotDimReduceCluster",
    signature(x = "vector"),
    function(x,
        dim1,
        dim2,
        size = 1,
        xlab = "Dimension_1",
        ylab = "Dimension_2",
        specificClusters = NULL,
        labelClusters = FALSE,
        groupBy = NULL,
        labelSize = 3.5) {

        g <- .plotDimReduceCluster(dim1 = dim1,
            dim2 = dim2,
            cluster = x,
            size = size,
            xlab = xlab,
            ylab = ylab,
            specificClusters = specificClusters,
            labelClusters = labelClusters,
            groupBy = groupBy,
            labelSize = labelSize)
        return(g)
    }
)


.plotDimReduceCluster <- function(dim1,
                                 dim2,
                                 cluster,
                                 size,
                                 xlab,
                                 ylab,
                                 specificClusters,
                                 labelClusters,
                                 groupBy,
                                 labelSize) {
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
      median1 <- stats::median(df.sub[, xlab])
      median2 <- stats::median(df.sub[, ylab])
      data.frame(median1 = median1, median2 = median2, x = x)
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
        mapping = ggplot2::aes(label = Cluster),
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


#' @title Feature Expression Violin Plot
#' @description Outputs a violin plot for feature expression data.
#' @param x Numeric matrix or a \linkS4class{SingleCellExperiment} object
#'  with the matrix located in the assay slot under \code{useAssay}. Rows
#'  represent features and columns represent cells.
#' @param features Character vector. Uses these genes for plotting.
#' @param useAssay A string specifying which \link{assay}
#'  slot to use if \code{x} is a
#'  \linkS4class{SingleCellExperiment} object. Default "counts".
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default "featureSubset".
#' @param celdaMod Celda object of class "celda_G" or "celda_CG". Used only if
#'  \code{x} is a matrix object.
#' @param exactMatch Logical. Whether an exact match or a partial match using
#'  \code{grep()} is used to look up the feature in the rownames of the counts
#'   matrix. Default \code{TRUE}.
#' @param plotDots Boolean. If \code{TRUE}, the
#'  expression of features will be plotted as points in addition to the violin
#'  curve. Default \code{TRUE}.
#' @param dotSize Numeric. Size of points if \code{plotDots = TRUE}.
#' Default \code{0.1}.
#' @param ... Ignored. Placeholder to prevent check warning.
#' @return Violin plot for each feature, grouped by celda cluster
#' @export
setGeneric("plotCeldaViolin", function(x, ...) {
    standardGeneric("plotCeldaViolin")})


#' @rdname plotCeldaViolin
#' @examples
#' data(sceCeldaCG)
#' plotCeldaViolin(x = sceCeldaCG, features = "Gene_1")
#' @export
setMethod("plotCeldaViolin",
    signature(x = "SingleCellExperiment"),
    function(x,
        features,
        useAssay = "counts",
        altExpName = "featureSubset",
        exactMatch = TRUE,
        plotDots = TRUE,
        dotSize = 0.1) {

        counts <- SummarizedExperiment::assay(x, i = useAssay)
        cluster <- celdaClusters(x, altExpName = altExpName)

        g <- .plotCeldaViolin(counts = counts,
            cluster = cluster,
            features = features,
            exactMatch = exactMatch,
            plotDots = plotDots,
            dotSize = dotSize)
        return(g)
    }
)


#' @rdname plotCeldaViolin
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' plotCeldaViolin(x = celdaCGSim$counts,
#'    celdaMod = celdaCGMod,
#'    features = "Gene_1")
#' @export
setMethod("plotCeldaViolin",
    signature(x = "matrix"),
    function(x,
        celdaMod,
        features,
        exactMatch = TRUE,
        plotDots = TRUE,
        dotSize = 0.1) {

        cluster <- celdaClusters(celdaMod)$z
        g <- .plotCeldaViolin(counts = x,
            cluster = cluster,
            features = features,
            exactMatch = exactMatch,
            plotDots = plotDots,
            dotSize = dotSize)
        return(g)
    }
)


.plotCeldaViolin <- function(counts,
    cluster,
    features,
    exactMatch = TRUE,
    plotDots = TRUE,
    dotSize = 0.1) {

    featuresIx <- retrieveFeatureIndex(features,
        counts,
        by = "rownames",
        exactMatch = exactMatch)
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
