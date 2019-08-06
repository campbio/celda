#' @title Plots heatmap based on Celda model
#' @description Renders a heatmap based on a matrix of counts where rows are
#'  features and columns are cells.
#' @param counts Numeric matrix. Normalized counts matrix where rows represent
#'  features and columns represent cells. .
#' @param z Numeric vector. Denotes cell population labels.
#' @param y Numeric vector. Denotes feature module labels.
#' @param featureIx Integer vector. Select features for display in heatmap. If
#'  NULL, no subsetting will be performed. Default NULL.
#' @param cellIx Integer vector. Select cells for display in heatmap. If NULL,
#'  no subsetting will be performed. Default NULL.
#' @param scaleRow Function. A function to scale each individual row. Set to
#'  NULL to disable. Occurs after normalization and log transformation. Defualt
#'  is 'scale' and thus will Z-score transform each row.
#' @param trim Numeric vector. Vector of length two that specifies the lower
#'  and upper bounds for the data. This threshold is applied after row scaling.
#'  Set to NULL to disable. Default c(-2,2).
#' @param clusterFeature Logical. Determines whether rows should be clustered.
#'  Default TRUE.
#' @param clusterCell Logical. Determines whether columns should be clustered.
#'  Default TRUE.
#' @param annotationCell Data frame. Additional annotations for each cell will
#'  be shown in the column color bars. The format of the data frame should be
#'  one row for each cell and one column for each annotation. Numeric variables
#'  will be displayed as continuous color bars and factors will be displayed as
#'  discrete color bars. Default NULL.
#' @param annotationFeature A data frame for the feature annotations (rows).
#' @param annotationColor List. Contains color scheme for all annotations. See
#'  `?pheatmap` for more details.
#' @param colorScheme "Character. One of ""divergent"" or ""sequential"". A
#'  ""divergent"" scheme is best for highlighting relative data (denoted by
#'  'colorSchemeCenter') such as gene expression data that has been normalized
#'  and centered. A ""sequential"" scheme is best for highlighting data that
#'  are ordered low to high such as raw counts or probabilities. Default
#'  "divergent".
#' @param colorSchemeSymmetric Logical. When the colorScheme is "divergent"
#'  and the data contains both positive and negative numbers, TRUE indicates
#'  that the color scheme should be symmetric from
#'  [-max(abs(data)), max(abs(data))]. For example, if the data ranges goes
#'  from -1.5 to 2, then setting this to TRUE will force the color scheme to
#'  range from -2 to 2. Default TRUE.
#' @param colorSchemeCenter Numeric. Indicates the center of a "divergent"
#'  colorScheme. Default 0.
#' @param col Color for the heatmap.
#' @param breaks Numeric vector. A sequence of numbers that covers the range
#'  of values in the normalized `counts`. Values in the normalized `matrix` are
#'  assigned to each bin in `breaks`. Each break is assigned to a unique color
#'  from `col`. If NULL, then breaks are calculated automatically. Default NULL.
#' @param legend Logical. Determines whether legend should be drawn. Default
#'  TRUE.
#' @param annotationLegend Logical. Whether legend for all annotations should
#'  be drawn. Default TRUE.
#' @param annotationNamesFeature Logical. Whether the names for features should
#'  be shown. Default TRUE.
#' @param annotationNamesCell Logical. Whether the names for cells should be
#'  shown. Default TRUE.
#' @param showNamesFeature Logical. Specifies if feature names should be shown.
#'  Default TRUE.
#' @param showNamesCell Logical. Specifies if cell names should be shown.
#'  Default FALSE.
#' @param hclustMethod Character. Specifies the method to use for the 'hclust'
#'  function. See `?hclust` for possible values. Default "ward.D2".
#' @param treeheightFeature Numeric. Width of the feature dendrogram. Set to 0
#'  to disable plotting of this dendrogram. Default: if clusterFeature == TRUE,
#'  then treeheightFeature = 50, else treeheightFeature = 0.
#' @param treeheightCell Numeric. Height of the cell dendrogram. Set to 0 to
#'  disable plotting of this dendrogram. Default: if clusterCell == TRUE, then
#'  treeheightCell = 50, else treeheightCell = 0.
#' @param silent Logical. Whether to plot the heatmap.
#' @param ... Other arguments to be passed to underlying pheatmap function.
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' plotHeatmap(celdaCGSim$counts,
#'     z = clusters(celdaCGMod)$z, y = clusters(celdaCGMod)$y)
#' @return list A list containing dendrogram information and the heatmap grob
#' @import graphics
#' @import grid
#' @export
plotHeatmap <- function(counts,
    z = NULL,
    y = NULL,
    rowGroupOrder = NULL,
    colGroupOrder = NULL, 
    scaleRow = scale,
    trim = c(-2, 2),
    featureIx = NULL,
    cellIx = NULL,
    clusterFeature = TRUE,
    clusterCell = TRUE,
    colorScheme = c("divergent", "sequential"),
    colorSchemeSymmetric = TRUE,
    colorSchemeCenter = 0,
    col = NULL,
    annotationCell = NULL,
    annotationFeature = NULL,
    annotationColor = NULL,
    breaks = NULL,
    legend = TRUE,
    annotationLegend = TRUE,
    annotationNamesFeature = TRUE,
    annotationNamesCell = TRUE,
    showNamesFeature = FALSE,
    showNamesCell = FALSE,
    hclustMethod = "ward.D2",
    treeheightFeature = ifelse(clusterFeature, 50, 0),
    treeheightCell = ifelse(clusterCell, 50, 0),
    silent = FALSE,
    ...) {

    # Check for same lengths for z and y group variables
    if (!is.null(z) & length(z) != ncol(counts)) {
        stop("Length of z must match number of columns in counts matrix")
    }

    if (!is.null(y) & length(y) != nrow(counts)) {
        stop("Length of y must match number of rows in counts matrix")
    }

    colorScheme <- match.arg(colorScheme)

    ## Create cell annotation
    if (!is.null(annotationCell) & !is.null(z)) {
        if (is.null(rownames(annotationCell))) {
            rownames(annotationCell) <- colnames(counts)
        } else {
            if (any(rownames(annotationCell) != colnames(counts))) {
                stop("Row names of 'annotationCell' are different than the
                    column names of 'counts'")
            }
        }
        annotationCell <- data.frame(cell = as.factor(z), annotationCell)
    } else if (is.null(annotationCell) & !is.null(z)) {
        annotationCell <- data.frame(cell = as.factor(z))
        rownames(annotationCell) <- colnames(counts)
    } else {
        annotationCell <- NA
    }

    # Set feature annotation
    if (!is.null(annotationFeature) & !is.null(y)) {
        if (is.null(rownames(annotationFeature))) {
            rownames(annotationFeature) <- rownames(counts)
        } else {
            if (any(rownames(annotationFeature) != rownames(counts))) {
                stop("Row names of 'annotationFeature' are different than the
                    row names of 'counts'")
            }
        }
        annotationFeature <- data.frame(module = as.factor(y),
            annotationFeature)
    } else if (is.null(annotationFeature) & !is.null(y)) {
        annotationFeature <- data.frame(module = as.factor(y))
        rownames(annotationFeature) <- rownames(counts)
    } else {
        annotationFeature <- NA
    }

    ## Select subsets of features/cells
    if (!is.null(featureIx)) {
        counts <- counts[featureIx, , drop = FALSE]
        if (length(annotationFeature) > 1 ||
                (length(annotationFeature) == 1 & !is.na(annotationFeature))) {
            annotationFeature <- annotationFeature[featureIx, , drop = FALSE]
        }
        if (!is.null(y)) {
            y <- y[featureIx]
        }
    }

    if (!is.null(cellIx)) {
        counts <- counts[, cellIx, drop = FALSE]
        if (length(annotationCell) > 1 ||
                (length(annotationCell) == 1 & !is.na(annotationCell))) {
            annotationCell <- annotationCell[cellIx, , drop = FALSE]
        }
        if (!is.null(z)) {
            z <- z[cellIx]
        }
    }

    ## Set annotation colors
    if (!is.null(z)) {
        K <- sort(unique(z))
        kCol <- distinctColors(length(K))
        names(kCol) <- K

        if (!is.null(annotationColor)) {
            if (!("cell" %in% names(annotationColor))) {
                annotationColor <- c(list(cell = kCol), annotationColor)
            }
        } else {
            annotationColor <- list(cell = kCol)
        }
    }

    if (!is.null(y)) {
        L <- sort(unique(y))
        lCol <- distinctColors(length(L))
        names(lCol) <- L
        if (!is.null(annotationColor)) {
            if (!("module" %in% names(annotationColor))) {
                annotationColor <- c(list(module = lCol), annotationColor)
            }
        } else {
            annotationColor <- list(module = lCol)
        }
    }

    # scale indivisual rows by scaleRow
    if (!is.null(scaleRow)) {
        if (is.function(scaleRow)) {
            cn <- colnames(counts)
            counts <- t(base::apply(counts, 1, scaleRow))
            colnames(counts) <- cn
        } else {
            stop("'scaleRow' needs to be of class 'function'")
        }
    }

    if (!is.null(trim)) {
        if (length(trim) != 2) {
            stop("'trim' should be a 2 element vector specifying the lower",
                " and upper boundaries")
        }
        trim <- sort(trim)
        counts[counts < trim[1]] <- trim[1]
        counts[counts > trim[2]] <- trim[2]
    }

    ## Set color scheme and breaks
    uBoundRange <- max(counts, na.rm = TRUE)
    lboundRange <- min(counts, na.rm = TRUE)

    if (colorScheme == "divergent") {
        if (colorSchemeSymmetric == TRUE) {
            uBoundRange <- max(abs(uBoundRange), abs(lboundRange))
            lboundRange <- -uBoundRange
        }
        if (is.null(col)) {
            col <- colorRampPalette(c("#1E90FF", "#FFFFFF", "#CD2626"),
                space = "Lab")(100)
        }
        colLen <- length(col)
        if (is.null(breaks)) {
            breaks <- c(seq(lboundRange, colorSchemeCenter,
                length.out = round(colLen / 2) + 1),
                seq(colorSchemeCenter + 1e-6, uBoundRange,
                    length.out = colLen - round(colLen / 2)))
        }
    } else {
        # Sequential color scheme
        if (is.null(col)) {
            col <- colorRampPalette(c("#FFFFFF", brewer.pal(n = 9,
                name = "Blues")))(100)
        }
        colLen <- length(col)
        if (is.null(breaks)) {
            breaks <- seq(lboundRange, uBoundRange, length.out = colLen)
        }
    }

    sp <- semiPheatmap(mat = counts,
        color = col,
        breaks = breaks,
        clusterCols = clusterCell,
        clusterRows = clusterFeature,
        annotationRow = annotationFeature,
        annotationCol = annotationCell,
        annotationColors = annotationColor,
        legend = legend,
        annotationLegend = annotationLegend,
        annotationNamesRow = annotationNamesFeature,
        annotationNamesCol = annotationNamesCell,
        showRownames = showNamesFeature,
        showColnames = showNamesCell,
        clusteringMethod = hclustMethod,
        treeHeightRow = treeheightFeature,
        treeHeightCol = treeheightCell,
        rowLabel = y,
        colLabel = z,
        rowGroupOrder = rowGroupOrder,
        colGroupOrder = colGroupOrder,
        silent = TRUE,
        ...)

    if (!isTRUE(silent)) {
        grid::grid.newpage()
        grid::grid.draw(sp$gtable)
    }

    invisible(sp)
}
