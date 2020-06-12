#' @title Heatmap for featureModules
#' @description Renders a heatmap for selected featureModules. Cells are
#'  ordered from those with the lowest probability of the module on the left to
#'  the highest probability on the right. If more than one module is used, then
#'  cells will be ordered by the probabilities of the first module only.
#'  Features are ordered from those with the highest probability in the module
#'  on the top to the lowest probability on the bottom.
#' @param x A numeric \link{matrix} of counts or a
#'  \linkS4class{SingleCellExperiment}
#'  with the matrix located in the assay slot under \code{useAssay}.
#'  Rows represent features and columns represent cells.
#' @param useAssay A string specifying which \link[SummarizedExperiment]{assay}
#'  slot to use if \code{x} is a
#'  \linkS4class{SingleCellExperiment} object. Default "counts".
#' @param celdaMod Celda object of class \link{celda_G} or \link{celda_CG}. Used
#'  only if \code{x} is a matrix object.
#' @param featureModule Integer Vector. The featureModule(s) to display.
#'  Multiple modules can be included in a vector.
#' @param topCells Integer. Number of cells with the highest and lowest
#'  probabilities for this module to include in the heatmap. For example, if
#'  `topCells` = 50, the 50 cells with the lowest probability and the 50 cells
#'  with the highest probability for that featureModule will be included. If
#'  NULL, all cells will be plotted. Default 100.
#' @param topFeatures Integer. Plot `topFeatures` with the highest probability
#'  in the featureModule. If NULL, plot all features in the module. Default
#'  NULL.
#' @param normalizedCounts Integer matrix. Rows represent features and columns
#'  represent cells. This matrix should correspond to the one provided for
#'  `counts`, but should be passed through. If NA, normalize `counts`.
#'  Default NA.
#'  `normalizeCounts(counts, "proportion", transformationFun=sqrt)`. Use of this
#'  parameter is particularly useful for plotting many moduleHeatmaps, where
#'  normalizing the counts matrix repeatedly would be too time consuming.
#' @param scaleRow Character. Which function to use to scale each individual
#'  row. Set to NULL to disable. Occurs after normalization and log
#'  transformation. For example, `scale` will Z-score transform each row.
#'  Default `scale`.
#' @param showFeaturenames Logical. Wheter feature names should be displayed.
#'  Default TRUE.
#' @return A list containing row and column dendrograms as well as a gtable for
#'  grob plotting
#' @importFrom methods .hasSlot
#' @export
setGeneric("moduleHeatmap", function(x, ...) {
    standardGeneric("moduleHeatmap")})


#' @rdname moduleHeatmap
#' @examples
#' data(sceCeldaCG)
#' moduleHeatmap(sceCeldaCG)
#' @export
setMethod("moduleHeatmap",
    signature(x = "SingleCellExperiment"),
    function(x,
        useAssay = "counts",
        featureModule = 1,
        topCells = 100,
        topFeatures = NULL,
        normalizedCounts = NA,
        scaleRow = scale,
        showFeaturenames = TRUE) {

        counts <- SummarizedExperiment::assay(x, i = useAssay)
        if (is.null(colnames(counts))) {
            stop("colnames(x) is NULL! Please assign column names to x and",
                " try again.")
        }

        if (is.null(rownames(counts))) {
            stop("rownames(x) is NULL! Please assign row names to x and",
                " try again.")
        }

        if (!(S4Vectors::metadata(x)$celda_parameters$model %in% c("celda_G",
            "celda_CG"))) {
            stop("metadata(x)$celda_parameters$model must be 'celda_G' or",
                " 'celda_CG'")
        }

        # factorize counts matrix
        factorizedMatrix <- factorizeMatrix(x, useAssay = useAssay)

        # take topRank
        if (!is.null(topFeatures) && (is.numeric(topFeatures)) |
                is.integer(topFeatures)) {
            topRanked <- topRank(
                matrix = factorizedMatrix$proportions$module,
                n = topFeatures
            )
        } else {
            topRanked <- topRank(
                matrix = factorizedMatrix$proportions$module,
                n = nrow(factorizedMatrix$proportions$module)
            )
        }

        # filter topRank using featureModule into featureIndices
        featureIndices <- lapply(
            featureModule,
            function(module) {
                topRanked$index[[module]]
            }
        )
        featureIndices <- unlist(featureIndices)

        # Determine cell order from factorizedMatrix$proportions$cell
        cellStates <- factorizedMatrix$proportions$cell
        cellStates <- cellStates[featureModule, , drop = FALSE]

        singleModule <- cellStates[1, ]
        singleModuleOrdered <- order(singleModule, decreasing = TRUE)

        if (!is.null(topCells)) {
            if (topCells * 2 < ncol(cellStates)) {
                cellIndices <- c(
                    utils::head(singleModuleOrdered, n = topCells),
                    utils::tail(singleModuleOrdered, n = topCells)
                )
            } else {
                cellIndices <- singleModuleOrdered
            }
        } else {
            cellIndices <- singleModuleOrdered
        }

        cellIndices <- rev(cellIndices)
        if (is.na(normalizedCounts)) {
            normCounts <- normalizeCounts(counts,
                normalize = "proportion",
                transformationFun = sqrt
            )
        } else {
            normCounts <- normalizedCounts
        }

        # filter counts based on featureIndices
        filteredNormCounts <-
            normCounts[featureIndices, cellIndices, drop = FALSE]

        filteredNormCounts <-
            filteredNormCounts[rowSums(filteredNormCounts > 0) > 0, ,
                drop = FALSE]

        geneIx <- match(rownames(filteredNormCounts), rownames(x))
        cellIx <- match(colnames(filteredNormCounts), colnames(x))
        zToPlot <- c()
        anno_cell_colors <- NULL
        if (S4Vectors::metadata(x)$celda_parameters$model == "celda_CG") {
            if ("celda_cell_cluster" %in%
                    colnames(SummarizedExperiment::colData(x))) {
                cell <-
                    distinctColors(length(unique(celdaClusters(x))))[
                        sort(unique(celdaClusters(x)[cellIx]))]
                names(cell) <- sort(unique(celdaClusters(x)[cellIx]))
                anno_cell_colors <- list(cell = cell)
                zToPlot <- celdaClusters(x)[cellIndices]
            }
        }

        plt <- plotHeatmap(
            filteredNormCounts,
            z = zToPlot,
            y = celdaModules(x)[geneIx],
            scaleRow = scaleRow,
            colorScheme = "divergent",
            showNamesFeature = showFeaturenames,
            clusterFeature = FALSE,
            clusterCell = FALSE,
            annotationColor = anno_cell_colors
        )
        invisible(plt)
    }
)


#' @rdname moduleHeatmap
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' moduleHeatmap(celdaCGSim$counts, celdaCGMod)
#' @export
setMethod("moduleHeatmap",
    signature(x = "matrix"),
    function(x,
        celdaMod,
        featureModule = 1,
        topCells = 100,
        topFeatures = NULL,
        normalizedCounts = NA,
        scaleRow = scale,
        showFeaturenames = TRUE) {

        counts <- x
        # Input checks
        if (is.null(counts) || !is.matrix(counts) & !is.data.frame(counts)) {
            stop("'counts' should be a numeric count matrix")
        }
        if (is.null(celdaMod) || !methods::is(celdaMod, "celda_G") &
                !methods::is(celdaMod, "celda_CG")) {
            stop("'celdaMod' should be an object of class celda_G or celda_CG")
        }
        compareCountMatrix(counts, celdaMod)

        # factorize counts matrix
        factorizedMatrix <- factorizeMatrix(x = counts, celdaMod = celdaMod)

        # take topRank
        if (!is.null(topFeatures) && (is.numeric(topFeatures)) |
                is.integer(topFeatures)) {
            topRanked <- topRank(
                matrix = factorizedMatrix$proportions$module,
                n = topFeatures
            )
        } else {
            topRanked <- topRank(
                matrix = factorizedMatrix$proportions$module,
                n = nrow(factorizedMatrix$proportions$module)
            )
        }

        # filter topRank using featureModule into featureIndices
        featureIndices <- lapply(
            featureModule,
            function(module) {
                topRanked$index[[module]]
            }
        )
        featureIndices <- unlist(featureIndices)

        # Determine cell order from factorizedMatrix$proportions$cell
        cellStates <- factorizedMatrix$proportions$cell
        cellStates <- cellStates[featureModule, , drop = FALSE]

        singleModule <- cellStates[1, ]
        singleModuleOrdered <- order(singleModule, decreasing = TRUE)

        if (!is.null(topCells)) {
            if (topCells * 2 < ncol(cellStates)) {
                cellIndices <- c(
                    utils::head(singleModuleOrdered, n = topCells),
                    utils::tail(singleModuleOrdered, n = topCells)
                )
            } else {
                cellIndices <- singleModuleOrdered
            }
        } else {
            cellIndices <- singleModuleOrdered
        }

        cellIndices <- rev(cellIndices)
        if (is.na(normalizedCounts)) {
            normCounts <- normalizeCounts(counts,
                normalize = "proportion",
                transformationFun = sqrt
            )
        } else {
            normCounts <- normalizedCounts
        }

        # filter counts based on featureIndices
        filteredNormCounts <-
            normCounts[featureIndices, cellIndices, drop = FALSE]

        filteredNormCounts <-
            filteredNormCounts[rowSums(filteredNormCounts > 0) > 0, ,
                drop = FALSE]

        geneIx <- match(rownames(filteredNormCounts),
            matrixNames(celdaMod)$row)
        cellIx <- match(colnames(filteredNormCounts),
            matrixNames(celdaMod)$column)
        zToPlot <- c()
        anno_cell_colors <- NULL
        if (class(celdaMod)[1] == "celda_CG") {
            if (methods::.hasSlot(celdaMod, "clusters")) {
                cell <-
                    distinctColors(length(unique(celdaClusters(celdaMod)$z)))[
                        sort(unique(celdaClusters(celdaMod)$z[cellIx]))
                        ]
                names(cell) <- sort(unique(celdaClusters(celdaMod)$z[cellIx]))
                anno_cell_colors <- list(cell = cell)
                zToPlot <- celdaClusters(celdaMod)$z[cellIndices]
            }
        }

        plt <- plotHeatmap(
            filteredNormCounts,
            z = zToPlot,
            y = celdaClusters(celdaMod)$y[geneIx],
            scaleRow = scaleRow,
            colorScheme = "divergent",
            showNamesFeature = showFeaturenames,
            clusterFeature = FALSE,
            clusterCell = FALSE,
            annotationColor = anno_cell_colors
        )
        invisible(plt)
    }
)
