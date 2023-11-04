#' @title Heatmap for featureModules
#' @description Renders a heatmap for selected \code{featureModule}. Cells are
#'  ordered from those with the lowest probability of the module on the left to
#'  the highest probability on the right. Features are ordered from those
#'  with the highest probability in the module
#'  on the top to the lowest probability on the bottom. Use of
#'  \link[multipanelfigure]{save_multi_panel_figure} is recommended for
#'  outputting figures in various formats.
#' @param x A numeric \link{matrix} of counts or a
#'  \linkS4class{SingleCellExperiment}
#'  with the matrix located in the assay slot under \code{useAssay}.
#'  Rows represent features and columns represent cells. Celda
#'  results must be present under \code{metadata(altExp(x, altExpName))}.
#' @param useAssay A string specifying which \link{assay}
#'  slot to use if \code{x} is a
#'  \linkS4class{SingleCellExperiment} object. Default "counts".
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default "featureSubset".
#' @param modules Integer Vector. The featureModule(s) to display.
#'  Multiple modules can be included in a vector. Default \code{NULL} which
#'  plots all module heatmaps.
#' @param featureModule Same as \code{modules}. Either can be used to specify
#' the modules to display.
#' @param col Passed to \link[ComplexHeatmap]{Heatmap}. Set color boundaries
#'  and colors.
#' @param topCells Integer. Number of cells with the highest and lowest
#'  probabilities for each module to include in the heatmap. For example, if
#'  \code{topCells = 50}, the 50 cells with the lowest probabilities and
#'  the 50 cells
#'  with the highest probabilities for each featureModule will be included. If
#'  NULL, all cells will be plotted. Default 100.
#' @param topFeatures Integer. Plot `topFeatures` features with the highest
#'  probabilities in the module heatmap for each featureModule. If \code{NULL},
#'  plot all features in the module. Default \code{NULL}.
#' @param normalizedCounts Integer matrix. Rows represent features and columns
#'  represent cells. If you have a normalized matrix result from
#'  \link{normalizeCounts}, you can pass through the result here to
#'  skip the normalization step in this function. Make sure the colnames and
#'  rownames match the object in x. This matrix should
#'  correspond to one generated from this count matrix
#'  \code{assay(altExp(x, altExpName), i = useAssay)}. If \code{NA},
#'  normalization will be carried out in the following form
#'  \code{normalizeCounts(assay(altExp(x, altExpName), i = useAssay),
#'  normalize = "proportion", transformationFun = sqrt)}.
#'  Use of this parameter is particularly useful for plotting many
#'  module heatmaps, where normalizing the counts matrix repeatedly would
#'  be too time consuming. Default NA.
#' @param normalize Character. Passed to \link{normalizeCounts} if
#'  \code{normalizedCounts} is \code{NA}.
#'  Divides counts by the library sizes for each cell. One of 'proportion',
#'  'cpm', 'median', or 'mean'. 'proportion' uses the total counts for each
#'  cell as the library size. 'cpm' divides the library size of each cell by
#'  one million to produce counts per million. 'median' divides the library
#'  size of each cell by the median library size across all cells. 'mean'
#'  divides the library size of each cell by the mean library size across all
#'  cells. Default "proportion".
#' @param transformationFun Function. Passed to \link{normalizeCounts} if
#'  \code{normalizedCounts} is \code{NA}. Applies a transformation such as
#'  \link{sqrt}, \link{log}, \link{log2}, \link{log10}, or \link{log1p}.
#'  If \code{NULL}, no transformation will be applied. Occurs after
#'  normalization. Default \link{sqrt}.
#' @param scaleRow Function. Which function to use to scale each individual
#'  row. Set to NULL to disable. Occurs after normalization and log
#'  transformation. For example, \link{scale} will Z-score transform each row.
#'  Default \link{scale}.
#' @param showFeatureNames Logical. Whether feature names should be displayed.
#'  Default TRUE.
#' @param displayName Character. The column name of
#'  \code{rowData(altExp(x, altExpName))} that specifies the display names for
#'  the features. Default \code{NULL}, which displays the row names. Only works
#'  if \code{showFeaturenames} is \code{TRUE} and \code{x} is a
#'  \linkS4class{SingleCellExperiment} object.
#' @param trim Numeric vector. Vector of length two that specifies the lower
#'  and upper bounds for plotting the data. This threshold is applied
#'  after row scaling. Set to NULL to disable. Default \code{c(-2,2)}.
#' @param rowFontSize Numeric. Font size for feature names. If \code{NULL},
#' then the size will automatically be determined. Default \code{NULL}.
#' @param showHeatmapLegend Passed to \link[ComplexHeatmap]{Heatmap}. Show
#'  legend for expression levels.
#' @param showTopAnnotationLegend Passed to
#'  \link[ComplexHeatmap]{HeatmapAnnotation}. Show legend for cell annotation.
#' @param showTopAnnotationName Passed to
#'  \link[ComplexHeatmap]{HeatmapAnnotation}. Show heatmap top annotation name.
#' @param topAnnotationHeight Passed to
#'  \link[ComplexHeatmap]{HeatmapAnnotation}. Column annotation height.
#'  \link[ComplexHeatmap]{rowAnnotation}. Show legend for module annotation.
#' @param showModuleLabel Show left side module labels.
#' @param moduleLabel The left side row titles for module heatmap. Must be
#'  vector of the same length as \code{featureModule}. Default "auto", which
#'  automatically pulls module labels from \code{x}.
#' @param moduleLabelSize Passed to \link{gpar}. The size of text (in points).
#' @param byrow Passed to \link{matrix}. logical. If \code{FALSE} (the default)
#'  the figure panel is filled by columns, otherwise the figure panel is filled
#'  by rows.
#' @param top Passed to \link[gridExtra]{marrangeGrob}. The title for each page.
#' @param unit Passed to \link[grid]{unit}. Single
#'  character object defining the unit of all dimensions defined.
#' @param ncol Integer. Number of columns of module heatmaps. If \code{NULL},
#' then this will be automatically calculated so that the number of columns
#' and rows will be approximately the same. Default \code{NULL}.
#' @param useRaster Boolean. Rasterizing will make the heatmap a single object
#' and reduced the memory of the plot and the size of a file. If \code{NULL},
#' then rasterization will be automatically determined by the underlying
#' \link[ComplexHeatmap]{Heatmap} function. Default \code{TRUE}.
#' @param returnAsList Boolean. If \code{TRUE}, then a list of plots will be
#' returned instead of a single multi-panel figure. These plots can be
#' displayed using the \link[grid]{grid.draw} function. Default \code{FALSE}.
#' @param ... Additional parameters passed to \link[ComplexHeatmap]{Heatmap}.
#' @return A \link[multipanelfigure]{multi_panel_figure} object if plotting
#'  more than one module heatmaps. Otherwise a
#'  \link[ComplexHeatmap]{HeatmapList} object is returned.
#' @importFrom methods .hasSlot
#' @importFrom multipanelfigure multi_panel_figure
#' @export
setGeneric("moduleHeatmap",
    function(x,
        useAssay = "counts",
        altExpName = "featureSubset",
        modules = NULL,
        featureModule = NULL,
        col = circlize::colorRamp2(c(-2, 0, 2),
            c("#1E90FF", "#FFFFFF", "#CD2626")),
        topCells = 100,
        topFeatures = NULL,
        normalizedCounts = NA,
        normalize = "proportion",
        transformationFun = sqrt,
        scaleRow = scale,
        showFeatureNames = TRUE,
        displayName = NULL,
        trim = c(-2, 2),
        rowFontSize = NULL,
        showHeatmapLegend = FALSE,
        showTopAnnotationLegend = FALSE,
        showTopAnnotationName = FALSE,
        topAnnotationHeight = 5,
        showModuleLabel = TRUE,
        moduleLabel = "auto",
        moduleLabelSize = NULL,
        byrow = TRUE,
        top = NA,
        unit = "mm",
        ncol = NULL,
        useRaster = TRUE,
        returnAsList = FALSE,
        ...) {
    standardGeneric("moduleHeatmap")})


#' @rdname moduleHeatmap
#' @examples
#' data(sceCeldaCG)
#' moduleHeatmap(sceCeldaCG, displayName = "rownames")
#' @export
setMethod("moduleHeatmap",
    signature(x = "SingleCellExperiment"),
    function(x,
        useAssay = "counts",
        altExpName = "featureSubset",
        modules = NULL,
        featureModule = NULL,
        col = circlize::colorRamp2(c(-2, 0, 2),
            c("#1E90FF", "#FFFFFF", "#CD2626")),
        topCells = 100,
        topFeatures = NULL,
        normalizedCounts = NA,
        normalize = "proportion",
        transformationFun = sqrt,
        scaleRow = scale,
        showFeatureNames = TRUE,
        displayName = NULL,
        trim = c(-2, 2),
        rowFontSize = NULL,
        showHeatmapLegend = FALSE,
        showTopAnnotationLegend = FALSE,
        showTopAnnotationName = FALSE,
        topAnnotationHeight = 5,
        showModuleLabel = TRUE,
        moduleLabel = "auto",
        moduleLabelSize = NULL,
        byrow = TRUE,
        top = NA,
        unit = "mm",
        ncol = NULL,
        useRaster = TRUE,
        returnAsList = FALSE,
        ...) {

        # 'modules' is an easier parameter name to remember so we include
        # support for both.
        if(!is.null(modules)) {
            featureModule <- modules
        }

        altExp <- SingleCellExperiment::altExp(x, altExpName)

        counts <- SummarizedExperiment::assay(altExp, i = useAssay)
        if (is.null(colnames(counts))) {
            stop("colnames(altExp(x, altExpName)) is NULL!",
                " Please assign column names to x and",
                " try again.")
        }

        if (is.null(rownames(counts))) {
            stop("rownames(altExp(x, altExpName)) is NULL!",
                " Please assign row names to x and",
                " try again.")
        }

        if (!(S4Vectors::metadata(altExp)$celda_parameters$model %in%
                c("celda_G", "celda_CG"))) {
            stop("metadata(altExp(x, altExpName))$",
                "celda_parameters$model must be 'celda_G' or",
                " 'celda_CG'")
        }

        if (is.null(featureModule)) {
            featureModule <- sort(unique(celdaModules(x)))
        }

        if (length(featureModule) == 1) {
            returnHeatmap <- TRUE
        } else {
            returnHeatmap <- FALSE
        }

        if (moduleLabel == "auto") {
            moduleLabel <- paste0("Module ", as.character(featureModule))
        } else if (length(moduleLabel) == 1 & length(featureModule) > 1) {
            moduleLabel <- rep(moduleLabel, length(featureModule))
        } else if (length(moduleLabel) != length(featureModule)) {
            stop("Invalid 'moduleLabel' length")
        }

        # factorize counts matrix
        factorizedMatrix <- factorizeMatrix(x,
            useAssay = useAssay,
            altExpName = altExpName,
            type = "proportion")
        allCellStates <- factorizedMatrix$proportions$cell

        if (is.na(normalizedCounts)) {
            normCounts <- normalizeCounts(counts,
                normalize = normalize,
                transformationFun = transformationFun)
        } else {
            normCounts <- normalizedCounts
        }

        # take topRank
        if (!is.null(topFeatures) && (is.numeric(topFeatures)) |
                is.integer(topFeatures)) {
            topRanked <- topRank(
                matrix = factorizedMatrix$proportions$module,
                n = topFeatures)
        } else {
            topRanked <- topRank(
                matrix = factorizedMatrix$proportions$module,
                n = nrow(factorizedMatrix$proportions$module))
        }

        # filter topRank using featureModule into featureIndices
        featureIndices <- lapply(
            featureModule,
            function(module) {
                topRanked$index[[module]]
            }
        )

        # Set up displayName variable if specified
        if (is.null(displayName)) {
            displayNames <- rownames(altExp)
        } else {
            displayNames <- SummarizedExperiment::rowData(altExp)[[
                displayName]]
        }

        z <- celdaClusters(x, altExpName = altExpName)
        y <- celdaModules(x, altExpName = altExpName)

        # Get max rowFontSize if multiple modules are selected
        if (is.null(rowFontSize)) {
          if (length(featureIndices) > 1 & !isTRUE(returnAsList)) {
            # If there is more than 1 module selected, then the miniumum size
            # size will be caculated for each module. This will ensure that
            # all modules will have the same rowFontSize and the module
            # heatmaps will have the same width.
            maxlen <- max(unlist(lapply(featureIndices, length)))
            maxlen <- maxlen * sqrt(length(featureIndices))
            rowFontSize <- rep(min(200 / maxlen, 20), length(featureIndices))
          } else {
            # If there is only one plot or each plot will be generated
            # separately and returned in a list, then the size of the labels,
            # will be caculated for each module separately.
            len <- unlist(lapply(featureIndices, length))
            rowFontSize <- pmin(200 / len, 20)
          }
        }

        plts <- vector("list", length = length(featureModule))
        for (i in seq(length(featureModule))) {
            plts[[i]] <- .plotModuleHeatmap(normCounts = normCounts,
                col = col,
                allCellStates = allCellStates,
                featureIndices = featureIndices[[i]],
                featureModule = featureModule[i],
                z = z,
                y = y,
                topCells = topCells,
                altExpName = altExpName,
                scaleRow = scaleRow,
                showFeatureNames = showFeatureNames,
                displayNames = displayNames[featureIndices[[i]]],
                trim = trim,
                rowFontSize = rowFontSize[i],
                showHeatmapLegend = showHeatmapLegend,
                showTopAnnotationLegend = showTopAnnotationLegend,
                showTopAnnotationName = showTopAnnotationName,
                topAnnotationHeight = topAnnotationHeight,
                showModuleLabel = showModuleLabel,
                moduleLabel = moduleLabel[i],
                moduleLabelSize = moduleLabelSize,
                useRaster = useRaster,
                unit = unit,
                ... = ...)
        }

        if (isTRUE(returnHeatmap)) {
            return(plts[[1]])
        } else {
            if (is.null(ncol)) {
                ncol <- floor(sqrt(length(plts)))
            }
            nrow <- ceiling(length(plts) / ncol)

            for (i in seq(length(plts))) {
                plts[[i]] <- grid::grid.grabExpr(
                    ComplexHeatmap::draw(plts[[i]]),
                    wrap.grobs = TRUE)
            }

            if (isTRUE(returnAsList)) {
                figure <- plts
            } else {
                figure <- marrangeGrob(plts,
                    ncol = ncol,
                    nrow = nrow,
                    layout_matrix = matrix(seq_len(nrow * ncol),
                        nrow = nrow,
                        ncol = ncol,
                        byrow = TRUE),
                    top = NA)
            }

            suppressWarnings(return(figure))
        }
    }
)


.plotModuleHeatmap <- function(normCounts,
    col,
    allCellStates,
    featureIndices,
    featureModule,
    z,
    y,
    topCells,
    altExpName,
    scaleRow,
    showFeatureNames,
    displayNames,
    trim,
    rowFontSize,
    showHeatmapLegend,
    showTopAnnotationLegend,
    showTopAnnotationName,
    topAnnotationHeight,
    showModuleLabel,
    moduleLabel,
    moduleLabelSize,
    useRaster,
    unit,
    ...) {

    # Determine cell order from factorizedMatrix$proportions$cell
    cellStates <- allCellStates[featureModule, , drop = TRUE]

    singleModuleOrdered <- order(cellStates, decreasing = TRUE)

    if (!is.null(topCells)) {
        if (topCells * 2 < ncol(allCellStates)) {
            cellIndices <- c(
                utils::head(singleModuleOrdered, n = topCells),
                utils::tail(singleModuleOrdered, n = topCells))
        } else {
            cellIndices <- singleModuleOrdered
        }
    } else {
        cellIndices <- singleModuleOrdered
    }

    cellIndices <- rev(cellIndices)

    # filter counts based on featureIndices
    filteredNormCounts <-
        normCounts[featureIndices, cellIndices, drop = FALSE]

    # Show/hide features with 0 counts in these cells in the module
    # filteredNormCounts <-
    #     filteredNormCounts[rowSums(filteredNormCounts > 0) > 0, ,
    #         drop = FALSE]

    geneIx <- match(rownames(filteredNormCounts), rownames(normCounts))
    cellIx <- match(colnames(filteredNormCounts), colnames(normCounts))

    zToPlot <- z[cellIx]

    uniquezToPlot <- sort(unique(zToPlot))
    ccols <- distinctColors(length(unique(z)))[uniquezToPlot]
    names(ccols) <- uniquezToPlot

    yToPlot <- y[geneIx]

    uniqueyToPlot <- sort(unique(yToPlot))
    rcols <- distinctColors(length(y))[uniqueyToPlot]
    names(rcols) <- uniqueyToPlot

    # scale indivisual rows by scaleRow
    if (!is.null(scaleRow)) {
        if (is.function(scaleRow)) {
            cn <- colnames(filteredNormCounts)
            filteredNormCounts <- t(base::apply(filteredNormCounts,
                1, scaleRow))
            colnames(filteredNormCounts) <- cn
        } else {
            stop("'scaleRow' needs to be of class 'function'")
        }
        # If the standard deviation was 0 then the values will be NA
        # Replacing the NAs with zero will keep the row the middle color
        # rather than grey (default with ComplexHeatmap)
        filteredNormCounts[is.na(filteredNormCounts)] <- 0
    }

    if (!is.null(trim)) {
        if (length(trim) != 2) {
            stop(
                "'trim' should be a 2 element vector specifying the lower",
                " and upper boundaries"
            )
        }
        trim <- sort(trim)
        filteredNormCounts[filteredNormCounts < trim[1]] <- trim[1]
        filteredNormCounts[filteredNormCounts > trim[2]] <- trim[2]
    }

    if (isTRUE(showModuleLabel)) {
        plt <- ComplexHeatmap::Heatmap(matrix = filteredNormCounts,
            col = col,
            row_title = moduleLabel,
            row_title_gp = grid::gpar(fontsize = moduleLabelSize),
            show_column_names = FALSE,
            show_row_names = showFeatureNames,
            row_labels = displayNames,
            row_names_gp = grid::gpar(fontsize = rowFontSize),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            heatmap_legend_param = list(title = "Expression"),
            show_heatmap_legend = showHeatmapLegend,
            use_raster = useRaster,
            top_annotation = ComplexHeatmap::HeatmapAnnotation(
                cell = factor(zToPlot,
                    levels = stringr::str_sort(unique(zToPlot),
                        numeric = TRUE)),
                show_legend = showTopAnnotationLegend,
                show_annotation_name = showTopAnnotationName,
                col = list(cell = ccols),
                simple_anno_size = grid::unit(topAnnotationHeight, unit),
                simple_anno_size_adjust = TRUE),
            ...)
    } else {
        plt <- ComplexHeatmap::Heatmap(matrix = filteredNormCounts,
            col = col,
            show_column_names = FALSE,
            show_row_names = showFeatureNames,
            row_labels = displayNames,
            row_names_gp = grid::gpar(fontsize = rowFontSize),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            heatmap_legend_param = list(title = "Expression"),
            show_heatmap_legend = showHeatmapLegend,
            use_raster = useRaster,
            top_annotation = ComplexHeatmap::HeatmapAnnotation(
                cell = factor(zToPlot,
                    levels = stringr::str_sort(unique(zToPlot),
                        numeric = TRUE)),
                show_legend = showTopAnnotationLegend,
                show_annotation_name = showTopAnnotationName,
                col = list(cell = ccols),
                simple_anno_size = grid::unit(topAnnotationHeight, unit),
                simple_anno_size_adjust = TRUE),
            ...)
    }
    return(plt)
}
