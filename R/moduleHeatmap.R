#' @title Heatmap for featureModules
#' @description Renders a heatmap for selected featureModules. Cells are
#'  ordered from those with the lowest probability of the module on the left to
#'  the highest probability on the right. Features are ordered from those
#'  with the highest probability in the module
#'  on the top to the lowest probability on the bottom.
#' @param x A numeric \link{matrix} of counts or a
#'  \linkS4class{SingleCellExperiment}
#'  with the matrix located in the assay slot under \code{useAssay}.
#'  Rows represent features and columns represent cells.
#' @param useAssay A string specifying which \link[SummarizedExperiment]{assay}
#'  slot to use if \code{x} is a
#'  \linkS4class{SingleCellExperiment} object. Default "counts".
#' @param altExpName The name for the \link[SingleCellExperiment]{altExp} slot
#'  to use. Default "featureSubset".
#' @param celdaMod Celda object of class \link{celda_G} or \link{celda_CG}. Used
#'  only if \code{x} is a matrix object.
#' @param featureModule Integer Vector. The featureModule(s) to display.
#'  Multiple modules can be included in a vector. Default \code{NULL} which
#'  plots all module heatmaps.
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
#'  \code{normalizedCounts} is \code{NA}. Applys a transformation such as
#'  \link{sqrt}, \link{log}, \link{log2}, \link{log10}, or \link{log1p}.
#'  If NULL, no transformation will be applied. Occurs after normalization.
#'  Default \link{sqrt}.
#' @param scaleRow Function. Which function to use to scale each individual
#'  row. Set to NULL to disable. Occurs after normalization and log
#'  transformation. For example, \link{scale} will Z-score transform each row.
#'  Default \link{scale}.
#' @param showFeaturenames Logical. Wheter feature names should be displayed.
#'  Default TRUE.
#' @param trim Numeric vector. Vector of length two that specifies the lower
#'  and upper bounds for plotting the data. This threshold is applied
#'  after row scaling. Set to NULL to disable. Default c(-2,2).
#' @param rowFontSize Integer. Font size for genes.
#' @param showHeatmapLegend Passed to \link[ComplexHeatmap]{Heatmap}. Show
#'  legend for expression levels.
#' @param showTopAnnotationLegend Passed to
#'  \link[ComplexHeatmap]{HeatmapAnnotation}. Show legend for cell annotation.
#' @param showTopAnnotationName Passed to
#'  \link[ComplexHeatmap]{HeatmapAnnotation}. Show heatmap top annotation name.
#' @param showLeftAnnotationLegend Passed to
#' @param annotationHeight Passed to
#'  \link[ComplexHeatmap]{HeatmapAnnotation}. Column annotation height.
#'  \link[ComplexHeatmap]{rowAnnotation}. Show legend for module annotation.
#' @param showLeftAnnotationName Passed to
#'  \link[ComplexHeatmap]{rowAnnotation}. Show lheatmap left annotation name.
#' @param annotationWidth Passed to
#'  \link[ComplexHeatmap]{rowAnnotation}. Row annotation width.
#' @param width Passed to \link[multipanelfigure]{multi_panel_figure}. The
#'  width of the output figure.
#' @param height Passed to \link[multipanelfigure]{multi_panel_figure}. The
#'  height of the output figure.
#' @param unit Passed to \link[multipanelfigure]{multi_panel_figure}. Single
#'  character object defining the unit of all dimensions defined.
#' @param ModuleLabel Must be
#'  vector of the same length as \code{length(unique(celdaModules(x)))} or
#'  \code{length(unique(celdaClusters(x)$y))}. Set to \code{""} to disable.
#' @param labelJust Passed to \link[multipanelfigure]{fill_panel}.
#'  Justification for the label within the interpanel spacing grob to the
#'  top-left of the panel content grob.
#' @param ... Additional parameters passed to \link[ComplexHeatmap]{Heatmap}.
#' @return A \link[multipanelfigure]{multi_panel_figure} object.
#' @importFrom methods .hasSlot
#' @export
setGeneric("moduleHeatmap", function(x, ...) {
    standardGeneric("moduleHeatmap")})


#' @rdname moduleHeatmap
#' @examples
#' data(sceCeldaCG)
#' moduleHeatmap(sceCeldaCG, width = 250, height = 250)
#' @export
setMethod("moduleHeatmap",
    signature(x = "SingleCellExperiment"),
    function(x,
        useAssay = "counts",
        altExpName = "featureSubset",
        featureModule = NULL,
        col = circlize::colorRamp2(c(-2, 0, 2),
            c("#1E90FF", "#FFFFFF", "#CD2626")),
        topCells = 100,
        topFeatures = NULL,
        normalizedCounts = NA,
        normalize = "proportion",
        transformationFun = sqrt,
        scaleRow = scale,
        showFeaturenames = TRUE,
        trim = c(-2, 2),
        rowFontSize = 6,
        showHeatmapLegend = FALSE,
        showTopAnnotationLegend = FALSE,
        showTopAnnotationName = FALSE,
        annotationHeight = 1.5,
        showLeftAnnotationLegend = FALSE,
        showLeftAnnotationName = FALSE,
        annotationWidth = 1.5,
        width = "auto",
        height = "auto",
        unit = "mm",
        ModuleLabel = "auto",
        labelJust = c("right", "bottom"),
        ...) {

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

        if (is.null(ModuleLabel)) {
            ModuleLabel <- NULL
        } else if (ModuleLabel == "auto") {
            ModuleLabel <- as.character(featureModule)
        } else if (ModuleLabel == "") {
            ModuleLabel <- rep("", length = length(unique(celdaModules(x,
                altExpName = altExpName))))
        } else if (length(ModuleLabel) != length(unique(celdaModules(x,
            altExpName = altExpName)))) {
            stop("Invalid 'ModuleLabel' length!")
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

        z <- celdaClusters(x, altExpName = altExpName)
        y <- celdaModules(x, altExpName = altExpName)

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
                showFeaturenames = showFeaturenames,
                trim = trim,
                rowFontSize = rowFontSize,
                showHeatmapLegend = showHeatmapLegend,
                showTopAnnotationLegend = showTopAnnotationLegend,
                showTopAnnotationName = showTopAnnotationName,
                annotationHeight = annotationHeight,
                showLeftAnnotationLegend = showLeftAnnotationLegend,
                showLeftAnnotationName = showLeftAnnotationName,
                annotationWidth = annotationWidth,
                unit = unit,
                ... = ...)
        }


        ncol <- floor(sqrt(length(plts)))
        nrow <- ceiling(length(plts) / ncol)

        for (i in seq(length(plts))) {
            plts[[i]] <- grid::grid.grabExpr(ComplexHeatmap::draw(plts[[i]]),
                wrap.grobs = TRUE)
        }

        figure <- multipanelfigure::multi_panel_figure(columns = ncol,
            rows = nrow,
            width = width,
            height = height,
            unit = unit)

        for (i in seq(length(plts))) {
            if (!is.null(ModuleLabel)) {
                figure <- suppressMessages(multipanelfigure::fill_panel(figure,
                    plts[[i]], label = ModuleLabel[i], label_just = labelJust))
            } else {
                figure <- suppressMessages(multipanelfigure::fill_panel(figure,
                    plts[[i]], label_just = labelJust))
            }
        }
        suppressWarnings(return(figure))
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
    showFeaturenames,
    trim,
    rowFontSize,
    showHeatmapLegend,
    showTopAnnotationLegend,
    showTopAnnotationName,
    annotationHeight,
    showLeftAnnotationLegend,
    showLeftAnnotationName,
    annotationWidth,
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

    filteredNormCounts <-
        filteredNormCounts[rowSums(filteredNormCounts > 0) > 0, ,
            drop = FALSE]

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

    plt <- ComplexHeatmap::Heatmap(matrix = filteredNormCounts,
        col = col,
        show_column_names = FALSE,
        show_row_names = showFeaturenames,
        row_names_gp = grid::gpar(fontsize = rowFontSize),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        heatmap_legend_param = list(title = "Expression"),
        show_heatmap_legend = showHeatmapLegend,
        top_annotation = ComplexHeatmap::HeatmapAnnotation(
            cell = factor(zToPlot,
                levels = stringr::str_sort(unique(zToPlot),
                    numeric = TRUE)),
            show_legend = showTopAnnotationLegend,
            show_annotation_name = showTopAnnotationName,
            col = list(cell = ccols),
            simple_anno_size = grid::unit(annotationHeight, unit)),
        left_annotation = ComplexHeatmap::rowAnnotation(
            module = factor(yToPlot,
                levels = stringr::str_sort(unique(yToPlot),
                    numeric = TRUE)),
            show_legend = showLeftAnnotationLegend,
            show_annotation_name = showLeftAnnotationName,
            col = list(module = rcols),
            simple_anno_size = grid::unit(annotationWidth, unit)),
        ...)
    return(plt)
}
