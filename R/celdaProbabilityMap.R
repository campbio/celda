#' @title Probability map for a celda model
#' @description Renders probability and relative expression heatmaps to
#'  visualize the relationship between features and cell populations (or cell
#'  populations and samples).
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object
#'  returned by \link{celda_C}, \link{celda_G}, or \link{celda_CG}.
#' @param useAssay A string specifying which \link{assay}
#'  slot to use. Default "counts".
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default "featureSubset".
#' @param level Character. One of "cellPopulation" or "Sample".
#'  "cellPopulation" will display the absolute probabilities and relative
#'  normalized expression of each module in each cell population.
#'  \strong{\code{level = "cellPopulation"} only works for celda_CG \code{sce}
#'  objects}. "sample" will display the absolute probabilities and relative
#'  normalized abundance of each cell population in each sample. Default
#'  "cellPopulation".
#' @param ncols The number of colors (>1) to be in the color palette of
#'  the absolute probability heatmap.
#' @param col2 Passed to \code{col} argument of \link[ComplexHeatmap]{Heatmap}.
#'  Set color boundaries and colors for the relative expression heatmap.
#' @param title1 Passed to \code{column_title} argument of
#'  \link[ComplexHeatmap]{Heatmap}. Figure title for the absolute probability
#'  heatmap.
#' @param title2 Passed to \code{column_title} argument of
#'  \link[ComplexHeatmap]{Heatmap}. Figure title for the relative expression
#'  heatmap.
#' @param showColumnNames Passed to \code{show_column_names} argument of
#'  \link[ComplexHeatmap]{Heatmap}. Show column names.
#' @param showRowNames Passed to \code{show_row_names} argument of
#'  \link[ComplexHeatmap]{Heatmap}. Show row names.
#' @param rowNamesgp Passed to \code{row_names_gp} argument of
#'  \link[ComplexHeatmap]{Heatmap}. Set row name font.
#' @param colNamesgp Passed to \code{column_names_gp} argument of
#'  \link[ComplexHeatmap]{Heatmap}. Set column name font.
#' @param clusterRows Passed to \code{cluster_rows} argument of
#'  \link[ComplexHeatmap]{Heatmap}. Cluster rows.
#' @param clusterColumns Passed to \code{cluster_columns} argument of
#'  \link[ComplexHeatmap]{Heatmap}. Cluster columns.
#' @param showHeatmapLegend Passed to \code{show_heatmap_legend} argument of
#'  \link[ComplexHeatmap]{Heatmap}. Show heatmap legend.
#' @param heatmapLegendParam Passed to \code{heatmap_legend_param} argument of
#'  \link[ComplexHeatmap]{Heatmap}. Heatmap legend parameters.
#' @param ... Additional parameters passed to \link[ComplexHeatmap]{Heatmap}.
#' @seealso \link{celda_C} for clustering cells. \link{celda_CG} for
#'  clustering features and cells
#' @return A \link[ComplexHeatmap]{HeatmapList} object containing 2
#'  \link[ComplexHeatmap]{Heatmap-class} objects
#' @export
setGeneric("celdaProbabilityMap",
    function(sce, ...) {
        standardGeneric("celdaProbabilityMap")
    })


#' @rdname celdaProbabilityMap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @examples
#' data(sceCeldaCG)
#' celdaProbabilityMap(sceCeldaCG)
#' @export
setMethod("celdaProbabilityMap", signature(sce = "SingleCellExperiment"),
    function(sce,
        useAssay = "counts",
        altExpName = "featureSubset",
        level = c("cellPopulation", "sample"),
        ncols = 100,
        col2 = circlize::colorRamp2(c(-2, 0, 2),
            c("#1E90FF", "#FFFFFF", "#CD2626")),
        title1 = "Absolute probability",
        title2 = "Relative expression",
        showColumnNames = TRUE,
        showRowNames = TRUE,
        rowNamesgp = grid::gpar(fontsize = 8),
        colNamesgp = grid::gpar(fontsize = 12),
        clusterRows = FALSE,
        clusterColumns = FALSE,
        showHeatmapLegend = TRUE,
        heatmapLegendParam = list(title = NULL,
            legend_height = grid::unit(6, "cm"))) {

        altExp <- SingleCellExperiment::altExp(sce, altExpName)
        level <- match.arg(level)
        if (celdaModel(sce, altExpName = altExpName) == "celda_C") {
            if (level == "cellPopulation") {
                warning("'level' has been set to 'sample'")
            }
            pm <- .celdaProbabilityMapC(sce = altExp,
                useAssay = useAssay,
                level = "sample",
                ncols = ncols,
                col2 = col2,
                title1 = title1,
                title2 = title2,
                showColumnNames = showColumnNames,
                showRowNames = showRowNames,
                rowNamesgp = rowNamesgp,
                colNamesgp = colNamesgp,
                clusterRows = clusterRows,
                clusterColumns = clusterColumns,
                showHeatmapLegend = showHeatmapLegend,
                heatmapLegendParam = heatmapLegendParam)
        } else if (celdaModel(sce, altExpName = altExpName) == "celda_CG") {
            pm <- .celdaProbabilityMapCG(sce = altExp,
                useAssay = useAssay,
                level = level,
                ncols = ncols,
                col2 = col2,
                title1 = title1,
                title2 = title2,
                showColumnNames = showColumnNames,
                showRowNames = showRowNames,
                rowNamesgp = rowNamesgp,
                colNamesgp = colNamesgp,
                clusterRows = clusterRows,
                clusterColumns = clusterColumns,
                showHeatmapLegend = showHeatmapLegend,
                heatmapLegendParam = heatmapLegendParam)
        } else {
            stop("S4Vectors::metadata(altExp(sce,",
                " altExpName))$celda_parameters$model must be",
                " one of 'celda_C', or 'celda_CG'!")
        }
        return(pm)
    }
)


.celdaProbabilityMapC <- function(sce,
    useAssay,
    level,
    ncols,
    col2,
    title1,
    title2,
    showColumnNames,
    showRowNames,
    rowNamesgp,
    colNamesgp,
    clusterRows,
    clusterColumns,
    showHeatmapLegend,
    heatmapLegendParam) {

    counts <- SummarizedExperiment::assay(sce, i = useAssay)
    counts <- .processCounts(counts)

    zInclude <- which(tabulate(SummarizedExperiment::colData(
        sce)$celda_cell_cluster,
        S4Vectors::metadata(sce)$celda_parameters$K) > 0)

    factorized <- .factorizeMatrixCelda_C(sce, useAssay = useAssay,
        type = "proportion")

    samp <- factorized$proportions$sample[zInclude, , drop = FALSE]
    col1 <- grDevices::colorRampPalette(c("white",
        "blue",
        "midnightblue",
        "springgreen4",
        "yellowgreen",
        "yellow",
        "orange",
        "red"))(100)
    breaks <- seq(0, 1, length.out = length(col1))

    g1 <- ComplexHeatmap::Heatmap(matrix = samp,
        col = circlize::colorRamp2(breaks, col1),
        column_title = title1,
        show_column_names = showColumnNames,
        show_row_names = showRowNames,
        row_names_gp = rowNamesgp,
        column_names_gp = colNamesgp,
        cluster_rows = clusterRows,
        cluster_columns = clusterColumns,
        show_heatmap_legend = showHeatmapLegend,
        heatmap_legend_param = heatmapLegendParam)

    if (ncol(samp) > 1) {
        sampNorm <- normalizeCounts(samp,
            normalize = "proportion",
            transformationFun = sqrt,
            scaleFun = base::scale)

        g2 <- ComplexHeatmap::Heatmap(matrix = sampNorm,
            col = col2,
            column_title = title2,
            show_column_names = showColumnNames,
            show_row_names = showRowNames,
            row_names_gp = rowNamesgp,
            column_names_gp = colNamesgp,
            cluster_rows = clusterRows,
            cluster_columns = clusterColumns,
            show_heatmap_legend = showHeatmapLegend,
            heatmap_legend_param = heatmapLegendParam)
        return(g1 + g2)
    } else {
        return(g1)
    }
}


.celdaProbabilityMapCG <- function(sce,
    useAssay,
    level,
    ncols,
    col2,
    title1,
    title2,
    showColumnNames,
    showRowNames,
    rowNamesgp,
    colNamesgp,
    clusterRows,
    clusterColumns,
    showHeatmapLegend,
    heatmapLegendParam) {

    counts <- SummarizedExperiment::assay(sce, i = useAssay)
    counts <- .processCounts(counts)

    factorized <- .factorizeMatrixCelda_CG(sce, useAssay,
        type = c("counts", "proportion"))
    zInclude <- which(tabulate(SummarizedExperiment::colData(
        sce)$celda_cell_cluster,
        S4Vectors::metadata(sce)$celda_parameters$K) > 0)
    yInclude <- which(tabulate(SummarizedExperiment::rowData(
        sce)$celda_feature_module,
        S4Vectors::metadata(sce)$celda_parameters$L) > 0)

    if (level == "cellPopulation") {
        pop <- factorized$proportions$cellPopulation[yInclude,
            zInclude,
            drop = FALSE]
        popNorm <- normalizeCounts(pop,
            normalize = "proportion",
            transformationFun = sqrt,
            scaleFun = base::scale)

        percentile9 <- round(stats::quantile(pop, .9), digits = 2) * 100
        cols11 <- grDevices::colorRampPalette(c("white",
            RColorBrewer::brewer.pal(n = 9, name = "Blues")))(percentile9)
        cols12 <- grDevices::colorRampPalette(c("midnightblue",
            c("springgreen4", "Yellowgreen", "Yellow", "Orange",
                "Red")))(ncols - percentile9)
        col1 <- c(cols11, cols12)
        breaks <- seq(0, 1, length.out = length(col1))

        g1 <- ComplexHeatmap::Heatmap(matrix = pop,
            col = circlize::colorRamp2(breaks, col1),
            column_title = title1,
            show_column_names = showColumnNames,
            show_row_names = showRowNames,
            row_names_gp = rowNamesgp,
            column_names_gp = colNamesgp,
            cluster_rows = clusterRows,
            cluster_columns = clusterColumns,
            show_heatmap_legend = showHeatmapLegend,
            heatmap_legend_param = heatmapLegendParam)
        g2 <- ComplexHeatmap::Heatmap(matrix = popNorm,
            col = col2,
            column_title = title2,
            show_column_names = showColumnNames,
            show_row_names = showRowNames,
            row_names_gp = rowNamesgp,
            column_names_gp = colNamesgp,
            cluster_rows = clusterRows,
            cluster_columns = clusterColumns,
            show_heatmap_legend = showHeatmapLegend,
            heatmap_legend_param = heatmapLegendParam)
        return(g1 + g2)
    } else {
        samp <- factorized$proportions$sample
        col1 <- grDevices::colorRampPalette(c(
            "white",
            "blue",
            "#08306B",
            "#006D2C",
            "yellowgreen",
            "yellow",
            "orange",
            "red"
        ))(100)
        breaks <- seq(0, 1, length.out = length(col1))

        g1 <- ComplexHeatmap::Heatmap(matrix = samp,
            col = circlize::colorRamp2(breaks, col1),
            column_title = title1,
            show_column_names = showColumnNames,
            show_row_names = showRowNames,
            row_names_gp = rowNamesgp,
            column_names_gp = colNamesgp,
            cluster_rows = clusterRows,
            cluster_columns = clusterColumns,
            show_heatmap_legend = showHeatmapLegend,
            heatmap_legend_param = heatmapLegendParam)

        if (ncol(samp) > 1) {
            sampNorm <- normalizeCounts(factorized$counts$sample,
                normalize = "proportion",
                transformationFun = sqrt,
                scaleFun = base::scale)
            g2 <- ComplexHeatmap::Heatmap(matrix = sampNorm,
                col = col2,
                column_title = title2,
                show_column_names = showColumnNames,
                show_row_names = showRowNames,
                row_names_gp = rowNamesgp,
                column_names_gp = colNamesgp,
                cluster_rows = clusterRows,
                cluster_columns = clusterColumns,
                show_heatmap_legend = showHeatmapLegend,
                heatmap_legend_param = heatmapLegendParam)
            return(g1 + g2)
        } else {
            return(g1 + g2)
        }
    }
}
