#' @title Probability map for a celda model
#' @description Renders probability and relative expression heatmaps to
#'  visualize the relationship between features and cell populations (or cell
#'  populations and samples).
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object
#'  returned by \link{celda_C}, \link{celda_G}, or \link{celda_CG}.
#' @param useAssay A string specifying which \link[SummarizedExperiment]{assay}
#'  slot to use. Default "counts".
#' @param level Character. One of "cellPopulation" or "Sample".
#'  "cellPopulation" will display the absolute probabilities and relative
#'  normalized expression of each module in each cell population.
#'  \strong{\code{level = "cellPopulation"} only works for celda_CG \code{sce}
#'  objects}. "sample" will display the absolute probabilities and relative
#'  normalized abundance of each cell population in each sample. Default
#'  "cellPopulation".
#' @param ... Additional parameters.
#' @seealso \link{celda_C} for clustering cells. \link{celda_CG} for
#'  clustering features and cells
#' @return A grob containing the specified plots
#' @export
setGeneric("celdaProbabilityMap",
    function(sce, ...) {
        standardGeneric("celdaProbabilityMap")
    })


#' @rdname celdaProbabilityMap
#' @importFrom gridExtra grid.arrange
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @examples
#' data(sceCeldaCG)
#' celdaProbabilityMap(sceCeldaCG)
#' @export
setMethod("celdaProbabilityMap", signature(sce = "SingleCellExperiment"),
    function(sce, useAssay = "counts", level = c("cellPopulation", "sample")) {
        level <- match.arg(level)
        if (celdaModel(sce) == "celda_C") {
            if (level == "cellPopulation") {
                warning("'level' has been set to 'sample'")
            }
            pm <- .celdaProbabilityMapC(sce = sce, useAssay = useAssay,
                level = "sample")
        } else if (celdaModel(sce) == "celda_CG") {
            pm <- .celdaProbabilityMapCG(sce = sce, useAssay = useAssay,
                level = level)
        } else {
            stop("S4Vectors::metadata(sce)$celda_parameters$model must be",
                " one of 'celda_C', or 'celda_CG'!")
        }
        return(pm)
    }
)


.celdaProbabilityMapC <- function(sce, useAssay, level) {
    counts <- SummarizedExperiment::assay(sce, i = useAssay)
    counts <- .processCounts(counts)

    zInclude <- which(tabulate(celdaClusters(sce),
        S4Vectors::metadata(sce)$celda_parameters$K) > 0)

    factorized <- factorizeMatrix(x = sce, useAssay = useAssay)

    samp <- factorized$proportions$sample[zInclude, , drop = FALSE]
    col <- grDevices::colorRampPalette(c("white",
        "blue",
        "#08306B",
        "#006D2C",
        "yellowgreen",
        "yellow",
        "orange",
        "red"))(100)
    breaks <- seq(0, 1, length.out = length(col))
    g1 <- plotHeatmap(samp,
        colorScheme = "sequential",
        scaleRow = NULL,
        clusterCell = FALSE,
        clusterFeature = FALSE,
        showNamesCell = TRUE,
        showNamesFeature = TRUE,
        breaks = breaks,
        col = col,
        main = "Absolute Probability",
        silent = TRUE)

    if (ncol(samp) > 1) {
        sampNorm <- normalizeCounts(samp,
            normalize = "proportion",
            transformationFun = sqrt,
            scaleFun = base::scale)
        g2 <- plotHeatmap(sampNorm,
            colorScheme = "divergent",
            clusterCell = FALSE,
            clusterFeature = FALSE,
            showNamesCell = TRUE,
            showNamesFeature = TRUE,
            main = "Relative Abundance",
            silent = TRUE)
        return(gridExtra::grid.arrange(g1$gtable, g2$gtable, ncol = 2))
    } else {
        return(gridExtra::grid.arrange(g1$gtable))
    }
}


.celdaProbabilityMapCG <- function(sce, useAssay, level) {
    counts <- SummarizedExperiment::assay(sce, i = useAssay)
    counts <- .processCounts(counts)

    factorized <- factorizeMatrix(sce, useAssay)
    zInclude <- which(tabulate(celdaClusters(sce),
        S4Vectors::metadata(sce)$celda_parameters$K) > 0)
    yInclude <- which(tabulate(celdaModules(sce),
        S4Vectors::metadata(sce)$celda_parameters$L) > 0)

    if (level == "cellPopulation") {
        pop <- factorized$proportions$cellPopulation[yInclude,
            zInclude,
            drop = FALSE
            ]
        popNorm <- normalizeCounts(pop,
            normalize = "proportion",
            transformationFun = sqrt,
            scaleFun = base::scale
        )

        percentile9 <- round(stats::quantile(pop, .9), digits = 2) * 100
        col1 <- grDevices::colorRampPalette(c(
            "#FFFFFF",
            RColorBrewer::brewer.pal(n = 9, name = "Blues")
        ))(percentile9)
        col2 <- grDevices::colorRampPalette(c(
            "#08306B",
            c(
                "#006D2C", "Yellowgreen", "Yellow", "Orange",
                "Red"
            )
        ))(100 - percentile9)
        col <- c(col1, col2)
        breaks <- seq(0, 1, length.out = length(col))

        g1 <- plotHeatmap(pop,
            colorScheme = "sequential",
            scaleRow = NULL,
            clusterCell = FALSE,
            clusterFeature = FALSE,
            showNamesCell = TRUE,
            showNamesFeature = TRUE,
            breaks = breaks,
            col = col,
            main = "Absolute Probability",
            silent = TRUE
        )
        g2 <- plotHeatmap(popNorm,
            colorScheme = "divergent",
            clusterCell = FALSE,
            clusterFeature = FALSE,
            showNamesCell = TRUE,
            showNamesFeature = TRUE,
            main = "Relative Expression",
            silent = TRUE
        )
        gridExtra::grid.arrange(g1$gtable, g2$gtable, ncol = 2)
    } else {
        samp <- factorized$proportions$sample
        col <- grDevices::colorRampPalette(c(
            "white",
            "blue",
            "#08306B",
            "#006D2C",
            "yellowgreen",
            "yellow",
            "orange",
            "red"
        ))(100)
        breaks <- seq(0, 1, length.out = length(col))
        g1 <- plotHeatmap(samp,
            colorScheme = "sequential",
            scaleRow = NULL,
            clusterCell = FALSE,
            clusterFeature = FALSE,
            showNamesCell = TRUE,
            showNamesFeature = TRUE,
            breaks = breaks,
            col = col,
            main = "Absolute Probability",
            silent = TRUE
        )

        if (ncol(samp) > 1) {
            sampNorm <- normalizeCounts(factorized$counts$sample,
                normalize = "proportion",
                transformationFun = sqrt,
                scaleFun = base::scale
            )
            g2 <- plotHeatmap(sampNorm,
                colorScheme = "divergent",
                clusterCell = FALSE,
                clusterFeature = FALSE,
                showNamesCell = TRUE,
                showNamesFeature = TRUE,
                main = "Relative Abundance",
                silent = TRUE
            )
            gridExtra::grid.arrange(g1$gtable, g2$gtable, ncol = 2)
        } else {
            gridExtra::grid.arrange(g1$gtable)
        }
    }
}

