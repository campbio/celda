#' @title Phenotype enrichment plot
#' @description Plot the enrichment of phenotypes within cell clusters or
#'  \emph{vice versa} in a given celda model. Generates a mosaic plot or
#'  a balloon plot.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells.
#' @param celda.mod Celda model returned by the \link{celda_C},
#'  \code{\link{celda_CG}} or \code{\link{subsetCeldaList}} functions. Must be
#'  a \linkS4class{celda_C} or \linkS4class{celda_CG} object.
#' @param phenoLabels Character vector of phenotype labels for the cells. Must
#'  have the same length as the column number of count matrix \code{counts}.
#' @param phenotype The phenotypes to inspect. If specified, cells will be
#'  categorized into two groups based on whether their phenotype is in this
#'  phenotype group or not. By default, no phenotype grouping is performed.
#' @param cellCluster The cell clusters to inspect. If specified, cells will be
#'  categorized into two groups based on whether their assigned cell cluster is
#'  in this cell cluster group or not. By default, no cell cluster grouping is
#'  performed.
#' @param perspective The comparing perspective. One of "phenotype" or
#'  "cellcluster". Indicates whether to look at the enrichment of phenotypes in
#'  cell clusters or the enrichment of cell clusters in phenotypes. Default is
#'  "phenotype", which looks at the enrichment of phenotype in cell clusters.
#' @param plot The plot of visualizing the result. One of "mosaic" or "balloon".
#'  Uses the \code{\link[vcd]{mosaic}} function for plotting mosaic plot and the
#'  \code{\link[ggpubr]{ggballoonplot}} function for plotting the balloon plot.
#' @param ... Additional parameters passed to the \code{\link[vcd]{mosaic}} or
#'  \code{\link[ggpubr]{ggballoonplot}} function.
#' @return NULL if \code{plot = "mosaic"} and a ggplot object if
#'  \code{plot = "balloon"}
#' @examples
#' # Use sample labels as a surrogate for phenotype labels
#' plotPhenoEnrich(celda.CG.sim$counts,
#'     celda.CG.mod,
#'     celda.CG.sim$sample.label)
#'
#' # Look at the enrichment of "Sample_1" in cell cluster 1 in balloon plot
#' plotPhenoEnrich(counts = celda.CG.sim$counts,
#'     celda.mod = celda.CG.mod,
#'     phenoLabels = celda.CG.sim$sample.label,
#'     phenotype = "Sample_1",
#'     cellCluster = 1,
#'     plot = "balloon")
#'
#' # Look at the enrichment of cell cluster 1 & 2 in
#' # phenotype "Sample_3" & "Sample_4" in mosaic plot
#' plotPhenoEnrich(counts = celda.CG.sim$counts,
#'     celda.mod = celda.CG.mod,
#'     phenoLabels = celda.CG.sim$sample.label,
#'     phenotype = c("Sample_3", "Sample_4"),
#'     cellCluster = c(1, 2),
#'     perspective = "cellcluster")
#' @export
plotPhenoEnrich <- function(counts,
    celda.mod,
    phenoLabels,
    phenotype = NULL,
    cellCluster = NULL,
    perspective = c("phenotype", "cellcluster"),
    plot = c("mosaic", "balloon"),
    ...) {

    if (!(is(celdaModel, "celda_CG") || is(celdaModel, "celda_C"))) {
        stop("celda.mod must be of class 'celda_C' or 'celda_CG'")
    }

    if (length(phenoLabels) != ncol(counts)) {
        stop("phenoLabels must be the same length as the column number of",
            " counts")
    }

    if (!is.null(phenotype)) {
        if (!(all(phenotype %in% phenoLabels))) {
            stop("phenotype ",
                paste(phenotype, collapse = " "),
                " must exist in phenoLabels")
        }
    }

    if (!is.null(cellCluster)) {
        if (!(all(cellCluster %in% unique(celda::clusters(celda.mod)$z)))) {
            stop("Cell cluster ",
                paste(cellCluster, collapse = " "),
                " must exist in celda.mod@clusters$z")
        }
    }

    plot <- match.arg(plot)
    perspective <- match.arg(perspective)

    phenoLabels <- factor(phenoLabels, levels = unique(phenoLabels))

    if (is.null(phenotype) && is.null(cellCluster)) {
        conti <- table(phenoLabels,
            celda::clusters(celda.mod)$z,
            dnn = c("Phenotype", "Cell cluster"))
    } else if (!is.null(phenotype) && is.null(cellCluster)) {
        conti <- table(phenoLabels %in% phenotype,
            celda::clusters(celda.mod)$z,
            dnn = c(paste0("In phenotype ", paste(phenotype, collapse = " ")),
                "Cell cluster"))
    } else if (is.null(phenotype) && !is.null(cellCluster)) {
        conti <- table(phenoLabels,
            celda::clusters(celda.mod)$z %in% cellCluster,
            dnn = c("Phenotype",
                paste0("In cell cluster ",
                    paste(cellCluster,
                        collapse = " "))))
    } else {
        conti <- table(phenoLabels %in% phenotype,
            celda::clusters(celda.mod)$z %in% cellCluster,
            dnn = c(paste0("In phenotype ", paste(phenotype, collapse = " ")),
                paste0("In cell cluster ", paste(cellCluster,
                    collapse = " "))))
    }

    if (perspective == "cellcluster") {
        conti <- t(conti)
    }

    if (plot == "mosaic") {
        vcd::mosaic(t(conti),
            shade = TRUE,
            legend = TRUE,
            direction = "v",
            ...)
    } else if (plot == "balloon") {
        g <- ggpubr::ggballoonplot(as.data.frame(t(conti)),
            fill = "value",
            ...) +
            ggplot2::scale_fill_viridis_c(option = "C")
        return(g)
    }
}
