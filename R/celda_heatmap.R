#' @title Plot celda Heatmap
#' @description Render a stylable heatmap of count data based on celda
#'  clustering results.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object
#'  returned by \link{celda_C}, \link{celda_G}, or \link{celda_CG}.
#' @param useAssay A string specifying which \link[SummarizedExperiment]{assay}
#'  slot to use. Default "counts".
#' @param altExpName The name for the \link[SingleCellExperiment]{altExp} slot
#'  to use. Default "featureSubset".
#' @param featureIx Integer vector. Select features for display in heatmap. If
#'  NULL, no subsetting will be performed. Default NULL. \strong{Only used for
#'  \code{sce} containing celda_C model result returned by \link{celda_C}.}
#' @param nfeatures Integer. Maximum number of features to select for each
#'  gene module. Default 25. \strong{Only used for \code{sce} containing
#'  celda_CG or celda_G model results returned by \link{celda_CG} or
#'  \link{celda_G}.}
#' @param ... Additional parameters passed to \link{plotHeatmap}.
#' @seealso `celdaTsne()` for generating 2-dimensional tSNE coordinates
#' @return list A list containing dendrogram information and the heatmap grob
#' @export
setGeneric("celdaHeatmap",
    function(sce, ...) {
        standardGeneric("celdaHeatmap")
    })


#' @rdname celdaHeatmap
#' @examples
#' data(sceCeldaCG)
#' celdaHeatmap(sceCeldaCG)
#' @export
setMethod("celdaHeatmap", signature(sce = "SingleCellExperiment"),
    function(sce, useAssay = "counts", altExpName = "featureSubset",
        featureIx = NULL, nfeatures = 25, ...) {

        aleExp <- SingleCellExperiment::altExp(sce, altExpName)

        if (celdaModel(sce) == "celda_C") {
            g <- .celdaHeatmapCelda_C(sce = sce,
                useAssay = useAssay,
                altExpName = altExpName,
                featureIx = featureIx,
                ...)
            return(g)
        } else if (celdaModel(sce) == "celda_CG") {
            g <- .celdaHeatmapCelda_CG(sce = sce,
                useAssay = useAssay,
                altExpName = altExpName,
                nfeatures = nfeatures,
                ...)
            return(g)
        } else if (celdaModel(sce) == "celda_G") {
            g <- .celdaHeatmapCelda_G(sce = sce,
                useAssay = useAssay,
                altExpName = altExpName,
                nfeatures = nfeatures,
                ...)
            return(g)
        } else {
            stop("S4Vectors::metadata(altExp(sce, altExpName))$",
                "celda_parameters$model must be",
                " one of 'celda_C', 'celda_G', or 'celda_CG'")
        }
    }
)


.celdaHeatmapCelda_C <- function(sce,
    useAssay, altExpName, featureIx, ...) {

    counts <- SummarizedExperiment::assay(sce, i = useAssay)
    counts <- .processCounts(counts)
    norm <- normalizeCounts(counts,
        normalize = "proportion",
        transformationFun = sqrt)

    if (is.null(featureIx)) {
        return(plotHeatmap(norm,
            z = celdaClusters(sce, altExpName = altExpName), ...))
    }

    return(plotHeatmap(norm[featureIx, ],
        z = celdaClusters(sce, altExpName = altExpName), ...))
}


.celdaHeatmapCelda_CG <- function(sce, useAssay, altExpName, nfeatures, ...) {
    counts <- SummarizedExperiment::assay(sce, i = useAssay)
    counts <- .processCounts(counts)
    fm <- factorizeMatrix(x = sce, useAssay = useAssay,
        altExpName = altExpName, type = "proportion")
    top <- topRank(fm$proportions$module, n = nfeatures)
    ix <- unlist(top$index)
    rn <- unlist(top$names)
    norm <- normalizeCounts(counts,
        normalize = "proportion",
        transformationFun = sqrt)
    plt <- plotHeatmap(norm[rn, ],
        z = celdaClusters(sce, altExpName = altExpName),
        y = celdaModules(sce)[ix],
        ...)
    invisible(plt)
}


.celdaHeatmapCelda_G <- function(sce, useAssay, altExpName, nfeatures, ...) {
    counts <- SummarizedExperiment::assay(sce, i = useAssay)
    counts <- .processCounts(counts)
    fm <- factorizeMatrix(x = sce, useAssay = useAssay,
        altExpName = altExpName, type = "proportion")
    top <- topRank(fm$proportions$module, n = nfeatures)
    ix <- unlist(top$index)
    rn <- unlist(top$names)
    norm <- normalizeCounts(counts,
        normalize = "proportion",
        transformationFun = sqrt)
    plt <- plotHeatmap(norm[rn, ], y = celdaModules(sce)[ix], ...)
    invisible(plt)
}
