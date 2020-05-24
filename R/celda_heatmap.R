#' @title Plot celda Heatmap
#' @description Render a stylable heatmap of count data based on celda
#'  clustering results.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object
#'  returned by \link{celda_C}, \link{celda_G}, or \link{celda_CG}.
#' @param useAssay A string specifying which \link[SummarizedExperiment]{assay}
#'  slot to use. Default "counts".
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
    function(sce, useAssay = "counts", featureIx = NULL, nfeatures = 25, ...) {
        if (celdaModel(sce) == "celda_C") {
            g <- .celdaHeatmapCelda_C(sce = sce,
                useAssay = useAssay,
                featureIx = featureIx,
                ...)
            return(g)
        } else if (celdaModel(sce) == "celda_CG") {
            g <- .celdaHeatmapCelda_CG(sce = sce,
                useAssay = useAssay,
                nfeatures = nfeatures,
                ...)
            return(g)
        } else if (celdaModel(sce) == "celda_G") {
            g <- .celdaHeatmapCelda_G(sce = sce,
                useAssay = useAssay,
                nfeatures = nfeatures,
                ...)
            return(g)
        } else {
            stop("S4Vectors::metadata(sce)$celda_parameters$model must be",
                " one of 'celda_C', 'celda_G', or 'celda_CG'")
        }
    }
)


.celdaHeatmapCelda_C <- function(sce,
    useAssay, featureIx, ...) {

    counts <- SummarizedExperiment::assay(sce, i = useAssay)
    norm <- normalizeCounts(counts,
        normalize = "proportion",
        transformationFun = sqrt)

    if (is.null(featureIx)) {
        return(plotHeatmap(norm,
            z = celdaClusters(sce), ...))
    }

    return(plotHeatmap(norm[featureIx, ],
        z = celdaClusters(sce), ...))
}


.celdaHeatmapCelda_CG <- function(sce, useAssay, nfeatures, ...) {
    counts <- SummarizedExperiment::assay(sce, i = useAssay)
    fm <- factorizeMatrix(x = sce, useAssay = useAssay, type = "proportion")
    top <- topRank(fm$proportions$module, n = nfeatures)
    ix <- unlist(top$index)
    norm <- normalizeCounts(counts,
        normalize = "proportion",
        transformationFun = sqrt)
    plt <- plotHeatmap(norm[ix, ],
        z = celdaClusters(sce),
        y = celdaModules(sce)[ix],
        ...
    )
    invisible(plt)
}


.celdaHeatmapCelda_G <- function(sce, useAssay, nfeatures, ...) {
    counts <- SummarizedExperiment::assay(sce, i = useAssay)
    fm <- factorizeMatrix(x = sce, useAssay = useAssay, type = "proportion")
    top <- celda::topRank(fm$proportions$module, n = nfeatures)
    ix <- unlist(top$index)
    norm <- normalizeCounts(counts,
        normalize = "proportion",
        transformationFun = sqrt
    )
    plt <- plotHeatmap(norm[ix, ], y = celdaModules(sce)[ix], ...)
    invisible(plt)
}
