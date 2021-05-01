#' @title Plot celda Heatmap
#' @description Render a stylable heatmap of count data based on celda
#'  clustering results.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object
#'  returned by \link{celda_C}, \link{celda_G}, or \link{celda_CG}.
#' @param useAssay A string specifying which \link{assay}
#'  slot to use. Default "counts".
#' @param altExpName The name for the \link{altExp} slot
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

        counts <- SummarizedExperiment::assay(sce, i = useAssay)
        counts <- .processCounts(counts)
        model <- celdaModel(sce, altExpName = altExpName)

        if (model == "celda_C") {
            z <- celdaClusters(sce, altExpName = altExpName)
            g <- .celdaHeatmapCelda_C(
                counts = counts,
                z = z,
                featureIx = featureIx,
                ...)
        } else if (model == "celda_CG") {
            fm <- factorizeMatrix(x = sce, useAssay = useAssay,
                altExpName = altExpName, type = "proportion")
            z <- celdaClusters(sce, altExpName = altExpName)
            y <- celdaModules(sce, altExpName = altExpName)
            g <- .celdaHeatmapCelda_CG(
                counts = counts,
                fm = fm,
                z = z,
                y = y,
                nfeatures = nfeatures,
                ...)
        } else if (model == "celda_G") {
            fm <- factorizeMatrix(x = sce, useAssay = useAssay,
                altExpName = altExpName, type = "proportion")
            y <- celdaModules(sce, altExpName = altExpName)
            g <- .celdaHeatmapCelda_G(counts,
                fm,
                y,
                nfeatures = nfeatures,
                ...)
        } else {
            stop("S4Vectors::metadata(altExp(sce, altExpName))$",
                "celda_parameters$model must be",
                " one of 'celda_C', 'celda_G', or 'celda_CG'")
        }
        return(g)
    }
)


.celdaHeatmapCelda_C <- function(
    counts,
    z,
    featureIx, ...) {

    norm <- normalizeCounts(counts,
        normalize = "proportion",
        transformationFun = sqrt)

    if (is.null(featureIx)) {
        return(plotHeatmap(norm, z, ...))
    }
    return(plotHeatmap(norm[featureIx, ], z = z, ...))
}


.celdaHeatmapCelda_CG <- function(
    counts = counts,
    fm = fm,
    z = z,
    y = y,
    nfeatures,
    ...) {

    top <- topRank(fm$proportions$module, n = nfeatures)
    ix <- unlist(top$index)
    rn <- unlist(top$names)
    norm <- normalizeCounts(counts,
        normalize = "proportion",
        transformationFun = sqrt)
    plt <- plotHeatmap(norm[rn, ],
        z = z,
        y = y[ix],
        ...)
    return(plt)
}


.celdaHeatmapCelda_G <- function(counts, fm, y, nfeatures, ...) {
    top <- topRank(fm$proportions$module, n = nfeatures)
    ix <- unlist(top$index)
    rn <- unlist(top$names)
    norm <- normalizeCounts(counts,
        normalize = "proportion",
        transformationFun = sqrt)
    plt <- plotHeatmap(norm[rn, ], y = y[ix], ...)
    return(plt)
}
