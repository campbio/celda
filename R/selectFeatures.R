#' @title Simple feature selection by feature counts
#' @description A simple heuristic feature selection procedure.
#'  Select features with at least \code{minCount} counts
#'  in at least \code{minCell} cells. A \linkS4class{SingleCellExperiment}
#'  object with subset features will be stored in the
#'  \link{altExp} slot with name \code{altExpName}.
#'  The name of the \code{assay} slot in \link{altExp}
#'  will be the same as \code{useAssay}.
#' @param x A numeric \link{matrix} of counts or a
#'  \linkS4class{SingleCellExperiment}
#'  with the matrix located in the assay slot under \code{useAssay}.
#'  Rows represent features and columns represent cells.
#' @param minCount Minimum number of counts required for feature selection.
#' @param minCell Minimum number of cells required for feature selection.
#' @param useAssay A string specifying the name of the
#'  \link{assay} slot to use. Default "counts".
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default "featureSubset".
#' @param ... Ignored. Placeholder to prevent check warning.
#' @return A \linkS4class{SingleCellExperiment} object with a
#'  \code{altExpName} \link{altExp} slot. Function
#'  parameter settings are stored in the \link{metadata}
#'  \code{"select_features"} slot.
#' @export
setGeneric("selectFeatures", function(x, ...) {
    standardGeneric("selectFeatures")})


#' @rdname selectFeatures
#' @examples
#' data(sceCeldaCG)
#' sce <- selectFeatures(sceCeldaCG)
#' @export
setMethod("selectFeatures",
    signature(x = "SingleCellExperiment"),
    function(x,
        minCount = 3,
        minCell = 3,
        useAssay = "counts",
        altExpName = "featureSubset") {

        assay <- SummarizedExperiment::assay(x, i = useAssay)
        sceSubset <- x[Matrix::rowSums(assay >= minCount) >= minCell, ]
        SingleCellExperiment::altExp(x, altExpName) <- sceSubset

        S4Vectors::metadata(x)[["select_features"]] <- list(
            xClass = "SingleCellExperiment",
            minCount = minCount,
            minCell = minCell,
            useAssay = useAssay,
            altExpName = altExpName)
        return(x)
    }
)


#' @rdname selectFeatures
#' @examples
#' data(celdaCGSim)
#' sce <- selectFeatures(celdaCGSim$counts)
#' @export
setMethod("selectFeatures",
    signature(x = "matrix"),
    function(x,
        minCount = 3,
        minCell = 3,
        useAssay = "counts",
        altExpName = "featureSubset") {

        ls <- list()
        ls[[useAssay]] <- x
        sce <- SingleCellExperiment::SingleCellExperiment(assays = ls)
        sceSubset <- sce[rowSums(x >= minCount) >= minCell, ]
        SingleCellExperiment::altExp(sce, altExpName) <- sceSubset
        S4Vectors::metadata(sce)[["select_features"]] <- list(
            xClass = "matrix",
            minCount = minCount,
            minCell = minCell,
            useAssay = useAssay,
            altExpName = altExpName)
        return(sce)
    }
)
