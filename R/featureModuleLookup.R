#' @title Obtain the gene module of a gene of interest
#' @description This function will output the corresponding feature module for
#'  a specified vector of genes from a celda_CG or celda_G celdaModel.
#'  \code{feature} must match the rownames of \code{sce}.
#' @param sce A \linkS4class{SingleCellExperiment} object returned by
#'  \link{celda_G}, or \link{celda_CG}, with the matrix
#'  located in the \code{useAssay} assay slot.
#'  Rows represent features and columns represent cells.
#' @param feature Character vector. Identify feature modules for the specified
#'  feature names. \code{feature} must match the rownames of \code{sce}.
#' @param exactMatch Logical. Whether to look for exactMatch of the gene name
#'  within counts matrix. Default \code{TRUE}.
#' @return List. Each entry corresponds to the feature module determined for
#' the provided features.
#' @export
setGeneric("featureModuleLookup",
    function(sce, ...) {standardGeneric("featureModuleLookup")})


#' @examples
#' data(sceCeldaCG)
#' module <- featureModuleLookup(sce = sceCeldaCG,
#'     feature = c("Gene_1", "Gene_XXX"))
#' @export
#' @rdname featureModuleLookup
setMethod("featureModuleLookup", signature(sce = "SingleCellExperiment"),
    function(sce,
        feature,
        exactMatch = TRUE) {

        if (celdaModel(sce) == "celda_CG") {
            featureList <- .featureModuleLookupCG(sce = sce, feature = feature,
                exactMatch = exactMatch)
        } else if (celdaModel(sce) == "celda_G") {
            featureList <- .featureModuleLookupG(sce = sce, feature = feature,
                exactMatch = exactMatch)
        } else {
            stop("S4Vectors::metadata(sce)$celda_parameters$model must be",
                " one of 'celda_G', or 'celda_CG'")
        }
        return(featureList)
    }
)


.featureModuleLookupCG <- function(sce,
    feature,
    exactMatch) {

    list <- list()
    if (!isTRUE(exactMatch)) {
        featureGrep <- c()
        for (x in seq(length(feature))) {
            featureGrep <- c(featureGrep, rownames(sce)[grep(
                feature[x],
                rownames(sce)
            )])
        }
        feature <- featureGrep
    }
    for (x in seq(length(feature))) {
        if (feature[x] %in% rownames(sce)) {
            list[x] <- celdaModules(sce)[which(rownames(sce) ==
                    feature[x])]
        } else {
            list[x] <- paste0(
                "No feature was identified matching '",
                feature[x],
                "'."
            )
        }
    }
    names(list) <- feature
    return(list)
}


.featureModuleLookupG <- function(sce, feature, exactMatch) {
    if (!isTRUE(exactMatch)) {
        feature <- unlist(lapply(
            seq(length(feature)),
            function(x) {
                rownames(sce)[grep(feature[x], rownames(sce))]
            }
        ))
    }

    featList <- lapply(
        seq(length(feature)),
        function(x) {
            if (feature[x] %in% rownames(sce)) {
                return(celdaModules(sce)[which(rownames(sce) ==
                        feature[x])])
            } else {
                return(paste0(
                    "No feature was identified matching '",
                    feature[x],
                    "'."
                ))
            }
        }
    )
    names(featList) <- feature
    return(featList)
}
