#' @title Obtain the gene module of a gene of interest
#' @description This function will output the corresponding feature module for
#'  a specified vector of genes from a celda_CG or celda_G \code{celdaModel}.
#'  \code{features} must match the rownames of \code{sce}.
#' @param sce A \linkS4class{SingleCellExperiment} object returned by
#'  \link{celda_G}, or \link{celda_CG}, with the matrix
#'  located in the \code{useAssay} assay slot.
#'  Rows represent features and columns represent cells.
#' @param features Character vector. Identify feature modules for the specified
#'  feature names. \code{feature} must match the rownames of \code{sce}.
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default "featureSubset".
#' @param exactMatch Logical. Whether to look for exactMatch of the gene name
#'  within counts matrix. Default \code{TRUE}.
#' @param by Character. Where to search for \code{features} in the sce object.
#' If set to \code{"rownames"} then the features will be searched for among
#' rownames(sce). This can also be set to one of the \code{colnames} of
#' rowData(sce). Default \code{"rownames"}.
#' @return Numeric vector containing the module numbers for each feature. If
#' the feature was not found, then an \code{NA} value will be returned in that
#' position. If no features were found, then an error will be given.
#' @export
setGeneric("featureModuleLookup",
    function(sce,
        features,
        altExpName = "featureSubset",
        exactMatch = TRUE,
        by = "rownames") {

        standardGeneric("featureModuleLookup")})


#' @examples
#' data(sceCeldaCG)
#' module <- featureModuleLookup(sce = sceCeldaCG,
#'     features = c("Gene_1", "Gene_XXX"))
#' @export
#' @rdname featureModuleLookup
setMethod("featureModuleLookup", signature(sce = "SingleCellExperiment"),
    function(sce,
        features,
        altExpName = "featureSubset",
        exactMatch = TRUE,
        by = "rownames") {

        modules <- as.numeric(celdaModules(sce, altExpName = altExpName))

        if (celdaModel(sce, altExpName = altExpName) %in%
                c("celda_CG", "celda_G")) {
          altExp <- SingleCellExperiment::altExp(sce, altExpName)
          featureIndex <- retrieveFeatureIndex(features, x = altExp,
                exactMatch = exactMatch, by = by)
          featureModules <- modules[featureIndex]
          names(featureModules) <- features
        } else {
            stop("S4Vectors::metadata(altExp(sce, altExpName))$",
                "celda_parameters$model must be",
                " one of 'celda_G', or 'celda_CG'")
        }
        return(featureModules)
    }
)
