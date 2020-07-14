#' @title Get or set the cell cluster labels from a celda
#'  \linkS4class{SingleCellExperiment} object or celda model
#'  object.
#' @description Return or set the cell cluster labels determined
#'  by \link{celda_C} or \link{celda_CG} models.
#' @param x Can be one of
#'  \itemize{
#'  \item A \linkS4class{SingleCellExperiment} object returned by
#'  \link{celda_C}, or \link{celda_CG}, with the matrix
#'  located in the \code{useAssay} assay slot. The
#'  a \link[SingleCellExperiment]{altExp} slot with name \code{altExpName} will
#'  be used. Rows represent features and columns represent cells.
#'  \item Celda model object.}
#' @param altExpName The name for the \link[SingleCellExperiment]{altExp} slot
#'  to use. Default "featureSubset".
#' @param value Character vector of cell cluster labels for replacements. Works
#'  only if \code{x} is a \linkS4class{SingleCellExperiment} object.
#' @return One of
#' \itemize{
#'  \item Character vector if \code{x} is a
#'  \linkS4class{SingleCellExperiment} object.
#'  Contains cell cluster labels for each cell in x.
#'  \item List if \code{x} is a celda model object. Contains cell cluster
#'  labels (for celda_C and celdaCG
#'  Models) and/or feature module labels (for celda_G and celdaCG Models).}
#' @export
setGeneric("celdaClusters",
    function(x, ...) {
        standardGeneric("celdaClusters")
    })


#' @rdname celdaClusters
#' @examples
#' data(sceCeldaCG)
#' celdaClusters(sceCeldaCG)
#' @export
setMethod("celdaClusters",
    signature(x = "SingleCellExperiment"),
    function(x, altExpName = "featureSubset") {
        altExp <- SingleCellExperiment::altExp(x, altExpName)
        return(SummarizedExperiment::colData(altExp)$celda_cell_cluster)
    })


#' @examples
#' data(celdaCGMod)
#' celdaClusters(celdaCGMod)
#' @rdname celdaClusters
#' @export
setMethod("celdaClusters",
    signature(x = "celdaModel"),
    function(x) {
        return(x@clusters)
    }
)


#' @rdname celdaClusters
#' @export
setGeneric("celdaClusters<-",
    function(x, altExpName, value) standardGeneric("celdaClusters<-")
)


#' @rdname celdaClusters
#' @export
setReplaceMethod("celdaClusters", signature(x = "SingleCellExperiment"),
    function(x, altExpName = "featureSubset", value) {
        altExp <- SingleCellExperiment::altExp(x, altExpName)
        SummarizedExperiment::colData(altExp)$celda_cell_cluster <- value
        SingleCellExperiment::altExp(x, altExpName) <- altExp
        return(x)
    })


#' @title Get or set the feature module labels from a celda
#'  \linkS4class{SingleCellExperiment} object.
#' @description Return or set the feature module cluster labels determined
#'  by \link{celda_G} or \link{celda_CG} models.
#' @param sce A \linkS4class{SingleCellExperiment} object returned by
#'  \link{celda_G}, or \link{celda_CG}, with the matrix
#'  located in the \code{useAssay} assay slot.
#'  Rows represent features and columns represent cells.
#' @param altExpName The name for the \link[SingleCellExperiment]{altExp} slot
#'  to use. Default "featureSubset".
#' @param value Character vector of feature module labels for replacements.
#'  Works only if \code{x} is a \linkS4class{SingleCellExperiment} object.
#' @return Character vector. Contains feature module labels for each
#'  feature in x.
#' @export
setGeneric("celdaModules",
    function(sce, altExpName) {
        standardGeneric("celdaModules")
    })


#' @rdname celdaModules
#' @examples
#' data(sceCeldaCG)
#' celdaModules(sceCeldaCG)
#' @export
setMethod("celdaModules",
    signature(sce = "SingleCellExperiment"),
    function(sce, altExpName = "featureSubset") {
        altExp <- SingleCellExperiment::altExp(sce, altExpName)
        return(SummarizedExperiment::rowData(altExp)$celda_feature_module)
    })


#' @rdname celdaModules
#' @export
setGeneric("celdaModules<-",
    function(sce, altExpName, value) standardGeneric("celdaModules<-")
)


#' @rdname celdaModules
#' @export
setReplaceMethod("celdaModules", signature(sce = "SingleCellExperiment"),
    function(sce, altExpName = "featureSubset", value) {
        altExp <- SingleCellExperiment::altExp(sce, altExpName)
        SummarizedExperiment::rowData(sce)$celda_feature_module <- value
        SingleCellExperiment::altExp(sce, altExpName) <- altExp
        return(sce)
    })


#' @title Get or set sample labels from a celda
#'  \linkS4class{SingleCellExperiment}  object
#' @description Return or set the sample labels for the cells in \code{sce}.
#' @param x Can be one of
#'  \itemize{
#'  \item A \linkS4class{SingleCellExperiment} object returned by
#'  \link{celda_C}, \link{celda_G}, or \link{celda_CG}, with the matrix
#'  located in the \code{useAssay} assay slot.
#'  Rows represent features and columns represent cells.
#'  \item A celda model object.}
#' @param altExpName The name for the \link[SingleCellExperiment]{altExp} slot
#'  to use. Default "featureSubset".
#' @param value Character vector of sample labels for replacements. Works
#'  only is \code{x} is a \linkS4class{SingleCellExperiment} object.
#' @return Character vector. Contains the sample labels provided at model
#'  creation, or those automatically generated by celda.
#' @export
setGeneric("sampleLabel",
    function(x, ...) {
        standardGeneric("sampleLabel")
    })


#' @rdname sampleLabel
#' @examples
#' data(sceCeldaCG)
#' sampleLabel(sceCeldaCG)
#' @export
setMethod("sampleLabel",
    signature(x = "SingleCellExperiment"),
    function(x, altExpName = "featureSubset") {
        altExp <- SingleCellExperiment::altExp(x, altExpName)
        return(SummarizedExperiment::colData(altExp)$celda_sample_label)
    })


#' @rdname sampleLabel
#' @export
setGeneric("sampleLabel<-",
    function(x, altExpName, value) standardGeneric("sampleLabel<-")
)
#' @rdname sampleLabel
#' @export
setReplaceMethod("sampleLabel", signature(x = "SingleCellExperiment"),
    function(x, altExpName = "featureSubset", value) {
        altExp <- SingleCellExperiment::altExp(x, altExpName)
        SummarizedExperiment::colData(altExp)$celda_sample_label <- value
        SingleCellExperiment::altExp(x, altExpName) <- altExp
        return(x)
    })


#' @examples
#' data(celdaCGMod)
#' sampleLabel(celdaCGMod)
#' @rdname sampleLabel
#' @export
setMethod("sampleLabel",
    signature(x = "celdaModel"),
    function(x) {
        x@sampleLabel
    }
)


#' @title Get parameter values provided for celdaModel creation
#' @description Retrieves the K/L, model priors (e.g. alpha, beta),
#'  and count matrix checksum parameters provided during the creation of the
#'  provided celdaModel.
#' @param celdaMod celdaModel. Options available in `celda::availableModels`.
#' @return List. Contains the model-specific parameters for the provided celda
#'  model object depending on its class.
#' @export
setGeneric(
    "params",
    function(celdaMod) {
        standardGeneric("params")
    }
)


#' @rdname params
#' @examples
#' data(celdaCGMod)
#' params(celdaCGMod)
#' @export
setMethod("params",
    signature = c(celdaMod = "celdaModel"),
    function(celdaMod) {
        celdaMod@params
    }
)


#' @title Get feature, cell and sample names from a celdaModel
#' @description Retrieves the row, column, and sample names used to generate
#'  a celdaModel.
#' @param celdaMod celdaModel. Options available in `celda::availableModels`.
#' @return List. Contains row, column, and sample character vectors
#'  corresponding to the values provided when the celdaModel was generated.
#' @export
setGeneric(
    "matrixNames",
    function(celdaMod) {
        standardGeneric("matrixNames")
    }
)


#' @rdname matrixNames
#' @examples
#' data(celdaCGMod)
#' matrixNames(celdaCGMod)
#' @export
setMethod("matrixNames",
    signature = c(celdaMod = "celdaModel"),
    function(celdaMod) {
        celdaMod@names
    }
)


#' @title Get run parameters from a celda model
#'  \code{SingleCellExperiment} or \code{celdaList} object
#' @description Returns details on the clustering parameters and model
#'  priors from the celdaList object when it was created.
#' @param x An object of class \linkS4class{SingleCellExperiment} or class
#'  \code{celdaList}.
#' @param altExpName The name for the \link[SingleCellExperiment]{altExp} slot
#'  to use. Default "featureSubset".
#' @return Data Frame. Contains details on the various K/L parameters, chain
#'  parameters, seed, and final log-likelihoods derived for each model in the
#'  provided celdaList.
#' @export
setGeneric("runParams",
    function(x, ...) {
        standardGeneric("runParams")
    }
)


#' @examples
#' data(sceCeldaCGGridSearch)
#' runParams(sceCeldaCGGridSearch)
#' @rdname runParams
#' @export
setMethod("runParams",
    signature(x = "SingleCellExperiment"),
    function(x, altExpName = "featureSubset") {
        altExp <- SingleCellExperiment::altExp(x, altExpName)
        return(altExp@metadata$celda_grid_search@runParams)
    }
)


#' @examples
#' data(celdaCGGridSearchRes)
#' runParams(celdaCGGridSearchRes)
#' @rdname runParams
#' @export
setMethod("runParams",
    signature(x = "celdaList"),
    function(x) {
        return(x@runParams)
    }
)


#' @title Get final celdaModels from a celda model \code{SCE} or celdaList
#'  object
#' @description Returns all celda models generated during a
#'  \link{celdaGridSearch} run.
#' @param x An object of class \linkS4class{SingleCellExperiment} or
#'  \code{celdaList}.
#' @param altExpName The name for the \link[SingleCellExperiment]{altExp} slot
#'  to use. Default "featureSubset".
#' @return List. Contains one celdaModel object for each of the parameters
#'  specified in \code{runParams(x)}.
#' @export
setGeneric(
    "resList",
    function(x, ...) {
        standardGeneric("resList")
    }
)


#' @examples
#' data(sceCeldaCGGridSearch)
#' celdaCGGridModels <- resList(sceCeldaCGGridSearch)
#' @rdname resList
#' @export
setMethod("resList",
    signature(x = "SingleCellExperiment"),
    function(x, altExpName = "featureSubset") {
        altExp <- SingleCellExperiment::altExp(x, altExpName)
        return(altExp@metadata$celda_grid_search@resList)
    }
)


#' @examples
#' data(celdaCGGridSearchRes)
#' celdaCGGridModels <- resList(celdaCGGridSearchRes)
#' @rdname resList
#' @export
setMethod("resList",
    signature(x = "celdaList"),
    function(x) {
        return(x@resList)
    }
)


#' @title Get celda model from a celda
#'  \link[SingleCellExperiment]{SingleCellExperiment} object
#' @description Return the celda model for \code{sce} returned by
#'  \link{celda_C}, \link{celda_G} or \link{celda_CG}.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object
#'  returned by \link{celda_C}, \link{celda_G}, or \link{celda_CG}.
#' @param altExpName The name for the \link[SingleCellExperiment]{altExp} slot
#'  to use. Default "featureSubset".
#' @return Character. The celda model. Can be one of "celda_C", "celda_G", or
#'  "celda_CG".
#' @examples
#' data(sceCeldaCG)
#' celdaModel(sceCeldaCG)
#' @export
setGeneric("celdaModel",
    function(sce, ...) {
        standardGeneric("celdaModel")
    })
#' @rdname celdaModel
#' @export
setMethod("celdaModel",
    signature(sce = "SingleCellExperiment"),
    function(sce, altExpName = "featureSubset") {

        if (!altExpName %in% SingleCellExperiment::altExpNames(sce)) {
            stop(altExpName, " not in 'altExpNames(sce)'. Run ",
                "selectFeatures(sce) first!")
        }

        altExp <- SingleCellExperiment::altExp(sce, altExpName)

        tryCatch(
            if (S4Vectors::metadata(altExp)$celda_parameters$model %in%
                    c("celda_C", "celda_G", "celda_CG")) {
                return(S4Vectors::metadata(altExp)$celda_parameters$model)
            } else {
                stop("S4Vectors::metadata(altExp(sce,",
                    " altExpName))$celda_parameters$model must be",
                    " one of 'celda_C', 'celda_G', or 'celda_CG'")
            },
            error = function(e) {
                message("S4Vectors::metadata(altExp(sce,",
                    " altExpName))$celda_parameters$model must",
                    " exist! Try running celda model (celda_C, celda_CG, or",
                    " celda_G) first.")
                stop(e)
            })
    })


#' @title Get perplexity for every model in a celdaList
#' @description Returns perplexity for each model in a celdaList as calculated
#'  by `perplexity().`
#' @param celdaList An object of class celdaList.
#' @return List. Contains one celdaModel object for each of the parameters
#'  specified in the `runParams()` of the provided celda list.
#' @examples
#' data(celdaCGGridSearchRes)
#' celdaCGGridModelPerplexities <- celdaPerplexity(celdaCGGridSearchRes)
#' @export
setGeneric(
    "celdaPerplexity",
    function(celdaList) {
        standardGeneric("celdaPerplexity")
    }
)
#' @title Get perplexity for every model in a celdaList
#' @description Returns perplexity for each model in a celdaList as calculated
#'  by `perplexity().`
#' @param celdaList An object of class celdaList.
#' @return List. Contains one celdaModel object for each of the parameters
#'  specified in the `runParams()` of the provided celda list.
#' @examples
#' data(celdaCGGridSearchRes)
#' celdaCGGridModelPerplexities <- celdaPerplexity(celdaCGGridSearchRes)
#' @export
setMethod("celdaPerplexity",
    signature = c(celdaList = "celdaList"),
    function(celdaList) {
        celdaList@perplexity
    }
)


#' @title Get the MD5 hash of the count matrix from the celdaList
#' @description Returns the MD5 hash of the count matrix used to generate the
#'  celdaList.
#' @param celdaList An object of class celdaList.
#' @return A character string of length 32 containing the MD5 digest of
#'  the count matrix.
#' @examples
#' data(celdaCGGridSearchRes)
#' countChecksum <- countChecksum(celdaCGGridSearchRes)
#' @export
setGeneric(
    "countChecksum",
    function(celdaList) {
        standardGeneric("countChecksum")
    }
)
#' @title Get the MD5 hash of the count matrix from the celdaList
#' @description Returns the MD5 hash of the count matrix used to generate the
#'  celdaList.
#' @param celdaList An object of class celdaList.
#' @return A character string of length 32 containing the MD5 digest of
#'  the count matrix.
#' @examples
#' data(celdaCGGridSearchRes)
#' countChecksum <- countChecksum(celdaCGGridSearchRes)
#' @export
setMethod("countChecksum",
    signature = c(celdaList = "celdaList"),
    function(celdaList) {
        celdaList@countChecksum
    }
)
