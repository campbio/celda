setClass(
  "celdaModel",
  representation(
    params = "list",
    # K, L, model priors, checksum
    names = "list",
    completeLogLik = "numeric",
    finalLogLik = "numeric",
    clusters = "list"
  )
) # z and or y


#' @title Get parameter values provided for celdaModel creation
#' @description Retrieves the K/L, model priors (e.g. alpha, beta),
#'  and count matrix checksum parameters provided during the creation of the
#'  provided celdaModel.
#' @param celdaMod celdaModel. Options available in `celda::availableModels`.
#' @return List. Contains the model-specific parameters for the provided celda
#'  model object depending on its class.
#' @examples
#' data(celdaCGMod)
#' params(celdaCGMod)
#' @export
setGeneric(
  "params",
  function(celdaMod) {
    standardGeneric("params")
  }
)
#' @title Get parameter values provided for celdaModel creation
#' @description Retrieves the K/L, model priors (e.g. alpha, beta),
#'  and count matrix checksum parameters provided during the creation of the
#'  provided celdaModel.
#' @param celdaMod celdaModel. Options available in `celda::availableModels`.
#' @return List. Contains the model-specific parameters for the provided celda
#'  model object depending on its class.
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
#' @examples
#' data(celdaCGMod)
#' matrixNames(celdaCGMod)
#' @export
setGeneric(
  "matrixNames",
  function(celdaMod) {
    standardGeneric("matrixNames")
  }
)
#' @title Get feature, cell and sample names from a celdaModel
#' @description Retrieves the row, column, and sample names used to generate a
#'  celdaModel.
#' @param celdaMod celdaModel. Options available in `celda::availableModels`.
#' @return List. Contains row, column, and sample character vectors
#'  corresponding to the values provided when the celdaModel was generated.
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


#' @title Get log-likelihood history
#' @description Retrieves the complete log-likelihood from all iterations of
#'  Gibbs sampling used to generate a celdaModel.
#' @param celdaMod celdaModel. Options available in `celda::availableModels`.
#' @return Numeric. The log-likelihood at each step of Gibbs sampling used to
#'  generate the model.
#' @examples
#' data(celdaCGMod)
#' logLikelihoodHistory(celdaCGMod)
#' @export
setGeneric(
  "logLikelihoodHistory",
  function(celdaMod) {
    standardGeneric("logLikelihoodHistory")
  }
)
#' @title Get log-likelihood history
#' @description Retrieves the complete log-likelihood from all iterations of
#'  Gibbs sampling used to generate a celdaModel.
#' @param celdaMod celdaModel. Options available in `celda::availableModels`.
#' @return Numeric. The log-likelihood at each step of Gibbs sampling used to
#'  generate the model.
#' @examples
#' data(celdaCGMod)
#' logLikelihoodHistory(celdaCGMod)
#' @export
setMethod("logLikelihoodHistory",
  signature = c(celdaMod = "celdaModel"),
  function(celdaMod) {
    celdaMod@completeLogLik
  }
)


#' @title Get the log-likelihood
#' @description Retrieves the final log-likelihood from all iterations of Gibbs
#'  sampling used to generate a celdaModel.
#' @return Numeric. The log-likelihood at the final step of Gibbs sampling used
#'  to generate the model.
#' @param celdaMod A celdaModel object of class celda_C, celda_G, or celda_CG.
#' @examples
#' data(celdaCGMod)
#' bestLogLikelihood(celdaCGMod)
#' @export
setGeneric(
  "bestLogLikelihood",
  function(celdaMod) {
    standardGeneric("bestLogLikelihood")
  }
)
#' @title Get the log-likelihood
#' @description Retrieves the final log-likelihood from all iterations of Gibbs
#'  sampling used to generate a celdaModel.
#' @param celdaMod A celdaModel object of class celda_C, celda_G, or celda_CG.
#' @return Numeric. The log-likelihood at the final step of Gibbs sampling used
#'  to generate the model.
#' @examples
#' data(celdaCGMod)
#' bestLogLikelihood(celdaCGMod)
#' @export
setMethod("bestLogLikelihood",
  signature = c(celdaMod = "celdaModel"),
  function(celdaMod) {
    celdaMod@finalLogLik
  }
)


#' @title Get or set the cell cluster labels from a celda
#'  \link[SingleCellExperiment]{SingleCellExperiment} object
#' @description Return or set the cell cluster labels determined
#'  by \link{celda_C} or \link{celda_CG} models.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object
#'  returned by \link{celda_C}, or \link{celda_CG}.
#' @return List. Contains cell cluster z (for celda_C and celda_CG Models)
#'  and/or gene module y (for celda_G and celda_CG Models).
#' @examples
#' data(sceCelda_CG)
#' clusters(sceCelda_CG)
#' @export
setGeneric("clusters",
    function(sce) {
        standardGeneric("clusters")
    })
#' @rdname clusters
#' @export
setMethod("clusters",
    signature(sce = "SingleCellExperiment"),
    function(sce) {
        return(SummarizedExperiment::colData(sce)$celda_cell_cluster)
    })


#' @rdname clusters
#' @export
setGeneric("clusters<-",
    function(sce, value) standardGeneric("clusters<-")
)
#' @rdname clusters
#' @export
setReplaceMethod("clusters", signature(sce = "SingleCellExperiment"),
    function(sce, value) {
        SummarizedExperiment::colData(sce)$celda_cell_cluster <- value
        return(sce)
    })


#' @title Get or set the feature module labels from a celda
#'  \link[SingleCellExperiment]{SingleCellExperiment} object
#' @description Return or set the feature module cluster labels determined
#'  by \link{celda_G} or \link{celda_CG} models.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object
#'  returned by \link{celda_G} or \link{celda_CG}.
#' @return List. Contains feature module y (for celda_G and celda_CG Models).
#' @examples
#' data(sceCelda_CG)
#' modules(sceCelda_CG)
#' @export
setGeneric("modules",
    function(sce) {
        standardGeneric("modules")
    })
#' @rdname modules
#' @export
setMethod("modules",
    signature(sce = "SingleCellExperiment"),
    function(sce) {
        return(SummarizedExperiment::rowData(sce)$celda_feature_module)
    })


#' @rdname modules
#' @export
setGeneric("modules<-",
    function(sce, value) standardGeneric("modules<-")
)
#' @rdname modules
#' @export
setReplaceMethod("modules", signature(sce = "SingleCellExperiment"),
    function(sce, value) {
        SummarizedExperiment::rowData(sce)$celda_feature_module <- value
        return(sce)
    })


setClass("celda_C",
  representation(sampleLabel = "factor"),
  contains = "celdaModel"
)


#' @title Get or set sample labels from a celda
#'  \link[SingleCellExperiment]{SingleCellExperiment} object
#' @description Return or set the sample labels for the cells in \code{sce}.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object
#'  returned by \link{celda_C}, \link{celda_G}, or \link{celda_CG}.
#' @return Character. Contains the sample labels provided at model creation,
#'  or those automatically generated by celda.
#' @examples
#' data(sceCelda_CG)
#' sampleLabel(sceCelda_CG)
#' @export
setGeneric("sampleLabel",
    function(sce) {
        standardGeneric("sampleLabel")
    })
#' @rdname sampleLabel
#' @export
setMethod("sampleLabel",
    signature(sce = "SingleCellExperiment"),
    function(sce) {
        return(SummarizedExperiment::colData(sce)$sample_label)
    })


#' @rdname sampleLabel
#' @export
setGeneric("sampleLabel<-",
    function(sce, value) standardGeneric("sampleLabel<-")
)
#' @rdname sampleLabel
#' @export
setReplaceMethod("sampleLabel", signature(sce = "SingleCellExperiment"),
    function(sce, value) {
        SummarizedExperiment::colData(sce)$sample_label <- value
        return(sce)
    })


#' @title Get celda model from a celda
#'  \link[SingleCellExperiment]{SingleCellExperiment} object
#' @description Return the celda model for \code{sce} returned by
#'  \link{celda_C}, \link{celda_G} or \link{celda_CG}.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object
#'  returned by \link{celda_C}, \link{celda_G}, or \link{celda_CG}.
#' @return Character. The celda model. Can be one of "celda_C", "celda_G", or
#'  "celda_CG".
#' @examples
#' data(sceCelda_CG)
#' celdaModel(sceCelda_CG)
#' @export
setGeneric("celdaModel",
    function(sce) {
        standardGeneric("celdaModel")
    })
#' @rdname celdaModel
#' @export
setMethod("celdaModel",
    signature(sce = "SingleCellExperiment"),
    function(sce) {
        tryCatch(
            if (S4Vectors::metadata(sce)$celda_parameters$model %in%
                    c("celda_C", "celda_G", "celda_CG")) {
                return(S4Vectors::metadata(sce)$celda_parameters$model)
            } else {
                stop("S4Vectors::metadata(sce)$celda_parameters$model must be",
                    " one of 'celda_C', 'celda_G', or 'celda_CG'")
            },
            error = function(e) {
                message("S4Vectors::metadata(sce)$celda_parameters$model must",
                    " exist! Try running celda model (celda_C, celda_CG, or",
                    " celda_G) first.")
                stop(e)
            })
    })
setClass("celda_G", contains = "celdaModel")

setClass("celda_CG", contains = c("celda_C", "celda_G"))

setClass(
  "celdaList",
  representation(
    runParams = "data.frame",
    resList = "list",
    countChecksum = "character",
    perplexity = "matrix"
  )
)


#' @title Get run parameters provided to `celdaGridSearch()`
#' @description Returns details on the clustering parameters, and model priors
#'  provided to `celdaGridSearch()` when the provided celdaList was
#'  created.
#' @param celdaList An object of class celdaList.
#' @return Data Frame. Contains details on the various K/L parameters, chain
#'  parameters, and final log-likelihoods derived for each model in the
#'  provided celdaList.
#' @examples
#' data(celdaCGGridSearchRes)
#' runParams(celdaCGGridSearchRes)
#' @export
setGeneric(
  "runParams",
  function(celdaList) {
    standardGeneric("runParams")
  }
)
#' @title Get run parameters provided to `celdaGridSearch()`
#' @description Returns details on the clustering parameters, and model priors
#'  provided to `celdaGridSearch()` when the provided celdaList was
#'  created.
#' @param celdaList An object of class celdaList.
#' @return Data Frame. Contains details on the various K/L parameters, chain
#'  parameters, and final log-likelihoods derived for each model in the
#'  provided celdaList.
#' @examples
#' data(celdaCGGridSearchRes)
#' runParams(celdaCGGridSearchRes)
#' @export
setMethod("runParams",
  signature = c(celdaList = "celdaList"),
  function(celdaList) {
    celdaList@runParams
  }
)


#' @title Get final celdaModels from a celdaList
#' @description Returns all models generated during a `celdaGridSearch()` run.
#' @param celdaList An object of class celdaList.
#' @return List. Contains one celdaModel object for each of the parameters
#'  specified in the `runParams()` of the provided celda list.
#' @examples
#' data(celdaCGGridSearchRes)
#' celdaCGGridModels <- resList(celdaCGGridSearchRes)
#' @export
setGeneric(
  "resList",
  function(celdaList) {
    standardGeneric("resList")
  }
)
#' @title Get final celdaModels from a celdaList
#' @description Returns all models generated during a `celdaGridSearch()` run.
#' @param celdaList An object of class celdaList.
#' @return List. Contains one celdaModel object for each of the parameters
#'  specified in the `runParams()` of the provided celda list.
#' @examples
#' data(celdaCGGridSearchRes)
#' celdaCGGridModels <- resList(celdaCGGridSearchRes)
#' @export
setMethod("resList",
  signature = c(celdaList = "celdaList"),
  function(celdaList) {
    celdaList@resList
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


#' @title Append two celdaList objects
#' @description Returns a single celdaList representing the combination of two
#'  provided celdaList objects.
#' @return A celdaList object. This object contains all resList entries and
#'  runParam records from both lists.
#' @param list1 A celda_list object
#' @param list2 A celda_list object to be joined with list_1
#' @examples
#' data(celdaCGGridSearchRes)
#' appendedList <- appendCeldaList(
#'   celdaCGGridSearchRes,
#'   celdaCGGridSearchRes
#' )
#' @importFrom methods new
#' @export
appendCeldaList <- function(list1, list2) {
  if (!is.element("celdaList", class(list1)) |
    !is.element("celdaList", class(list2))) {
    stop("Both parameters to appendCeldaList must be of class celdaList.")
  }
  if (!(countChecksum(list1) == countChecksum(list2))) {
    warning(
      "Provided lists have different countChecksums and may have",
      " been generated from different count matrices. Using checksum",
      " from first list..."
    )
  }
  newList <- methods::new(
    "celdaList",
    runParams = rbind(runParams(list1), runParams(list2)),
    resList = c(resList(list1), resList(list2)),
    countChecksum = countChecksum(list1),
    perplexity = matrix(nrow = 0, ncol = 0)
  )
  return(newList)
}


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

###############################################################################
# Generics
###############################################################################


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
#' @examples
#' data(sceCelda_CG)
#' celdaHeatmap(sceCelda_CG)
#' @return list A list containing dendrogram information and the heatmap grob
#' @export
setGeneric("celdaHeatmap",
    function(sce, ...) {
        standardGeneric("celdaHeatmap")
    })


#' @export
#' @rdname celdaHeatmap
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
    })


#' @title Calculate the Log-likelihood of a celda model
#' @description Calculate the log-likelihood for user-provided cell population
#'  and feature module cluster
#'  assignments on the count matrix, per the desired celda model.
#' @param sce A \linkS4class{SingleCellExperiment} object returned by
#'  \link{celda_C}, \link{celda_G}, or \link{celda_CG}, with the matrix
#'  located in the \code{useAssay} assay slot.
#'  Rows represent features and columns represent cells.
#' @param useAssay A string specifying which \link[SummarizedExperiment]{assay}
#'  slot to use. Default "counts".
#' @return The log-likelihood of the cluster assignment for the
#'  provided \linkS4class{SingleCellExperiment}.
#' @seealso `celda_C()` for clustering cells
#' @examples
#' data(sceCelda_C, sceCelda_CG)
#' loglikC <- logLikelihood(sceCelda_C)
#' loglikCG <- logLikelihood(sceCelda_CG)
#' @export
setGeneric("logLikelihood",
    function(sce, ...) {
        standardGeneric("logLikelihood")
    })


#' @rdname logLikelihood
#' @export
setMethod("logLikelihood", signature(sce = "SingleCellExperiment"),
    function(sce, useAssay = "counts") {
        if (celdaModel(sce) == "celda_C") {
            ll <- .logLikelihoodcelda_C(sce = sce, useAssay = useAssay)
            return(ll)
        } else if (celdaModel(sce) == "celda_CG") {
            ll <- .logLikelihoodcelda_CG(sce = sce, useAssay = useAssay)
            return(ll)
        } else if (celdaModel(sce) == "celda_G") {
            ll <- .logLikelihoodcelda_G(sce = sce, useAssay = useAssay)
            return(ll)
        } else {
            stop("S4Vectors::metadata(sce)$celda_parameters$model must be",
                " one of 'celda_C', 'celda_G', or 'celda_CG'")
        }
    })


#' @title Get the conditional probabilities of cell in subpopulations from celda
#'  model
#' @description Calculate the conditional probability of each cell belonging to
#'  each subpopulation given all other cell cluster assignments and/or
#'  each feature belonging to each module given all other feature cluster
#'  assignments in a celda model.
#' @param sce A \linkS4class{SingleCellExperiment} object returned by
#'  \link{celda_C}, \link{celda_G}, or \link{celda_CG}, with the matrix
#'  located in the \code{useAssay} assay slot.
#'  Rows represent features and columns represent cells.
#' @param useAssay A string specifying which \link[SummarizedExperiment]{assay}
#'  slot to use. Default "counts".
#' @param log Logical. If \code{FALSE}, then the normalized conditional
#'  probabilities will be returned. If \code{TRUE}, then the unnormalized log
#'  probabilities will be returned. Default \code{FALSE}.
#' @examples
#' data(sceCelda_CG)
#' clusterProb <- clusterProbability(sceCelda_CG, log = TRUE)
#' @return A list containging a matrix for the conditional cell subpopulation
#'  cluster and/or feature module probabilities.
#' @export
setGeneric("clusterProbability",
    function(sce, ...) {
        standardGeneric("clusterProbability")
    })


#' @seealso `celda_C()` for clustering cells
#' @examples
#' data(sceCelda_C)
#' clusterProb <- clusterProbability(sceCelda_C)
#' @rdname clusterProbability
#' @export
setMethod("clusterProbability", signature(sce = "SingleCellExperiment"),
    function(sce, useAssay = "counts", log = FALSE) {

        if (celdaModel(sce) == "celda_C") {
            cp <- .clusterProbabilityCeldaC(sce = sce,
                useAssay = useAssay,
                log = log)
            return(cp)
        } else if (celdaModel(sce) == "celda_CG") {
            cp <- .clusterProbabilityCeldaCG(sce = sce,
                useAssay = useAssay,
                log = log)
            return(cp)
        } else if (celdaModel(sce) == "celda_G") {
            cp <- .clusterProbabilityCeldaG(sce = sce,
                useAssay = useAssay,
                log = log)
            return(cp)
        } else {
            Stop()
        }
    })


#' @title Calculate the perplexity of a celda model
#' @description Perplexity is a statistical measure of how well a probability
#'  model can predict new data. Lower perplexity indicates a better model.
#' @param sce A \linkS4class{SingleCellExperiment} object returned by
#'  \link{celda_C}, \link{celda_G}, or \link{celda_CG}, with the matrix
#'  located in the \code{useAssay} assay slot.
#'  Rows represent features and columns represent cells.
#' @param useAssay A string specifying which \link[SummarizedExperiment]{assay}
#'  slot to use. Default "counts".
#' @return Numeric. The perplexity for the provided
#'  \linkS4class{SingleCellExperiment}.
#' @examples
#' data(sceCelda_CG)
#' clusterProb <- clusterProbability(sceCelda_CG)
#' @export
setGeneric("perplexity",
    function(sce, ...) {
        standardGeneric("perplexity")
    })


#' @importFrom matrixStats logSumExp
#' @rdname perplexity
#' @export
setMethod("perplexity", signature(sce = "SingleCellExperiment"),
    function(sce, useAssay = "counts") {

        if (celdaModel(sce) == "celda_C") {
            p <- .perplexityCelda_C(sce = sce, useAssay = useAssay)
            return(p)
        } else if (celdaModel(sce) == "celda_CG") {
            p <- .perplexityCelda_CG(sce = sce, useAssay = useAssay)
            return(p)
        } else if (celdaModel(sce) == "celda_G") {
            p <- .perplexityCelda_G(sce = sce, useAssay = useAssay)
            return(p)
        } else {
            stop("S4Vectors::metadata(sce)$celda_parameters$model must be",
                " one of 'celda_C', 'celda_G', or 'celda_CG'")
        }
    })


#' @title Simulate count data from the celda generative models.
#' @description This function generates a \linkS4class{SingleCellExperiment}
#'  containing a simulated counts matrix in the \code{"counts"} assay slot, as
#'  well as various parameters used in the simulation which can be
#'  useful for running celda and are stored in \code{metadata} slot. The user
#'  must provide the desired model (one of celda_C, celda_G, celda_CG) as well
#'  as any desired tuning parameters for those model's simulation functions
#'  as detailed below.
#' @param model Character. Options available in \code{celda::availableModels}.
#'  Can be one of \code{"celda_CG"}, \code{"celda_C"}, or \code{"celda_G"}.
#'  Default \code{"celda_CG"}.
#' @param S Integer. Number of samples to simulate. Default 5. Only used if
#'  \code{model} is one of \code{"celda_CG"} or \code{"celda_C"}.
#' @param CRange Integer vector. A vector of length 2 that specifies the lower
#'  and upper bounds of the number of cells to be generated in each sample.
#'  Default c(50, 100). Only used if
#'  \code{model} is one of \code{"celda_CG"} or \code{"celda_C"}.
#' @param NRange Integer vector. A vector of length 2 that specifies the lower
#'  and upper bounds of the number of counts generated for each cell. Default
#'  c(500, 1000).
#' @param C Integer. Number of cells to simulate. Default 100. Only used if
#'  \code{model} is \code{"celda_G"}.
#' @param G Integer. The total number of features to be simulated. Default 100.
#' @param K Integer. Number of cell populations. Default 5. Only used if
#'  \code{model} is one of \code{"celda_CG"} or \code{"celda_C"}.
#' @param L Integer. Number of feature modules. Default 10. Only used if
#'  \code{model} is one of \code{"celda_CG"} or \code{"celda_G"}.
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount
#'  to each cell population in each sample. Default 1. Only used if
#'  \code{model} is one of \code{"celda_CG"} or \code{"celda_C"}.
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to
#'  each feature module in each cell population. Default 1.
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to
#'  the number of features in each module. Default 5. Only used if
#'  \code{model} is one of \code{"celda_CG"} or \code{"celda_G"}.
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to
#'  each feature in each module. Default 1. Only used if
#'  \code{model} is one of \code{"celda_CG"} or \code{"celda_G"}.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  simulated count matrix stored in the "counts" assay slot. Function
#'  parameter settings are stored in the \link[S4Vectors]{metadata} slot. For
#'  \code{"celda_CG"} and \code{"celda_C"} models,
#'  columns \code{sample_label} and \code{celda_cell_cluster} in
#'  \link[SummarizedExperiment]{colData} contain simulated sample labels and
#'  cell population clusters. For \code{"celda_CG"} and \code{"celda_G"}
#'  models, column \code{celda_feature_module} in
#'  \link[SummarizedExperiment]{rowData} contains simulated gene modules.
#' @examples
#' sce <- simulateCells()
#' @export
simulateCells <- function(
    model = c("celda_CG", "celda_C", "celda_G"),
    S = 5,
    CRange = c(50, 100),
    NRange = c(500, 1000),
    C = 100,
    G = 100,
    K = 5,
    L = 10,
    alpha = 1,
    beta = 1,
    gamma = 5,
    delta = 1,
    seed = 12345) {

    model <- match.arg(model)

    if (model == "celda_C") {
        sce <- .simulateCellsMaincelda_C(model = model,
            S = S,
            CRange = CRange,
            NRange = NRange,
            G = G,
            K = K,
            alpha = alpha,
            beta = beta,
            seed = seed)
    } else if (model == "celda_CG") {
        sce <- .simulateCellsMaincelda_CG(
            model = model,
            S = S,
            CRange = CRange,
            NRange = NRange,
            G = G,
            K = K,
            L = L,
            alpha = alpha,
            beta = beta,
            gamma = gamma,
            delta = delta,
            seed = seed)
    } else if (model == "celda_G") {
        sce <- .simulateCellsMaincelda_G(
            model = model,
            C = C,
            L = L,
            NRange = NRange,
            G = G,
            beta = beta,
            delta = delta,
            gamma = gamma,
            seed = seed)
    } else {
        stop("'model' must be one of 'celda_C', 'celda_G', or 'celda_CG'")
    }

    return(sce)
}


#' @title Generate factorized matrices showing each feature's influence on cell
#'  / gene clustering
#' @description Generates factorized matrices showing the contribution of each
#'  feature in each cell population or each cell population in each sample.
#' @param sce A \linkS4class{SingleCellExperiment} object returned by
#'  \link{celda_C}, \link{celda_G}, or \link{celda_CG}, with the matrix
#'  located in the \code{useAssay} assay slot.
#'  Rows represent features and columns represent cells.
#' @param useAssay A string specifying which \link[SummarizedExperiment]{assay}
#'  slot to use. Default "counts".
#' @param type Character vector. A vector containing one or more of "counts",
#'  "proportion", or "posterior". "counts" returns the raw number of counts for
#'  each factorized matrix. "proportions" returns the normalized probabilities
#'  for each factorized matrix, which are calculated by dividing the raw counts
#'  in each factorized matrix by the total counts in each column. "posterior"
#'  returns the posterior estimates which include the addition of the Dirichlet
#'  concentration parameter (essentially as a pseudocount). Default
#'  \code{"counts"}.
#' @examples
#' data(sceCelda_CG)
#' factorizedMatrices <- factorizeMatrix(sceCelda_CG, type = "posterior")
#' @return A list with elements for `counts`, `proportions`, or `posterior`
#'  probabilities. Each element will be a list containing factorized matrices
#'  for `module` and `sample`.
#' @export
setGeneric("factorizeMatrix",
    function(sce, ...) {standardGeneric("factorizeMatrix")})


#' @examples
#' data(sceCelda_C)
#' factorizedMatrices <- factorizeMatrix(sceCelda_C, type = "posterior")
#' @seealso `celda_C()` for clustering cells
#' @rdname factorizeMatrix
#' @export
setMethod("factorizeMatrix", signature(sce = "SingleCellExperiment"),
    function(sce,
        useAssay = "counts",
        type = c("counts", "proportion", "posterior")) {

        if (celdaModel(sce) == "celda_C") {
            res <- .factorizeMatrixCelda_C(sce = sce, useAssay = useAssay,
                type = type)
        } else if (celdaModel(sce) == "celda_CG") {
            res <- .factorizeMatrixCelda_CG(sce = sce, useAssay = useAssay,
                type = type)
        } else if (celdaModel(sce) == "celda_G") {
            res <- .factorizeMatrixCelda_G(sce = sce, useAssay = useAssay,
                type = type)
        } else {
            stop("S4Vectors::metadata(sce)$celda_parameters$model must be",
                " one of 'celda_C', 'celda_G', or 'celda_CG'")
        }
        return(res)
    })


#' @title Renders probability and relative expression heatmaps to visualize the
#'  relationship between feature modules and cell populations.
#' @description It is often useful to visualize to what degree each feature
#' influences each cell cluster. This can also be useful for identifying
#' features which may be redundant or unassociated with cell clustering.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda object of class "celda_C" or "celda_CG".
#' @param ... Additional parameters.
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' celdaProbabilityMap(celdaCGSim$counts, celdaCGMod)
#' @return A grob containing the specified plots
#' @export
setGeneric("celdaProbabilityMap",
  signature = "celdaMod",
  function(counts, celdaMod, ...) {
    standardGeneric("celdaProbabilityMap")
  }
)


#' @title Embeds cells in two dimensions using tSNE based on celda_CG results.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda object of class `celda_CG`.
#' @param maxCells Integer. Maximum number of cells to plot. Cells will be
#'  randomly subsampled if ncol(counts) > maxCells. Larger numbers of cells
#'  requires more memory. Default \code{25000}.
#' @param minClusterSize Integer. Do not subsample cell clusters below this
#'  threshold. Default \code{100}.
#' @param initialDims integer. The number of dimensions that should be retained
#'  in the initial PCA step. Default \code{20}.
#' @param modules Integer vector. Determines which features modules to use for
#'  tSNE. If NULL, all modules will be used. Default NULL.
#' @param perplexity Numeric. Perplexity parameter for tSNE. Default \code{20}.
#' @param maxIter Integer. Maximum number of iterations in tSNE generation.
#'  Default \code{2500}.
#' @param ... Additional parameters.
#' @return Numeric Matrix of dimension `ncol(counts)` x 2, colums representing
#'  the "X" and "Y" coordinates in the data's t-SNE represetation.
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' tsneRes <- celdaTsne(celdaCGSim$counts, celdaCGMod)
#' @export
setGeneric("celdaTsne",
  signature = "celdaMod",
  function(counts,
           celdaMod,
           maxCells = 25000,
           minClusterSize = 100,
           initialDims = 20,
           modules = NULL,
           perplexity = 20,
           maxIter = 2500,
           ...) {
    # counts = processCounts(counts)
    # compareCountMatrix(counts, celdaMod)
    standardGeneric("celdaTsne")
  }
)


#' @title Embeds cells in two dimensions using umap.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda object of class `celda_CG`.
#' @param maxCells Integer. Maximum number of cells to plot. Cells will be
#'  randomly subsampled if ncol(counts) > maxCells. Larger numbers of cells
#'  requires more memory. Default 25000.
#' @param minClusterSize Integer. Do not subsample cell clusters below this
#'  threshold. Default 100.
#' @param modules Integer vector. Determines which features modules to use for
#'  tSNE. If NULL, all modules will be used. Default NULL.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#'  the UMAP algorithm.
#' @param ... Additional parameters to `uwot::umap`
#' @return A two column matrix of UMAP coordinates#' @examples
#' data(celdaCGSim, celdaCGMod)
#' umapRes <- celdaUmap(celdaCGSim$counts, celdaCGMod)
#' @export
setGeneric("celdaUmap",
  signature = "celdaMod",
  function(counts,
           celdaMod,
           maxCells = NULL,
           minClusterSize = 100,
           modules = NULL,
           seed = 12345,
           ...) {
    standardGeneric("celdaUmap")
  }
)


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
#' data(sceCelda_CG)
#' module <- featureModuleLookup(sce = sceCelda_CG,
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
            list[x] <- modules(sce)[which(rownames(sce) ==
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


#' @title Uniform Manifold Approximation and Projection (UMAP) dimension
#'  reduction for celda \code{sce} object
#' @description Embeds cells in two dimensions using \link[uwot]{umap} based on
#'  a celda model. For celda_C \code{sce} objects, PCA on the normalized counts
#'  is used to reduce the number of features before applying UMAP. For celda_CG
#'  \code{sce} object, UMAP is run on module probabilities to reduce the number
#'  of features instead of using PCA. Module probabilities are square-root
#'  transformed before applying UMAP.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object
#'  returned by \link{celda_C}, \link{celda_G}, or \link{celda_CG}.
#' @param useAssay A string specifying which \link[SummarizedExperiment]{assay}
#'  slot to use. Default "counts".
#' @param maxCells Integer. Maximum number of cells to plot. Cells will be
#'  randomly subsampled if \code{ncol(sce) > maxCells}. Larger numbers of cells
#'  requires more memory. If NULL, no subsampling will be performed.
#'  Default NULL.
#' @param minClusterSize Integer. Do not subsample cell clusters below this
#'  threshold. Default 100.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @param nNeighbors The size of local neighborhood used for
#'   manifold approximation. Larger values result in more global
#'   views of the manifold, while smaller values result in more
#'   local data being preserved. Default 30.
#'   See \link[uwot]{umap} for more information.
#' @param minDist The effective minimum distance between embedded points.
#'   Smaller values will result in a more clustered/clumped
#'   embedding where nearby points on the manifold are drawn
#'   closer together, while larger values will result on a more
#'   even dispersal of points. Default 0.75.
#'   See \link[uwot]{umap} for more information.
#' @param spread The effective scale of embedded points. In combination with
#'  \code{min_dist}, this determines how clustered/clumped the
#'   embedded points are. Default 1. See \link[uwot]{umap} for more information.
#' @param pca Logical. Whether to perform
#' dimensionality reduction with PCA before UMAP. Only works for celda_C
#'  \code{sce} objects.
#' @param initialDims Integer. Number of dimensions from PCA to use as
#' input in UMAP. Default 50. Only works for celda_C \code{sce} objects.
#' @param cores Number of threads to use. Default 1.
#' @param ... Additional parameters to pass to \link[uwot]{umap}.
#' @examples
#' data(sceCelda_CG)
#' umapRes <- celdaUmap(sceCelda_CG)
#' @return \code{sce} with UMAP coordinates
#'  (columns "celda_UMAP1" & "celda_UMAP2") added to
#'  \code{\link[SummarizedExperiment]{colData}(sce)}.
#' @export
setGeneric("celdaUmap",
    function(sce, ...) {
        standardGeneric("celdaUmap")
    })


#' @rdname celdaUmap
#' @export
setMethod("celdaUmap", signature(sce = "SingleCellExperiment"),
    function(sce,
        useAssay = "counts",
        maxCells = NULL,
        minClusterSize = 100,
        modules = NULL,
        seed = 12345,
        nNeighbors = 30,
        minDist = 0.75,
        spread = 1,
        pca = TRUE,
        initialDims = 50,
        cores = 1,
        ...) {

        if (is.null(seed)) {
            res <- .celdaUmap(sce = sce,
                useAssay = useAssay,
                maxCells = maxCells,
                minClusterSize = minClusterSize,
                modules = modules,
                seed = seed,
                nNeighbors = nNeighbors,
                minDist = minDist,
                spread = spread,
                pca = pca,
                initialDims = initialDims,
                cores = cores,
                ...)
        } else {
            with_seed(seed,
                res <- .celdaUmap(sce = sce,
                    useAssay = useAssay,
                    maxCells = maxCells,
                    minClusterSize = minClusterSize,
                    modules = modules,
                    seed = seed,
                    nNeighbors = nNeighbors,
                    minDist = minDist,
                    spread = spread,
                    pca = pca,
                    initialDims = initialDims,
                    cores = cores,
                    ...))
        }

        SummarizedExperiment::colData(sce)["celda_UMAP1"] <- res$UMAP1
        SummarizedExperiment::colData(sce)["celda_UMAP2"] <- res$UMAP2
        return(sce)
    })


.celdaUMAP <- function(sce,
    useAssay,
    maxCells,
    minClusterSize,
    modules,
    seed,
    nNeighbors,
    minDist,
    spread,
    pca,
    initialDims,
    cores,
    ...) {

    celdaMod <- celdaModel(sce)

    if (celdaMod == "celda_C") {
        res <- .celdaUmapC(sce = sce,
            useAssay = useAssay,
            maxCells = maxCells,
            minClusterSize = minClusterSize,
            nNeighbors = nNeighbors,
            minDist = minDist,
            spread = spread,
            pca = pca,
            initialDims = initialDims,
            cores = cores,
            ...)
    } else if (celdaMod == "celda_CG") {
        res <- .celdaUmapCG(sce = sce,
            useAssay = useAssay,
            maxCells = maxCells,
            minClusterSize = minClusterSize,
            modules = modules,
            seed = seed,
            nNeighbors = nNeighbors,
            minDist = minDist,
            spread = spread,
            cores = cores,
            ...)
    } else if (celdaMod == "celda_G") {
        res <- .celdaUmapG(sce = sce,
            useAssay = useAssay,
            maxCells = maxCells,
            minClusterSize = minClusterSize,
            modules = modules,
            seed = seed,
            nNeighbors = nNeighbors,
            minDist = minDist,
            spread = spread,
            cores = cores,
            ...)
    } else {
        stop("S4Vectors::metadata(sce)$celda_parameters$model must be",
            " one of 'celda_C', 'celda_G', or 'celda_CG'")
    }
    return(res)

}


.celdaUmapC <- function(sce,
    useAssay,
    maxCells,
    minClusterSize,
    nNeighbors,
    minDist,
    spread,
    pca,
    initialDims,
    cores,
    ...) {

    preparedCountInfo <- .prepareCountsForDimReductionCeldaC(sce,
        useAssay,
        maxCells,
        minClusterSize)
    umapRes <- .calculateUmap(preparedCountInfo$norm,
        nNeighbors = nNeighbors,
        minDist = minDist,
        spread = spread,
        pca = pca,
        initialDims = initialDims,
        cores = cores,
        ...
    )

    final <- matrix(NA, nrow = ncol(sce), ncol = 2)
    final[preparedCountInfo$cellIx, ] <- umapRes
    rownames(final) <- colnames(sce)
    colnames(final) <- c("UMAP1", "UMAP2")
    return(final)
}


.celdaUmapCG <- function(sce,
    useAssay,
    maxCells,
    minClusterSize,
    modules,
    seed,
    nNeighbors,
    minDist,
    spread,
    cores,
    ...) {

    preparedCountInfo <- .prepareCountsForDimReductionCeldaCG(sce,
        useAssay,
        maxCells,
        minClusterSize,
        modules)
    umapRes <- .calculateUmap(preparedCountInfo$norm,
        nNeighbors = nNeighbors,
        minDist = minDist,
        spread = spread,
        cores = cores,
        ...)

    final <- matrix(NA, nrow = ncol(sce), ncol = 2)
    final[preparedCountInfo$cellIx, ] <- umapRes
    rownames(final) <- colnames(sce)
    colnames(final) <- c("UMAP1", "UMAP2")
    return(final)
}


.celdaUmapG <- function(sce = sce,
    useAssay = useAssay,
    maxCells = maxCells,
    minClusterSize = minClusterSize,
    modules = modules,
    seed = seed,
    nNeighbors = nNeighbors,
    minDist = minDist,
    spread = spread,
    cores = cores,
    ...) {

    preparedCountInfo <- .prepareCountsForDimReductionCeldaCG(sce = sce,
        useAssay = useAssay,
        maxCells = maxCells,
        minClusterSize = minClusterSize,
        modules = modules)
    umapRes <- .calculateUmap(preparedCountInfo$norm,
        nNeighbors = nNeighbors,
        minDist = minDist,
        spread = spread,
        cores = cores,
        ...)

    final <- matrix(NA, nrow = ncol(sce), ncol = 2)
    final[preparedCountInfo$cellIx, ] <- umapRes
    rownames(final) <- colnames(sce)
    colnames(final) <- c("UMAP1", "UMAP2")
    return(final)
}


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
#' @examples
#' data(sceCelda_CG)
#' celdaProbabilityMap(sceCelda_CG)
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
        } else if (celdaModel(sce) == "celda_G") {
            pm <- .celdaProbabilityMapG(sce = sce, useAssay = useAssay,
                level = level)
        } else {
            stop("S4Vectors::metadata(sce)$celda_parameters$model must be",
                " one of 'celda_C', 'celda_G', or 'celda_CG'")
        }
        return(pm)
    })
