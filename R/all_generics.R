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


#' @title Get or set the cluster labels from a celda
#'  \link[SingleCellExperiment]{SingleCellExperiment} object
#' @description Return or set the cell (and/or gene) cluster labels determined
#'  by celda model.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object
#'  returned by \link{celda_C}, \link{celda_G}, or \link{celda_CG}.
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
        return(SummarizedExperiment::colData(sce)$cell_cluster)
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
        SummarizedExperiment::colData(sce)$cell_cluster <- value
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
' @examples
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
#'  celda_CG or celda_G model results returned by \link{celda_C} or
#'  \link{celda_CG}.}
#' @param ... Additional parameters passed to \link{plotHeatmap}.
#' @examples
#' data(sceCelda_CG)
#' celdaHeatmap(sceCelda_CG)
#' @return list A list containing dendrogram information and the heatmap grob
#' @export
setGeneric("celdaHeatmap",
    function(sce, ...) {
        standardGeneric("celdaHeatmap")
    })


#' @seealso `celda_C()` for clustering cells and `celdaTsne()` for generating
#'  2-dimensional coordinates
#' @examples
#' data(sceCelda_C)
#' celdaHeatmap(sceCelda_C)
#' @return list A list containing dendrograms and the heatmap grob
#' @export
#' @rdname celdaHeatmap
setMethod("celdaHeatmap", signature(sce = "SingleCellExperiment"),
    function(sce, useAssay = "counts", featureIx = NULL, nfeatures = 25, ...) {

        if (celdaModel(sce) == "celda_C") {
            counts <- SummarizedExperiment::assay(sce, i = useAssay)
            norm <- normalizeCounts(counts,
                normalize = "proportion",
                transformationFun = sqrt)

            if (is.null(featureIx)) {
                return(plotHeatmap(norm,
                    z = clusters(sce), ...))
            }

            return(plotHeatmap(norm[featureIx, ],
                z = clusters(sce), ...))
        } else if (celdaModel(sce) %in% c("celda_G", "celda_CG")) {

            # fm <- factorizeMatrix(counts, celdaMod, type = "proportion")
            # top <- celda::topRank(fm$proportions$module, n = nfeatures)
            # ix <- unlist(top$index)
            # norm <- normalizeCounts(counts,
            #     normalize = "proportion",
            #     transformationFun = sqrt)
            # plotHeatmap(norm[ix, ],
            #     z = clusters(celdaMod)$z,
            #     y = clusters(celdaMod)$y[ix],
            #     ...)
        }
    })


#' @title Calculate the Log-likelihood of celda model
#' @description Calculate the log-likelihood for a user-provided cluster
#'  assignment on the count matrix, per the desired celda model.
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
        if (S4Vectors::metadata(sce)$celda_parameters$model == "celda_C") {
            ll <- .logLikelihoodcelda_C(sce = sce, useAssay = useAssay)
            return(ll)
        }

        stop("metadata(sce)$celda_parameters$model must be 'celda_C'!")
    })


#' @title Get the conditional probabilities of cell in subpopulations from celda
#'  model
#' @description Calculate the conditional probability of each cell belonging to
#'  each subpopulation given all other cell cluster assignments in a celda
#'  model.
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
#'  cluster probabilities.
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

        if (S4Vectors::metadata(sce)$celda_parameters$model == "celda_C") {
            counts <- SummarizedExperiment::assay(sce, i = useAssay)

            z <- clusters(sce)
            sampleLabel <- sampleLabel(sce)
            s <- as.integer(sampleLabel)

            K <- S4Vectors::metadata(sce)$celda_parameters$K
            alpha <- S4Vectors::metadata(sce)$celda_parameters$alpha
            beta <- S4Vectors::metadata(sce)$celda_parameters$beta

            p <- .cCDecomposeCounts(counts, s, z, K)

            nextZ <- .cCCalcGibbsProbZ(counts = counts,
                mCPByS = p$mCPByS,
                nGByCP = p$nGByCP,
                nByC = p$nByC,
                nCP = p$nCP,
                z = z,
                s = s,
                K = K,
                nG = p$nG,
                nM = p$nM,
                alpha = alpha,
                beta = beta,
                doSample = FALSE)
            zProb <- t(nextZ$probs)

            if (!isTRUE(log)) {
                zProb <- .normalizeLogProbs(zProb)
            }

            return(list(zProbability = zProb))
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

        if (S4Vectors::metadata(sce)$celda_parameters$model == "celda_C") {
            counts <- SummarizedExperiment::assay(sce, i = useAssay)


            counts <- .processCounts(counts)

            factorized <- factorizeMatrix(sce = sce,
                useAssay = useAssay,
                type = "posterior")
            theta <- log(factorized$posterior$sample)
            phi <- log(factorized$posterior$module)
            s <- as.integer(sampleLabel(sce))

            inner.log.prob <- eigenMatMultInt(phi, Counts) + theta[, s]
            logPx <- sum(apply(inner.log.prob, 2, matrixStats::logSumExp))

            perplexity <- exp(- (logPx / sum(Counts)))
            return(perplexity)
        } else {
            stop("metadata(sce)$celda_parameters$model is not 'celda_C',",
                " 'celda_G', or 'celda_CG'!")
        }
    })


#' @title Simulate count data from the celda generative models.
#' @description This function generates a list containing a simulated counts
#'  matrix, as well as various parameters used in the simulation which can be
#'  useful for running celda. The user must provide the desired model
#'  (one of celda_C, celda_G, celda_CG) as well as any desired tuning
#'  parameters for those model's simulation functions as detailed below.
#' @param model Character. Options available in `celda::availableModels`.
#' @param ... Additional parameters.
#' @return List. Contains the simulated counts matrix, derived cell cluster
#'  assignments, the provided parameters, and estimated Dirichlet distribution
#'  parameters for the model.
#' @examples
#' data(celdaCGSim)
#' dim(celdaCGSim$counts)
#' @export
simulateCells <- function(model, ...) {
    do.call(paste0(".simulateCellsMain", model), args = list(...))
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
#'  \code{c("counts", "proportion", "posterior")}.
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

        if (S4Vectors::metadata(sce)$celda_parameters$model == "celda_C") {
            res <- .factorizeMatrixCelda_C(sce = sce, useAssay = useAssay,
                type = type)
            return(res)
        }
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
#' a specified list of genes from a celdaModel.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Model of class "celda_G" or "celda_CG".
#' @param feature Character vector. Identify feature modules for the specified
#'  feature names.
#' @param exactMatch Logical. Whether to look for exactMatch of the gene name
#'  within counts matrix. Default \code{TRUE}.
#' @return List. Each entry corresponds to the feature module determined for
#' the provided features.
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' featureModuleLookup(
#'   counts = celdaCGSim$counts,
#'   celdaMod = celdaCGMod, "Gene_1"
#' )
#' @export
setGeneric("featureModuleLookup",
  signature = "celdaMod",
  function(counts, celdaMod, feature, exactMatch = TRUE) {
    standardGeneric("featureModuleLookup")
  }
)
