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
setGeneric("params",
    function(celdaMod) {
        standardGeneric("params")
    })
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
    })


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
setGeneric("matrixNames",
    function(celdaMod) {
        standardGeneric("matrixNames")
    })
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
    })


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
setGeneric("logLikelihoodHistory",
    function(celdaMod) {
        standardGeneric("logLikelihoodHistory")
    })
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
    })


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
setGeneric("bestLogLikelihood",
    function(celdaMod) {
        standardGeneric("bestLogLikelihood")
    })
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
    })


#' @title Get clustering outcomes from a celdaModel
#' @description Returns the z / y results corresponding to the cell / gene
#'  cluster labels determined by the provided celdaModel.
#' @param celdaMod celdaModel. Options available in `celda::availableModels`.
#' @return List. Contains z (for celda_C and celdaCGModels) and/or y
#'  (for celda_G and celdaCGModels)
#' @examples
#' data(celdaCGMod)
#' clusters(celdaCGMod)
#' @export
setGeneric("clusters",
    function(celdaMod) {
        standardGeneric("clusters")
    })
#' @title Get clustering outcomes from a celdaModel
#' @description Returns the z / y results corresponding to the cell / gene
#'  cluster labels determined by the provided celdaModel.
#' @param celdaMod celdaModel. Options available in `celda::availableModels`.
#' @return List. Contains z (for celda_C and celdaCGModels) and/or y
#'  (for celda_G and celdaCGModels)
#' @examples
#' data(celdaCGMod)
#' clusters(celdaCGMod)
#' @export
setMethod("clusters",
    signature = c(celdaMod = "celdaModel"),
    function(celdaMod) {
        return(celdaMod@clusters)
    })


setClass("celda_C",
    representation(sampleLabel = "factor"),
    contains = "celdaModel")


#' @title Get sampleLabels from a celdaModel
#' @description Returns the sampleLabels for the count matrix provided for
#'  generation of a given celdaModel.
#' @param celdaMod celdaModel. Options available in `celda::availableModels`.
#' @return Character. Contains the sampleLabels provided at model creation time,
#'  or those automatically generated by celda.
#' @examples
#' data(celdaCGMod)
#' sampleLabel(celdaCGMod)
#' @export
setGeneric("sampleLabel",
    function(celdaMod) {
        standardGeneric("sampleLabel")
    })
#' @title Get sampleLabels from a celdaModel
#' @description Returns the sampleLabels for the count matrix provided for
#'  generation of a given celdaModel.
#' @param celdaMod celdaModel. Options available in `celda::availableModels`.
#' @return Character. Contains the sampleLabels provided at model creation time,
#'  or those automatically generated by celda.
#' @examples
#' data(celdaCGMod)
#' sampleLabel(celdaCGMod)
#' @export
setMethod("sampleLabel",
    signature = c(celdaMod = "celdaModel"),
    function(celdaMod) {
        celdaMod@sampleLabel
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
#' @param celdaMod An object of class celdaList.
#' @return Data Frame. Contains details on the various K/L parameters, chain
#'  parameters, and final log-likelihoods derived for each model in the provided
#'  celdaList.
#' @examples
#' data(celdaCGGridSearchRes)
#' runParams(celdaCGGridSearchRes)
#' @export
setGeneric("runParams",
    function(celdaMod) {
        standardGeneric("runParams")
    })
#' @title Get run parameters provided to `celdaGridSearch()`
#' @description Returns details on the clustering parameters, and model priors
#'  provided to `celdaGridSearch()` when the provided celdaList was
#'  created.
#' @param celdaMod An object of class celdaList.
#' @return Data Frame. Contains details on the various K/L parameters, chain
#'  parameters, and final log-likelihoods derived for each model in the provided
#'  celdaList.
#' @examples
#' data(celdaCGGridSearchRes)
#' runParams(celdaCGGridSearchRes)
#' @export
setMethod("runParams",
    signature = c(celdaMod = "celdaList"),
    function(celdaMod) {
        celdaMod@runParams
    })


#' @title Get final celdaModels from a celdaList
#' @description Returns all models generated during a `celdaGridSearch()` run.
#' @param celdaMod An object of class celdaList.
#' @return List. Contains one celdaModel object for each of the parameters
#'  specified in the `runParams()` of the provided celda list.
#' @examples
#' data(celdaCGGridSearchRes)
#' celdaCGGridModels <- resList(celdaCGGridSearchRes)
#' @export
setGeneric("resList",
    function(celdaMod) {
        standardGeneric("resList")
    })
#' @title Get final celdaModels from a celdaList
#' @description Returns all models generated during a `celdaGridSearch()` run.
#' @param celdaMod An object of class celdaList.
#' @return List. Contains one celdaModel object for each of the parameters
#'  specified in the `runParams()` of the provided celda list.
#' @examples
#' data(celdaCGGridSearchRes)
#' celdaCGGridModels <- resList(celdaCGGridSearchRes)
#' @export
setMethod("resList",
    signature = c(celdaMod = "celdaList"),
    function(celdaMod) {
        celdaMod@resList
    })


#' @title Get perplexity for every model in a celdaList
#' @description Returns perplexity for each model in a celdaList as calculated
#'  by `perplexity().`
#' @param celdaMod A celdaModel object of class "celda_C", "celda_G", or
#'  "celda_CG".
#' @return List. Contains one celdaModel object for each of the parameters
#'  specified in the `runParams()` of the provided celda list.
#' @examples
#' data(celdaCGGridSearchRes)
#' celdaCGGridModelPerplexities <- celdaPerplexity(celdaCGGridSearchRes)
#' @export
setGeneric("celdaPerplexity",
    function(celdaMod) {
        standardGeneric("celdaPerplexity")
    })
#' @title Get perplexity for every model in a celdaList
#' @description Returns perplexity for each model in a celdaList as calculated
#'  by `perplexity().`
#' @param celdaMod A celdaModel object of class "celda_C", "celda_G", or
#'  "celda_CG".
#' @return List. Contains one celdaModel object for each of the parameters
#'  specified in the `runParams()` of the provided celda list.
#' @examples
#' data(celdaCGGridSearchRes)
#' celdaCGGridModelPerplexities <- celdaPerplexity(celdaCGGridSearchRes)
#' @export
setMethod("celdaPerplexity",
    signature = c(celdaMod = "celdaList"),
    function(celdaMod) {
        celdaMod@perplexity
    })


#' @title Append two celdaList objects
#' @description Returns a single celdaList representing the combination of two
#'  provided celdaList objects.
#' @return A celdaList object. This object contains all resList entries and
#'  runParam records from both lists.
#' @param list1 A celda_list object
#' @param list2 A celda_list object to be joined with list_1
#' @examples
#' data(celdaCGGridSearchRes)
#' appendedList <- appendCeldaList(celdaCGGridSearchRes,
#'     celdaCGGridSearchRes)
#' @importFrom methods new
#' @export
appendCeldaList <- function(list1, list2) {
    if (!is.element("celdaList", class(list1)) |
            !is.element("celdaList", class(list2))) {
        stop("Both parameters to appendCeldaList must be of class celdaList.")
    }
    if (!(list1@countChecksum == list2@countChecksum)) {
        warning("Provided lists have different countChecksums and may have",
            " been generated from different count matrices. Using checksum",
            " from first list...")
    }
    newList <- methods::new(
        "celdaList",
        runParams = rbind(list1@runParams, list2@runParams),
        resList = c(list1@resList, list2@resList),
        countChecksum = list1@countChecksum,
        perplexity = matrix(nrow = 0, ncol = 0))
    return(newList)
}

################################################################################
# Generics
################################################################################


#' @title Plot celda Heatmap
#' @description Render a stylable heatmap of count data based on celda
#'  clustering results.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod A celdaModel object of class "celda_C", "celda_G", or
#'  "celda_CG".
#' @param featureIx Integer vector. Select features for display in heatmap. If
#'  NULL, no subsetting will be performed. Default NULL.
#' @param ... Additional parameters.
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' celdaHeatmap(celdaCGSim$counts, celdaCGMod)
#' @return list A list containing dendrogram information and the heatmap grob
#' @export
setGeneric("celdaHeatmap",
    signature = "celdaMod",
    function(counts, celdaMod, featureIx, ...) {
        standardGeneric("celdaHeatmap")
    })

#' @title Calculate LogLikelihood
#' @description Calculate a log-likelihood for a user-provided cluster
#'  assignment and count matrix, per the desired celdaModel.
#' @param counts The counts matrix used to generate the provided cluster
#'  assignments.
#' @param model celdaModel. Options available in `celda::availableModels`.
#' @param ... Additional parameters.
#' @return The log-likelihood of the provided cluster assignment for the
#'  provided counts matrix.
#' @examples
#' data(celdaCGSim)
#' loglik <- logLikelihood(celdaCGSim$counts,
#'     model = "celda_CG",
#'     sampleLabel = celdaCGSim$sampleLabel,
#'     z = celdaCGSim$z, y = celdaCGSim$y,
#'     K = celdaCGSim$K, L = celdaCGSim$L,
#'     alpha = celdaCGSim$alpha, beta = celdaCGSim$beta,
#'     gamma = celdaCGSim$gamma, delta = celdaCGSim$delta
#' )
#' @export
#'
#'
logLikelihood <- function(counts, model, ...) {
    do.call(paste0("logLikelihood", model),
        args = list(counts = counts, ...))
}


#' @title Get cluster probability
#' @description Get the probability of the cluster assignments generated during
#'  a celda run.
#' @param counts Integer matrix. Rows represent features and columns represent
#' cells. This matrix should be the same as the one used to generate `celdaMod`.
#' @param celdaMod celdaModel. Options available in `celda::availableModels`.
#' @param log Logical. If FALSE, then the normalized conditional probabilities
#'  will be returned. If TRUE, then the unnormalized log probabilities will be
#'  returned. Default FALSE.
#' @param ... Additional parameters.
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' clusterProb <- clusterProbability(celdaCGSim$counts, celdaCGMod)
#' @return A numeric vector of the cluster assignment probabilties
#' @export
setGeneric("clusterProbability",
    signature = "celdaMod",
    function(counts, celdaMod, log = FALSE, ...) {
        standardGeneric("clusterProbability")
    })


#' @title Calculate the perplexity from a single celdaModel
#' @description Perplexity can be seen as a measure of how well a provided set
#'  of cluster assignments fit the data being clustered.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod celdaModel. Options available in `celda::availableModels`.
#' @param newCounts A newCounts matrix used to calculate perplexity. If NULL,
#'  perplexity will be calculated for the 'counts' matrix. Default NULL.
#' @return Numeric. The perplexity for the provided count data and model.
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' perplexity <- perplexity(celdaCGSim$counts, celdaCGMod)
#' @export
setGeneric("perplexity",
    signature = "celdaMod",
    function(counts, celdaMod, newCounts = NULL) {
        standardGeneric("perplexity")
    })


#' @title Simulate count data from the celda generative models.
#' @description This function generates a list containing a simulated counts
#'  matrix, as well as various parameters used in the simulation which can be
#'  useful for running celda. The user must provide the desired model
#'  (one of celda_C, celda_G, celda_CG) as well as any desired tuning parameters
#'  for those model's simulation functions as detailed below.
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
    do.call(paste0("simulateCells", model), args = list(...))
}


#' @title Generate factorized matrices showing each feature's influence on cell
#'  / gene clustering
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda object of class "celda_C", "celda_G", or "celda_CG".
#' @param type A character vector containing one or more of "counts",
#'  "proportions", or "posterior". "counts" returns the raw number of counts for
#'  each entry in each matrix. "proportions" returns the counts matrix where
#'  each vector is normalized to a probability distribution. "posterior" returns
#'  the posterior estimates which include the addition of the Dirichlet
#'  concentration parameter (essentially as a pseudocount).
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' factorizedMatrices <- factorizeMatrix(
#'     celdaCGSim$counts, celdaCGMod,
#'     "posterior"
#' )
#' @return A list of lists of the types of factorized matrices specified
#' @export
setGeneric("factorizeMatrix",
    signature = "celdaMod",
    function(counts,
        celdaMod,
        type = c("counts", "proportion", "posterior")) {
        standardGeneric("factorizeMatrix")
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
    })


#' @title Embeds cells in two dimensions using tSNE based on celda_CG results.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda object of class `celda_CG`.
#' @param maxCells Integer. Maximum number of cells to plot. Cells will be
#'  randomly subsampled if ncol(counts) > maxCells. Larger numbers of cells
#'  requires more memory. Default 25000.
#' @param minClusterSize Integer. Do not subsample cell clusters below this
#'  threshold. Default 100.
#' @param initialDims integer. The number of dimensions that should be retained
#'  in the initial PCA step. Default 20.
#' @param modules Integer vector. Determines which features modules to use for
#'  tSNE. If NULL, all modules will be used. Default NULL.
#' @param perplexity Numeric. Perplexity parameter for tSNE. Default 20.
#' @param maxIter Integer. Maximum number of iterations in tSNE generation.
#'  Default 2500.
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
    })


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
#' @param initialDims Integer. PCA will be used to reduce the dimentionality
#'  of the dataset. The top 'initialDims' principal components will be used
#'  for umap. Default 20.
#' @param modules Integer vector. Determines which features modules to use for
#'  tSNE. If NULL, all modules will be used. Default NULL.
#' @param umapConfig An object of class "umapConfig" specifying parameters to
#'  the UMAP algorithm.
#' @param ... Additional parameters.
#' @return Numeric Matrix of dimension `ncol(counts)` x 2, colums representing
#'  the "X" and "Y" coordinates in the data's t-SNE represetation.
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' tsneRes <- celdaUmap(celdaCGSim$counts, celdaCGMod)
#' @importFrom umap umap.defaults
#' @export
setGeneric("celdaUmap",
    signature = "celdaMod",
    function(counts,
        celdaMod,
        maxCells = 25000,
        minClusterSize = 100,
        initialDims = 20,
        modules = NULL,
        umapConfig = umap::umap.defaults) {
        standardGeneric("celdaUmap")
    })


#' @title Obtain the gene module of a gene of interest
#' @description This function will output the corresponding feature module for a
#'  specified list of genes from a celdaModel.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Model of class "celda_G" or "celda_CG".
#' @param feature Character vector. Identify feature modules for the specified
#'  feature names.
#' @param exactMatch Logical. Whether to look for exactMatch of the gene name
#'  within counts matrix. Default TRUE.
#' @return List. Each entry corresponds to the feature module determined for the
#'  provided features
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' featureModuleLookup(
#'     counts = celdaCGSim$counts,
#'     celdaMod = celdaCGMod, "Gene_1")
#' @export
setGeneric("featureModuleLookup",
    signature = "celdaMod",
    function(counts, celdaMod, feature, exactMatch = TRUE) {
        standardGeneric("featureModuleLookup")
    })
