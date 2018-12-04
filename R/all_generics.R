setClass("celdaModel", 
         representation(completeLogLik = "numeric", 
                        finalLogLik = "numeric",
                        seed = "numeric", 
                        count.checksum = "character",
                        names = "list",
                        clustering = "list",    # K, L, z, y, etc.
                        modelPriors = "list"))  # alpha, beta, delta, etc.

#' @title Get log-likelihood history
#' @description Retrieves the complete log-likelihood from all iterations of Gibbs sampling used to generate a celda model.
#' 
#' @return Numeric. The log-likelihood at each step of Gibbs sampling used to generate the model.
#' @examples
#' completeLogLik(celda.CG.mod)
#' @export
setGeneric("completeLogLik",
           function(celda.mod){ standardGeneric("completeLogLik") })
setMethod("completeLogLik",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod){  celda.mod@completeLogLik  })

# @title Set log-likelihood history on a celda model
# @description Set the complete log-likelihood from all iterations of Gibbs sampling used to generate a celda model.
# @return A celdaModel object. The provided input object, with the updated completeLogLik property.
# @examples
# completeLogLik(celda.CG.mod) = c(0.00, 0.01, 0.02)  # Lenth must match num.iter
#@export
# setGeneric("completeLogLik<-",
#            function(celda.mod, value){ standardGeneric("completeLogLik<-") })
# setReplaceMethod("completeLogLik", "celdaModel",
#                  function(celda.mod, value){
#                    celda.mod@completeLogLik = value
#                    celda.mod
#                  })


#' @title Get final log-likelihood 
#' @description Retrieves the final log-likelihood from all iterations of Gibbs sampling used to generate a celda model.
#' 
#' @return Numeric. The log-likelihood at the final step of Gibbs sampling used to generate the model.
#' @examples
#' finalLogLik(celda.CG.mod)
#' @export
setGeneric("finalLogLik",
           function(celda.mod){ standardGeneric("finalLogLik") })
setMethod("finalLogLik",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod){  celda.mod@finalLogLik  })

# @title Set log-likelihood history on a celda model
# @description Set the final log-likelihood of Gibbs sampling used to generate a celda model.
# @return A celdaModel object. The provided input object, with the updated completeLogLik property.
# @examples
# finalLogLik(celda.CG.mod) = 0.97 
#@export
# setGeneric("finalLogLik<-",
#            function(celda.mod, value){ standardGeneric("finalLogLik<-") })
# setReplaceMethod("finalLogLik", "celdaModel",
#                  function(celda.mod, value){
#                    celda.mod@finalLogLik = value
#                    celda.mod
#                  })


#' @title Get seed used to generate model
#' @description Retrieves the random seed used to generate a celda model.
#' @return Numeric. The random seed used to generate the provided celda model.
#' @examples
#' initialSeed(celda.CG.mod)
#' @export
setGeneric("initialSeed",
           function(celda.mod){ standardGeneric("initialSeed") })
setMethod("initialSeed",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod){  celda.mod@seed  })
#@export
# setGeneric("initialSeed<-",
#            function(celda.mod, value){ standardGeneric("initialSeed<-") })
# setReplaceMethod("initialSeed", "celdaModel",
#                  function(celda.mod, value){
#                    celda.mod@seed = value
#                    celda.mod
#                  })


#' @title Get count matrix checksum for comparison
#' @description Retrieves the MD5 checksum of the count matrix used to generate the provided celda mdoel.
#' @return Character. The MD5 hash of the count matrix used to generate the provided celda model.
#' @examples
#' countChecksum(celda.CG.mod)
#' @export
setGeneric("countChecksum",
           function(celda.mod){ standardGeneric("countChecksum") })
setMethod("countChecksum",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod){  celda.mod@count.checksum  })
#@export
# setGeneric("countChecksum<-",
#            function(celda.mod, value){ standardGeneric("countChecksum<-") })
# setReplaceMethod("countChecksum", "celdaModel",
#                  function(celda.mod, value){
#                    celda.mod@count.checksum = value
#                    celda.mod
#                  })

#' @title Get feature, cell and sample names from a celda model
#' @description Retrieves the row, column, and sample names used to generate a celda model.
#' @return List. Contains row, column, and sample character vectors corresponding to the values provided when the celda model was generated.
#' @examples
#' matrixNames(celda.CG.mod)
#' @export
setGeneric("matrixNames",
           function(celda.mod){ standardGeneric("matrixNames") })
setMethod("matrixNames",
           signature=c(celda.mod="celdaModel"),
#@export
# setGeneric("matrixNames<-",
#            function(celda.mod, value){ standardGeneric("matrixNames<-") })
          function(celda.mod){  celda.mod@names  })
# setReplaceMethod("matrixNames", "celdaModel",
#                  function(celda.mod, value){
#                    celda.mod@names = value
#                    celda.mod
#                  })


#' @title Get clustering parameters and outcomes from a celda model.
#' @description Returns the K/L parameters provided for modeling, as well as the corresponding z/y results.
#' @return List. Contains K, z (for celda_C and celda_CG models), and/or L, y (for celda_G and celda_CG models.) 
#' @examples
#' clustering(celda.CG.mod)
#' @export
setGeneric("clustering",
           function(celda.mod){ standardGeneric("clustering")})
setMethod("clustering", signature=c(celda.mod="celdaModel"),
          function(celda.mod){
            return(celda.mod@clustering)
          })
#@export
# setGeneric("clustering<-",
#            function(celda.mod, value){ standardGeneric("clustering<-") })


#' @title Get model prior parameters from a celda model.
#' @description Returns the model priors (e.g. alpha, beta) provided at model creation for a given celda model.
#' @return List. Contains alpha, beta (for celda_C and celda_CG models), or delta, gamma (for celda_G and celda_CG models).
#' @examples
#' modelPriors(celda.CG.mod)
#' @export
setGeneric("modelPriors",
           function(celda.mod){ standardGeneric("modelPriors")})
setMethod("modelPriors", signature=c(celda.mod="celdaModel"),
          function(celda.mod){
            return(celda.mod@modelPriors)
          })
#@export
# setGeneric("modelPriors<-",
#            function(celda.mod, value){ standardGeneric("modelPriors<-") })


setClass("celda_C",
         representation(sample.label = "factor"),
         contains = "celdaModel")
# setReplaceMethod("clustering", "celda_C",
#                  function(celda.mod, value){
#                    lapply(names(value),
#                           function(key) {
#                             if (!is.element(key, c("K", "z"))) {
#                               stop(paste0("Invalid parameter provided: ", key))
#                             }
#                             celda.mod@clustering[[key]] = value[[key]]
#                           })
#                  })
# setReplaceMethod("modelPriors", "celda_C",
#                  function(celda.mod, value){
#                    lapply(names(value),
#                           function(key) {
#                             if (!is.element(key, c("alpha", "beta"))) {
#                               stop(paste0("Invalid parameter provided: ", key))
#                             }
#                             celda.mod@modelPriors[[key]] = value[[key]]
#                           })
#                  })


#' @title Get sample labels
#' @description Returns the sample labels for the count matrix provided for generation of a given celda model.
#' @return Character. Contains the sample labels provided at model creation time, or those automatically generated by celda.
#' @examples
#' modelPriors(celda.CG.mod)
#' @export
setGeneric("sampleLabel",
           function(celda.mod){ standardGeneric("sampleLabel") })
setMethod("sampleLabel",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod){  celda.mod@sample.label  })
#@export
# setGeneric("sampleLabel<-",
#            function(celda.mod, value){ standardGeneric("sampleLabel<-") })
# setReplaceMethod("sampleLabel", "celdaModel",
#                  function(celda.mod, value){
#                    celda.mod@sample.label = value
#                    celda.mod
#                  })


setClass("celda_G",
         contains = "celdaModel")
# setReplaceMethod("clustering", "celda_G",
#                  function(celda.mod, value){
#                    lapply(names(value),
#                           function(key) {
#                             if (!is.element(key, c("L", "y"))) {
#                               stop(paste0("Invalid parameter provided: ", key))
#                             }
#                             celda.mod@clustering[[key]] = value[[key]]
#                           })
#                  })
# setReplaceMethod("modelPriors", "celda_G",
#                  function(celda.mod, value){
#                    lapply(names(value),
#                           function(key) {
#                             if (!is.element(key, c("beta", "delta", "gamma"))) {
#                               stop(paste0("Invalid parameter provided: ", key))
#                             }
#                             celda.mod@modelPriors[[key]] = value[[key]]
#                           })
#                  })


setClass("celda_CG",
         contains = c("celda_C", "celda_G"))
# setReplaceMethod("clustering", "celda_CG",
#                  function(celda.mod, value){
#                    lapply(names(value),
#                           function(key) {
#                             if (!is.element(key, c("K", "L", "y", "z"))) {
#                               stop(paste0("Invalid parameter provided: ", key))
#                             }
#                             celda.mod@clustering[[key]] = value[[key]]
#                           })
#                  })
# setReplaceMethod("modelPriors", "celda_CG",
#                  function(celda.mod, value){
#                    lapply(names(value),
#                           function(key) {
#                             if (!is.element(key, c("alpha", "beta", "delta", "gamma"))) {
#                               stop(paste0("Invalid parameter provided: ", key))
#                             }
#                             celda.mod@modelPriors[[key]] = value[[key]]
#                           })
#                  })




setClass("celdaList",
         representation(run.params = "data.frame",
                        res.list = "list",
                        count.checksum = "character",
                        perplexity = "matrix"))

#' @title Get run parameters provided to `celdaGridSearch()`
#' @description Returns details on the clustering parameters, model priors, and seeds provided to `celdaGridSearch()` when the provided celdaList was created.
#' @return Data Frame. Contains details on the various K/L parameters, chain parameters, and final log-likelihoods derived for each model in the provided celdaList.
#' @examples
#' runParams(celda.CG.grid.search.res)
#' @export
setGeneric("runParams",
           function(celda.mod){ standardGeneric("runParams") })
setMethod("runParams",
           signature=c(celda.mod="celdaList"),
           function(celda.mod){  celda.mod@run.params  })
#@export
# setGeneric("runParams<-",
#            function(celda.mod, value){ standardGeneric("runParams<-") })
# setReplaceMethod("runParams", "celdaModel",
#                  function(celda.mod, value){
#                    celda.mod@run.params = value
#                    celda.mod
#                  })

#' @title Get final celda models from a celdaList
#' @description Returns all models generated during a `celdaGridSearch()` run.
#' @return List. Contains one celdaModel object for each of the parameters specified in the `runParams()` of the provided celda list.
#' @examples
#' celda.CG.grid.models = resList(celda.CG.grid.search.res)
#' @export
setGeneric("resList",
           function(celda.mod){ standardGeneric("resList") })
setMethod("resList",
           signature=c(celda.mod="celdaList"),
           function(celda.mod){  celda.mod@res.list  })
#@export
#setGeneric("resList<-",
#           function(celda.mod, value){ standardGeneric("resList<-") })
# setReplaceMethod("resList", "celdaModel",
#                  function(celda.mod, value){
#                    celda.mod@resList = value
#                    celda.mod
#                  })

setMethod("countChecksum",
           signature=c(celda.mod="celdaList"),
           function(celda.mod){  celda.mod@count.checksum  })


#' @title Get perplexity for every model in a celdaList
#' @description Returns perplexity for each model in a celdaList as calculated by `perplexity().`
#' @return List. Contains one celdaModel object for each of the parameters specified in the `runParams()` of the provided celda list.
#' @examples
#' celda.CG.grid.model.perplexities = getPerplexity(celda.CG.grid.search.res)
#' @export
setGeneric("getPerplexity",
           function(celda.mod){ standardGeneric("getPerplexity") })
setMethod("getPerplexity",
           signature=c(celda.mod="celdaList"),
           function(celda.mod){  celda.mod@perplexity  })
# setGeneric("setPerplexity<-",
#            function(celda.mod, value){ standardGeneric("setPerplexity<-") })
# setReplaceMethod("setPerplexity", "celdaModel",
#                  function(celda.mod, value){
#                    celda.mod@perplexity = value
#                    celda.mod
#                  })


################################################################################
# Generics
################################################################################


#' Render a stylable heatmap of count data based on celda clustering results.
#'
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`. 
#' @param celda.mod Celda object of class "celda_C", "celda_G", or "celda_CG".
#' @param ... Additional parameters.
#' @examples 
#' celdaHeatmap(celda.CG.sim$counts, celda.CG.mod)
#' @return list A list containing dendrogram information and the heatmap grob
#' @export 
setGeneric("celdaHeatmap", 
           signature="celda.mod",
           function(counts, celda.mod, feature.ix, ...) {
             standardGeneric("celdaHeatmap")
           })


#' Calculate a log-likelihood for a user-provided cluster assignment and count matrix, per the desired celda model. 
#' 
#' @param counts The counts matrix used to generate the provided cluster assignments.
#' @param model Celda model. Options available in `celda::available.models`.
#' @param ... Additional parameters.
#' @return The log-likelihood of the provided cluster assignment for the provided counts matrix.
#' @examples
#' loglik = logLikelihood(celda.CG.sim$counts, model="celda_CG", 
#'                        sample.label=celda.CG.sim$sample.label,
#'                        z=celda.CG.sim$z, y=celda.CG.sim$y,
#'                        K=celda.CG.sim$K, L=celda.CG.sim$L,
#'                        alpha=celda.CG.sim$alpha, beta=celda.CG.sim$beta,
#'                        gamma=celda.CG.sim$gamma, delta=celda.CG.sim$delta)
#' @export
#' 
#' 
logLikelihood = function(counts, model, ...) {
  do.call(paste0("logLikelihood.", model), 
          args=list(counts=counts, ...))
}


#' Get the probability of the cluster assignments generated during a celda run.
#'
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @param log Logical. If FALSE, then the normalized conditional probabilities will be returned. If TRUE, then the unnormalized log probabilities will be returned. Default FALSE.  
#' @examples
#' cluster.prob = clusterProbability(celda.CG.sim$counts, celda.CG.mod)
#' @return A numeric vector of the cluster assignment probabilties
#' @export
setGeneric("clusterProbability", 
           signature="celda.mod",
           function(counts, celda.mod, log=FALSE, modules=NULL, ...) {
             standardGeneric("clusterProbability")
           })


#' Calculate the perplexity from a single celda model
#' 
#' Perplexity can be seen as a measure of how well a provided set of 
#' cluster assignments fit the data being clustered.
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class "celda_C", "celda_G" or "celda_CG".
#' @param new.counts A new counts matrix used to calculate perplexity. If NULL, perplexity will be calculated for the 'counts' matrix. Default NULL.
#' @return Numeric. The perplexity for the provided count data and model.
#' @examples
#' perplexity = perplexity(celda.CG.sim$counts, celda.CG.mod)
#' @export
setGeneric("perplexity",
           signature="celda.mod",
           function(counts, celda.mod, new.counts=NULL) {
             standardGeneric("perplexity")
           })


#' Simulate count data from the celda generative models.
#' 
#' This function generates a list containing a simulated counts matrix, as well as various parameters
#' used in the simulation which can be useful for running celda. The user must provide the desired model
#' (one of celda_C, celda_G, celda_CG) as well as any desired tuning parameters for those model's simulation
#' functions as detailed below.
#' 
#' @param model Character. Options available in `celda::available.models`.
#' @param ... Additional parameters.
#' @return List. Contains the simulated counts matrix, derived cell cluster assignments, the provided parameters, and estimated Dirichlet distribution parameters for the model.
#' @examples
#' dim(celda.CG.sim$counts)
#' @export
simulateCells = function(model, ...) {
  do.call(paste0("simulateCells.", model), args=list(...))
}


#' Generate factorized matrices showing each feature's influence on cell / gene clustering
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class "celda_C", "celda_G", or "celda_CG".
#' @param type A character vector containing one or more of "counts", "proportions", or "posterior". "counts" returns the raw number of counts for each entry in each matrix. "proportions" returns the counts matrix where each vector is normalized to a probability distribution. "posterior" returns the posterior estimates which include the addition of the Dirichlet concentration parameter (essentially as a pseudocount).
#' @examples 
#' factorized.matrices = factorizeMatrix(celda.CG.sim$counts, celda.CG.mod, 
#'                                       "posterior")
#' @return A list of lists of the types of factorized matrices specified
#' @export
setGeneric("factorizeMatrix",
           signature = "celda.mod",
           function(counts, celda.mod,  
                    type=c("counts", "proportion", "posterior")){ 
             standardGeneric("factorizeMatrix") 
           })

#' Renders probability and relative expression heatmaps to visualize the relationship between feature modules and cell populations.
#' 
#' It is often useful to visualize to what degree each feature influences each 
#' cell cluster. This can also be useful for identifying features which may
#' be redundant or unassociated with cell clustering.
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class "celda_C" or "celda_CG".
#' @param ... Additional parameters.
#' @examples
#' celdaProbabilityMap(celda.CG.sim$counts, celda.CG.mod)
#' @return A grob containing the specified plots
#' @export
setGeneric("celdaProbabilityMap",
           signature="celda.mod",
           function(counts, celda.mod, ...) {
             standardGeneric("celdaProbabilityMap")
           })


#' Embeds cells in two dimensions using tSNE based on celda_CG results.
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @param ... Additional parameters.
#' @return Numeric Matrix of dimension `ncol(counts)` x 2, colums representing the "X" and "Y" coordinates in the data's t-SNE represetation.
#' @examples 
#' tsne.res = celdaTsne(celda.CG.sim$counts, celda.CG.mod)
#' @export
setGeneric("celdaTsne",
           signature = "celda.mod",
           function(counts, celda.mod, max.cells=25000, min.cluster.size=100,
                    initial.dims=20, modules=NULL, perplexity=20, max.iter=2500, 
                    seed=12345, ...) {
             # counts = processCounts(counts)
             # compareCountMatrix(counts, celda.mod)
             standardGeneric("celdaTsne")
           })


#' Obtain the gene module of a gene of interest
#' 
#' This function will output the corresponding feature module for a specified list of genes from a celda model.
#'  
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Model of class "celda_G" or "celda_CG".
#' @param feature Character vector. Identify feature modules for the specified feature names. 
#' @param exact.match Logical. Whether to look for exact match of the gene name within counts matrix. Default TRUE. 
#' @return List. Each entry corresponds to the feature module determined for the provided features
#' @examples 
#' featureModuleLookup(counts = celda.CG.sim$counts, 
#'                     celda.mod = celda.CG.mod, "Gene_1")
#' @export
setGeneric("featureModuleLookup",
           signature = "celda.mod",
           function(counts, celda.mod, feature, exact.match = TRUE){
             standardGeneric("featureModuleLookup")
           })
