setClass("celdaModel", 
         representation(completeLogLik = "numeric", 
                        finalLogLik = "numeric",
                        seed = "numeric", 
                        count.checksum = "character",
                        names = "list"))

setGeneric("completeLogLik",
           function(celda.mod){ standardGeneric("completeLogLik") })
setGeneric("completeLogLik<-",
           function(celda.mod, value){ standardGeneric("completeLogLik<-") })
setMethod("completeLogLik",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod){  celda.mod@completeLogLik  })
setReplaceMethod("completeLogLik", "celdaModel",
                 function(celda.mod, value){
                   celda.mod@completeLogLik = value
                   celda.mod
                 })

setGeneric("finalLogLik",
           function(celda.mod){ standardGeneric("finalLogLik") })
setGeneric("finalLogLik<-",
           function(celda.mod, value){ standardGeneric("finalLogLik<-") })
setMethod("finalLogLik",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod){  celda.mod@finalLogLik  })
setReplaceMethod("finalLogLik", "celdaModel",
                 function(celda.mod, value){
                   celda.mod@finalLogLik = value
                   celda.mod
                 })

setGeneric("seed",
           function(celda.mod){ standardGeneric("seed") })
setGeneric("seed<-",
           function(celda.mod, value){ standardGeneric("seed<-") })
setMethod("seed",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod){  celda.mod@seed  })
setReplaceMethod("seed", "celdaModel",
                 function(celda.mod, value){
                   celda.mod@seed = value
                   celda.mod
                 })

setGeneric("countChecksum",
           function(celda.mod){ standardGeneric("countChecksum") })
setGeneric("countChecksum<-",
           function(celda.mod, value){ standardGeneric("countChecksum<-") })
setMethod("countChecksum",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod){  celda.mod@count.checksum  })
setReplaceMethod("countChecksum", "celdaModel",
                 function(celda.mod, value){
                   celda.mod@count.checksum = value
                   celda.mod
                 })

setGeneric("matrixNames",
           function(celda.mod){ standardGeneric("matrixNames") })
setGeneric("matrixNames<-",
           function(celda.mod, value){ standardGeneric("matrixNames<-") })
setMethod("matrixNames",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod){  celda.mod@names  })
setReplaceMethod("matrixNames", "celdaModel",
                 function(celda.mod, value){
                   celda.mod@names = value
                   celda.mod
                 })


setClass("celdaList",
         representation(runParams = "data.frame",
                        res.list = "list",
                        count.checksum = "numeric"))


setClass("celda_C",
         representation(z = "numeric",
                        K = "numeric",
                        alpha = "numeric",
                        beta = "numeric",
                        sample.label = "factor"),
         contains = "celdaModel")

setGeneric("z",
           function(celda.mod){ standardGeneric("z") })
setGeneric("z<-",
           function(celda.mod, value){ standardGeneric("z<-") })
setMethod("z",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod){  celda.mod@z  })
setReplaceMethod("z", "celdaModel",
                 function(celda.mod, value){
                   celda.mod@z = value
                   celda.mod
                 })
setGeneric("K",
           function(celda.mod){ standardGeneric("K") })
setGeneric("K<-",
           function(celda.mod, value){ standardGeneric("K<-") })
setMethod("K",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod){  celda.mod@K  })
setReplaceMethod("K", "celdaModel",
                 function(celda.mod, value){
                   celda.mod@K = value
                   celda.mod
                 })
setGeneric("alpha",
           function(celda.mod){ standardGeneric("alpha") })
setGeneric("alpha<-",
           function(celda.mod, value){ standardGeneric("alpha<-") })
setMethod("alpha",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod){  celda.mod@alpha  })
setReplaceMethod("alpha", "celdaModel",
                 function(celda.mod, value){
                   celda.mod@alpha = value
                   celda.mod
                 })
# TODO Overwrites base::beta()...
setGeneric("beta",
           function(celda.mod){ standardGeneric("beta") })
setGeneric("beta<-",
           function(celda.mod, value){ standardGeneric("beta<-") })
setMethod("beta",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod){  celda.mod@beta  })
setReplaceMethod("beta", "celdaModel",
                 function(celda.mod, value){
                   celda.mod@beta = value
                   celda.mod
                 })

setGeneric("sample.label",
           function(celda.mod){ standardGeneric("sample.label") })
setGeneric("sample.label<-",
           function(celda.mod, value){ standardGeneric("sample.label<-") })
setMethod("sample.label",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod){  celda.mod@sample.label  })
setReplaceMethod("sample.label", "celdaModel",
                 function(celda.mod, value){
                   celda.mod@sample.label = value
                   celda.mod
                 })


setClass("celda_G",
         representation(y = "numeric",
                        L = "numeric",
                        delta = "numeric",
                        gamma = "numeric"),
         contains = "celdaModel")

setGeneric("y",
           function(celda.mod){ standardGeneric("y") })
setGeneric("y<-",
           function(celda.mod, value){ standardGeneric("y<-") })
setMethod("y",
           signature=c(celda.mod="celda_G"),
           function(celda.mod){  celda.mod@y  })
setReplaceMethod("y", "celdaModel",
                 function(celda.mod, value){
                   celda.mod@y = value
                   celda.mod
                 })
setGeneric("L",
           function(celda.mod){ standardGeneric("L") })
setGeneric("L<-",
           function(celda.mod, value){ standardGeneric("L<-") })
setMethod("L",
           signature=c(celda.mod="celda_G"),
           function(celda.mod){  celda.mod@L  })
setReplaceMethod("L", "celdaModel",
                 function(celda.mod, value){
                   celda.mod@L = value
                   celda.mod
                 })
setGeneric("delta",
           function(celda.mod){ standardGeneric("delta") })
setGeneric("delta<-",
           function(celda.mod, value){ standardGeneric("delta<-") })
setMethod("delta",
           signature=c(celda.mod="celda_G"),
           function(celda.mod){  celda.mod@delta  })
setReplaceMethod("delta", "celdaModel",
                 function(celda.mod, value){
                   celda.mod@delta = value
                   celda.mod
                 })
# TODO Overwrites base::gamma()...
# setGeneric("gamma",
#            function(celda.mod){ standardGeneric("gamma") })
# setGeneric("gamma<-",
#            function(celda.mod, value){ standardGeneric("gamma<-") })
# setMethod("gamma",
#            signature=c(celda.mod="celda_G"),
#            function(celda.mod){  celda.mod@gamma  })
# setReplaceMethod("gamma", "celdaModel",
#                  function(celda.mod, value){
#                    celda.mod@gamma = value
#                    celda.mod
#                  })

setClass("celda_CG",
         contains = c("celda_C", "celda_G"))





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
setGeneric("logLikelihood", 
           signature="model",
           function(counts, model, ...) {
             standardGeneric("logLikelihood")
           })


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
           function(counts, celda.mod, log=FALSE) {
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
                    initial.dims=20, perplexity=20, max.iter=2500, seed=12345, 
                    ...) {
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
