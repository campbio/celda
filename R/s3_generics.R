################################################################################
# Generics
################################################################################

#' Get run parameters for a celda run.
#'
#' @param celda.list Object of class "celda_list". An object containing celda models returned from `celdaGridSearch()`.
#' @examples 
#' celda.sim = simulateCells("celda_CG")
#' celda.mod = celda_CG(celda.sim$counts, K=celda.sim$K, L=celda.sim$L,
#'                      nchains=1, max.iter=1)
#' runParams(celda.mod)
#' @return A data.frame containing the run parameters used to generate the provided celda_list.
#' @export
runParams = function(celda.list) {
  return(celda.list$run.params)
}


#' Get the complete log likelihood for a given celda model.
#'
#' @param celda.mod Celda object of class "celda_C", "celda_G", or "celda_CG".
#' @return Numeric Vector. The log-likelihood of the model's cluster assignments during each iteration.
#' @examples 
#' celda.sim = simulateCells("celda_CG")
#' celda.mod = celda_CG(celda.sim$counts, K=celda.sim$K, L=celda.sim$L,
#'                      nchains=1, max.iter=1)
#' complete.loglik = completeLogLikelihood(celda.mod)
#' @export
completeLogLikelihood = function(celda.mod) {
  return(celda.mod$completeLogLik)
}


#' Get the log likelihood from the final iteration of Gibbs sampling
#' for a given celda model.
#'
#' @param celda.mod Celda object of class "celda_C", "celda_G", or "celda_CG".
#' @examples
#' celda.sim = simulateCells("celda_CG")
#' celda.mod = celda_CG(celda.sim$counts, K=celda.sim$K, L=celda.sim$L,
#'                      nchains=1, max.iter=1)
#' finalLogLikelihood(celda.mod = celda.mod)
#' @return The final log-likelihood determined by Gibbs sampling for this model
#' @export
finalLogLikelihood = function(celda.mod) {
  return(celda.mod$finalLogLik)
}


#' Get the probability of the cluster assignments generated during a celda run.
#'
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @param log Logical. If FALSE, then the normalized conditional probabilities will be returned. If TRUE, then the unnormalized log probabilities will be returned. Default FALSE.  
#' @examples
#' celda.sim = simulateCells("celda_CG")
#' celda.mod = celda_CG(celda.sim$counts, K=celda.sim$K, L=celda.sim$L,
#'                      nchains=1, max.iter=1)
#' cluster.prob = clusterProbability(celda.sim$counts, celda.mod)
#' @return A numeric vector of the cluster assignment probabilties
#' @export
clusterProbability = function(counts, celda.mod, log=FALSE) {
  UseMethod("clusterProbability", celda.mod)
}


#' Render a stylable heatmap of count data based on celda clustering results.
#'
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`. 
#' @param celda.mod Celda object of class "celda_C", "celda_G", or "celda_CG".
#' @param ... Additional parameters.
#' @examples 
#' celda.sim = simulateCells("celda_CG")
#' celda.mod = celda_CG(celda.sim$counts, K=celda.sim$K, L=celda.sim$L,
#'                      nchains=1, max.iter=1)
#' celdaHeatmap(celda.sim$counts, celda.mod)
#' @return list A list containing dendrogram information and the heatmap grob
#' @export 
celdaHeatmap <- function(counts, celda.mod, ...) {
  UseMethod("celdaHeatmap", celda.mod)
}


#' Calculate a log-likelihood for a user-provided cluster assignment and count matrix, per the desired celda model. 
#' 
#' @param counts The counts matrix used to generate the provided cluster assignments.
#' @param model Celda model. Options available in `celda::available.models`.
#' @param ... Additional parameters.
#' @return The log-likelihood of the provided cluster assignment for the provided counts matrix.
#' @examples
#' celda.sim = simulateCells(model="celda_CG")
#' loglik = logLikelihood(celda.sim$counts, model="celda_CG", 
#'                        sample.label=celda.sim$sample.label,
#'                        z=celda.sim$z, y=celda.sim$y,
#'                        K=celda.sim$K, L=celda.sim$L,
#'                        alpha=celda.sim$alpha, beta=celda.sim$beta,
#'                        gamma=celda.sim$gamma, delta=celda.sim$delta)
#' @export
logLikelihood = function(counts, model, ...) {
  class(counts) = c(model)
  do.call(paste("logLikelihood.", model, sep=""),
          list(counts, ...))
}


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
#' celda.sim = simulateCells(model = "celda_CG")
#' dim(celda.sim$counts)
#' @export
simulateCells = function(model, ...) {
  class(model) = c(class(model), model)
  UseMethod("simulateCells", model)
}


#' Generate factorized matrices showing each feature's influence on cell / gene clustering
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class "celda_C", "celda_G", or "celda_CG".
#' @param type A character vector containing one or more of "counts", "proportions", or "posterior". "counts" returns the raw number of counts for each entry in each matrix. "proportions" returns the counts matrix where each vector is normalized to a probability distribution. "posterior" returns the posterior estimates which include the addition of the Dirichlet concentration parameter (essentially as a pseudocount).
#' @examples 
#' celda.sim = simulateCells("celda_CG")
#' celda.mod = celda_CG(celda.sim$counts, K=celda.sim$K, L=celda.sim$L,
#'                      nchains=1, max.iter=1)
#' factorized.matrices = factorizeMatrix(celda.sim$counts, celda.mod, "posterior")
#' @return A list of lists of the types of factorized matrices specified
#' @export
factorizeMatrix = function(counts, celda.mod, type) {
  
  UseMethod("factorizeMatrix", celda.mod)
}

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
#' celda.sim = simulateCells("celda_CG")
#' celda.mod = celda_CG(celda.sim$counts, K=celda.sim$K, L=celda.sim$L,
#'                      nchains=1, max.iter=1)
#' celdaProbabilityMap(celda.sim$counts, celda.mod)
#' @return A grob containing the specified plots
#' @export
celdaProbabilityMap = function(counts, celda.mod, ...) {
  UseMethod("celdaProbabilityMap", celda.mod)
}

#' Embeds cells in two dimensions using tSNE based on celda_CG results.
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @param ... Additional parameters.
#' @return Numeric Matrix of dimension `ncol(counts)` x 2, colums representing the "X" and "Y" coordinates in the data's t-SNE represetation.
#' @examples 
#' celda.sim = simulateCells("celda_CG")
#' celda.mod = celda_CG(celda.sim$counts, K=celda.sim$K, L=celda.sim$L,
#'                      nchains=1, max.iter=1)
#' tsne.res = celdaTsne(celda.sim$counts, celda.mod)
#' @export
celdaTsne = function(counts, celda.mod, ...) {
  counts = processCounts(counts)
  compareCountMatrix(counts, celda.mod)
  if (!isTRUE(class(celda.mod) %in% c("celda_CG","celda_C","celda_G"))) {
    stop("celda.mod argument is not of class celda_C, celda_G or celda_CG")
  }
  UseMethod("celdaTsne", celda.mod)
}

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
#' celda.sim = simulateCells("celda_CG")
#' celda.mod = celda_CG(celda.sim$counts, K=celda.sim$K, L=celda.sim$L,
#'                      nchains=1, max.iter=1)
#' featureModuleLookup(counts = celda.sim$counts, celda.mod = celda.mod, "Gene_1")
#' @export
featureModuleLookup = function(counts, celda.mod, feature, exact.match = TRUE){
  class(celda.mod) = c(class(celda.mod), celda.mod)
  UseMethod("featureModuleLookup", celda.mod)
}
