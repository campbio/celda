################################################################################
<<<<<<< HEAD
# S3 Methods                                                                   #
################################################################################
# Below are getters for the various types of celda models.                     #
# Concrete implementations of these functions are in their corresponding model #
# files (e.g. getZ.celda_C is in celda_C.R).                                   #
#                                                                              #
# TODO:                                                                        #
#        * Collapse ROxygen documentation into single page for these functions #
#        * Consider moving model specific implementations to their             #
#          corresponding files                                                 #
#        * Can reduce redundancy for celda_C / celda_G getters by renaming the #
#          fields on their respective return objects to match.                 #  
################################################################################



################################################################################
=======
>>>>>>> upstream-devel
# Generics
################################################################################

#' Get run parameters for a celda run.
#'
#' @param celda.list Object of class "celda_list". An object containing celda models returned from `celdaGridSearch()`.
<<<<<<< HEAD
=======
#' @examples 
#' celda.sim = simulateCells("celda_CG")
#' celda.mod = celda_CG(celda.sim$counts, K=celda.sim$K, L=celda.sim$L,
#'                      nchains=1, max.iter=1)
#' runParams(celda.mod)
#' @return A data.frame containing the run parameters used to generate the provided celda_list.
>>>>>>> upstream-devel
#' @export
runParams = function(celda.list) {
  return(celda.list$run.params)
}


<<<<<<< HEAD
#' Get the random seed for a given celda model.
#'
#' @param celda.mod Celda object of class "celda_C", "celda_G", or "celda_CG".
#' @export
seed = function(celda.mod) {
  return(celda.mod$seed)
}


#' Get the complete log likelihood for a given celda model.
#'
#' @param celda.mod Celda object of class "celda_C", "celda_G", or "celda_CG".
=======
#' Get the complete log likelihood for a given celda model.
#'
#' @param celda.mod Celda object of class "celda_C", "celda_G", or "celda_CG".
#' @return Numeric Vector. The log-likelihood of the model's cluster assignments during each iteration.
#' @examples 
#' celda.sim = simulateCells("celda_CG")
#' celda.mod = celda_CG(celda.sim$counts, K=celda.sim$K, L=celda.sim$L,
#'                      nchains=1, max.iter=1)
#' complete.loglik = completeLogLikelihood(celda.mod)
>>>>>>> upstream-devel
#' @export
completeLogLikelihood = function(celda.mod) {
  return(celda.mod$completeLogLik)
}


#' Get the log likelihood from the final iteration of Gibbs sampling
#' for a given celda model.
#'
#' @param celda.mod Celda object of class "celda_C", "celda_G", or "celda_CG".
<<<<<<< HEAD
=======
#' @examples
#' celda.sim = simulateCells("celda_CG")
#' celda.mod = celda_CG(celda.sim$counts, K=celda.sim$K, L=celda.sim$L,
#'                      nchains=1, max.iter=1)
#' finalLogLikelihood(celda.mod = celda.mod)
#' @return The final log-likelihood determined by Gibbs sampling for this model
>>>>>>> upstream-devel
#' @export
finalLogLikelihood = function(celda.mod) {
  return(celda.mod$finalLogLik)
}


<<<<<<< HEAD
#' Get the final gene / cell / gene & cell cluster assignments generated during
#' a celda run, dependent on the model provided.
#'
#' @param celda.mod Celda object of class "celda_C", "celda_G", or "celda_CG".
#' @export
finalClusterAssignment = function(celda.mod) {
  UseMethod("finalClusterAssignment", celda.mod)
}


#' Get the probability of the cluster assignments generated during a celda run.
#'
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param log Logical. If FALSE, then the normalized conditional probabilities will be returned. If TRUE, then the unnormalized log probabilities will be returned. Default FALSE.  
#' @export
clusterProbability = function(celda.mod, counts, log=FALSE) {
=======
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
>>>>>>> upstream-devel
  UseMethod("clusterProbability", celda.mod)
}


<<<<<<< HEAD
#' Get the K value used for each chain in a celda run.
#'
#' @param celda.mod Celda object of class "celda_C", "celda_G", or "celda_CG".
#' @export
getK = function(celda.mod) {
  UseMethod("getK", celda.mod)
}


#' Get the L value used for each chain in a celda run.
#'
#' @param celda.mod Celda object of class "celda_C", "celda_G", or "celda_CG".
#' @export
getL = function(celda.mod) {
  UseMethod("getL", celda.mod)
}


=======
>>>>>>> upstream-devel
#' Render a stylable heatmap of count data based on celda clustering results.
#'
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`. 
#' @param celda.mod Celda object of class "celda_C", "celda_G", or "celda_CG".
#' @param ... Additional parameters.
<<<<<<< HEAD
=======
#' @examples 
#' celda.sim = simulateCells("celda_CG")
#' celda.mod = celda_CG(celda.sim$counts, K=celda.sim$K, L=celda.sim$L,
#'                      nchains=1, max.iter=1)
#' celdaHeatmap(celda.sim$counts, celda.mod)
#' @return list A list containing dendrogram information and the heatmap grob
>>>>>>> upstream-devel
#' @export 
celdaHeatmap <- function(counts, celda.mod, ...) {
  UseMethod("celdaHeatmap", celda.mod)
}


#' Calculate a log-likelihood for a user-provided cluster assignment and count matrix, per the desired celda model. 
#' 
<<<<<<< HEAD
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @param counts The counts matrix used to generate the provided cluster assignments.
#' @return The log-likelihood of the provided cluster assignment for the provided counts matrix.
#' @param ... Additional parameters.
#' @export
calculateLoglikFromVariables <- function(counts, celda.mod, ...) {
  class(counts) = c(celda.mod, class(counts))
  do.call(paste("calculateLoglikFromVariables.", celda.mod, sep=""),
=======
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
>>>>>>> upstream-devel
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
<<<<<<< HEAD
=======
#' @return List. Contains the simulated counts matrix, derived cell cluster assignments, the provided parameters, and estimated Dirichlet distribution parameters for the model.
#' @examples
#' celda.sim = simulateCells(model = "celda_CG")
#' dim(celda.sim$counts)
>>>>>>> upstream-devel
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
<<<<<<< HEAD
=======
#' @examples 
#' celda.sim = simulateCells("celda_CG")
#' celda.mod = celda_CG(celda.sim$counts, K=celda.sim$K, L=celda.sim$L,
#'                      nchains=1, max.iter=1)
#' factorized.matrices = factorizeMatrix(celda.sim$counts, celda.mod, "posterior")
#' @return A list of lists of the types of factorized matrices specified
>>>>>>> upstream-devel
#' @export
factorizeMatrix = function(counts, celda.mod, type) {
  
  UseMethod("factorizeMatrix", celda.mod)
}

<<<<<<< HEAD
#' Generate factorized matrices showing each feature's influence on cell / gene clustering
=======
#' Renders probability and relative expression heatmaps to visualize the relationship between feature modules and cell populations.
#' 
#' It is often useful to visualize to what degree each feature influences each 
#' cell cluster. This can also be useful for identifying features which may
#' be redundant or unassociated with cell clustering.
>>>>>>> upstream-devel
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class "celda_C" or "celda_CG".
#' @param ... Additional parameters.
<<<<<<< HEAD
#' @export
celdaProbabilityMap = function(counts, celda.mod, ...) {
  
  UseMethod("celdaProbabilityMap", celda.mod)
}
#' Runs tSNE via Rtsne based on the CELDA model and specified cell states.
=======
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
>>>>>>> upstream-devel
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @param ... Additional parameters.
<<<<<<< HEAD
=======
#' @return Numeric Matrix of dimension `ncol(counts)` x 2, colums representing the "X" and "Y" coordinates in the data's t-SNE represetation.
#' @examples 
#' celda.sim = simulateCells("celda_CG")
#' celda.mod = celda_CG(celda.sim$counts, K=celda.sim$K, L=celda.sim$L,
#'                      nchains=1, max.iter=1)
#' tsne.res = celdaTsne(celda.sim$counts, celda.mod)
>>>>>>> upstream-devel
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
<<<<<<< HEAD
#' @export
featureModuleLookup = function(counts, celda.mod, feature){
  class(celda.mod) = c(class(celda.mod), celda.mod)
  UseMethod("featureModuleLookup", celda.mod)
}
=======
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
>>>>>>> upstream-devel
