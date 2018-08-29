################################################################################
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
# Generics
################################################################################

#' Get run parameters for a celda run.
#'
#' @param celda.list Object of class "celda_list". An object containing celda models returned from `celdaGridSearch()`.
#' @export
runParams = function(celda.list) {
  return(celda.list$run.params)
}


#' Get the complete log likelihood for a given celda model.
#'
#' @param celda.mod Celda object of class "celda_C", "celda_G", or "celda_CG".
#' @export
completeLogLikelihood = function(celda.mod) {
  return(celda.mod$completeLogLik)
}


#' Get the log likelihood from the final iteration of Gibbs sampling
#' for a given celda model.
#'
#' @param celda.mod Celda object of class "celda_C", "celda_G", or "celda_CG".
#' @export
finalLogLikelihood = function(celda.mod) {
  return(celda.mod$finalLogLik)
}


#' Get the probability of the cluster assignments generated during a celda run.
#'
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param log Logical. If FALSE, then the normalized conditional probabilities will be returned. If TRUE, then the unnormalized log probabilities will be returned. Default FALSE.  
#' @export
clusterProbability = function(celda.mod, counts, log=FALSE) {
  UseMethod("clusterProbability", celda.mod)
}


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


#' Render a stylable heatmap of count data based on celda clustering results.
#'
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`. 
#' @param celda.mod Celda object of class "celda_C", "celda_G", or "celda_CG".
#' @param ... Additional parameters.
#' @export 
celdaHeatmap <- function(counts, celda.mod, ...) {
  UseMethod("celdaHeatmap", celda.mod)
}


#' Calculate a log-likelihood for a user-provided cluster assignment and count matrix, per the desired celda model. 
#' 
#' @param counts The counts matrix used to generate the provided cluster assignments.
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @return The log-likelihood of the provided cluster assignment for the provided counts matrix.
#' @param ... Additional parameters.
#' @export
calculateLoglikFromVariables <- function(counts, celda.mod, model, ...) {
  class(counts) = c(model)
  do.call(paste("calculateLoglikFromVariables.", model, sep=""),
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
#' @export
factorizeMatrix = function(counts, celda.mod, type) {
  
  UseMethod("factorizeMatrix", celda.mod)
}

#' Generate factorized matrices showing each feature's influence on cell / gene clustering
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class "celda_C" or "celda_CG".
#' @param ... Additional parameters.
#' @export
celdaProbabilityMap = function(counts, celda.mod, ...) {
  
  UseMethod("celdaProbabilityMap", celda.mod)
}
#' Runs tSNE via Rtsne based on the CELDA model and specified cell states.
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @param ... Additional parameters.
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
#' @export
featureModuleLookup = function(counts, celda.mod, feature){
  class(celda.mod) = c(class(celda.mod), celda.mod)
  UseMethod("featureModuleLookup", celda.mod)
}
