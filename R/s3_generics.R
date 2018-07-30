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
#' @param celda.list A celda_list object, as returned from celda()
#' @export
runParams = function(celda.list) {
  return(celda.list$run.params)
}


#' Get the random seed for a given celda model.
#'
#' @param celda.mod A celda model object (of class "celda_C", "celda_G" or "celda_CG")
#' @export
seed = function(celda.mod) {
  return(celda.mod$seed)
}


#' Get the complete log likelihood for a given celda model.
#'
#' @param celda.mod A celda model object (of class "celda_C", "celda_G" or "celda_CG")
#' @export
completeLogLikelihood = function(celda.mod) {
  return(celda.mod$completeLogLik)
}


#' Get the log likelihood from the final iteration of Gibbs sampling
#' for a given celda model.
#'
#' @param celda.mod A celda model object (of class "celda_C", "celda_G" or "celda_CG")
#' @export
finalLogLikelihood = function(celda.mod) {
  return(celda.mod$finalLogLik)
}


#' Get the final gene / cell / gene & cell cluster assignments generated during
#' a celda run, dependent on the model provided.
#'
#' @param celda.mod A celda model object (of class "celda_C", "celda_G" or "celda_CG")
#' @export
finalClusterAssignment = function(celda.mod) {
  UseMethod("finalClusterAssignment", celda.mod)
}


#' Get the probability of the cluster assignments generated during a celda run.
#'
#' @param celda.mod A celda model object
#' @param counts The count matrix used to generate the model
#' @param log If FALSE, then the normalized conditional probabilities will be returned. If TRUE, then the unnormalized log probabilities will be returned.  
#' @export
clusterProbability = function(celda.mod, counts, log=FALSE) {
  UseMethod("clusterProbability", celda.mod)
}


#' Get the K value used for each chain in a celda run.
#'
#' @param celda.mod A celda model object (of class "celda_C", "celda_G" or "celda_CG")
#' @export
getK = function(celda.mod) {
  UseMethod("getK", celda.mod)
}


#' Get the L value used for each chain in a celda run.
#'
#' @param celda.mod A celda model object (of class "celda_C", "celda_G" or "celda_CG")
#' @export
getL = function(celda.mod) {
  UseMethod("getL", celda.mod)
}


#' Render a stylable heatmap of count data based on celda clustering results.
#'
#' @param celda.mod A celda model object (of class "celda_C", "celda_G" or "celda_CG")
#' @param counts the counts matrix 
#' @param ... extra parameters passed onto celda_heatmap
#' @export 
celdaHeatmap <- function(celda.mod, counts, ...) {
  UseMethod("celdaHeatmap", celda.mod)
}


#' Calculate a log-likelihood for a user-provided cluster assignment and count matrix, per the desired celda model. 
#' 
#' @param model Model to use for calculating log-likelihood of assignments; one of ("celda_C", "celda_CG", "celda_G")
#' @param counts The counts matrix used to generate the provided cluster assignments
#' @return The log-likelihood of the provided cluster assignment for the provided counts matrix.
#' @param ... extra parameters passed onto calculateLoglikFromVariables
#' @export
calculateLoglikFromVariables <- function(model, counts, ...) {
  class(counts) = c(model, class(counts))
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
#' @param model The celda generative model to use (one of celda_C, celda_G, celda_CG)
#' @param ... Parameters to pass to underlying generative model simulation
#' @export
simulateCells = function(model, ...) {
  class(model) = c(class(model), model)
  UseMethod("simulateCells", model)
}


#' Generate factorized matrices showing each feature's influence on cell / gene clustering
#' 
#' @param counts A numeric count matrix
#' @param celda.mod An object from a celda_list's res.list property
#' @param type A character vector containing one or more of "counts", "proportions", or "posterior". "counts" returns the raw number of counts for each entry in each matrix. "proportions" returns the counts matrix where each vector is normalized to a probability distribution. "posterior" returns the posterior estimates which include the addition of the Dirichlet concentration parameter (essentially as a pseudocount).
#' @param validate.counts Whether to verify that the counts matrix provided was used to generate the results in celda.mod. Defaults to TRUE.
#' @export
factorizeMatrix = function(counts, celda.mod, type, validate.counts=TRUE) {
  counts = processCounts(counts)  # Ensure counts are integer and have corresponding storage mode
  
  UseMethod("factorizeMatrix", celda.mod)
}


#' Runs tSNE via Rtsne based on the CELDA model and specified cell states.
#' 
#' @param counts Counts matrix, should have cell name for column name and gene name for row name.
#' @param celda.mod Celda model to use for tsne. class "celda_C","celda_G" or "celda_CG".
#' @param ... Other arguments to be passed to model-specific tSNE functions. Use methods("celdaTsne") for list of 'celdaTsne' functions.
#' @export
celdaTsne = function(counts, celda.mod, ...) {
  compareCountMatrix(counts, celda.mod)
  if (!isTRUE(class(celda.mod) %in% c("celda_CG","celda_C","celda_G"))) {
    stop("celda.mod argument is not of class celda_C, celda_G or celda_CG")
  }
  UseMethod("celdaTsne", celda.mod)
}
