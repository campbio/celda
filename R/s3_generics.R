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


#' Get the complete history of gene / cell / gene & cell cluster assignments 
#' generated during a celda run, dependent on the model provided.
#'
#' @param celda.mod A celda model object (of class "celda_C", "celda_G" or "celda_CG")
#' @export
completeClusterHistory = function(celda.mod) {
  UseMethod("completeClusterHistory", celda.mod)
}


#' Get the probability of the cluster assignments generated during a celda run.
#'
#' @param celda.mod A celda model object (of class "celda_C", "celda_G" or "celda_CG")
#' @export
clusterProbabilities = function(celda.mod) {
  UseMethod("clusterProbabilities", celda.mod)
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
celda_heatmap <- function(celda.mod, counts, ...) {
  UseMethod("celda_heatmap", celda.mod)
}


#' Visualize various performance metrics as a function of K / L to aid parameter choice.
#' 
#' @param celda.list A celda_list object as returned from *celda()*
#' @param method Which performance metric to visualize. One of ("perplexity", "harmonic", "loglik"). "perplexity" calculates the inverse of the geometric mean of the log likelihoods from each iteration of Gibbs sampling. "harmonic" calculates the marginal likelihood has the harmonic mean of the likelihoods. "loglik" plots the highest log-likelihood during Gibbs iteration.
#' @param title Title for the visualize_model_performance
#' @param log Set log to TRUE to visualize the log(perplexity) of Celda_CG objects. Does not work for "harmonic" metric
#' @return A ggplot object containing the requested plot(s)
#' @export
visualize_model_performance <- function(celda.list, method, title, log = FALSE) {
  # Dispatch on the list's content type
  UseMethod("visualize_model_performance", celda.list$res.list[[1]])
}


#' Calculate a log-likelihood for a user-provided cluster assignment and count matrix, per the desired celda model. 
#' 
#' @param model Model to use for calculating log-likelihood of assignments; one of ("celda_C", "celda_CG", "celda_G")
#' @return The log-likelihood of the provided cluster assignment for the provided counts matrix.
#' @param ... extra parameters passed onto calculate_loglik_from_variables
#' @export
calculate_loglik_from_variables <- function(model, ...) {
  # Dispatch on the specified model
  if (model == "celda_C") calculate_loglik_from_variables.celda_C(...)
  else if (model == "celda_G") calculate_loglik_from_variables.celda_G(...)
  else if (model == "celda_CG") calculate_loglik_from_variables.celda_CG(...)
  else stop("Invalid model specified.")
}


#' Simulate count data from the celda generative models.
#' 
#' This function generates a list containing a simulated counts matrix, as well as various parameters
#' used in the simulation which can be useful for running celda. The user must provide the desired model
#' (one of celda_C, celda_G, celda_CG) as well as any desired tuning parameters for those model's simulation
#' functions as detailed below.
#' 
#' 
#' @param S Total number of samples (celda_C, celda_CG)
#' @param C The number of cells (celda_G)
#' @param C.Range Vector of length 2 given the range (min,max) of number of cells for each sample to be randomly generated from the uniform distribution (celda_C, celda_CG)
#' @param N.Range Vector of length 2 given the range (min,max) of number of counts for each cell to be randomly generated from the uniform distribution (all model types)
#' @param G Total number of Genes to be simulated (celda_C, celda_CG)
#' @param K An integer or range of integers indicating the desired number of cell clusters (celda_C, celda_CG)
#' @param L The number of transcriptional states (celda_G, celda_CG)
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution (celda_C, celda_CG)
#' @param beta Non-zero concentration parameter for gene Dirichlet distribution (all models)
#' @param delta The Dirichlet distribution parameter for Psi; adds a pseudocount to each gene within each transcriptional state (celda_G, celda_CG)
#' @param gamma The Dirichlet distribution parameter for Psi; adds a pseudocount to each gene within each transcriptional state (celda_G, celda_CG)
#' @param seed Parameter to set.seed() for random number generation (all models)
#' 
#' @export
simulateCells = function(model, ...) {
  switch(model, celda_C = simulateCells.celda_C(...),
         celda_G = simulateCells.celda_G(...),
         celda_CG = simulateCells.celda_CG(...),
         stop("Invalid model specified"))
}
