################################################################################
# S3 Methods                                                                   #
################################################################################
# Below are getters for the various types of celda models. They expect output  #
# in the format provided by the celda() wrapper function in celda.R, *NOT* as  # 
# provided in the individual model functions (e.g. celda_C()).                 #
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
#' @param celda.res A celda_list object, as returned from celda()
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
#' @param z A numeric vector of cluster assignments for cell. Resolved automatically from celda object when available.
#' @param y A numeric vector of cluster assignments for gene. Resolved automatically from celda object when available.
#' @param scale.log specify the transformation type of the matrix for (semi-)heatmap, can be "log","row"(z-acore by row),"col"(z-score by column), etc. #To be completed
#' @param scale.row specify the transformation type of the matrix for (semi-)heatmap, can be "log","row"(z-acore by row),"col"(z-score by column), etc. #To be completed
#' @param z.trim two element vector to specify the lower and upper cutoff of the z-score normalization result by default it is set to NULL so no trimming will be done.
#' @param scale_fun specify the function for scaling 
#' @param cluster.row boolean values determining if rows should be clustered
#' @param cluster.column boolean values determining if columns should be clustered
#' @export 
celda_heatmap <- function(celda.mod, counts, ...) {
  UseMethod("celda_heatmap", celda.mod)
}
 

################################################################################
# celda_C                                                                      #
################################################################################
#' @export
finalClusterAssignment.celda_C = function(celda.mod) {
  return(celda.mod$z)
}


#' @export
completeClusterHistory.celda_C = function(celda.mod) {
  return(celda.mod$complete.z)
}


#' @export
clusterProbabilities.celda_C = function(celda.mod) {
  return(celda.mod$z.probability)
}


#' @export
getK.celda_C = function(celda.mod) {
  return(celda.mod$K)
}


#' @export
getL.celda_C = function(celda.mod) { return(NA) }


#' @export
celda_heatmap.celda_C = function(celda.mod, counts, ...) {
  render_celda_heatmap(counts, z=celda.mod$z, ...)
}



################################################################################
# celda_G                                                                      #
################################################################################
#' @export
finalClusterAssignment.celda_G = function(celda.mod) {
  return(celda.mod$y)
}


#' @export
completeClusterHistory.celda_G = function(celda.mod) {
  return(celda.mod$complete.y)
}


#' @export
clusterProbabilities.celda_G = function(celda.mod) {
  return(celda.mod$y.probability)
}


#' @export
getK.celda_G = function(celda.mod) { return(NA) }


#' @export
getL.celda_G = function(celda.mod) {
  return(celda.mod$L)
}


#' @export
celda_heatmap.celda_G = function(celda.mod, counts, ...) {
  render_celda_heatmap(counts, y=celda.mod$y, ...)
}



################################################################################
# celda_CG                                                                     #
################################################################################
#' @export
finalClusterAssignment.celda_CG = function(celda.mod) {
  return(list(z=celda.mod$z, y=celda.mod$y))
}


#' @export
completeClusterHistory.celda_CG = function(celda.mod) {
  return(list(complete.z=celda.mod$complete.z, complete.y=celda.mod$complete.y))
}


#' @export
clusterProbabilities.celda_CG = function(celda.mod) {
  return(list(z.prob=celda.mod$z.prob, y.prob=celda.mod$y.prob))
}


#' @export
getK.celda_CG = function(celda.mod) {
  return(celda.mod$K)
}


#' @export
getL.celda_CG = function(celda.mod) {
  return(celda.mod$L)
}


#' @export
celda_heatmap.celda_CG = function(celda.mod, counts, ...) {
  render_celda_heatmap(counts, z=celda.mod$z, y=celda.mod$y, ...)
}
