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
#' @param celda.res A result object returned from celda()
#' @export
runParams = function(celda.res) {
  return(celda.res$run.params)
}


#' Get seeds for all chains from a celda run.
#'
#' @param celda.res A result object returned from celda()
#' @export
chainSeeds = function(celda.res) {
  seeds = sapply(celda.res$res.list,
                 function(chain) { chain$seed })
  return(seeds)
}


#' Get the complete log likelihood for each iteration of Gibbs sampling
#' for each chain generated during a celda run.
#'
#' @param celda.res A result object returned from celda()
#' @export
completeLogLikelihood = function(celda.res) {
  complete.log.lik = sapply(celda.res$res.list, 
                            function(chain) { chain$completeLogLik })
  return(complete.log.lik)
}


#' Get the log likelihood from the final iteration of Gibbs sampling
#' for each chain generated during a celda run.
#'
#' @param celda.res A result object returned from celda()
#' @export
finalLogLikelihood = function(celda.res) {
  final.log.lik = sapply(celda.res$res.list,
                         function(chain) { chain$finalLogLik })
  return(final.log.lik)
}


#' Get the final gene / cell / gene & cell cluster assignments generated during
#' a celda run, dependent on the model provided.
#'
#' @param celda.res A result object returned from celda()
#' @export
finalClusterAssignment = function(celda.res) {
  UseMethod("finalClusterAssignment", celda.res)
}


#' Get the complete history of gene / cell / gene & cell cluster assignments 
#' generated during a celda run, dependent on the model provided.
#'
#' @param celda.res A result object returned from celda()
#' @export
completeClusterHistory = function(celda.res) {
  UseMethod("completeClusterHistory", celda.res)
}


#' Get the probability of the cluster assignments generated during a celda run.
#'
#' @param celda.res A result object returned from celda()
#' @export
clusterProbabilities = function(celda.res) {
  UseMethod("clusterProbabilities", celda.res)
}


#' Get the K value used for each chain in a celda run.
#'
#' @param celda.res A result object returned from celda()
#' @export
chainKs = function(celda.res) {
  UseMethod("chainKs", celda.res)
}


#' Get the L value used for each chain in a celda run.
#'
#' @param celda.res A result object returned from celda()
#' @export
chainLs = function(celda.res) {
  UseMethod("chainLs", celda.res)
}


#' Render a stylable heatmap of count data based on celda clustering results.
#'
#' @param celda.res A result object returned from celda()
#' plot the heatmap of the counts data
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
celda_heatmap <- function(celda.res, counts, ...) {
  UseMethod("celda_heatmap", celda.res)
}
 

################################################################################
# celda_C                                                                      #
################################################################################
#' @export
finalClusterAssignment.celda_C = function(celda.res) {
  final.assignments = sapply(celda.res$res.list, 
                             function(chain) { chain$z })
  return(final.assignments)
}


#' @export
completeClusterHistory.celda_C = function(celda.res) {
  complete.history = sapply(celda.res$res.list, 
                            function(chain) { chain$complete.z })
  return(complete.history)
}


#' @export
clusterProbabilities.celda_C = function(celda.res) {
  cluster.probs = sapply(celda.res$res.list,
                         function(chain) { chain$z.probability })
  return(cluster.probs)
}


#' @export
chainKs.celda_C = function(celda.res) {
  cluster.probs = sapply(celda.res$res.list,
                         function(chain) { chain$K })
  return(cluster.probs)
}


#' @export
chainLs.celda_C = function(celda.res) { return(NA) }


#' @export
celda_heatmap.celda_C = function(celda.res, counts, ...) {
  render_celda_heatmap(counts, z=celda.res$z, ...)
}



################################################################################
# celda_G                                                                      #
################################################################################
#' @export
finalClusterAssignment.celda_G = function(celda.res) {
  final.assignments = sapply(celda.res$res.list, 
                             function(chain) { chain$y})
  return(final.assignments)
}


#' @export
completeClusterHistory.celda_G = function(celda.res) {
  complete.history = sapply(celda.res$res.list, 
                            function(chain) { chain$complete.y })
  return(complete.history)
}


#' @export
clusterProbabilities.celda_G = function(celda.res) {
  cluster.probs = sapply(celda.res$res.list,
                         function(chain) { chain$y.probability })
  return(cluster.probs)
}


#' @export
chainKs.celda_G = function(celda.res) { return(NA) }


#' @export
chainLs.celda_G = function(celda.res) {
  cluster.probs = sapply(celda.res$res.list,
                         function(chain) { chain$L })
  return(cluster.probs)
}


#' @export
celda_heatmap.celda_G = function(celda.res, counts, ...) {
  render_celda_heatmap(counts, y=celda.res$y, ...)
}



################################################################################
# celda_CG                                                                      #
################################################################################
#' @export
finalClusterAssignment.celda_CG = function(celda.res) {
  z.assignments = sapply(celda.res$res.list, function(chain) { chain$z })
  y.assignments = sapply(celda.res$res.list, function(chain) { chain$y })
  return(list(z.assignments=z.assignments, y.assignments=y.assignments))
}


#' @export
completeClusterHistory.celda_CG = function(celda.res) {
  complete.z = sapply(celda.res$res.list, function(chain) { chain$complete.z })
  complete.y = sapply(celda.res$res.list, function(chain) { chain$complete.y })
  return(list(complete.z=complete.z, complete.y=complete.y))
}


#' @export
clusterProbabilities.celda_CG = function(celda.res) {
  z.prob = sapply(celda.res$res.list, function(chain) { chain$z.prob })
  y.prob = sapply(celda.res$res.list, function(chain) { chain$y.prob })
  return(list(z.prob=z.prob, y.prob=y.prob))
}


#' @export
chainKs.celda_CG = function(celda.res) {
  cluster.probs = sapply(celda.res$res.list,
                         function(chain) { chain$K })
  return(cluster.probs)
}


#' @export
chainLs.celda_CG = function(celda.res) {
  cluster.probs = sapply(celda.res$res.list,
                         function(chain) { chain$L })
  return(cluster.probs)
}


#' @export
celda_heatmap.celda_CG = function(celda.res, counts, ...) {
  render_celda_heatmap(counts, z=celda.res$z, y=celda.res$y, ...)
}
