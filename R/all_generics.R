setClass("celdaModel", 
         representation(params = "list", # K, L, model priors, seed, checksum
                        names = "list",
                        completeLogLik = "numeric", 
                        finalLogLik = "numeric",
                        clusters = "list"))  # z and or y


#' @title Get parameter values provided for celda model creation
#' @description Retrieves the K/L, model priors (e.g. alpha, beta), random seed, and count matrix checksum parameters provided during the creation of the provided celda model.
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @return List. Contains the model-specific parameters for the provided celda model object depending on its class.
#' @examples
#' params(celda.CG.mod)
#' @export
setGeneric("params",
           function(celda.mod){ standardGeneric("params") })
#' @title Get parameter values provided for celda model creation
#' @description Retrieves the K/L, model priors (e.g. alpha, beta), random seed, and count matrix checksum parameters provided during the creation of the provided celda model.
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @return List. Contains the model-specific parameters for the provided celda model object depending on its class.
#' @examples
#' params(celda.CG.mod)
#' @export
setMethod("params",
          signature=c(celda.mod="celdaModel"),
          function(celda.mod){  celda.mod@params  })


#' @title Get feature, cell and sample names from a celda model
#' @description Retrieves the row, column, and sample names used to generate a celda model.
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @return List. Contains row, column, and sample character vectors corresponding to the values provided when the celda model was generated.
#' @examples
#' matrixNames(celda.CG.mod)
#' @export
setGeneric("matrixNames",
           function(celda.mod){ standardGeneric("matrixNames") })
#' @title Get feature, cell and sample names from a celda model
#' @description Retrieves the row, column, and sample names used to generate a celda model.
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @return List. Contains row, column, and sample character vectors corresponding to the values provided when the celda model was generated.
#' @examples
#' matrixNames(celda.CG.mod)
#' @export
setMethod("matrixNames",
          signature=c(celda.mod="celdaModel"),
          function(celda.mod){  celda.mod@names  })


#' @title Get log-likelihood history
#' @description Retrieves the complete log-likelihood from all iterations of Gibbs sampling used to generate a celda model.
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @return Numeric. The log-likelihood at each step of Gibbs sampling used to generate the model.
#' @examples
#' logLikelihoodHistory(celda.CG.mod)
#' @export
setGeneric("logLikelihoodHistory",
           function(celda.mod){ standardGeneric("logLikelihoodHistory") })
#' @title Get log-likelihood history
#' @description Retrieves the complete log-likelihood from all iterations of Gibbs sampling used to generate a celda model.
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @return Numeric. The log-likelihood at each step of Gibbs sampling used to generate the model.
#' @examples
#' logLikelihoodHistory(celda.CG.mod)
#' @export
setMethod("logLikelihoodHistory",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod){  celda.mod@completeLogLik  })


#' @title Get the log-likelihood 
#' @description Retrieves the final log-likelihood from all iterations of Gibbs sampling used to generate a celda model.
#' @return Numeric. The log-likelihood at the final step of Gibbs sampling used to generate the model.
#' @param celda.mod A celda model object of class celda_C, celda_G, or celda_CG.
#' @examples
#' bestLogLikelihood(celda.CG.mod)
#' @export
setGeneric("bestLogLikelihood",
           function(celda.mod){ standardGeneric("bestLogLikelihood") })
#' @title Get the log-likelihood 
#' @description Retrieves the final log-likelihood from all iterations of Gibbs sampling used to generate a celda model.
#' @param celda.mod A celda model object of class celda_C, celda_G, or celda_CG.
#' @return Numeric. The log-likelihood at the final step of Gibbs sampling used to generate the model.
#' @examples
#' bestLogLikelihood(celda.CG.mod)
#' @export
setMethod("bestLogLikelihood",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod){  celda.mod@finalLogLik  })


#' @title Get clustering outcomes from a celda model
#' @description Returns the z / y results corresponding to the cell / gene cluster labels determined by the provided celda model.
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @return List. Contains z (for celda_C and celda_CG models) and/or y (for celda_G and celda_CG models)
#' @examples
#' clusters(celda.CG.mod)
#' @export
setGeneric("clusters",
           function(celda.mod){ standardGeneric("clusters")})
#' @title Get clustering outcomes from a celda model
#' @description Returns the z / y results corresponding to the cell / gene cluster labels determined by the provided celda model.
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @return List. Contains z (for celda_C and celda_CG models) and/or y (for celda_G and celda_CG models)
#' @examples
#' clusters(celda.CG.mod)
#' @export
setMethod("clusters", signature=c(celda.mod="celdaModel"),
          function(celda.mod){
            return(celda.mod@clusters)
          })


setClass("celda_C",
         representation(sample.label = "factor"),
         contains = "celdaModel")


#' @title Get sample labels from a celda model
#' @description Returns the sample labels for the count matrix provided for generation of a given celda model.
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @return Character. Contains the sample labels provided at model creation time, or those automatically generated by celda.
#' @examples
#' sampleLabel(celda.CG.mod)
#' @export
setGeneric("sampleLabel",
           function(celda.mod){ standardGeneric("sampleLabel") })
#' @title Get sample labels from a celda model
#' @description Returns the sample labels for the count matrix provided for generation of a given celda model.
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @return Character. Contains the sample labels provided at model creation time, or those automatically generated by celda.
#' @examples
#' sampleLabel(celda.CG.mod)
#' @export
setMethod("sampleLabel",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod){  celda.mod@sample.label  })


setClass("celda_G",
         contains = "celdaModel")

setClass("celda_CG",
         contains = c("celda_C", "celda_G"))

setClass("celdaList",
         representation(run.params = "data.frame",
                        res.list = "list",
                        count.checksum = "character",
                        perplexity = "matrix"))


#' @title Get run parameters provided to `celdaGridSearch()`
#' @description Returns details on the clustering parameters, model priors, and seeds provided to `celdaGridSearch()` when the provided celdaList was created.
#' @param celda.mod An object of class celdaList.
#' @return Data Frame. Contains details on the various K/L parameters, chain parameters, and final log-likelihoods derived for each model in the provided celdaList.
#' @examples
#' runParams(celda.CG.grid.search.res)
#' @export
setGeneric("runParams",
           function(celda.mod){ standardGeneric("runParams") })
#' @title Get run parameters provided to `celdaGridSearch()`
#' @description Returns details on the clustering parameters, model priors, and seeds provided to `celdaGridSearch()` when the provided celdaList was created.
#' @param celda.mod An object of class celdaList.
#' @return Data Frame. Contains details on the various K/L parameters, chain parameters, and final log-likelihoods derived for each model in the provided celdaList.
#' @examples
#' runParams(celda.CG.grid.search.res)
#' @export
setMethod("runParams",
           signature=c(celda.mod="celdaList"),
           function(celda.mod){  celda.mod@run.params  })


#' @title Get final celda models from a celdaList
#' @description Returns all models generated during a `celdaGridSearch()` run.
#' @param celda.mod An object of class celdaList.
#' @return List. Contains one celdaModel object for each of the parameters specified in the `runParams()` of the provided celda list.
#' @examples
#' celda.CG.grid.models = resList(celda.CG.grid.search.res)
#' @export
setGeneric("resList",
           function(celda.mod){ standardGeneric("resList") })
#' @title Get final celda models from a celdaList
#' @description Returns all models generated during a `celdaGridSearch()` run.
#' @param celda.mod An object of class celdaList.
#' @return List. Contains one celdaModel object for each of the parameters specified in the `runParams()` of the provided celda list.
#' @examples
#' celda.CG.grid.models = resList(celda.CG.grid.search.res)
#' @export
setMethod("resList",
           signature=c(celda.mod="celdaList"),
           function(celda.mod){  celda.mod@res.list  })


#' @title Get perplexity for every model in a celdaList
#' @description Returns perplexity for each model in a celdaList as calculated by `perplexity().`
#' @param celda.mod A celda model object of class "celda_C", "celda_G", or "celda_CG".
#' @return List. Contains one celdaModel object for each of the parameters specified in the `runParams()` of the provided celda list.
#' @examples
#' celda.CG.grid.model.perplexities = celdaPerplexity(celda.CG.grid.search.res)
#' @export
setGeneric("celdaPerplexity",
           function(celda.mod){ standardGeneric("celdaPerplexity") })
#' @title Get perplexity for every model in a celdaList
#' @description Returns perplexity for each model in a celdaList as calculated by `perplexity().`
#' @param celda.mod A celda model object of class "celda_C", "celda_G", or "celda_CG".
#' @return List. Contains one celdaModel object for each of the parameters specified in the `runParams()` of the provided celda list.
#' @examples
#' celda.CG.grid.model.perplexities = celdaPerplexity(celda.CG.grid.search.res)
#' @export
setMethod("celdaPerplexity",
           signature=c(celda.mod="celdaList"),
           function(celda.mod){  celda.mod@perplexity  })


#' @title Append two celdaList objects
#' @description Returns a single celdaList representing the combination of two provided celdaList objects.
#' @return A celdaList object. This object contains all resList entries and runParam records from both lists.
#' @param list1 A celda_list object
#' @param list2 A celda_list object to be joined with list_1
#' @examples
#' appended.list = appendCeldaList(celda.CG.grid.search.res, celda.CG.grid.search.res)
#' @export
appendCeldaList = function(list1, list2) {
  if (!is.element("celdaList", class(list1)) | !is.element("celdaList", class(list2))) {
    stop("Both parameters to appendCeldaList must be of class celdaList.")
  }
  if (!(list1@count.checksum == list2@count.checksum)) {
    warning("Provided lists have different count.checksums and may have been generated from different count matrices. Using checksum from first list...")
  }
  newList = methods::new("celdaList",
                         run.params = rbind(list1@run.params, list2@run.params),
                         res.list = c(list1@res.list, list2@res.list),
                         count.checksum = list1@count.checksum,
                         perplexity = matrix(nrow=0, ncol=0))
  return(newList)
}

################################################################################
# Generics
################################################################################


#' Render a stylable heatmap of count data based on celda clustering results.
#'
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`. 
#' @param celda.mod A celda model object of class "celda_C", "celda_G", or "celda_CG".
#' @param feature.ix Integer vector. Select features for display in heatmap. If NULL, no subsetting will be performed. Default NULL.
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
#' @param ... Additional parameters.
#' @examples
#' cluster.prob = clusterProbability(celda.CG.sim$counts, celda.CG.mod)
#' @return A numeric vector of the cluster assignment probabilties
#' @export
setGeneric("clusterProbability", 
           signature="celda.mod",
           function(counts, celda.mod, log=FALSE, ...) {
             standardGeneric("clusterProbability")
           })


#' Calculate the perplexity from a single celda model
#' 
#' Perplexity can be seen as a measure of how well a provided set of 
#' cluster assignments fit the data being clustered.
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda model. Options available in `celda::available.models`.
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
#' @param S Integer. Number of samples to simulate. Default 5. celda_C and celda_CG simulations only.
#' @param C Integer. Number of cells to simulate. Default 100. celda_G simulations only.
#' @param C.Range Vector of length 2 given the range (min,max) of number of cells for each sample to be randomly generated from the uniform distribution. Default c(50, 100). celda_C and celda_CG simulations only.
#' @param N.Range Integer vector. A vector of length 2 that specifies the lower and upper bounds of the number of counts generated for each cell. Default c(500, 1000). 
#' @param G Integer. The total number of features to be simulated. Default 100.
#' @param K Integer. Number of cell populations. Default 5. celda_C and celda_CG simulations only.
#' @param L Integer. Number of feature modules. Default 10.
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount to each cell population in each sample. Default 1. celda_C and celda_CG simulations only.
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to each feature in each cell population. Default 1. 
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to each feature in each module. Default 1. celda_G and celda_CG simulations only. 
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to the number of features in each module. Default 5. celda_G and celda_CG simulations only. 
#' @param seed Integer. Passed to `set.seed()`. Default 12345. If NULL, no calls to `set.seed()` are made.
#' 
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
#' @param celda.mod Celda object of class `celda_CG`. 
#' @param max.cells Integer. Maximum number of cells to plot. Cells will be randomly subsampled if ncol(counts) > max.cells. Larger numbers of cells requires more memory. Default 25000.
#' @param min.cluster.size Integer. Do not subsample cell clusters below this threshold. Default 100.
#' @param initial.dims integer. The number of dimensions that should be retained in the initial PCA step. Default 20.
#' @param modules Integer vector. Determines which features modules to use for tSNE. If NULL, all modules will be used. Default NULL.
#' @param perplexity Numeric. Perplexity parameter for tSNE. Default 20.
#' @param max.iter Integer. Maximum number of iterations in tSNE generation. Default 2500.
#' @param seed Integer. Passed to `set.seed()`. Default 12345. If NULL, no calls to `set.seed()` are made.
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


#' Embeds cells in two dimensions using umap.
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class `celda_CG`. 
#' @param max.cells Integer. Maximum number of cells to plot. Cells will be randomly subsampled if ncol(counts) > max.cells. Larger numbers of cells requires more memory. Default 25000.
#' @param min.cluster.size Integer. Do not subsample cell clusters below this threshold. Default 100. 
#' @param modules Integer vector. Determines which features modules to use for tSNE. If NULL, all modules will be used. Default NULL.
#' @param umap.config An object of class "umap.config" specifying parameters to the UMAP algorithm.
#' @param ... Additional parameters.
#' @return Numeric Matrix of dimension `ncol(counts)` x 2, colums representing the "X" and "Y" coordinates in the data's t-SNE represetation.
#' @examples 
#' tsne.res = celdaUmap(celda.CG.sim$counts, celda.CG.mod)
#' @export
setGeneric("celdaUmap",
           signature = "celda.mod",
           function(counts, celda.mod, max.cells=25000, min.cluster.size=100,
                   modules=NULL, umap.config=umap::umap.defaults) {
             standardGeneric("celdaUmap")
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
