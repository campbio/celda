################################################################################
# S3 Methods for celda_list objects                                            #
################################################################################


#' Select models from a list of results generated in a celda run.
#' 
#' Convenience function for picking out specific models from a celda_list object.
#' Models can be selected by various parameters, most importantly the K/L parameters (number of cell
#'  clusters / number of feature clusters). 
#' 
#' @param celda.list Object of class "celda_list". An object containing celda models returned from `celdaGridSearch()`.
#' @param K The K parameter for the desired model in the results list. Matches all K by default. Accepts ranges.
#' @param L The L parameter for the desired model in the results list. Matches all L by default. Accepts ranges.
#' @param chain The desired chain(s) for the specified model, for the specified K/L. Matches all chains by default. Accepts ranges.
#' @param index The index of the desired model in the run.parameters in the provided celda_list. Overrides all other parameters if provided. Defaults to NULL.
#' @return A celda model object matching the provided parameters, or a list of celda model objects if multiple models were matched (of class "celda_C", "celda_G", "celda_CG" accordingly), or NA if one is not found.
#' @examples
#' celda.mods = celda_CG(celda::pbmc_select, K.to.test=c(10,15), 
#'                       L=c(50, 100), max.iter=2, nchains=1)
#' desired.mod = filterCeldaList(celda.mods, K=10, L=100)
#' @export
filterCeldaList = function(celda.list, K=c(), L=c(), chain=c(), index=NULL) {
  validateGetModelParams(celda.list, K, L, chain) 
  
  # If user provides index / range of indices to select, override all else
  if (!is.null(index)) {
    return(celda.list$res.list[index])
  }
   
  # Ensure we have accurate run.params, in case the res.list or run.params
  # was modified
  if (isTRUE(validateRunParams(celda.list))) {
    filtered.run.params = celda.list$run.params
  } else {
    filtered.run.params = newRunParamsFromResList(celda.list)
  }
  
  # Filter the run params to find matching indices.
  # Is there a more concise way to do this ..?
  if (length(K) > 0) {
    filtered.run.params = filtered.run.params[filtered.run.params$K %in% K, ]
  }
  if (length(L) > 0) {
    filtered.run.params = filtered.run.params[filtered.run.params$L %in% L, ]
  }
  if (length(chain) > 0) {
    filtered.run.params = filtered.run.params[filtered.run.params$chain %in% chain, ]
  }
  return(celda.list$res.list[filtered.run.params$index])
}


#' Select the best model from a celda_list by final log-likelihood.
#' 
#' This function returns the celda model from a celda_list object containing
#' the maxiumim final log-likelihood. K, L, or combination of K and L must be 
#' provided for celda_list objects containing celda_C, celda_G, and celda_CG 
#' models, respectively.
#' 
#' @param celda.list Object of class "celda_list". An object containing celda models returned from `celdaGridSearch()`.
#' @param K Limit search for best model to models with this number of cell clusters.
#' @param L Limit search for best model to models with this number of feature clusters.
#' @return The celda model object with the highest finalLogLik attribute, meeting any K/L criteria provided
#' @examples
#' celda.mods = celda_CG(celda::pbmc_select, K.to.test=c(5,10), 
#'                       L=c(10,20,30), max.iter=2, nchains=1)
#' best.mod.k5.l20 = selectBestModel(celda.mods, K=5, L=20)
#' @export
selectBestModel = function(celda.list, K=c(), L=c()) {
  if (class(celda.list)[1] != "celda_list") {
    stop("Provided object is not of class celda_list")
  }
 
  matching.models = filterCeldaList(celda.list, K=K, L=L)
  
  logliks = unlist(sapply(matching.models, function(mod) { mod[["finalLogLik"]] }))
  max.idx = which(logliks == max(logliks, na.rm=TRUE))
  
  return(matching.models[[max.idx[1]]])
}


validateGetModelParams = function(celda.list, K, L, chain) {
  if (class(celda.list)[1] != "celda_list") stop("First argument to getModel() should be an object of class 'celda_list'")
  
  if ((is.null(K) | is.null(L)) & celda.list$content.type == "celda_CG") {
    stop("Both K and L parameters needed for subsetting celda_CG result lists")
  }

  if (is.null(K) & celda.list$content.type == "celda_C") {
    stop("K parameter needed when subsetting celda_C result lists")
  }

  if (is.null(L) & celda.list$content.type == "celda_G") {
    stop("L parameter needed when subsetting celda_G result lists")
  }
  
  if (!is.null(K)){
    if(!(K %in% celda.list$run.params$K)) {
      stop("Provided K was not profiled in the provided celda_list object")
    } 
  }
  
  if (!is.null(L)){ 
    if(!(L %in% celda.list$run.params$L)) {
      stop("Provided L was not profiled in the provided celda_list object")
    }
  }
}
