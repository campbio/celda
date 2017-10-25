################################################################################
# S3 Methods for celda_list objects                                            #
################################################################################

#' Run the celda Bayesian hierarchical model on a matrix of counts.
#' TODO: - If no chain provided, automatically choose best chain
#'       - Smarter subsetting of the run.params to DRY this function up
#' 
#' @param celda.list A celda_list object returned from celda()
#' @param K The K parameter for the desired model in the results list
#' @param L The L parameter for the desired model in the results list
#' @param chain The desired chain for the specified model
#' @param best Method for choosing best chain automatically. Options are c("perplexity", "loglik"). See documentation for chooseBestModel for details. Overrides chain parameter if provided.
#' @return A celda model object matching the provided parameters (of class "celda_C", "celda_G", "celda_CG" accordingly), or NA if one is not found.
#' @export
getModel = function(celda.list, K=NULL, L=NULL, chain=1, best=NULL) {
  validate_get_model_params(celda.list, K, L, chain, best)  # Sanity check params
  
  requested.chain = NA
  run.params = celda.list$run.params
  
  if (celda.list$content.type == "celda_CG") {
    if (!is.null(best)) {
      matching.chain.idx = run.params[run.params$K == K & run.params$L == L, "index"]
      requested.chain = chooseBestChain(celda.list$res.list[matching.chain.idx], best)
    } else {
      requested.chain.idx = run.params[run.params$K == K & run.params$L == L & run.params$chain == chain,
                                       "index"]
      requested.chain = celda.list$res.list[[requested.chain.idx]]
    }
  }
  
  
  if (celda.list$content.type == "celda_C") {
    if (!is.null(best)) {
      matching.chain.idx = run.params[run.params$K == K, "index"]
      requested.chain = chooseBestChain(celda.list$res.list[matching.chain.idx], best)
    } else {
      requested.chain.idx = run.params[run.params$K == K & run.params$chain == chain, "index"]
      requested.chain = celda.list$res.list[[requested.chain.idx]]
    }
  }
  
  
  if (celda.list$content.type == "celda_G") {
    if (!is.null(best)) {
      matching.chain.idx = run.params[run.params$L == L, "index"]
      requested.chain = chooseBestChain(celda.list$res.list[matching.chain.idx], best)
    } else { 
      requested.chain.idx = run.params[run.params$L == L & run.params$chain == chain, "index"]
      requested.chain = celda.list$res.list[[requested.chain.idx]]
    }
  }
  
  # Ensure that the chain we grabbed actually has the requested K/L.
  # This should only happen if the user alters the celda_list's run.params dataframe.
  if ((!is.null(K) & K != requested.chain$K) | !is.null(L) & L != requested.chain$L) {
    requested.chain = search_res_list(celda.list, K=K, L=L)
  }
  
  return(requested.chain)
}


validate_get_model_params = function(celda.list, K, L, chain, best) {
  if (class(celda.list) != "celda_list") stop("First argument to getModel() should be an object of class 'celda_list'")
  
  if ((is.null(K) | is.null(L)) & celda.list$content.type == "celda_CG") {
    stop("Both K and L parameters needed for subsetting celda_CG result lists")
  }

  if (is.null(K) & celda.list$content.type == "celda_C") {
    stop("K parameter needed when subsetting celda_C result lists")
  }

  if (is.null(L) & celda.list$content.type == "celda_G") {
    stop("L parameter needed when subsetting celda_G result lists")
  }
  
  if (!is.null(K) & !(K %in% celda.list$run.params$K)) {
    stop("Provided K was not profiled in the provided celda_list object")
  }
  
  if (!is.null(L) & !(L %in% celda.list$run.params$L)) {
    stop("Provided L was not profiled in the provided celda_list object")
  }
}


chooseBestChain = function(celda.mods, method="perplexity") {
  # We want to get the *most negative* perplexity, as opposed to the *least* negative
  # for the other metrics...
  if (method == "perplexity"){
    metrics = lapply(celda.mods, function(mod) { calculatePerplexity(mod$completeLogLik) })
    metrics = methods::new("mpfr", unlist(metrics))
    best = which(metrics == min(metrics))
    return(celda.mods[[best]])
  } 
  
  else if (method == "loglik"){
    metrics = lapply(celda.mods, function(mod) { max(mod$completeLogLik) })
    metrics = unlist(metrics)
  }  else {
    stop("Invalid method specified.")
  }
  best = which(metrics == max(metrics))
  if (length(best) > 1) best = best[1]  # Choose first chain if there's a tie
  return(celda.mods[[best]])
}


# Search through a celda_list's res.list model-by-model for one with the corresponding
# K/L.
search_res_list = function(celda_list, K=NULL, L=NULL) {
  requested.chain = NULL
  for (model in celda_list$res.list) {
    requested.chain = model
    if (K != model$K) next
    if (L != model$L) next
    break
  }
  if (is.null(requested.chain)) {
    stop("K/L parameter(s) requested did not appear for any model in the celda_list. Did you modify the run.params?")
  }
  return(requested.chain)
}
