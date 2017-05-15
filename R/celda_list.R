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
#' @param best Method for choosing best chain automatically. Options are c("perplexity", "harmonic", "loglik"). See documentation for chooseBestModel for details. Overrides chain parameter if provided.
#' @return A celda model object matching the provided parameters (of class "celda_C", "celda_G", "celda_CG" accordingly), or NA if one is not found.
#' @export
getModel = function(celda.list, K=NULL, L=NULL, chain=1, best=NULL) {
  validate_get_model_params(celda.list, K, L, chain, best)  # Sanity check params
  
  requested.chain = NA
  run.params = celda.list$run.params
  
  if (celda.list$content.type == "celda_CG") {
    if (!is.null(best)) {
      matching.chain.idx = run.params[run.params$K == K & run.params$L == L, "index"]
      requested.chain = chooseBestChain(celda.list[matching.chain.idx])
    } else {
      requested.chain.idx = run.params[run.params$K == K & run.params$L == L & run.params$chain == chain,
                                       "index"]
      requested.chain = celda.list$res.list[[requested.chain.idx]]
    }
  }
  
  
  if (celda.list$content.type == "celda_C") {
    if (!is.null(best)) {
      matching.chain.idx = run.params[run.params$K == K, "index"]
      requested.chain = chooseBestChain(celda.list[matching.chain.idx])
    } else {
      requested.chain.idx = run.params[run.params$K == K & run.params$chain == chain, "index"]
      requested.chain = celda.list$res.list[[requested.chain.idx]]
    }
  }
  
  
  if (celda.list$content.type == "celda_G") {
    if (!is.null(best)) {
      matching.chain.idx = run.params[run.params$L == L, "index"]
      requested.chain = chooseBestChain(celda.list[matching.chain.idx])
    } else { 
      requested.chain.idx = run.params[run.params$L == L & run.params$chain == chain, "index"]
      requested.chain = celda.list$res.list[[requested.chain.idx]]
    }
  }
  
  return (requested.chain)
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
}


#' Determine the best chain among a set of celda_* objects with
#' otherwise uniform K/L choices.
#' @param celda.mods A list of celda class objects (celda_C, celda_CG, celda_G)
#' @param method How to choose the best chain. Choices are c("perplexity", "harmonic", "loglik"). Defaults to perplexity. "perplexity" calculates each chain's perplexity as the inverse of the geometric mean, per the original LDA description. "harmonic" calculates each chain's marginal likelihood as the harmonic mean of each iteration of Gibbs sampling's log likelihoods. "loglik" chooses the chain which reached the maximal log likelihood during Gibbs sampling.
chooseBestChain = function(celda.mods, method="perplexity") {
  # We want to get the *most negative* perplexity, as opposed to the *least* negative
  # for the other metrics...
  if (method == "perplexity"){
    metrics = lapply(celda.mods, function(mod) { calculate_perplexity(mod$completeLogLik) })
    best = which(metrics == min(metrics))
    return(celda.mods[[best]])
  } 
  
  else if (method == "harmonic"){
    metrics = lapply(celda.mods, function(mod) { calculate_perplexity(mod$completeLogLik) })
  } 
  else if (method == "loglik"){
    metrics = lapply(celda.mods, function(mod) { max(mod$completeLogLik) })
  }  else {
    stop("Invalid method specified.")
  }
  best = which(metrics == max(metrics))
  return(celda.mods[[best]])
  
}



