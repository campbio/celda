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
#' @return A celda model object matching the provided parameters (of class "celda_C", "celda_G", "celda_CG" accordingly), or NA if one is not found.
#' @export
getModel = function(celda.list, K=NULL, L=NULL, chain=1) {
  if (class(celda.list) != "celda_list") stop("First argument to getModel() should be an object of class 'celda_list'")
  
  requested.chain = NA
  run.params = celda.list$run.params
  
  if (celda.list$list.contents == "celda_CG") {
    if (is.null(K) | is.null(L)) stop("Both K and L parameters needed for subsetting celda_CG result lists")
    requested.chain.idx = run.params[run.params$K == K & run.params$L == L & run.params$chain == chain,
                                     "index"]
    requested.chain = celda.list$res.list[[requested.chain.idx]]
  }
  
  
  if (celda.list$list.contents == "celda_C") {
    if (is.null(K)) stop("K parameter needed when subsetting celda_C result lists")
    requested.chain.idx = run.params[run.params$K == K & run.params$chain == chain, "index"]
    requested.chain = celda.list$res.list[[requested.chain.idx]]
  }
  
  
  if (celda.list$list.contents == "celda_G") {
    if (is.null(L)) stop("L parameter needed when subsetting celda_C result lists")
    requested.chain.idx = run.params[run.params$L == L & run.params$chain == chain, "index"]
    requested.chain = celda.list$res.list[[requested.chain.idx]]
  }
  
  return (requested.chain)
}
