#' available models
#' @export
available_models = c("celda_C", "celda_G", "celda_CG")


#' Run the celda Bayesian hierarchical model on a matrix of counts.
#' 
#' Yields assigments of genes/cells to clusters, depending on the provided model type.
#' 
#' @param counts A count matrix.
#' @param model Celda model. Options available in `celda::available.models`.
#' @param max.iter Integer. Maximum number of iterations of Gibbs sampling to perform. Default 200. 
#' @param nchains Integer. Number of random cluster initializations. Default 1. 
#' @param best.only Logical. Whether to return only the chain with the highest log likelihood per combination of parameters. Default TRUE. 
#' @param cores Integer. The number of cores to use for parallel Gibbs sampling. Default 1.
#' @param seed Integer. Passed to set.seed(). Default 12345.  
#' @param verbose Logical. Whether to print log messages during celda chain execution. Default TRUE. 
#' @param logfile.prefix Character. Prefix for log files from worker threads and main process. Default "Celda". 
#' @return Object of class "celda_list", which contains results for all model parameter combinations and summaries of the run parameters
#' @examples
#' celda.sim = simulateCells(model="celda_CG")
#' celda.mods = celdaGridSearch(celda.sim$counts, model="celda_CG", 
#'                              sample.label=celda.sim$sample.label,
#'                              K.to.test=2:4, L.to.test=9:11, max.iter=2, nchains=1)
#' @import foreach
#' @export
celdaGridSearch = function(counts, model, params,
				    max.iter=200, cores=1, best.only=TRUE,
				    seed=12345, verbose=TRUE, nchains=3,
                    logfile.prefix="Celda") {
 
  ## Check parameters
  model.params = as.list(formals(model))
  if(!all(names(params) %in% names(model.params))) {
    bad.params = setdiff(names(params), names(model.params))
    stop(paste0("The following elements in 'params' are not arguments of '", model, "': ", paste(bad.params, collapse=",")))
  }
  
  model.params.required = setdiff(names(model.params[model.params == ""]), "counts")
  if(!all(model.params.required %in% names(params))) {
    missing.params = setdiff(model.params.required, names(params))
    stop(paste0("The following arguments are not in 'params' but are required for '", model, "': ", paste(missing.params, collapse=",")))
  }
  if(any(c("z.init", "y.init") %in% names(params))) {
    stop("Setting initialization parameters such as 'z.init' and 'y.init' in 'params' is not currently supported.")
  }
  
  # Set up parameter combinations for each individual chain
  run.params = expand.grid(c(chain=list(1:nchains), params))
  run.params = cbind(index = 1:nrow(run.params), run.params)

  # Pre-generate a set of random seeds to be used for each chain
  all.seeds = seed:(seed + nchains - 1)

  logMessages("--------------------------------------------------------------------", logfile=NULL, append=FALSE, verbose=verbose)  
  logMessages("Starting celdaGridSearch with", model, logfile=NULL, append=TRUE, verbose=verbose)
  logMessages("Number of cores:", cores, logfile=NULL, append=TRUE, verbose=verbose)  
  logMessages("--------------------------------------------------------------------", logfile=NULL, append=TRUE, verbose=verbose)  
  start.time = Sys.time()

  # An MD5 checksum of the count matrix. Passed to models so
  # later on, we can check on celda_* model objects which
  # count matrix was used.
  counts = processCounts(counts)
  count.checksum = digest::digest(counts, algo="md5")

  ## Use DoParallel to loop through each combination of parameters
  cl = parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)   
  res.list = foreach(i = 1:nrow(run.params), .export=model, .combine = c, .multicombine=TRUE) %dopar% {
    
    ## Set up chain parameter list
    current.run = c(run.params[i,])
    chain.params = list()
    for(j in names(params)) {
      chain.params[[j]] = current.run[[j]]
    }
    chain.params$counts = counts
    chain.params$seed = all.seeds[ifelse(i %% nchains == 0, nchains, i %% nchains)]
    chain.params$max.iter = max.iter
    chain.params$nchain = 1
    chain.params$count.checksum = count.checksum
    chain.params$verbose = verbose
    chain.params$logfile = paste0(logfile.prefix, "_", paste(paste(colnames(run.params), run.params[i,], sep="-"), collapse="_"),  "_Seed-", chain.params$seed, "_log.txt")                                
    
    ## Run model
    res = do.call(model, chain.params)
    return(list(res))
  }
  parallel::stopCluster(cl)  
  
  logliks = sapply(res.list, function(mod) { mod[["finalLogLik"]] })
  run.params = cbind(run.params, log_likelihood=logliks)
    
  celda.res = list(run.params=run.params, params=params, res.list=res.list, 
                   content.type=model, count.checksum=count.checksum)
  class(celda.res) = c("celda_list", model)
  
  if (isTRUE(best.only)) {
    celda.res = selectBestModel(celda.res) 
  }
  
  end.time = Sys.time()
  logMessages("--------------------------------------------------------------------", logfile=NULL, append=TRUE, verbose=verbose)  
  logMessages("Completed celdaGridSearch. Total time:", format(difftime(end.time, start.time)), logfile=NULL, append=TRUE, verbose=verbose)
  logMessages("--------------------------------------------------------------------", logfile=NULL, append=TRUE, verbose=verbose)  

  return(celda.res)
}

    

#' Deprecation warning for old grid search function
#' 
#' @param ... Additional parameters.
#' @export
celda = function(...) {
  warning("Warning: The celda() wrapper function has been deprecated. Please see celdaGridSearch().")
}