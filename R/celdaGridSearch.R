#' available models
#' @export
available_models = c("celda_C", "celda_G", "celda_CG")


#' Run the celda Bayesian hierarchical model on a matrix of counts.
#' 
#' Yields assigments of genes/cells to clusters, depending on the provided model type.
#' 
#' @param counts A count matrix.
#' @param model Celda model. Options available in `celda::available.models`.
#' @param params.test List. A list denoting the combinations of parameters to run in a celda model. For example, "list(K=5:10, L=15:20)" will run all combinations of K from 5 to 10 and L from 15 to 20 in model 'celda_CG'.
#' @param params.fixed List. A list denoting additional parameters to use in each celda model. Default NULL.
#' @param max.iter Integer. Maximum number of iterations of Gibbs sampling to perform. Default 200. 
#' @param nchains Integer. Number of random cluster initializations. Default 1. 
#' @param cores Integer. The number of cores to use for parallel Gibbs sampling. Default 1.
#' @param best.only Logical. Whether to return only the chain with the highest log likelihood per combination of parameters or return all chains. Default TRUE. 
#' @param seed Integer. Passed to set.seed(). Default 12345.  
#' @param verbose Logical. Whether to print log messages during celda chain execution. Default TRUE. 
#' @param logfile.prefix Character. Prefix for log files from worker threads and main process. Default "Celda". 
#' @return Object of class "celda_list", which contains results for all model parameter combinations and summaries of the run parameters
#' @examples
#' ## Simulate a small dataset with 5 cell clusters and 10 feature modules
#' celda.sim = simulateCells(model="celda_CG", K=5, L=10)
#'
#' ## Run various combinations of parameters with 'celdaGridSearch'
#' celda.gs = celdaGridSearch(celda.sim$counts, model="celda_CG", 
#'                              params.test=list(K=4:6, L=9:11),
#'								params.fixed=list(sample.label=celda.sim$sample.label))                              
#' @import foreach
#' @export
celdaGridSearch = function(counts, model, params.test, params.fixed=NULL,
				    max.iter=200, nchains=3, cores=1,
				    best.only=TRUE, seed=12345, verbose=TRUE, 
                    logfile.prefix="Celda") {
 
  ## Check parameters
  model.params = as.list(formals(model))
  if(!all(names(params.test) %in% names(model.params))) {
    bad.params = setdiff(names(params.test), names(model.params))
    stop(paste0("The following elements in 'params.test' are not arguments of '", model, "': ", paste(bad.params, collapse=",")))
  }
  if(!is.null(params.fixed) && !all(names(params.fixed) %in% names(model.params))) {
    bad.params = setdiff(names(params.fixed), names(model.params))
    stop(paste0("The following elements in 'params.fixed' are not arguments of '", model, "': ", paste(bad.params, collapse=",")))
  }
  
  model.params.required = setdiff(names(model.params[model.params == ""]), "counts")
  if(!all(model.params.required %in% c(names(params.test), names(params.fixed)))) {
    missing.params = setdiff(model.params.required, c(names(params.test), names(params.fixed)))
    stop(paste0("The following arguments are not in 'params.test' or 'params.fixed' but are required for '", model, "': ", paste(missing.params, collapse=",")))
  }
  if(any(c("z.init", "y.init", "sample.label") %in% names(params.test))) {
    stop("Setting parameters such as 'z.init', 'y.init', and 'sample.label' in 'params.test' is not currently supported.")
  }
  
  # Set up parameter combinations for each individual chain
  run.params = expand.grid(c(chain=list(1:nchains), params.test))
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
    for(j in names(params.test)) {
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
    res = do.call(model, c(chain.params, params.fixed))
    return(list(res))
  }
  parallel::stopCluster(cl)  
  
  logliks = sapply(res.list, function(mod) { mod[["finalLogLik"]] })
  run.params = cbind(run.params, log_likelihood=logliks)
    
  celda.res = list(run.params=run.params, res.list=res.list, 
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


################################################################################
# Methods for manipulating celda_list objects                                  #
################################################################################
#' Select a subset of models from a 'celda_list' object generated by celdaGridSearch.
#' 
#' Convenience function for picking out specific models from a celda_list object.
#' Models can be selected by various parameters, most importantly the K/L parameters (number of cell
#'  clusters / number of feature clusters). 
#' 
#' @param celda.list Object of class "celda_list". An object containing celda models returned from `celdaGridSearch()`.
#' @return A new 'celda_list' object containing all models matching the provided criteria in 'params'. If entry in the list matches, the results for the matching model will be returned.
#' @examples
#' celda.sim = simulateCells(model="celda_CG", K=5, L=10)
#' celda.gs = celdaGridSearch(celda.sim$counts, model="celda_CG", params.test=list(K=4:6, L=9:11), params.fixed=list(sample.label=celda.sim$sample.label))
#' res.K5.L10 = subsetCeldaList(celda.gs, params=list(K=5, L=10))
#' @export
subsetCeldaList = function(celda.list, params) {
  if (!isTRUE(class(celda.list)[1] == "celda_list")) stop("celda.list parameter was not of class celda_list.")

  ## Check for bad parameter names
  if(!all(names(params) %in% colnames(celda.list$run.params))) {
	bad.params = setdiff(names(params), colnames(celda.list$run.params))
	stop(paste0("The following elements in 'params' are not columns in celda.list$run.params: ", paste(bad.params, collapse=",")))
  }
  
  ## Subset 'run.params' based on items in 'params'
  new.run.params = celda.list$run.params
  for(i in names(params)) {
	new.run.params = subset(new.run.params, new.run.params[,i] %in% params[[i]])
	
	if(nrow(new.run.params) == 0) {
	  stop("No runs matched the criteria given in 'params'. Check 'celda.list$run.params' for complete list of parameters used to generate 'celda.list'.")
	}
  }
  
  ## Get index of selected models, subset celda.list, and return
  ix = match(new.run.params$index, celda.list$run.params$index)
  if(length(ix) == 1) {
	return(celda.list$res.list[[ix]])
  } else {
	celda.list$run.params = as.data.frame(new.run.params)
	celda.list$res.list = celda.list$res.list[ix]
	return(celda.list)
  }
}


#' Select models with the best log likelihood from a 'celda_list' object gererated by celdaGridSearch.
#' 
#' This function returns the celda model from a celda_list object containing
#' the maxiumim final log-likelihood. K, L, or combination of K and L must be 
#' provided for celda_list objects containing celda_C, celda_G, and celda_CG 
#' models, respectively.
#' 
#' @param celda.list Object of class "celda_list". An object containing celda models returned from `celdaGridSearch()`.
#' @return The celda model object with the highest finalLogLik attribute, meeting any K/L criteria provided
#' @examples
#' ## Run various combinations of parameters with 'celdaGridSearch'
#' celda.gs = celdaGridSearch(celda.sim$counts, model="celda_CG", params.test=list(K=4:6, L=9:11), params.fixed=list(sample.label=celda.sim$sample.label), best.only=FALSE)
#'
#' ## Returns same result as running celdaGridSearch with "best.only = TRUE"
#' celda.gs.best = selectBestModel(celda.gs)  
#' @export
selectBestModel = function(celda.list) {
  if (!isTRUE(class(celda.list)[1] == "celda_list")) stop("celda.list parameter was not of class celda_list.")
 
  group = setdiff(colnames(celda.list$run.params), c("index", "chain", "log_likelihood"))
  new.run.params = celda.list$run.params %>% dplyr::group_by(.dots=group) %>% dplyr::slice(which.max(log_likelihood))
  
  ix = match(new.run.params$index, celda.list$run.params$index)
  if(nrow(new.run.params) == 1) {
    return(celda.list$res.list[[ix]])
  } else {
    celda.list$run.params = as.data.frame(new.run.params)
    celda.list$res.list = celda.list$res.list[ix]
    return(celda.list)
  }
}

    

#' Deprecation warning for old grid search function
#' 
#' @param ... Additional parameters.
#' @export
celda = function(...) {
  warning("Warning: The celda() wrapper function has been deprecated. Please see celdaGridSearch().")
}