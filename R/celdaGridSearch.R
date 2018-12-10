#' available models
#' @export
available_models = c("celda_C", "celda_G", "celda_CG")


#' @title Run Celda in parallel with multiple parameters
#' @description Run Celda with different combinations of parameters and multiple chains in parallel. The variable `available_models` contains the potential models that can be utilized. Different parameters to be tested should be stored in a list and passed to the argument `params.test`. Fixed parameters to be used in all models, such as `sample.label`, can be passed as a list to the argument `params.fixed`. When `verbose=TRUE`, output from each chain will be sent to a log file but not be displayed in stdout.
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. 
#' @param model Celda model. Options available in `celda::available.models`.
#' @param params.test List. A list denoting the combinations of parameters to run in a celda model. For example, `list(K=5:10, L=15:20)` will run all combinations of K from 5 to 10 and L from 15 to 20 in model `celda_CG()`.
#' @param params.fixed List. A list denoting additional parameters to use in each celda model. Default NULL.
#' @param max.iter Integer. Maximum number of iterations of sampling to perform. Default 200. 
#' @param nchains Integer. Number of random cluster initializations. Default 3. 
#' @param cores Integer. The number of cores to use for parallel estimation of chains. Default 1.
#' @param best.only Logical. Whether to return only the chain with the highest log likelihood per combination of parameters or return all chains. Default TRUE. 
#' @param seed Integer. Passed to `set.seed()`. Default 12345.  
#' @param verbose Logical. Whether to print log messages during celda chain execution. Default TRUE. 
#' @param logfile.prefix Character. Prefix for log files from worker threads and main process. Default "Celda". 
#' @return Object of class `celda_list`, which contains results for all model parameter combinations and summaries of the run parameters
#' @seealso `celda_G()` for feature clustering, `celda_C()` for clustering of cells, and `celda_CG()` for simultaneous clustering of features and cells. `subsetCeldaList()` can subset the `celda_list` object. `selectBestModel()` can get the best model for each combination of parameters. 
#' @examples
#' ## Run various combinations of parameters with 'celdaGridSearch'
#' cgs = celdaGridSearch(celda.CG.sim$counts, model="celda_CG", 
#'                       params.test=list(K=4:6, L=9:11), 
#'                       params.fixed=list(sample.label=celda.CG.sim$sample.label),
#'                       best.only=TRUE, nchains=1)
#' @import foreach
#' @export
celdaGridSearch = function(counts, model, params.test, params.fixed=NULL,
                           max.iter=200, nchains=3, cores=1,
                           best.only=TRUE, seed=12345, verbose=TRUE, 
                           logfile.prefix="Celda") {
  
  ## Check parameters
  validateCounts(counts)
  
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
  if(any(c("nchains") %in% names(params.test))){
    warning("Parameter 'nchains' should not be used within the params.test list")
    params.test[["nchains"]] <- NULL
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
  i = NULL  # Setting visible binding for R CMD CHECK
  #res.list = foreach(i = 1:nrow(run.params), .export=model, .combine = c, .multicombine=TRUE) %dopar% {
  #res.list = foreach(i = 1:nrow(run.params), .export=model, .combine = c, .multicombine=TRUE) %dopar% {
  res.list = sapply(1:nrow(run.params), function(i) {
    
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
  })
  parallel::stopCluster(cl)  
  
  logliks = sapply(res.list, function(mod) { mod@finalLogLik })
  run.params = cbind(run.params, log_likelihood=logliks)
  
  celda.res = methods::new("celdaList", run.params=run.params, res.list=res.list, 
                           count.checksum=count.checksum)
  
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
#' @title Subset celda_list object from celdaGridSearch
#' @description Select a subset of models from a `celda_list` object generated by `celdaGridSearch()` that match the criteria in the argument `params`.
#' 
#' @param celda.list celda.list Object of class `celda_list`. An object containing celda models returned from `celdaGridSearch`.
#' @param params List. List of parameters used to subset celda.list.
#' @return A new `celda_list` object containing all models matching the provided criteria in `params`. If only one item in the `celda_list` matches the given criteria, the matching model will be returned directly instead of a `celda_list` object.
#' @seealso `celdaGridSearch()` can run Celda with multiple parameters and chains in parallel. `selectBestModel()` can get the best model for each combination of parameters. 
#' @examples
#' res.K5.L10 = subsetCeldaList(celda.CG.grid.search.res, params=list(K=5, L=10))
#' @export
subsetCeldaList = function(celda.list, params) {
  if (!methods::is(celda.list, "celdaList")) stop("celda.list parameter was not of class celdaList.")
  
  ## Check for bad parameter names
  if(!all(names(params) %in% colnames(celda.list@run.params))) {
    bad.params = setdiff(names(params), colnames(celda.list@run.params))
    stop(paste0("The following elements in 'params' are not columns in run.params(celdaList)", paste(bad.params, collapse=",")))
  }
  
  ## Subset 'run.params' based on items in 'params'
  new.run.params = celda.list@run.params
  for(i in names(params)) {
    new.run.params = subset(new.run.params, new.run.params[,i] %in% params[[i]])
    
    if(nrow(new.run.params) == 0) {
      stop("No runs matched the criteria given in 'params'. Check 'run.params(celda.list)' for complete list of parameters used to generate 'celda.list'.")
    }
  }
  
  ## Get index of selected models, subset celda.list, and return
  ix = match(new.run.params$index, celda.list@run.params$index)
  if(length(ix) == 1) {
    return(celda.list@res.list[[ix]])
  } else {
    celda.list@run.params = as.data.frame(new.run.params)
    celda.list@res.list = celda.list@res.list[ix]
    return(celda.list)
  }
}

#' @title Select best chain within each combination of parameters
#' @description Select the chain with the best log likelihood for each combination of tested parameters from a `celda_list` object gererated by `celdaGridSearch()`.
#' 
#' @param celda.list Object of class `celda_list`. An object containing celda models returned from `celdaGridSearch()`.
#' @return A new `celda_list` object containing one model with the best log likelihood for each set of parameters. If only one set of parameters is in the `celda_list`, the best model will be returned directly instead of a `celda_list` object.
#' @seealso `celdaGridSearch()` can run Celda with multiple parameters and chains in parallel. `subsetCeldaList()` can subset the `celda_list` object.
#' @examples
#' ## Returns same result as running celdaGridSearch with "best.only = TRUE"
#' cgs.best = selectBestModel(celda.CG.grid.search.res)
#' @import data.table
#' @export
selectBestModel = function(celda.list) {
  if (!methods::is(celda.list, "celdaList")) stop("celda.list parameter was not of class celdaList.")
  
  log_likelihood = NULL
  group = setdiff(colnames(celda.list@run.params), c("index", "chain", "log_likelihood"))
  dt = data.table::as.data.table(celda.list@run.params)
  new.run.params = as.data.frame(dt[,.SD[which.max(log_likelihood)], by=group])
  new.run.params = new.run.params[,colnames(celda.list@run.params)]
  
  ix = match(new.run.params$index, celda.list@run.params$index)
  if(nrow(new.run.params) == 1) {
    return(celda.list@res.list[[ix]])
  } else {
    celda.list@run.params = as.data.frame(new.run.params)
    celda.list@res.list = celda.list@res.list[ix]
    return(celda.list)
  }
}



#' Celda models
#'
#' List of available Celda models with correpsonding descriptions.
#' 
#' @export
#' @return None
celda = function() {
  cat("celda_C: Clusters the columns of a count matrix containing single-cell data into K subpopulations.\n") 
  cat("celda_G: Clusters the rows of a count matrix containing single-cell data into L modules.\n")
  cat("celda_CG: Clusters the rows and columns of a count matrix containing single-cell data into L modules and K subpopulations, respectively.\n")
  cat("celdaGridSearch: Run Celda with different combinations of parameters and multiple chains in parallel.\n")  
}
