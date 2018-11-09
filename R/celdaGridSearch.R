#' available models
#' @export
available_models = c("celda_C", "celda_G", "celda_CG")


<<<<<<< HEAD
#' Run the celda Bayesian hierarchical model on a matrix of counts.
#'
#' Yields assigments of genes/cells to clusters, depending on the provided model type.
#'
#' @param counts A count matrix
#' @param model Which celda sub-model to run. Options include "celda_C" (cell clustering), "celda_G" (gene clustering), "celda_CG" (gene and cell clustering)
#' @param sample.label A numeric vector, character vector, or factor indicating the originating sample for each cell (column) in the count matrix. By default, every cell will be assumed to be from an independent sample.
#' @param K Number of desired cell subpopulation clusters. Required for celda_C and celda_CG models.
#' @param L Number of desired gene clusters. Required for celda_G and celda_CG models.
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution (celda_C / celda_CG only)
#' @param beta Non-zero concentration parameter for gene Dirichlet distribution
#' @param delta The Dirichlet distribution parameter for Eta; adds a gene pseudocount to the numbers of genes each state (celda_G / celda_CG only)
#' @param gamma The Dirichlet distribution parameter for Psi; adds a pseudocount to each gene within each transcriptional state (celda_G / celda_CG only)
#' @param max.iter The maximum number of iterations
#' @param z.init Initial values of z. If NULL, z will be randomly sampled. Default NULL. (celda_C / celda_CG only)
#' @param y.init Initial values of y. If NULL, y will be randomly sampled. Default NULL. (celda_G / celda_CG only)
#' @param stop.iter Number of iterations without improvement in the log likelihood to stop the Gibbs sampler. Default 10.
#' @param split.on.iter On every 'split.on.iter' iteration, a heuristic will be applied to determine if a gene/cell cluster should be reassigned and another gene/cell cluster should be split into two clusters. Default 10.
#' @param nchains The number of chains of Gibbs sampling to run for every combination of K/L parameters. Defaults to 1 (celda_G / celda_CG only)
#' @param bestChainsOnly Return only the best chain (by final log-likelihood) per K/L combination.
#' @param cores The number of cores to use for parallell Gibbs sampling. Defaults to 1.
#' @param seed The base seed for random number generation
#' @param verbose Print log messages during celda chain execution
#' @param logfile_prefix Prefix for log files from worker threads and main process.
#' @return Object of class "celda_list", which contains results for all model parameter combinations and summaries of the run parameters
=======
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
>>>>>>> s4
#' @import foreach
#' @export
celdaGridSearch = function(counts, model, sample.label=NULL, K=NULL, L=NULL, alpha=1, beta=1, 
                 delta=1, gamma=1, max.iter=200, z.init=NULL, y.init=NULL,
                 stop.iter=10, split.on.iter=10, nchains=1,
                 bestChainsOnly=TRUE, cores=1, seed=12345, verbose=FALSE,
                 logfile_prefix="Celda") {

  validateArgs(counts, model, sample.label, nchains, cores, seed, K=K, L=L)
  params.list = buildParamList(counts, model, sample.label, alpha, beta, delta,
                               gamma, max.iter, z.init, y.init, stop.iter, split.on.iter,
                               nchains, cores, seed)

  # Redirect stderr from the worker threads if user asks for verbose
  if(!is.null(logfile_prefix)) {
    logfile = paste0(logfile_prefix, "_main_log.txt")
  } else {
    logfile = NULL
  }
  if (isTRUE(verbose)) logMessages(date(), "... Starting ", model, logfile=logfile, append=FALSE)
  params.list$logfile = logfile
  cl = if (verbose) parallel::makeCluster(cores, outfile=logfile) else parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  # Details for each model parameter / chain combination
  run.params = expand.grid(plyr::compact(list(chain=1:nchains, K=K, L=L)))
  run.params$index = as.numeric(rownames(run.params))

  # Pre-generate a set of random seeds to be used for each chain
  all.seeds = seed:(seed + nchains - 1)

  # An MD5 checksum of the count matrix. Passed to models so
  # later on, we can check on celda_* model objects which
  # count matrix was used.
  counts = processCounts(counts)
  count.checksum = digest::digest(counts, algo="md5")
<<<<<<< HEAD
  params.list$count.checksum = count.checksum
  params.list$nchains = 1
   
  res.list = foreach(i = 1:nrow(run.params), .export=model, .combine = c, .multicombine=TRUE) %dopar% {
    chain.params = append(params.list,
                          as.list(dplyr::select(run.params[i,],
                                                dplyr::matches("K|L"))))
=======

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
>>>>>>> s4
    chain.params$seed = all.seeds[ifelse(i %% nchains == 0, nchains, i %% nchains)]

    if (isTRUE(verbose)) {
      ## Generate a unique log file name based on given prefix and parameters
      chain.params$logfile = paste0(logfile_prefix, "_",
                                    paste(paste(colnames(run.params), run.params[i,], sep="-"), collapse="_"),  "_Seed-", chain.params$seed, "_log.txt")
      res = do.call(model, chain.params)
    } else {
      chain.params$logfile = NULL
      res = suppressMessages(do.call(model, chain.params))
    }
    return(list(res))
<<<<<<< HEAD
  }
  parallel::stopCluster(cl)
  celda.res = list(run.params=run.params, res.list=res.list,
                   content.type=model, count.checksum=count.checksum)
  class(celda.res) = "celda_list"

  if (isTRUE(bestChainsOnly)) {
    new.run.params = unique(dplyr::select(run.params, -index, -chain))
    new.run.params$index = 1:nrow(new.run.params)
    best.chains = apply(new.run.params, 1,
                        function(params) {
                          k = if ("K" %in% names(params)) params[["K"]] else NULL
                          l = if ("L" %in% names(params)) params[["L"]] else NULL
                          getBestModel(celda.res, k, l)
                        })
    celda.res$run.params = new.run.params
    celda.res$res.list = best.chains
=======
  })
  parallel::stopCluster(cl)  
  
  logliks = sapply(res.list, function(mod) { mod@finalLogLik })
  run.params = cbind(run.params, log_likelihood=logliks)
    
  celda.res = methods::new("celdaList", run.params=run.params, res.list=res.list, 
                           count.checksum=count.checksum)
  
  if (isTRUE(best.only)) {
    celda.res = selectBestModel(celda.res) 
>>>>>>> s4
  }

  if (isTRUE(verbose)) logMessages(date(), "... Completed ", model, logfile=logfile, append=TRUE)
  return(celda.res)
}


<<<<<<< HEAD
# Build a list of parameters tailored to the specific celda model being run,
# validating the provided parameters along the way
buildParamList = function(counts, model, sample.label, alpha, beta, delta,
                          gamma, max.iter, z.init, y.init, stop.iter, split.on.iter,
                          nchains, cores, seed) {

  params.list = list(counts=counts,
                     max.iter=max.iter,
                     stop.iter=stop.iter,
                     split.on.iter=split.on.iter)

  if (model %in% c("celda_C", "celda_CG")) {
    params.list$alpha = alpha
    params.list$beta = beta
    params.list$z.init=z.init
    params.list$sample.label=sample.label
  }
  if (model %in% c("celda_G", "celda_CG")) {
    params.list$beta = beta
    params.list$delta = delta
    params.list$gamma = gamma
    params.list$y.init = y.init
=======
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
>>>>>>> s4
  }

  return(params.list)
}

<<<<<<< HEAD

# Sanity check arguments to celda() to ensure a smooth run.
# See parameter descriptions from celda() documentation.
validateArgs = function(counts, model, sample.label,
                         nchains, cores, seed, K=NULL, L=NULL) { #, ...) {
  model_args = names(formals(model))
  if ("K" %in% model_args) {
    if (is.null(K)) {
      stop("Must provide a K parameter when running a celda_C or celda_CG model")
    } else if (is.numeric(K) && K <= 1) {
      stop("K parameter must be > 1")
    }

=======
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
>>>>>>> s4
  }
  if ("L" %in% model_args) {
    if (is.null(L)) {
      stop("Must provide a L parameter when running a celda_G or celda_CG model")
    } else if (is.numeric(L) && L <= 1) {
      stop("L parameter must be > 1")
    }
  }

  validateCounts(counts, K, L)

  if (!(model %in% available_models)) stop("Unavailable model specified")

  if (!is.null(sample.label)) {
    if (!(class(sample.label) %in% c("numeric", "character", "factor"))) {
      stop("Invalid sample.label; parameter should be either a numeric vector, character vector, or factor")
    }

    if (ncol(counts) != length(sample.label)) stop("Length of sample.label does not match number of columns (cells) in counts matrix")
  }

  if (!is.numeric(nchains) | length(nchains) > 1 | nchains == 0) stop("Invalid nchains specified")

  if (!is.numeric(cores) | length(cores) > 1 | cores == 0) stop("Invalid cores specified")
  if (!is.numeric(seed) | length(seed) > 1) stop("Invalid seed specified")
}


<<<<<<< HEAD
# Perform some simple checks on the counts matrix, to ensure celda won't choke.
# See parameter descriptions from celda() documentation.
validateCounts = function(counts, K, L) {
  # counts has to be a matrix...
  if (class(counts) != "matrix") stop("counts argument must be of class 'matrix'")

  # And each row/column of the count matrix must have at least one count
  count.row.sum = rowSums(counts)
  count.col.sum = colSums(counts)

  if (sum(count.row.sum == 0) >= 1 | sum(count.col.sum == 0) >= 1) {
    stop("Each row and column of the count matrix must have at least one count")
  }

  # Ensure that number of genes / cells is never more than
  # the number of requested clusters for each
  if (!is.null(L) && nrow(counts) < L) {
    stop("Number of genes (rows) in count matrix must be >= L")
  }
  if (!is.null(K) && ncol(counts) < K) {
    stop("Number of cells (columns) in count matrix must be >= K")
  }
=======
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
>>>>>>> s4
}
