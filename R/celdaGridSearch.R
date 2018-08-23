#' available models
#' @export
available_models = c("celda_C", "celda_G", "celda_CG")


#' Run the celda Bayesian hierarchical model on a matrix of counts.
#' 
#' Yields assigments of genes/cells to clusters, depending on the provided model type.
#' 
#' @param counts A count matrix.
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @param sample.label Vector or factor. Denotes the sample label for each cell (column) in the count matrix.
#' @param K.to.test Integer vector. List of K's to evaluate, where each K is the number of cell populations. 
#' @param L Integer vector. List of L's to evaluate, where each L is the number of cell populations. 
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount to each cell population in each sample. Default 1. 
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to each feature module in each cell population. Default 1. 
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to each feature in each module. Default 1. 
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to the number of features in each module. Default 1. 
#' @param max.iter Integer. Maximum number of iterations of Gibbs sampling to perform. Default 20. 
#' @param z.init Integer vector. Sets initial starting values of z. If NULL, starting values for each cell will be randomly sampled from 1:K. Default NULL.
#' @param y.init Integer vector. Sets initial starting values of y. If NULL, starting values for each feature will be randomly sampled from 1:L. Default NULL.
#' @param stop.iter Integer. Number of iterations without improvement in the log likelihood to stop inference. Default 20.
#' @param split.on.iter Integer. On every `split.on.iter` iteration, a heuristic will be applied to determine if a cell population or feature module should be reassigned and another cell population or feature module should be split into two clusters. To disable splitting, set to -1. Default 20.
#' @param nchains Integer. Number of random cluster initializations. Default 1. 
#' @param bestChainsOnly Logical. Whether to return only the best chain (by final log-likelihood) per K/L combination. Default TRUE. 
#' @param cores Integer. The number of cores to use for parallel Gibbs sampling. Default 1.
#' @param seed Integer. Passed to set.seed(). Default 12345.  
#' @param verbose Logical. Whether to print log messages during celda chain execution. Default FALSE. 
#' @param logfile.prefix Character. Prefix for log files from worker threads and main process. Default "Celda". 
#' @return Object of class "celda_list", which contains results for all model parameter combinations and summaries of the run parameters
#' @import foreach
#' @export
celdaGridSearch = function(counts, celda.mod, sample.label=NULL, K.to.test=NULL, L=NULL, alpha=1, beta=1, 
                 delta=1, gamma=1, max.iter=20, z.init=NULL, y.init=NULL,
                 stop.iter=10, split.on.iter=10, nchains=1, 
                 bestChainsOnly=TRUE, cores=1, seed=12345, verbose=FALSE, 
                 logfile.prefix="Celda") {
 
  validateArgs(counts, celda.mod, sample.label, nchains, cores, seed, K.to.test=K.to.test, L=L)
  params.list = buildParamList(counts, celda.mod, sample.label, alpha, beta, delta,
                               gamma, max.iter, z.init, y.init, stop.iter, split.on.iter,
                               nchains, cores, seed)
  
  # Redirect stderr from the worker threads if user asks for verbose
  if(!is.null(logfile.prefix)) {
    logfile = paste0(logfile.prefix, "_main_log.txt")
  } else {
    logfile = NULL
  }  
  if (isTRUE(verbose)) logMessages(date(), "... Starting ", celda.mod, logfile=logfile, append=FALSE)
  params.list$logfile = logfile
  cl = if (verbose) parallel::makeCluster(cores, outfile=logfile) else parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  
  # Details for each model parameter / chain combination 
  run.params = expand.grid(plyr::compact(list(chain=1:nchains, K.to.test=K.to.test, L=L)))
  run.params$index = as.numeric(rownames(run.params))
  
  # Pre-generate a set of random seeds to be used for each chain
  all.seeds = seed:(seed + nchains - 1)
  
  # An MD5 checksum of the count matrix. Passed to models so
  # later on, we can check on celda_* model objects which
  # count matrix was used.
  counts = processCounts(counts)
  count.checksum = digest::digest(counts, algo="md5")
  params.list$count.checksum = count.checksum
  params.list$nchains = 1
   
  res.list = foreach(i = 1:nrow(run.params), .export=celda.mod, .combine = c, .multicombine=TRUE) %dopar% {
    chain.params = append(params.list,
                          as.list(dplyr::select(run.params[i,],
                                                dplyr::matches("K.to.test|L"))))
    chain.params$seed = all.seeds[ifelse(i %% nchains == 0, nchains, i %% nchains)]
    
    if (isTRUE(verbose)) {
      ## Generate a unique log file name based on given prefix and parameters
      chain.params$logfile = paste0(logfile.prefix, "_",  
                                    paste(paste(colnames(run.params), run.params[i,], sep="-"), collapse="_"),  "_Seed-", chain.params$seed, "_log.txt")
      res = do.call(celda.mod, chain.params)
    } else {
      chain.params$logfile = NULL
      res = suppressMessages(do.call(celda.mod, chain.params))
    }
    return(list(res))
  }
  parallel::stopCluster(cl)
  celda.res = list(run.params=run.params, res.list=res.list, 
                   content.type=celda.mod, count.checksum=count.checksum)
  class(celda.res) = c("celda_list", celda.mod)
  
  if (isTRUE(bestChainsOnly)) {
    new.run.params = unique(dplyr::select(run.params, -index, -chain))
    new.run.params$index = 1:nrow(new.run.params)
    best.chains = apply(new.run.params, 1,
                        function(params) {
                          k = if ("K.to.test" %in% names(params)) params[["K.to.test"]] else NULL
                          l = if ("L" %in% names(params)) params[["L"]] else NULL
                          selectBestModel(celda.res, k, l)
                        })
    celda.res$run.params = new.run.params
    celda.res$res.list = best.chains
  }
  
  if (isTRUE(verbose)) logMessages(date(), "... Completed ", celda.mod, logfile=logfile, append=TRUE)
  return(celda.res)
}


# Build a list of parameters tailored to the specific celda model being run,
# validating the provided parameters along the way
buildParamList = function(counts, celda.mod, sample.label, alpha, beta, delta,
                          gamma, max.iter, z.init, y.init, stop.iter, split.on.iter,
                          nchains, cores, seed) {
  
  params.list = list(counts=counts,
                     max.iter=max.iter,
                     stop.iter=stop.iter,
                     split.on.iter=split.on.iter)
  
  if (celda.mod %in% c("celda_C", "celda_CG")) {
    params.list$alpha = alpha
    params.list$beta = beta
    params.list$z.init=z.init
    params.list$sample.label=sample.label
  } 
  if (celda.mod %in% c("celda_G", "celda_CG")) {
    params.list$beta = beta
    params.list$delta = delta
    params.list$gamma = gamma
    params.list$y.init = y.init
  }
  
  return(params.list)
}


# Sanity check arguments to celda() to ensure a smooth run.
# See parameter descriptions from celda() documentation.
validateArgs = function(counts, celda.mod, sample.label, 
                         nchains, cores, seed, K.to.test=NULL, L=NULL) { #, ...) {
  model_args = names(formals(celda.mod))
  if ("K.to.test" %in% model_args) {
    if (is.null(K.to.test)) { 
      stop("Must provide a K.to.test parameter when running a celda_C or celda_CG model")
    } else if (is.numeric(K.to.test) && K.to.test <= 1) {
      stop("K.to.test parameter must be > 1")
    }
    
  }
  if ("L" %in% model_args) {
    if (is.null(L)) {
      stop("Must provide a L parameter when running a celda_G or celda_CG model")
    } else if (is.numeric(L) && L <= 1) {
      stop("L parameter must be > 1")
    }
  }
  
  validateCounts(counts, K.to.test, L)
  
  if (!(celda.mod %in% available_models)) stop("Unavailable model specified")
      
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
    
    
# Perform some simple checks on the counts matrix, to ensure celda won't choke.
# See parameter descriptions from celda() documentation.
validateCounts = function(counts, K.to.test, L) {
  # counts has to be a matrix...
  if (class(counts) != "matrix") stop("counts argument must be of class 'matrix'")
  
  # And each row/column of the count matrix must have at least one count
  count.row.sum = rowSums(counts)
  count.col.sum = colSums(counts)
  
  if (sum(count.row.sum == 0) > 1 | sum(count.col.sum == 0) > 1) {
    stop("Each row and column of the count matrix must have at least one count")
  }
  
  # Ensure that number of genes / cells is never more than
  # the number of requested clusters for each
  if (!is.null(L) && nrow(counts) < L) {
    stop("Number of genes (rows) in count matrix must be >= L")
  }
  if (!is.null(K.to.test) && nrow(counts) < K.to.test) {
    stop("Number of cells (counts) in count matrix must be >= K")
  }
}


#' Deprecation warning for old grid search function
#' 
#' @param ... Additional parameters.
#' @export
celda = function(...) {
  warning("Warning: The celda() wrapper function has been deprecated. Please see celdaGridSearch().")
}