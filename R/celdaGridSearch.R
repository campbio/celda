#' available models
#' @export
available_models = c("celda_C", "celda_G", "celda_CG")


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
  params.list$count.checksum = count.checksum
  params.list$nchains = 1
   
  res.list = foreach(i = 1:nrow(run.params), .export=model, .combine = c, .multicombine=TRUE) %dopar% {
    chain.params = append(params.list,
                          as.list(dplyr::select(run.params[i,],
                                                dplyr::matches("K|L"))))
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
  }
  parallel::stopCluster(cl)
  celda.res = list(run.params=run.params, res.list=res.list, 
                   content.type=model, count.checksum=count.checksum)
  class(celda.res) = c("celda_list", model)
  
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
  }
  
  if (isTRUE(verbose)) logMessages(date(), "... Completed ", model, logfile=logfile, append=TRUE)
  return(celda.res)
}


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
  }
  
  return(params.list)
}


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
    
    
# Perform some simple checks on the counts matrix, to ensure celda won't choke.
# See parameter descriptions from celda() documentation.
validateCounts = function(counts, K, L) {
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
  if (!is.null(K) && nrow(counts) < K) {
    stop("Number of cells (counts) in count matrix must be >= K")
  }
}
