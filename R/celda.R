#' available models
#' @export
available_models = c("celda_C", "celda_G", "celda_CG")


#' Run the celda Bayesian hierarchical model on a matrix of counts.
#' 
#' Yields assigments of genes/cells to clusters, depending on the provided model type.
#' 
#' @param counts A count matrix
#' @param model Which celda sub-model to run. Options include "celda_C" (cell clustering), "celda_G" (gene clustering), "celda_CG" (gene and cell clustering)
#' @param sample.label A numeric vector indicating the originating sample for each cell (column) in the count matrix. By default, every cell will be assumed to be from an independent sample.
#' @param nchains The number of chains of Gibbs sampling to run for every combination of K/L parameters. Defaults to 1
#' @param cores The number of cores to use for parallell Gibbs sampling. Defaults to 1.
#' @param seed The base seed for random number generation
#' @param verbose Print log messages during celda chain execution
#' @param logfile_prefix Prefix for log files from worker threads and main process. 
#' @param ... Model specific parameters
#' @return Object of class "celda_list", which contains results for all model parameter combinations and summaries of the run parameters
#' @import foreach
#' @export
celda = function(counts, model, sample.label=NULL, nchains=1, cores=1, seed=12345, verbose=TRUE, logfile_prefix="Celda", ...) {
  message("Starting celda...")
  validate_args(counts, model, sample.label, nchains, cores, seed, ...)
  
  # Redirect stderr from the worker threads if user asks for verbose
  logfile = paste0(logfile_prefix, "_main_log.txt")
  cl = if (verbose) parallel::makeCluster(cores, outfile=logfile) else parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  
  # Details for each model parameter / chain combination 
  runs = expand.grid(chain=1:nchains, ...)
  runs$index = as.numeric(rownames(runs))
  if (verbose) print(runs)
  
  # Pre-generate a set of random seeds to be used for each chain
  all.seeds = seed:(seed + nchains - 1)
  
  # An MD5 checksum of the count matrix. Passed to models so
  # later on, we can check on celda_* model objects which
  # count matrix was used.
  count.checksum = digest::digest(counts, algo="md5")
    
  
  res.list = foreach(i = 1:nrow(runs), .export=model, .combine = c, .multicombine=TRUE) %dopar% {
    chain.seed = all.seeds[ifelse(i %% nchains == 0, nchains, i %% nchains)]
    
    if (verbose) {
      ## Generate a unique log file name based on given prefix and parameters
      logfile = paste0(logfile_prefix, "_", paste(paste(colnames(runs), runs[i,], sep="-"), collapse="_"), "_log.txt")
      res = do.call(model, c(list(counts=counts, sample.label=sample.label, count.checksum=count.checksum, seed=chain.seed, logfile=logfile), c(runs[i,-1])))
    } else {
      res = suppressMessages(do.call(model, c(list(counts=counts, sample.label=sample.label, count.checksum=count.checksum, seed=chain.seed, logfile=NULL), c(runs[i,-1]))))
    }
    return(list(res))
  }  
  parallel::stopCluster(cl)
  
  
  celda.res = list(run.params=runs, res.list=res.list, 
                   content.type=model, count.checksum=count.checksum)
  class(celda.res) = "celda_list"
  return(celda.res)
}


# Sanity check arguments to celda() to ensure a smooth run.
# See parameter descriptions from celda() documentation.
validate_args = function(counts, model, sample.label, 
                         nchains, cores, seed, K=NULL, L=NULL, ...) {
  if (model %in% c("celda_C", "celda_CG") && is.null(K)) {
    stop("Must provide a K parameter when running a celda_C or celda_CG model")
  }
  if (model %in% c("celda_G", "celda_CG") && is.null(L)) {
    stop("Must provide a L parameter when running a celda_G or celda_CG model")
  }
  
  validate_counts(counts, K, L)
  
  if (!(model %in% available_models)) stop("Unavailable model specified")
      
  if (!is.null(sample.label)) {
    if (!is.numeric(sample.label)) stop("Invalid sample.label; parameter should be a numeric vector")
    if (ncol(counts) != length(sample.label)) stop("Length of sample.label does not match number of columns (cells) in counts matrix") 
  }    
  
  if (!is.numeric(nchains) | length(nchains) > 1 | nchains == 0) stop("Invalid nchains specified")
  
  if (!is.numeric(cores) | length(cores) > 1 | cores == 0) stop("Invalid cores specified")
  if (!is.numeric(seed) | length(seed) > 1) stop("Invalid seed specified")
}
    
    
# Perform some simple checks on the counts matrix, to ensure celda won't choke.
# See parameter descriptions from celda() documentation.
validate_counts = function(counts, K, L) {
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
