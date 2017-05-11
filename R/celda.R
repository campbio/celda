available_models = c("celda_C", "celda_G", "celda_CG")


#' Run the celda Bayesian hierarchical model on a matrix of counts.
#' 
#' @param counts A count matrix
#' @param model Which celda sub-model to run. Options include "celda_C" (cell clustering), "celda_G" (gene clustering), "celda_CG" (gene and cell clustering)
#' @param sample.label A numeric vector indicating the sample for each cell (column) in the count matrix
#' @param nchains The number of chains of Gibbs sampling to run for every combination of K/L parameters
#' @param K An integer or range of integers indicating the desired number of cell clusters (for celda_C / celda_CG models)
#' @param L An integer or range of integers indicating the desired number of gene clusters (for celda_G / celda_CG models)
#' @param cores The number of cores to use to speed up Gibbs sampling
#' @param seed The base seed for random number generation. Each chain celda runs with have a seed index off of this one.
#' @param verbose Print messages during celda chain execution
#' @import foreach
#' @export
celda = function(counts, model, sample.label=NULL, nchains=1, cores=1, seed=12345, verbose=F, ...) {
  message("Starting celda...")
  validate_args(counts, model, sample.label, nchains, cores, seed)
  
  # Redirect stderr from the worker threads if user asks for verbose
  cl = if (verbose) parallel::makeCluster(cores, outfile="") else parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  
  runs = expand.grid(chain=1:nchains, ...)
  if (verbose) print(runs)
  
  all.seeds = seed:(seed + nchains - 1)
  
  res.list = foreach(i = 1:nrow(runs), .export=model, .combine = c, .multicombine=TRUE) %dopar% {
    chain.seed = all.seeds[ifelse(i %% nchains == 0, nchains, i %% nchains)]
    
    if (verbose) {
      res = do.call(model, c(list(counts=counts, sample.label=sample.label, seed=chain.seed, thread=i), c(runs[i,-1])))
    } else {
      res = suppressMessages(do.call(model, c(list(counts=counts, sample.label=sample.label, seed=chain.seed, thread=i), c(runs[i,-1]))))
    }
    return(list(res))
  }  
  parallel::stopCluster(cl)
  
  celda.res = list(run.params=runs, res.list=res.list)
  class(celda.res) = model
  return(celda.res)
}


#' Sanity check arguments to celda() to ensure a smooth run.
validate_args = function(counts, model, sample.label, nchains, cores, seed) {
  validate_counts(counts)
  
  if (!(model %in% available_models)) stop("Unavailable model specified")
      
  if (!is.null(sample.label)) {
    if (!is.numeric(sample.label)) stop("Invalid sample.label; parameter should be a numeric vector")
    if (ncol(counts) != length(sample.label)) stop("Length of sample.label does not match number of columns (cells) in counts matrix") 
  }    
  
  if (!is.numeric(nchains) | length(nchains) > 1 | nchains == 0) stop("Invalid nchains specified")
  
  if (!is.numeric(cores) | length(cores) > 1 | cores == 0) stop("Invalid cores specified")
  if (!is.numeric(seed) | length(seed) > 1) stop("Invalid seed specified")
}
    
    
#' Perform some simple checks on the counts matrix, to ensure celda won't choke.
validate_counts = function(counts) {
  # counts has to be a matrix...
  if (class(counts) != "matrix") stop("counts argument must be of class 'matrix'")
  
  # And each row/column of the count matrix must have at least one count
  count.row.sum = rowSums(counts)
  count.col.sum = colSums(counts)
  
  if (sum(count.row.sum == 0) > 1 | sum(count.col.sum == 0) > 1) {
    stop("Each row and column of the count matrix must have at least one count")
  }
}




