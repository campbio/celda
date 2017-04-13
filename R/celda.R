#' @import foreach
#' @export
celda = function(counts, model, sample.label=NULL, nchains=1, cores=1, ...) {
  cl = parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  
  runs = expand.grid(chain=1:nchains, ...)
  print(runs)  
  
  
  # Generate unique random seeds for each chain. 
  # Reusable between parameter combinations (e.g. unique combinations of K and L)
  seeds = sample(1:1000000, size=nchains)
  
  res.list = foreach(i = 1:nrow(runs), .export=model, .combine = c, .multicombine=TRUE) %dopar% {
    seed = seeds[ifelse(i %% nchains == 0, nchains, i %% nchains)]
    res = do.call(model, c(list(counts=counts, sample.label=sample.label, seed=seed), c(runs[i,-1])))
    all = list(runs[i,], res)
    return(list(res))
  }  
  
  return(list(run.params=runs, res.list=res.list))
}
