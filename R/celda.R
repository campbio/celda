#' @import foreach
#' @export
celda = function(counts, model, sample.label=NULL, nchains=1, cores=1, seed=12345, ...) {
  cl = parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  
  runs = expand.grid(chain=1:nchains, ...)
  print(runs)  # TODO Toggle based on log level
  
  all.seeds = seed:(seed + nchains - 1)
  
  res.list = foreach(i = 1:nrow(runs), .export=model, .combine = c, .multicombine=TRUE) %dopar% {
    chain.seed = all.seeds[ifelse(i %% nchains == 0, nchains, i %% nchains)]
    res = do.call(model, c(list(counts=counts, sample.label=sample.label, seed=chain.seed), c(runs[i,-1])))
    all = list(runs[i,], res)
    return(list(res))
  }  
  
  celda.res = list(run.params=runs, res.list=res.list)
  class(celda.res) = model
  return(celda.res)
}
