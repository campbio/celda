#' @import foreach
#' @export
celda = function(counts, model, sample.label=NULL, nchains=1, cores=1, ...) {
  cl = parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  
  runs = expand.grid(chain=1:nchains, ...)
  print(runs)  
  
  res.list = foreach(i = 1:nrow(runs), .export=model, .combine = c, .multicombine=TRUE) %dopar% {
    res = do.call(model, c(list(counts=counts, sample.label=sample.label), c(runs[i,-1])))
    all = list(runs[i,], res)
    return(list(res))
  }  
  
  return(list(run.params=runs, res.list=res.list))
}
