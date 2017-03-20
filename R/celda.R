#' @export
celda = function(counts, model, sample.label=NULL, nchains=1, cores=1, ...) {
  suppressPackageStartupMessages(require(foreach))
  suppressPackageStartupMessages(require(doParallel))
  
  cl<-makeCluster(cores)
  registerDoParallel(cl)
  
  runs =expand.grid(chain=1:nchains, ...)
  print(runs)  
  res.list = foreach(i = 1:nrow(runs), .export=model, .combine = c, .multicombine=TRUE) %dopar% {
    source("celda_CG.R")  ## quick hack to bring in auxillary functions, will need to eventually hanlde via ".export", not sure about how to handle in a package
    source("celda_functions.R")
    res = do.call(model, c(list(counts=counts, sample.label=sample.label), c(runs[i,-1])))
    all = list(runs[i,], res)
    return(list(res))
  }  
  
  return(res.list)
}
