

recursiveSplitCell = function(counts, s, K, K.subcluster, split.all=TRUE, alpha, beta, min.cell = 3, seed=12345) {
  current.K = 2
  overall.z = rep(1, ncol(counts))
  cluster.splits = matrix(NA, nrow=ncol(counts), ncol=K)
  cluster.split.flag = rep(TRUE, K)
  
  res = .celda_C(counts, K=K.subcluster, max.iter=20, split.on.iter=-1, split.on.last=FALSE, verbose=FALSE)
  overall.z = as.integer(as.factor(res@clusters$z))
  current.K = length(unique(overall.z)) + 1
  
  clust.split = vector("list", K)    
  while(current.K <= K) {
    z.ta = tabulate(overall.z, max(overall.z))
    z.to.split = which(z.ta > min.cell)
    
    for(i in z.to.split) {
      if(is.null(clust.split[[i]]) || is.na(clust.split[[i]])) {
        if(z.ta[i] <= K.subcluster) {
          clust.split[[i]] = seq.int(1, z.ta[i]) 
        } else {
          clustLabel = .celda_C(counts[,overall.z == i,drop=FALSE], K=K.subcluster, max.iter=20, split.on.iter=-1, split.on.last=FALSE, verbose=FALSE)
          clust.split[[i]] = as.integer(as.factor(clustLabel@clusters$z))
        }  
      } 
    }
    
    if(isTRUE(split.all)) {
      for(i in z.to.split) {
        split.ix = clust.split[[i]] > 1
        ix = overall.z == i
        new.z = overall.z[ix]
        
        new.z[split.ix] = current.K + clust.split[[i]][split.ix] - 2
        
        overall.z[ix] = new.z
        current.K = max(overall.z) + 1
        clust.split[[i]] = NA
      }   
    } else {
      ## Calculate likelihood for each cluster split
      temp.z = matrix(overall.z, nrow=ncol(counts), ncol=length(z.to.split))
      ll = vector(mode="numeric", length=ncol(temp.z))    
      for(i in 1:ncol(temp.z)) {
        split.ix = clust.split[[z.to.split[i]]] > 1
        ix = overall.z == z.to.split[i]
        new.z = temp.z[ix,i]
        new.z[split.ix] = current.K + clust.split[[z.to.split[i]]][split.ix] - 2
        temp.z[ix,i] = new.z
        temp.K = max(temp.z[,i])
        ll[i] = logLikelihood.celda_C(counts, "celda_C", s, temp.z[,i], temp.K, alpha, beta)	    	  
        
        ## Choose best split. Reset flags so the old cluster will be re-evaluated for splitting
        best.ix = which.max(ll)
        overall.z = temp.z[,best.ix]  
        clust.split[[z.to.split[best.ix]]] = NA
        
        current.K = max(overall.z, na.rm=TRUE) + 1  
      }
    }
  }
  
  return(overall.z)
}  

recursiveSplitCell_v2 = function(counts, sample.label=NULL, initial.K=5, max.K=25, alpha=1, beta=1, min.cell = 3, seed=12345, count.checksum=NULL) {
  counts = processCounts(counts)
  if(is.null(count.checksum)) {
    count.checksum = digest::digest(counts, algo="md5")
  }
  
  sample.label = processSampleLabels(sample.label, num.cells = ncol(counts))
  s = as.integer(sample.label)
  names = list(row=rownames(counts), column=colnames(counts), 
               sample=levels(sample.label))
  
  
  overall.z = rep(1, ncol(counts))
  cluster.splits = matrix(NA, nrow=ncol(counts), ncol=max.K)
  cluster.split.flag = rep(TRUE, max.K)
  
  model.initial = .celda_C(counts, sample.label=s, K=initial.K, alpha=alpha, beta=beta, verbose=TRUE)
  overall.z = as.integer(as.factor(model.initial@clusters$z))
  current.K = length(unique(overall.z)) + 1
  current.ll = model.initial@finalLogLik
  
  res.list = list(model.initial)
  clust.split = vector("list", max.K)    
  while(current.K <= max.K) {
    z.ta = tabulate(overall.z, current.K)
    z.to.split = which(z.ta > min.cell)
clust.split = vector("list", max.K)    
    for(i in z.to.split) {
      if(is.null(clust.split[[i]]) || is.na(clust.split[[i]])) {
        clustLabel = .celda_C(counts[,overall.z == i,drop=FALSE], K=2, split.on.iter=-1, split.on.last=FALSE, verbose=FALSE)
        if(length(unique(clustLabel@clusters$z)) == 2) {
          clust.split[[i]] = as.integer(as.factor(clustLabel@clusters$z))
        } else {
          z.to.split = setdiff(z.to.split, i) 
        }
      } 
    }
    
      ## Calculate likelihood for each cluster split
    temp.z = matrix(overall.z, nrow=ncol(counts), ncol=length(z.to.split))
    ll = vector(mode="numeric", length=ncol(temp.z))    
    for(i in 1:ncol(temp.z)) {
#      split.ix = clust.split[[z.to.split[i]]] > 1
      ix = overall.z == z.to.split[i]
#      new.z = temp.z[ix,i]
#      new.z[split.ix] = current.K + clust.split[[z.to.split[i]]][split.ix] - 2
      temp.z[ix,i] = ifelse(clust.split[[z.to.split[i]]] == 1, z.to.split[i], current.K)
#      temp.z[ix,i] = new.z
      temp.K = max(temp.z[,i])
      ll[i] = logLikelihood.celda_C(counts, "celda_C", s, temp.z[,i], temp.K, alpha, beta)	    	  
    }        
    ## Choose best split. Reset flags so the old cluster will be re-evaluated for splitting
    best.ix = which.max(ll)
    overall.z = temp.z[,best.ix]  
    clust.split[[z.to.split[best.ix]]] = NA
#print(table(overall.z))
#print(ll)
#print(best.ix)
    
    print(current.K)

  

    ## Peform reordering on final Z and Y assigments:
    ll = logLikelihood.celda_C(counts, "celda_C", s, overall.z, current.K, alpha, beta)

    temp.model = .celda_C(counts, sample.label=s, K=current.K, nchains=1, alpha=alpha, beta=beta, split.on.last=FALSE, verbose=TRUE, z.init=overall.z, seed=seed)
    overall.z = temp.model@clusters$z
    current.K = length(unique(overall.z))
    
    new.model = methods::new("celda_C", 
                             clusters=list(z=overall.z),
                             params=list(K=current.K,  alpha=alpha, beta=beta, 
                             seed=seed, count.checksum=count.checksum),
                             finalLogLik=ll,
                             sample.label=sample.label, 
                             names=names)
    res.list = c(res.list, list(new.model))
    current.K = max(overall.z, na.rm=TRUE) + 1  
  }
  
  ## Summarize paramters of different models
  logliks = sapply(res.list, function(mod) { mod@finalLogLik })
  runK = sapply(res.list, function(mod) { mod@params$K })
  run.params = data.frame(index=seq.int(1, length(res.list)), K=runK, log_likelihood=logliks, stringsAsFactors=FALSE)
  
  
  celda.res = methods::new("celdaList", run.params=run.params, res.list=res.list, 
                           count.checksum=count.checksum)
  
  return(celda.res)
}  

recursiveSplitCell_v3 = function(counts, sample.label=NULL, initial.K=5, max.K=25, alpha=1, beta=1, min.cell = 3, seed=12345) {
  counts = processCounts(counts)
  count.checksum = digest::digest(counts, algo="md5")
  
  sample.label = processSampleLabels(sample.label, num.cells = ncol(counts))
  s = as.integer(sample.label)
  names = list(row=rownames(counts), column=colnames(counts), 
               sample=levels(sample.label))
  
  
  model.initial = .celda_C(counts, sample.label=s, K=initial.K, alpha=alpha, beta=beta, verbose=TRUE)
  overall.z = as.integer(as.factor(model.initial@clusters$z))
  current.K = length(unique(overall.z)) + 1
  current.ll = model.initial@finalLogLik
  
  res.list = list(model.initial)
  clust.split = vector("list", max.K)  
  while(current.K <= max.K) {
    print(paste0("currentK: ", current.K))    
    z.ta = tabulate(overall.z, current.K)
    z.to.split = which(z.ta > min.cell)
    clust.split = vector("list", max.K)    
    best.z = overall.z
    best.ll = -Inf
    
    for(i in z.to.split) {
      clustLabel = .celda_C(counts[,overall.z == i,drop=FALSE], K=2, split.on.iter=-1, split.on.last=FALSE, verbose=FALSE)
      if(length(unique(clustLabel@clusters$z)) == 2) {
        clust.split[[i]] = as.integer(as.factor(clustLabel@clusters$z))

        ix = overall.z == i
        new.z = overall.z
        new.z[ix] = ifelse(clustLabel@clusters$z == 2, i, current.K)
        ll = logLikelihood.celda_C(counts, "celda_C", s, new.z, current.K, alpha, beta)
        
        if(ll > best.ll) {
          best.z = new.z
          best.ll = ll
        }
      } 
    }

    temp.model = .celda_C(counts, sample.label=s, K=current.K, nchains=1, alpha=alpha, beta=beta, split.on.iter=-1, split.on.last=FALSE, verbose=FALSE, z.init=best.z, seed=seed, reorder=FALSE)
    ## If the number of clusters is still "current.K", then keep the reordering, otherwise keep the previous configuration
    if(length(unique(temp.model@clusters$z)) == current.K) {
      overall.z = temp.model@clusters$z
    } else {
      overall.z = best.z
    }
    
    ## Peform reordering on final Z and Y assigments:
    ll = logLikelihood.celda_C(counts, "celda_C", s, overall.z, current.K, alpha, beta)
    
    new.model = methods::new("celda_C", 
                             clusters=list(z=overall.z),
                             params=list(K=current.K,  alpha=alpha, beta=beta, 
                                         seed=seed, count.checksum=count.checksum),
                             finalLogLik=ll,
                             sample.label=sample.label, 
                             names=names)
    res.list = c(res.list, list(new.model))
    current.K = current.K + 1  
  }
  
  ## Summarize paramters of different models
  logliks = sapply(res.list, function(mod) { mod@finalLogLik })
  runK = sapply(res.list, function(mod) { mod@params$K })
  run.params = data.frame(index=seq.int(1, length(res.list)), K=runK, log_likelihood=logliks, stringsAsFactors=FALSE)
  
  
  celda.res = methods::new("celdaList", run.params=run.params, res.list=res.list, 
                           count.checksum=count.checksum)
  
  return(celda.res)
}  





recursiveSplitModule = function(counts, initial.L=10, max.L=100, temp.K=100, K.breaks=5, z=NULL, beta=1, delta=1, gamma=1, min.feature=3, seed=12345, verbose=TRUE, logfile=NULL) {
  
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=FALSE, verbose=verbose)  
  logMessages("Starting celda_G with recursive module splitting.", logfile=logfile, append=TRUE, verbose=verbose)
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose)  
  
  counts = processCounts(counts)
  count.checksum = digest::digest(counts, algo="md5")   
  names = list(row=rownames(counts), column=colnames(counts))
  
  start.time = Sys.time()
  s = rep(1, ncol(counts))
  if(!is.null(z) | (!is.null(temp.K) && K.breaks > 0)) {
    if(is.null(z)) {
      z = initialize.splitZ(counts, K=temp.K, K.subcluster=K.breaks, min.cell = 3, seed=seed)  
    }      
    logMessages(date(), ".. Collapsed to", length(unique(z)), "temporary cell populations", append=TRUE, verbose=verbose, logfile=logfile)
    new.counts = colSumByGroup(counts, z, length(unique(z)))
  }
  
  logMessages(date(), ".. Initializing with", initial.L, "modules", append=TRUE, verbose=verbose, logfile=logfile)
  model.initial = .celda_G(new.counts, L=initial.L, max.iter=20, nchains=1, verbose=FALSE, seed=seed)

  ## Create decomposed counts for this initial model
  overall.y = model.initial@clusters$y
  current.L = length(unique(overall.y)) + 1  
  
  ## Decomposed counts for full count matrix
  p = cG.decomposeCounts(counts, overall.y, current.L)
  n.TS.by.C = p$n.TS.by.C
  n.by.TS = p$n.by.TS
  nG.by.TS = p$nG.by.TS  
  n.by.G = p$n.by.G
  nG = p$nG
  nM = p$nM

  ## Fix initial model to contain correct LogLik with full counts matrix
  model.initial@finalLogLik = cG.calcLL(n.TS.by.C=n.TS.by.C, n.by.TS=n.by.TS, n.by.G=n.by.G, nG.by.TS=nG.by.TS, nM=nM, nG=nG, L=current.L, beta=beta, delta=delta, gamma=gamma)
  model.initial@completeLogLik = model.initial@finalLogLik
  model.initial@params$count.checksum = count.checksum
  model.initial@names = names

  ## Perform splitting for y labels
  res.list = list(model.initial)
  while(current.L <= max.L) {
    
    y.ta = tabulate(overall.y, max.L)
    y.to.split = which(y.ta > min.feature)
    
    best.y = overall.y
    best.ll = -Inf
    previous.y = overall.y
    for(i in y.to.split) {
      clustLabel = .celda_G(new.counts[overall.y == i,,drop=FALSE], L=2, split.on.iter=-1, split.on.last=FALSE, nchains=1, verbose=FALSE, seed=seed)
      
      if(length(unique(clustLabel@clusters$y)) == 2) {
        ix = overall.y == i
        new.y = overall.y
        new.y[ix] = ifelse(clustLabel@clusters$y == 2, i, current.L)
        ll = logLikelihood.celda_G(new.counts, new.y, current.L, beta, delta, gamma)

        if(ll > best.ll) {
          best.y = new.y
          best.ll = ll
        }
      }    
    }  
    
    # Allow features to cluster further
    temp.model = .celda_G(new.counts, L=current.L, stop.iter=5, split.on.iter=-1, split.on.last=FALSE, nchains=1, verbose=FALSE, y.init=best.y, reorder=FALSE)
    overall.y = temp.model@clusters$y
  
    # Adjust decomposed count matrices
    p = cG.reDecomposeCounts(counts, overall.y, previous.y, n.TS.by.C, n.by.G, L = current.L)
    n.TS.by.C = p$n.TS.by.C
    n.by.TS = p$n.by.TS
    nG.by.TS = p$nG.by.TS 
    previous.y = overall.y

    ## Create the final model object with correct info on full counts matrix  
    temp.model@finalLogLik = cG.calcLL(n.TS.by.C=n.TS.by.C, n.by.TS=n.by.TS, n.by.G=n.by.G, nG.by.TS=nG.by.TS, nM=nM, nG=nG, L=current.L, beta=beta, delta=delta, gamma=gamma)
    temp.model@completeLogLik = temp.model@finalLogLik
    temp.model@params$count.checksum = count.checksum
    temp.model@names = names
    
    ## Add extra row/column for next round of L
    n.TS.by.C = rbind(n.TS.by.C, rep(0L, ncol(n.TS.by.C)))
    n.by.TS = c(n.by.TS, 0L)    
    nG.by.TS = c(nG.by.TS, 0L)
    
    ## Add new model to results list and increment L
    logMessages(date(), ".. Created module", current.L, "| logLik:", temp.model@finalLogLik, append=TRUE, verbose=verbose, logfile=NULL)
    res.list = c(res.list, list(temp.model))
    current.L = current.L + 1
  }
  
  ## Summarize paramters of different models
  logliks = sapply(res.list, function(mod) { mod@finalLogLik })
  runL = sapply(res.list, function(mod) { mod@params$L })
  run.params = data.frame(index=seq.int(1, length(res.list)), L=runL, log_likelihood=logliks, stringsAsFactors=FALSE)
  
  
  celda.res = methods::new("celdaList", run.params=run.params, res.list=res.list, 
                           count.checksum=count.checksum)
  
  end.time = Sys.time()
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose)  
  logMessages("Completed celda_G with recursive module splitting. Total time:", format(difftime(end.time, start.time)), logfile=logfile, append=TRUE, verbose=verbose)
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose)  
  
  return(celda.res)
} 

splitModule = function(counts, celda.mod, module, n=2) {
  ix = celda.mod@clusters$y == module
  new.L = max(celda.mod@clusters$y) + n - 1
  if(sum(ix) > 1) {
    temp.model = .celda_G(counts[ix,,drop=FALSE], L=n, split.on.iter=-1, split.on.last=FALSE, nchains=1, verbose=FALSE)  
    new.y = celda.mod@clusters$y
    new.y[ix] = ifelse(temp.model@clusters$y == 1, module, new.L)
    #model = .celda_G(counts, L=new.L, split.on.iter=-1, split.on.last=FALSE, nchains=1, verbose=TRUE, y.init=new.y, reorder=FALSE)
    new.ll = logLikelihood.celda_G(counts=counts, y=new.y, L=new.L, beta=celda.mod@params$beta, delta=celda.mod@params$delta, gamma=celda.mod@params$gamma)
    
    model = methods::new("celda_G", clusters=list(y=new.y),
                         params=list(L=new.L, beta=celda.mod@params$beta, delta=celda.mod@params$delta, gamma=celda.mod@params$gamma, count.checksum=celda.mod@params$count.checksum),
                         names=celda.mod@names, 
                         finalLogLik=new.ll)
    
  } else {
    message("No additional splitting was able to be performed.")
    model = celda.mod
  }
  
  return(model)
}

