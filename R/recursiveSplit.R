
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



recursiveSplitModule = function(counts, initial.L=10, max.L=100, temp.cell.cluster=1000, cell.cluster.breaks=10, beta=1, delta=1, gamma=1, min.feature=3, seed=12345) {
  
  s = rep(1, ncol(counts))
  if(!is.null(temp.cell.cluster) && temp.cell.cluster > 0) {
    new.z = recursive.splitZ(counts, s, K=temp.cell.cluster, K.subcluster=cell.cluster.breaks, split.all=TRUE, alpha=1, beta=1, min.cell = 3, seed=seed)
    counts = colSumByGroup(counts, new.z, length(unique(new.z)))
  }
  count.checksum = digest::digest(counts, algo="md5")    
  #  model.initial = .celda_G(counts, L=initial.L, split.on.iter=-1, split.on.last=FALSE, nchains=1, verbose=T)
  model.initial = .celda_G(counts, L=initial.L, nchains=3, verbose=T)
  
  ## Perform splitting for y labels
  current.L = initial.L + 1
  overall.y = model.initial@clusters$y
  current.ll = model.initial@finalLogLik
  
  y.history = matrix(NA, nrow=nrow(counts), ncol=L)
  y.history[,2] = overall.y
  
  ll.history = c()
  res.list = list(model.initial)
  
  while(current.L <= max.L) {
    print(paste0("currentL: ", current.L))
    y.ta = tabulate(overall.y, max.L)
    y.to.split = which(y.ta > min.feature)
    
    best.y = overall.y
    best.ll = -Inf
    for(i in y.to.split) {
      clustLabel = .celda_G(counts[overall.y == i,,drop=FALSE], L=2, split.on.iter=-1, split.on.last=FALSE, nchains=1, verbose=FALSE)
      #      clust.split[[i]] = as.integer(as.factor())
      if(length(unique(clustLabel@clusters$y)) == 2) {
        ix = overall.y == i
        temp.y = overall.y
        temp.y[ix] = ifelse(clustLabel@clusters$y == 2, i, current.L)
        ll = logLikelihood.celda_G(counts=counts, y=temp.y, L=current.L, beta=beta, delta=delta, gamma=gamma)
        if(ll > best.ll) {
          best.y = temp.y
          current.ll = ll
        }
      }    
    }  
    
    
    
    ## Calculate likelihood for each cluster split
    #    temp.y = matrix(overall.y, nrow=nrow(counts), ncol=length(y.to.split))
    #    ll = vector(mode="numeric", length=ncol(temp.y))    
    #    for(i in 1:ncol(temp.y)) {
    #      ix = overall.y == y.to.split[i]
    #      temp.y[ix,i] = ifelse(clust.split[[y.to.split[i]]] == 2, y.to.split[i], current.L)
    #      ll[i] = logLikelihood.celda_G(counts=counts, y=temp.y[,i], L=current.L, beta=beta, delta=delta, gamma=gamma)
    #    }  
    
    ## Choose best split. Reset flags so the old cluster will be re-evaluated for splitting
    #best.ix = which.max(ll)
    #temp.y = temp.y[,best.ix] 
    #temp.model = .celda_G(counts, L=current.L, stop.iter=5, split.on.iter=-1, split.on.last=FALSE, nchains=1, verbose=FALSE, y.init=temp.y[,best.ix])
    temp.model = .celda_G(counts, L=current.L, stop.iter=5, split.on.iter=-1, split.on.last=FALSE, nchains=1, verbose=T, y.init=best.y)
    #temp.model@count.checksum = cs
    overall.y = temp.model@clusters$y
    current.ll = temp.model@finalLogLik
    
    ## Get list of modules whose features have changed membership.
    ##The best split for these clusters need to b
    #    prev.y = y.history[,current.L-1]
    #    module.change = unique(prev.y[overall.y != prev.y])
    #    print(paste0("module change: ", module.change))
    #    for(i in module.change) {
    #      clust.split[[i]] = NA 
    #    }
    res.list = c(res.list, list(temp.model))
    #    ll.history = c(ll.history, temp.model@finalLogLik)
    #    clust.split = vector("list", L)
    #    y.history[,current.L] = overall.y
    current.L = current.L + 1
  }
  
  run.params = data.frame(index=seq.int(1, max.L-initial.L + 1), L=seq.int(initial.L, max.L))
  logliks = sapply(res.list, function(mod) { mod@finalLogLik })
  run.params = cbind(run.params, log_likelihood=logliks)
  
  
  celda.res = methods::new("celdaList", run.params=run.params, res.list=res.list, 
                           count.checksum=count.checksum)
  
  print("calculating perplexity")
  celda.res = resamplePerplexity(counts, celda.res, seed=12345)
  
  #return(list(y=y.history, ll=ll.history))
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

