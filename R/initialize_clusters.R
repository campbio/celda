
initialize.cluster = function(N, len, z = NULL, initial = NULL, fixed = NULL, seed=12345) {

  ## If initial values are given, then they will not be randomly initialized
  if(!is.null(initial)) {
    init.values = sort(unique(initial))
    if(length(unique(initial)) != N || length(initial) != len || !all(init.values %in% 1:N)) {
      stop("'initial' needs to be a vector of length 'len' containing N unique values.")
    }
    z = as.numeric(as.factor(initial))
  } else {
    z = rep(NA, len)
  } 
  
  ## Set any values that need to be fixed during sampling
  if(!is.null(fixed)) {
    fixed.values = sort(unique(fixed))
    if(length(fixed) != len || !all(fixed.values %in% 1:N)) {
      stop("'fixed' to be a vector of length 'len' where each entry is one of N unique values or NA.")
    }
    fixed.ix = !is.na(fixed)
    z[fixed.ix] = fixed[fixed.ix]
    z.not.used = setdiff(1:N, unique(fixed[fixed.ix]))
  } else {
    z.not.used = 1:N
    fixed.ix = rep(FALSE, len)
  }  

  ## Randomly sample remaining values
  if (!is.null(seed)) {
    set.seed(seed)
  }
  z.na = which(is.na(z))
  if(length(z.na) > 0) {
    z[z.na] = sample(z.not.used, length(z.na), replace=TRUE)
  }  
  
  ## Check to ensure each value is in the vector at least once
  missing = setdiff(1:N, z)
  for(i in missing) {
    ta = sort(table(z[!fixed.ix]), decreasing = TRUE)
    if(ta[1] == 1) {
      stop("'len' is not long enough to accomodate 'N' unique values")
    }
    ix = which(z == as.numeric(names(ta))[1] & !fixed.ix)
    z[sample(ix, 1)] = i
  }
  
  return(z)
}


recursive.splitZ = function(counts, s, K, K.subcluster, split.all=TRUE, alpha, beta, min.cell = 3, seed=12345) {
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
  	    if(current.K > K) {
  	      break()
  	    }
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


recursive.splitY = function(counts, L, beta, delta, gamma, z=NULL, K=NULL, K.subclusters=NULL, min.feature=3, max.cells=100, seed=12345) {

  ## Decrease number cells into cell populations to save time
  ## If z and K are supplied from previous cell clustering, then each cluster will be split into K.subclusters
  if(!is.null(z) & !is.null(K) & !is.null(K.subclusters)) {
    if(K * K.subclusters > max.cells) {
      K.subclusters = round(max.cells / K)
    }

	z.ta = tabulate(z, K)
	z.non.empty = which(z.ta > 0)
	temp.z = rep(0, length(z))
	current.top.z = 0
	for(i in z.non.empty) { 
	  ix = z == i
	  if(z.ta[i] <= K.subclusters) {
		temp.z[ix] = (current.top.z + 1):(current.top.z + z.ta[i])
	  } else {
		clustLabel = suppressMessages(.celda_C(counts[,z == i], K=K.subclusters, max.iter=5, stop.iter=1, algorithm="EM", nchains=1, split.on.iter=-1, split.on.last=FALSE, verbose=FALSE, initialize="random"))
		temp.z[ix] = clustLabel@clusters$z + current.top.z 
	  }
	  current.top.z = max(temp.z, na.rm=TRUE)
	}
  ## If z or K are not supplied from previous cell clustering, then each cluster will be split into 'max.cells' subclusters
  } else {
	if(ncol(counts) > max.cells) {
	  res = .celda_C(counts, K=max.cells, stop.iter = 1, split.on.iter=-1, split.on.last=FALSE, nchains=3, seed=seed, verbose=FALSE)
	  temp.z = res@clusters$z
	} else {
	  temp.z = 1:ncol(counts)
	}
  }
  
  ## Make collapsed count matrix based on z labels
  counts = colSumByGroup(counts, as.integer(as.factor(temp.z)), length(unique(temp.z)))
  
  ## Perform splitting for y labels
  current.L = 2
  overall.y = rep(1, nrow(counts))
  cluster.splits = matrix(NA, nrow=nrow(counts), ncol=L)
  cluster.split.flag = rep(TRUE, L)
    
  while(current.L <= L) {

    y.ta = tabulate(overall.y, L)
    y.to.split = which(y.ta > min.feature)

    ## Split each cluster <= current.L and number of features > min.feature
    for(i in y.to.split) {
      if(isTRUE(cluster.split.flag[i])) {
        ix = which(overall.y == i)
        res = suppressMessages(.celda_G(counts[ix,], L=2, max.iter = 5, split.on.iter=-1, split.on.last=FALSE, nchains=1, seed=seed, verbose=FALSE, initialize="random"))
        cluster.splits[cbind(ix, i)] = res@clusters$y
        cluster.split.flag[i] = FALSE
      }
    }

    ## Calculate likelihood for each cluster split
    temp.y = matrix(overall.y, nrow=nrow(counts), ncol=length(y.to.split))
	  ll = vector(mode="numeric", length=ncol(temp.y))    
	  for(i in 1:ncol(temp.y)) {
	    temp.y[cluster.splits[,i] == 2,i] = current.L
	    ll[i] = logLikelihood.celda_G(counts=counts, y=temp.y[,i], 
	                                  L=current.L, beta=beta, delta=delta, 
	                                  gamma=gamma)
	}  

    ## Choose best split. Reset flags so the old cluster will be re-evaluated for splitting
    best.ix = which.max(ll)
  	overall.y = temp.y[,best.ix]  
  	cluster.splits[,y.to.split[best.ix]] = NA
    cluster.split.flag[y.to.split[best.ix]] = TRUE
    	
	  current.L = current.L + 1
  }
  
  return(overall.y)
} 



recursive.splitY_v2 = function(counts, L, temp.cell.cluster=1000, cell.cluster.breaks=10, beta=1, delta=1, gamma=1, min.feature=3, seed=12345) {
  
  s = rep(1, ncol(counts))
  if(!is.null(temp.cell.cluster) && temp.cell.cluster > 0) {
    new.z = recursive.splitZ(counts, s, K=temp.cell.cluster, K.subcluster=cell.cluster.breaks, split.all=TRUE, alpha=1, beta=1, min.cell = 3, seed=seed)
    counts = colSumByGroup(counts, new.z, length(unique(new.z)))
  }
  
  model.initial = .celda_G(counts, L=2, split.on.iter=-1, split.on.last=FALSE, nchains=1, verbose=T)

  ## Perform splitting for y labels
  current.L = 3
  overall.y = model.initial@clusters$y
  
  y.history = matrix(NA, nrow=nrow(counts), ncol=L)
  y.history[,2] = overall.y
  
  ll.history = c()
  

  while(current.L <= L) {
print(paste0("currentL: ", current.L))
    y.ta = tabulate(overall.y, L)
    y.to.split = which(y.ta > min.feature)

    clust.split = vector("list", current.L)
    for(i in y.to.split) {
      if(is.null(clust.split[[i]]) || is.na(clust.split[[i]])) {
        if(y.ta[i] >= min.feature) {
          print(i)
          clustLabel = .celda_G(counts[overall.y == i,,drop=FALSE], L=2, split.on.iter=-1, split.on.last=FALSE, nchains=1, verbose=FALSE)
          clust.split[[i]] = as.integer(as.factor(clustLabel@clusters$y))
          if(length(unique(clust.split[[i]])) != 2) {
            y.to.split = setdiff(y.to.split, i)
          }
        }  
      } 
    }
    

    ## Calculate likelihood for each cluster split
    temp.y = matrix(overall.y, nrow=nrow(counts), ncol=length(y.to.split))
    ll = vector(mode="numeric", length=ncol(temp.y))    
    for(i in 1:ncol(temp.y)) {
      ix = overall.y == y.to.split[i]
      temp.y[ix,i] = ifelse(clust.split[[y.to.split[i]]] == 2, y.to.split[i], current.L)
      ll[i] = logLikelihood.celda_G(counts=counts, y=temp.y[,i], L=current.L, beta=beta, delta=delta, gamma=gamma)
    }  
    
    ## Choose best split. Reset flags so the old cluster will be re-evaluated for splitting
    best.ix = which.max(ll)
    #temp.y = temp.y[,best.ix] 
    temp.model = .celda_G(counts, L=current.L, stop.iter=5, split.on.iter=-1, split.on.last=FALSE, nchains=1, verbose=FALSE, y.init=temp.y[,best.ix])
    overall.y = temp.model@clusters$y
    
    ## Get list of modules whose features have changed membership.
    ##The best split for these clusters need to b
    prev.y = y.history[,current.L-1]
    module.change = unique(prev.y[overall.y != prev.y])
    print(paste0("module change: ", module.change))
    for(i in module.change) {
      clust.split[[i]] = NA 
    }
    ll.history = c(ll.history, temp.model@finalLogLik)
    clust.split = vector("list", L)
    y.history[,current.L] = overall.y
    current.L = current.L + 1
  }
  
  return(list(y=y.history, ll=ll.history))
} 

splitModule = function(counts, celda.mod, module, n=2) {
  ix = celda.mod@clusters$y == module
  new.L = max(celda.mod@clusters$y) + n - 1
  if(sum(ix) > 1) {
    temp.model = .celda_G(counts[ix,,drop=FALSE], L=n, split.on.iter=-1, split.on.last=FALSE, nchains=1, verbose=FALSE)  
    new.y = celda.mod@clusters$y
    new.y[ix] = ifelse(temp.model@clusters$y == 1, module, new.L)
    model = .celda_G(counts, L=new.L, split.on.iter=-1, split.on.last=FALSE, nchains=1, verbose=FALSE, y.init=new.y)
  } else {
    model = celda.mod
  }
  
  return(model)
}


