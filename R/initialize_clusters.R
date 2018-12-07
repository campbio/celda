
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
  set.seed(seed)
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


recursive.splitZ = function(counts, s, K, alpha, beta, min.cell = 3, seed=12345) {
  current.K = 2
  overall.z = rep(1, ncol(counts))
  cluster.splits = matrix(NA, nrow=ncol(counts), ncol=K)
  cluster.split.flag = rep(TRUE, K)
    
  while(current.K <= K) {

    z.ta = tabulate(overall.z, K)
    z.to.split = which(z.ta > min.cell)

    ## Split each cluster <= current.K and number of cells > min.cell
    for(i in z.to.split) {
      if(isTRUE(cluster.split.flag[i])) {
        ix = which(overall.z == i)
        res = suppressMessages(.celda_C(counts[,ix], K=2, stop.iter = 1, split.on.iter=-1, split.on.last=FALSE, nchains=1, seed=seed, verbose=FALSE, initialize="random"))
        cluster.splits[cbind(ix, i)] = res@clusters$z
        cluster.split.flag[i] = FALSE
      }
    }

    ## Calculate likelihood for each cluster split
    temp.z = matrix(overall.z, nrow=ncol(counts), ncol=length(z.to.split))
  	ll = rep(NA, ncol(temp.z))    
  	for(i in 1:ncol(temp.z)) {
  	  temp.z[cluster.splits[,i] == 2,i] = current.K 
  	  ll[i] = logLikelihood.celda_C(counts, "celda_C", s, temp.z[,i], current.K, alpha, beta)	  
  	}  

    ## Choose best split. Reset flags so the old cluster will be re-evaluated for splitting
    best.ix = which.max(ll)
  	overall.z = temp.z[,best.ix]  
  	cluster.splits[,z.to.split[best.ix]] = NA
    cluster.split.flag[z.to.split[best.ix]] = TRUE
    	
	  current.K = current.K + 1
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
	ll = rep(NA, ncol(temp.y))    
	for(i in 1:ncol(temp.y)) {
	  temp.y[cluster.splits[,i] == 2,i] = current.L
	  ll[i] = logLikelihood.celda_G(counts=counts, y=temp.y[,i], L=current.L, beta=beta, delta=delta, gamma=gamma)
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





