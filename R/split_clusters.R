# cC.calcLL = function(m.CP.by.S, n.G.by.CP, s, z, K, nS, nG, alpha, beta) 
cC.splitZ = function(counts, m.CP.by.S, n.G.by.CP, n.CP, s, z, K, nS, nG, alpha, beta, z.prob, max.clusters.to.try=10, min.cell=3) {

  ## Identify clusters to split
  z.ta = tabulate(z, K)
  z.to.split = which(z.ta >= min.cell)
  z.non.empty = which(z.ta > 0)
  
  if(length(z.to.split) == 0) {
    m = paste0(date(), " .... Cluster sizes too small. No additional splitting was performed.") 
    return(list(z=z, m.CP.by.S, n.G.by.CP, n.CP=n.CP, message=m))  
  }
  
  ## Loop through each split-able Z and perform split
  clust.split = vector("list", K)
  for(i in z.to.split) { 
    clustLabel = suppressMessages(.celda_C(counts[,z == i], K=2, max.iter=5, split.on.iter=-1, split.on.last=FALSE))
    clust.split[[i]] = clustLabel@clusters$z
  }

  ## Find second best assignment give current assignments for each cell
  z.prob[cbind(1:nrow(z.prob), z)] = NA
  z.second = apply(z.prob, 1, which.max)

  ## Set up initial variables
  z.split = matrix(NA, nrow=length(z), ncol=length(z.to.split) * max.clusters.to.try)
  z.split.ll = rep(NA, ncol=length(z.to.split) * max.clusters.to.try)  
  z.split.ll[1] = cC.calcLL(m.CP.by.S, n.G.by.CP, s, z, K, nS, nG, alpha, beta) 
  z.split[,1] = z

  ## Select worst clusters to test for reshuffling  
  previous.z = z
  ll.shuffle = rep(NA, K)
  for(i in z.non.empty) {
    ix = z == i
    new.z = z
    new.z[ix] = z.second[ix]
    
    p = cC.reDecomposeCounts(counts, s, new.z, previous.z, n.G.by.CP, K)
    n.G.by.CP = p$n.G.by.CP
    m.CP.by.S = p$m.CP.by.S
    ll.shuffle[i] = cC.calcLL(m.CP.by.S, n.G.by.CP, s, z, K, nS, nG, alpha, beta) 
    previous.z = new.z
  } 
  z.to.shuffle = utils::head(order(ll.shuffle, decreasing = TRUE, na.last=NA), n = max.clusters.to.try)

  
  pairs = c(NA, NA)
  split.ix = 2
  for(i in z.to.shuffle) {
  
    other.clusters = setdiff(z.to.split, i)
   
	for(j in other.clusters) {
	  new.z = z
	  
      ## Assign cluster i to the next most similar cluster (excluding cluster j) 
      ## as defined above by the correlation      
      ix.to.move = z == i
      new.z[ix.to.move] = z.second[ix.to.move]
            
      ## Split cluster j according to the clustering defined above
      ix.to.split = z == j
      new.z[ix.to.split] = ifelse(clust.split[[j]] == 1, j, i)

      p = cC.reDecomposeCounts(counts, s, new.z, previous.z, n.G.by.CP, K)
      n.G.by.CP = p$n.G.by.CP
      m.CP.by.S = p$m.CP.by.S
      
	  ## Calculate likelihood of split
	  z.split.ll[split.ix] = cC.calcLL(m.CP.by.S, n.G.by.CP, s, z, K, nS, nG, alpha, beta) 
	  z.split[,split.ix] = new.z
	  split.ix = split.ix + 1L
      previous.z = new.z
      
	  pairs = rbind(pairs, c(i, j))
	}  
  }

  select = which.max(z.split.ll) 

  if(select == 1) {
    m = paste0(date(), " .... No additional splitting was performed.") 
  } else {
    m = paste0(date(), " .... Cluster ", pairs[select,1], " was reassigned and cluster ", pairs[select,2], " was split in two.")
  } 

  p = cC.reDecomposeCounts(counts, s, z.split[,select], previous.z, n.G.by.CP, K)
  return(list(z=z.split[,select], m.CP.by.S=p$m.CP.by.S, n.G.by.CP=p$n.G.by.CP, n.CP=p$n.CP, message=m))
}


# cCG.calcLL = function(K, L, m.CP.by.S, n.TS.by.CP, n.by.G, n.by.TS, nG.by.TS, nS, nG, alpha, beta, delta, gamma) 
cCG.splitZ = function(counts, m.CP.by.S, n.TS.by.C, n.TS.by.CP, n.by.G, n.by.TS, nG.by.TS, n.CP, s, z, K, L, nS, nG, alpha, beta, delta, gamma, z.prob, max.clusters.to.try=10, min.cell=3) {

  ## Identify clusters to split
  z.ta = tabulate(z, K)
  z.to.split = which(z.ta >= min.cell)
  z.non.empty = which(z.ta > 0)
  
  if(length(z.to.split) == 0) {
    m = paste0(date(), " .... Cluster sizes too small. No additional splitting was performed.") 
    return(list(z=z, m.CP.by.S=m.CP.by.S, n.TS.by.CP=n.TS.by.CP, n.CP=n.CP, message=m))  
  }
  
  ## Loop through each split-able Z and perform split
  clust.split = vector("list", K)
  for(i in z.to.split) { 
    clustLabel = suppressMessages(.celda_C(counts[,z == i], K=2, max.iter=5, split.on.iter=-1, split.on.last=FALSE))
    clust.split[[i]] = clustLabel@clusters$z
  }

  ## Find second best assignment give current assignments for each cell
  z.prob[cbind(1:nrow(z.prob), z)] = NA
  z.second = apply(z.prob, 1, which.max)

  ## Set up initial variables
  z.split = matrix(NA, nrow=length(z), ncol=length(z.to.split) * max.clusters.to.try)
  z.split.ll = rep(NA, ncol=length(z.to.split) * max.clusters.to.try)  
  z.split.ll[1] = cCG.calcLL(K, L, m.CP.by.S, n.TS.by.CP, n.by.G, n.by.TS, nG.by.TS, nS, nG, alpha, beta, delta, gamma) 
  z.split[,1] = z

  ## Select worst clusters to test for reshuffling  
  previous.z = z
  ll.shuffle = rep(NA, K)
  for(i in z.non.empty) {
    ix = z == i
    new.z = z
    new.z[ix] = z.second[ix]
    
    p = cC.reDecomposeCounts(n.TS.by.C, s, new.z, previous.z, n.TS.by.CP, K)
    n.TS.by.CP = p$n.G.by.CP
    m.CP.by.S = p$m.CP.by.S
    ll.shuffle[i] = cCG.calcLL(K, L, m.CP.by.S, n.TS.by.CP, n.by.G, n.by.TS, nG.by.TS, nS, nG, alpha, beta, delta, gamma) 
    previous.z = new.z
  } 
  z.to.shuffle = utils::head(order(ll.shuffle, decreasing = TRUE, na.last=NA), n = max.clusters.to.try)

  
  pairs = c(NA, NA)
  split.ix = 2
  for(i in z.to.shuffle) {
  
    other.clusters = setdiff(z.to.split, i)
   
	for(j in other.clusters) {
	  new.z = z
	  
      ## Assign cluster i to the next most similar cluster (excluding cluster j) 
      ## as defined above by the correlation      
      ix.to.move = z == i
      new.z[ix.to.move] = z.second[ix.to.move]
            
      ## Split cluster j according to the clustering defined above
      ix.to.split = z == j
      new.z[ix.to.split] = ifelse(clust.split[[j]] == 1, j, i)

      p = cC.reDecomposeCounts(n.TS.by.C, s, new.z, previous.z, n.TS.by.CP, K)
      n.TS.by.CP = p$n.G.by.CP
      m.CP.by.S = p$m.CP.by.S
      
	  ## Calculate likelihood of split
	  z.split.ll[split.ix] = cCG.calcLL(K, L, m.CP.by.S, n.TS.by.CP, n.by.G, n.by.TS, nG.by.TS, nS, nG, alpha, beta, delta, gamma) 
	  z.split[,split.ix] = new.z
	  split.ix = split.ix + 1L
      previous.z = new.z
      
	  pairs = rbind(pairs, c(i, j))
	}  
  }

  select = which.max(z.split.ll) 

  if(select == 1) {
    m = paste0(date(), " .... No additional splitting was performed.") 
  } else {
    m = paste0(date(), " .... Cluster ", pairs[select,1], " was reassigned and cluster ", pairs[select,2], " was split in two.")
  } 

  p = cC.reDecomposeCounts(n.TS.by.C, s, z.split[,select], previous.z, n.TS.by.CP, K)
  return(list(z=z.split[,select], m.CP.by.S=p$m.CP.by.S, n.TS.by.CP=p$n.G.by.CP, n.CP=p$n.CP, message=m))
}



cCG.splitY = function(counts, y, m.CP.by.S, n.G.by.CP, n.TS.by.C, n.TS.by.CP, n.by.G, n.by.TS, nG.by.TS, n.CP, s, z, K, L, nS, nG, alpha, beta, delta, gamma, y.prob, max.clusters.to.try=10, K.subclusters=10, min.cell=3) {

  #########################
  ## First, the cell dimension of the original matrix will be reduced by splitting each z cluster into 'K.subclusters'.
  #########################
  
  ## This will not be as big as the original matrix (which can take a lot of time to process with large number of cells), but not as small as the 'n.G.by.CP' with current z assignments
  z.ta = tabulate(z, K)
  z.non.empty = which(z.ta > 0)
  temp.z = rep(0, length(z))
  current.top.z = 0
  for(i in z.non.empty) { 
    ix = z == i
    if(z.ta[i] <= K.subclusters) {
      temp.z[ix] = (current.top.z + 1):(current.top.z + z.ta[i])
    } else {
      clustLabel = suppressMessages(.celda_C(counts[,z == i], K=K.subclusters, max.iter=5, split.on.iter=-1, split.on.last=FALSE))
      temp.z[ix] = clustLabel@clusters$z + current.top.z 
    }
    current.top.z = max(temp.z, na.rm=TRUE)
  }
  
  ## Decompose counts according to new/temp z labels
  temp.n.G.by.CP = colSumByGroup(counts, group=temp.z, K=current.top.z)

  #########################
  ## Second, different y splits will be estimated and tested
  #########################
  
  ## Identify clusters to split
  y.ta = tabulate(y, L)
  y.to.split = which(y.ta >= min.cell)
  y.non.empty = which(y.ta > 0)
  
  if(length(y.to.split) == 0) {
    m = paste0(date(), " .... Cluster sizes too small. No additional splitting was performed.") 
    return(list(y=y, m.CP.by.S=m.CP.by.S, n.TS.by.CP=n.TS.by.CP, n.CP=n.CP, message=m))  
  }

  ## Loop through each split-able Z and perform split
  clust.split = vector("list", L)
  for(i in y.to.split) { 
    clustLabel = suppressMessages(.celda_G(temp.n.G.by.CP[y == i,], L=2, max.iter=5, split.on.iter=-1, split.on.last=FALSE))
    clust.split[[i]] = clustLabel@clusters$y
  }

  ## Find second best assignment give current assignments for each cell
  y.prob[cbind(1:nrow(y.prob), y)] = NA
  y.second = apply(y.prob, 1, which.max)

  ## Set up initial variables
  y.split = matrix(NA, nrow=length(y), ncol=length(y.to.split) * max.clusters.to.try)
  y.split.ll = rep(NA, ncol=length(y.to.split) * max.clusters.to.try)  
  y.split.ll[1] = cCG.calcLL(K, L, m.CP.by.S, n.TS.by.CP, n.by.G, n.by.TS, nG.by.TS, nS, nG, alpha, beta, delta, gamma) 
  y.split[,1] = y

  ## Select worst clusters to test for reshuffling  
  previous.y = y
  ll.shuffle = rep(NA, L)
  for(i in y.non.empty) {
    ix = y == i
    new.y = y
    new.y[ix] = y.second[ix]
    
    p = cG.reDecomposeCounts(n.G.by.CP, new.y, previous.y, n.TS.by.CP, n.by.G, L)
    n.TS.by.CP = p$n.TS.by.C
    n.by.TS = p$n.by.TS
    nG.by.TS = p$nG.by.TS
    
    ll.shuffle[i] = cCG.calcLL(K, L, m.CP.by.S, n.TS.by.CP, n.by.G, n.by.TS, nG.by.TS, nS, nG, alpha, beta, delta, gamma) 
    previous.y = new.y
  } 
  y.to.shuffle = utils::head(order(ll.shuffle, decreasing = TRUE, na.last=NA), n = max.clusters.to.try)
  
  pairs = c(NA, NA)
  split.ix = 2
  for(i in y.to.shuffle) {
  
    other.clusters = setdiff(y.to.split, i)
   
	for(j in other.clusters) {
	  new.y = y
	  
      ## Assign cluster i to the next most similar cluster (excluding cluster j) 
      ## as defined above by the correlation      
      ix.to.move = y == i
      new.y[ix.to.move] = y.second[ix.to.move]
            
      ## Split cluster j according to the clustering defined above
      ix.to.split = y == j
      new.y[ix.to.split] = ifelse(clust.split[[j]] == 1, j, i)

      p = cG.reDecomposeCounts(n.G.by.CP, new.y, previous.y, n.TS.by.CP, n.by.G, L)
      n.TS.by.CP = p$n.TS.by.C
      n.by.TS = p$n.by.TS
      nG.by.TS = p$nG.by.TS
      
	  ## Calculate likelihood of split
	  y.split.ll[split.ix] = cCG.calcLL(K, L, m.CP.by.S, n.TS.by.CP, n.by.G, n.by.TS, nG.by.TS, nS, nG, alpha, beta, delta, gamma) 
	  y.split[,split.ix] = new.y
	  split.ix = split.ix + 1L
      previous.y = new.y
      
	  pairs = rbind(pairs, c(i, j))
	}  
  }

  select = which.max(y.split.ll) 

  if(select == 1) {
    m = paste0(date(), " .... No additional splitting was performed.") 
  } else {
    m = paste0(date(), " .... Cluster ", pairs[select,1], " was reassigned and cluster ", pairs[select,2], " was split in two.")
  } 

  p = cG.reDecomposeCounts(n.G.by.CP, y.split[,select], previous.y, n.TS.by.CP, n.by.G, L)
  return(list(y=y.split[,select], n.TS.by.CP=p$n.TS.by.C, n.by.TS=p$n.by.TS, nG.by.TS=p$nG.by.TS, message=m))
}





#cG.calcLL = function(n.TS.by.C, n.by.TS, n.by.G, nG.by.TS, nM, nG, L, beta, delta, gamma) {
cG.splitY = function(counts, y, n.TS.by.C, n.by.TS, n.by.G, nG.by.TS, nM, nG, L, beta, delta, gamma, y.prob, min.feature=3, max.clusters.to.try=10) { 
  
  ## Identify clusters to split
  y.ta = table(factor(y, levels=1:L))
  y.to.split = which(y.ta >= min.feature)
  y.non.empty = which(y.ta > 0)

  if(length(y.to.split) == 0) {
    m = paste0(date(), " .... Cluster sizes too small. No additional splitting was performed.") 
    return(list(y=y, n.TS.by.C=n.TS.by.C, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, message=m))
  }

  ## Loop through each split-able y and find best split
  clust.split = vector("list", L)
  for(i in y.to.split) {
    clustLabel = suppressMessages(.celda_G(counts[y == i,], L=2, max.iter=5, split.on.iter=-1, split.on.last=FALSE))
    clust.split[[i]] = clustLabel@clusters$y
  }

  ## Find second best assignment give current assignments for each cell
  y.prob[cbind(1:nrow(y.prob), y)] = NA
  y.second = apply(y.prob, 1, which.max)

  ## Set up initial variables
  y.split = matrix(NA, nrow=length(y), ncol=length(y.to.split) * max.clusters.to.try)
  y.split.ll = rep(NA, ncol=length(y.to.split) * max.clusters.to.try)  
  y.split.ll[1] = cG.calcLL(n.TS.by.C, n.by.TS, n.by.G, nG.by.TS, nM, nG, L, beta, delta, gamma)
  y.split[,1] = y

  ## Select worst clusters to test for reshuffling  
  ll.shuffle = rep(NA, L)
  previous.y = y
  for(i in y.non.empty) {
    ix = y == i
    new.y = y
    new.y[ix] = y.second[ix]
    p = cG.reDecomposeCounts(counts, new.y, previous.y, n.TS.by.C, n.by.G, L)
    ll.shuffle[i] = cG.calcLL(p$n.TS.by.C, p$n.by.TS, n.by.G, p$nG.by.TS, nM, nG, L, beta, delta, gamma)
    previous.y = new.y
  }
  y.to.shuffle = utils::head(order(ll.shuffle, decreasing = TRUE, na.last=NA), n = max.clusters.to.try)
  
  pairs = c(NA, NA)
  split.ix = 2  
  for(i in y.to.shuffle) {

    other.clusters = setdiff(y.to.split, i)
    
    for(j in other.clusters) {
      new.y = y
      
      ## Assign cluster i to the next most similar cluster (excluding cluster j) 
      ## as defined above by the spearman correlation      
      ix.to.move = y == i
      new.y[ix.to.move] = y.second[ix.to.move]
            
      ## Split cluster j according to the clustering defined above
      ix.to.split = y == j
      new.y[ix.to.split] = ifelse(clust.split[[j]] == 1, j, i)
      
      ## Calculate likelihood of split
      p = cG.reDecomposeCounts(counts, new.y, previous.y, n.TS.by.C, n.by.G, L)
	  y.split.ll[split.ix] = cG.calcLL(p$n.TS.by.C, p$n.by.TS, n.by.G, p$nG.by.TS, nM, nG, L, beta, delta, gamma)
	  y.split[,split.ix] = new.y
	  split.ix = split.ix + 1L
      previous.y = new.y
      
      pairs = rbind(pairs, c(i, j))
    }
  }

  select = which.max(y.split.ll) 
  
  if(select == 1) {
    m = paste0(date(), " .... No additional splitting was performed.") 
  } else {
    m = paste0(date(), " .... Cluster ", pairs[select,1], " was reassigned and cluster ", pairs[select,2], " was split in two.")
  } 
 
  p = cG.reDecomposeCounts(counts, y.split[,select], previous.y, n.TS.by.C, n.by.G, L)
  return(list(y=y.split[,select], n.TS.by.C=p$n.TS.by.C, n.by.TS=p$n.by.TS, nG.by.TS=p$nG.by.TS, message=m))
}





