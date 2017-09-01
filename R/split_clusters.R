split.z = function(counts, z, empty.K, K, min.cell=3, LLFunction, ...) { 
  
  ## Normalize counts to fraction for cosine clustering
  counts.norm = normalizeCounts(counts, scale.factor=1)

  ## Identify other clusters to split
  z.ta = table(factor(z, levels=1:K))
  k.pass.min = which(z.ta >= min.cell)
  k.to.test = setdiff(k.pass.min, empty.K)

  if(length(k.to.test) == 0) {
    message(date(), " ... Cluster sizes too small. No additional splitting was performed.") 
    break
  }

  ## Set up variables for holding results
  z.split = matrix(z, ncol=length(k.to.test), nrow=length(z))
  k.split.ll = rep(NA, length(k.to.test))

  ## Loop through each cluster, split, and determine logLik
  for(i in 1:length(k.to.test)) {
    
    k.dist = cosineDist(counts.norm[,z == k.to.test[i]])
    clustLabel = clusterByHC(k.dist, min=min.cell)
	clustLabel.final = ifelse(clustLabel == 1, k.to.test[i], empty.K)
	
	## Assign new labels to test cluster    
	ix = (z == k.to.test[i])
	z.split[ix,i] = clustLabel.final

    ## Calculate likelihood of split
    params = c(list(counts=counts, z=z.split[,i], K=K), list(...))
    k.split.ll[i] = do.call(LLFunction, params)
  }
  k.to.test.select = sample.ll(k.split.ll)

  message(date(), " ... Cell cluster ", empty.K, " had ", z.ta[empty.K], " cells. Splitting Cluster ", k.to.test[k.to.test.select], " into two clusters.")
  return(z.split[,k.to.test.select])
}




split.each.z = function(counts, z, K, LLFunction, min=3, ...) { 
  ## Normalize counts to fraction for cosine clustering
  counts.norm = normalizeCounts(counts, scale.factor=1)
  
  ## Identify clusters to split
  z.ta = table(factor(z, levels=1:K))
  z.to.split = which(z.ta >= min)

  if(length(z.to.split) == 0) {
    message(date(), " ... Cluster sizes too small. No additional splitting was performed.") 
    break
  }
  
  ## For each cluster, determine which other cluster is most closely related
  counts.z.collapse = t(rowsum(t(counts), group=z, reorder=TRUE))
  counts.z.collapse.norm = normalizeCounts(counts.z.collapse, scale.factor=1)
  counts.z.cor = cosine(t(counts.z.collapse.norm))
  diag(counts.z.cor) = 0
  
  ## Loop through each split-able Z and perform split
  clust.split = list()
  for(i in 1:K) { 
    if(i %in% z.to.split) {
      d = cosineDist(counts.norm[,z == i])
      clustLabel = clusterByHC(d)
      clust.split = c(clust.split, list(clustLabel))
    } else {
      clust.split = c(clust.split, list(NA))
    }  
  }

  ## Set up variables for holding results with current iteration of z
  params = c(list(counts=counts, z=z, K=K), list(...))
  z.split.ll = do.call(LLFunction, params)
  z.split = z

  pairs = c(NA, NA, NA)
  for(i in 1:K) {
	for(j in setdiff(z.to.split, i)) {
	  new.z = z
	  
      ## Assign cluster i to the next most similar cluster (excluding cluster j) 
      ## as defined above by the spearman correlation      
      ix.to.move = z == i
      h = setdiff(order(counts.z.cor[,i], decreasing=TRUE), j)[1]
      new.z[ix.to.move] = h
            
      ## Split cluster j according to the clustering defined above
      ix.to.split = z == j
      new.z[ix.to.split] = ifelse(clust.split[[j]] == 1, j, i)

	  ## Calculate likelihood of split
	  params = c(list(counts=counts, z=new.z, K=K), list(...))
	  new.ll = do.call(LLFunction, params)

	  z.split.ll = c(z.split.ll, new.ll)
	  z.split = cbind(z.split, new.z)
   
	  pairs = rbind(pairs, c(i, j, h))
	}  
  }
  
  select = sample.ll(z.split.ll) 

  if(select == 1) {
    message(date(), " ... No additional splitting was performed.") 
  } else {
    message(date(), " ... Cluster ", pairs[select,1], " was moved to cluster ", pairs[select,3], " and cluster ", pairs[select,2], " was split in two.")
  } 
  
  return(z.split[,select])
}







split.y = function(counts, y, empty.L, L, min.gene=3, LLFunction, ...) { 

  ## Normalize counts to fraction for cosine clustering
  counts.norm = normalizeCounts(counts, scale.factor=1)
  
  ## Identify other clusters to split
  y.ta = table(factor(y, levels=1:L))
  l.pass.min = which(y.ta >= min.gene)
  l.to.test = setdiff(l.pass.min, empty.L)
  
  if(length(l.to.test) == 0) {
    message(date(), " ... Cluster sizes too small. No additional splitting was performed.") 
    break
  }

  ## Set up variables for holding results
  y.split = matrix(y, ncol=length(l.to.test), nrow=length(y))
  l.split.ll = rep(NA, length(l.to.test))
  
  ## Loop through each cluster, split, and determine logLik
  for(i in 1:length(l.to.test)) {
    
    ind = y == l.to.test[i]
    l.dist = spearmanDist(t(counts.norm[ind,]))

    clustLabel = clusterByHC(l.dist, min=min.gene)
	clustLabel.final = ifelse(clustLabel == 1, l.to.test[i], empty.L)

	## Assign new labels to test cluster    
	ix = (y == l.to.test[i])
	y.split[ix,i] = clustLabel.final

    ## Calculate likelihood of split
    params = c(list(counts=counts, y=y.split[,i], L=L), list(...))
    l.split.ll[i] = do.call(LLFunction, params)
  }
 
  l.to.test.select = sample.ll(l.split.ll)
  
  message(date(), " ... Gene cluster ", empty.L, " had ", y.ta[empty.L], " genes. Splitting Cluster ", l.to.test[l.to.test.select], " into two clusters.")
  return(y.split[,l.to.test.select])
}




split.each.y = function(counts, y, L, LLFunction, min=3, ...) { 
  ## Normalize counts to fraction for hierarchical clustering
  counts.norm = normalizeCounts(counts, scale.factor=1)
  
  ## Identify clusters to split
  y.ta = table(factor(y, levels=1:L))
  y.to.split = which(y.ta >= min)

  if(length(y.to.split) == 0) {
    message(date(), " ... Cluster sizes too small. No additional splitting was performed.") 
    break
  }

  ## For each cluster, determine which other cluster is most closely related
  counts.y.collapse = rowsum(counts, group=y, reorder=TRUE)
  counts.y.collapse.norm = normalizeCounts(counts.y.collapse, scale.factor=1)
  counts.y.cor = cor(t(counts.y.collapse.norm), method="spearman")
  diag(counts.y.cor) = 0
  
  ## Loop through each split-able y and find best split
  clust.split = list()
  for(i in 1:L) {
    if(i %in% y.to.split) {    
      d = spearmanDist(t(counts.norm[y == i,]))
      clustLabel = clusterByHC(d)
      clust.split = c(clust.split, list(clustLabel))
    } else {
      clust.split = c(clust.split, list(NA))
    }  
  }

  ## Set up variables for holding results with current iteration of y
  params = c(list(counts=counts, y=y, L=L), list(...))
  y.split.ll = do.call(LLFunction, params)
  y.split = y

  pairs = c(NA, NA, NA)
  for(i in 1:L) {
    for(j in setdiff(y.to.split, i)) {
      new.y = y
      
      ## Assign cluster i to the next most similar cluster (excluding cluster j) 
      ## as defined above by the spearman correlation      
      ix.to.move = y == i
      h = setdiff(order(counts.y.cor[,i], decreasing=TRUE), j)[1]
      new.y[ix.to.move] = h
            
      ## Split cluster j according to the clustering defined above
      ix.to.split = y == j
      new.y[ix.to.split] = ifelse(clust.split[[j]] == 1, j, i)
      
      ## Calculate likelihood of split
      params = c(list(counts=counts, y=new.y, L=L), list(...))
      new.ll = do.call(LLFunction, params)

      y.split.ll = c(y.split.ll, new.ll)
      y.split = cbind(y.split, new.y)
      
      pairs = rbind(pairs, c(i, j, h))
    }
  }
  
  select = sample.ll(y.split.ll) 
  
  if(select == 1) {
    message(date(), " ... No additional splitting was performed.") 
  } else {
    message(date(), " ... Cluster ", pairs[select,1], " was moved to cluster ", pairs[select,3], " and cluster ", pairs[select,2], " was split in two.")
  } 
  
  return(y.split[,select])
}





clusterByHC = function(d, min=3, method="ward.D") {
  label = cluster::pam(x = d, k=2)$clustering

  ## If PAM split is too small, perform secondary hclust procedure to split into roughly equal groups
  if(min(table(label)) < min) {
	d.hclust = hclust(d, method=method)

	## Get maximum sample size of each subcluster
	d.hclust.size = sapply(1:length(d.hclust$height), function(i) max(table(cutree(d.hclust, h=d.hclust$height[i]))))
  
	## Find the height of the dendrogram that best splits the samples in "roughly" half
	sample.size = round(length(d.hclust$order)/ 2)
	d.hclust.select = which.min(abs(d.hclust.size - sample.size))
	temp = cutree(d.hclust, h=d.hclust$height[d.hclust.select])
	
	## Set half of the samples as the largest cluster and the other samples to the other cluster
	label.max.cluster = which.max(table(temp))
	label = ifelse(temp == label.max.cluster, 1, 2)
  } 
  
  return(label)
}




