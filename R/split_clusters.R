split.z = function(counts, z, empty.K, K, min.cell=3, LLFunction, ...) { 
  
  ## Normalize counts to fraction for cosine clustering
  counts.norm = normalizeCounts(counts, scale.factor=1)

  ## Identify other clusters to split
  z.ta = table(factor(z, levels=1:K))
  k.pass.min = which(z.ta >= min.cell)
  k.to.test = setdiff(k.pass.min, empty.K)

  ## Set up variables for holding results
  z.split = matrix(z, ncol=K, nrow=length(z))
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
  
  message(date(), " ... Cluster ", empty.K, " had ", z.ta[empty.K], " cells. Splitting Cluster ", k.to.test[k.to.test.select], " into two clusters.")
  return(z.split[,k.to.test.select])
}




split.each.z = function(counts, z, K, min.cell=3, LLFunction, ...) { 
  ## Normalize counts to fraction for cosine clustering
  counts.norm = normalizeCounts(counts, scale.factor=1)

  ## Identify clusters to split
  z.ta = table(factor(z, levels=1:K))
  z.to.split = which(z.ta >= min.cell)

  ## Set up variables for holding results with current iteration of z
  params = c(list(counts=counts, z=z, K=K), list(...))
  z.split.ll = do.call(LLFunction, params)
  z.split = z

  ## Loop through each split-able Z and perform split
  clust.split = list()
  for(i in z.to.split) {
    d = cosineDist(counts.norm[,z == i])
    clustLabel = clusterByHC(d, min=min.cell)
    clust.split = c(clust.split, list(clustLabel))
  }

  pairs = c()
  for(i in 1:K) {
    z.to.test = setdiff(z.to.split, i)
    temp.z = z

    ## Randomly assign z labels of current cluster to other clusters
    ix = z == i
    temp.z[ix] = sample(1:K, sum(ix), replace=TRUE)

    for(j in z.to.test) {
      new.z = temp.z
      ix = z == j
      new.z[ix] = ifelse(clust.split[[j]] == 1, j, i)

      ## Calculate likelihood of split
      params = c(list(counts=counts, z=new.z, K=K), list(...))
      new.ll = do.call(LLFunction, params)

      z.split.ll = c(z.split.ll, new.ll)
      z.split = cbind(z.split, new.z)
      
      pairs = rbind(pairs, c(i, j))
    }
  }
  
  select = sample.ll(z.split.ll) 
  
  if(select == 1) {
    message(date(), " ... No additional splitting was performed.") 
  } else {
    message(date(), " ... Cluster ", pairs[select,1], " was randomly redistributed and cluster ", pairs[select,2], " was split in two.")
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
  
  ## Set up variables for holding results
  y.split = matrix(y, ncol=L, nrow=length(y))
  l.pass.min = which(y.ta >= min.gene)
  l.split.ll = rep(NA, length(l.pass.min))
  
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
  
  message(date(), " ... Cluster ", empty.L, " had ", y.ta[empty.L], " genes. Splitting Cluster ", l.to.test[l.to.test.select], " into two clusters.")
  return(y.split[,l.to.test.select])
}






clusterByHC = function(d, min=1, method="ward.D") {
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




