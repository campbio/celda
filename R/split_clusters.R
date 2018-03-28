split.each.z = function(counts, z, K, LLFunction, z.prob, min.cell=3, max.clusters.to.try = 10, ...) { 

  dot.args = list(...)
  
  ## Identify clusters to split
  z.ta = table(factor(z, levels=1:K))
  z.to.split = which(z.ta >= min.cell)
  z.non.empty = which(z.ta > 0)
  
  ## Loop through each split-able Z and perform split
  clust.split = list()
  for(i in 1:K) { 
    if(i %in% z.to.split) {
      clustLabel = suppressMessages(celda_C(counts[,z == i], K=2, max.iter=3, split.on.iter=-1))
      clust.split = c(clust.split, list(clustLabel$z))
    } else {
      clust.split = c(clust.split, list(NA))
    }  
  }
  
  ## Find second best assignment give current assignments for each cell
  z.prob[cbind(1:nrow(z.prob), z)] = NA
  z.second = apply(z.prob, 1, which.max)
  
  ll.shuffle = rep(NA, K)
  for(i in z.non.empty) {
    ix = z == i
    new.z = z
    new.z[ix] = z.second[ix]
    params = c(list(counts=counts, z=new.z, K=K), dot.args)
    ll.shuffle[i] = do.call(LLFunction, params)
  }
  z.to.shuffle = head(order(ll.shuffle, decreasing = TRUE, na.last=NA), n = max.clusters.to.try)
  
  if(length(z.to.split) == 0 | length(z.to.shuffle) == 0) {
    m = paste0(date(), " ... Cluster sizes too small. No additional splitting was performed.") 
    return(list(z=z,message=m))
  }
  
  ## Set up variables for holding results with current iteration of z
  params = c(list(counts=counts, z=z, K=K), dot.args)
  z.split.ll = do.call(LLFunction, params)
  z.split = z

  pairs = c(NA, NA)
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

	  ## Calculate likelihood of split
	  params = c(list(counts=counts, z=new.z, K=K), dot.args)
	  new.ll = do.call(LLFunction, params)

	  z.split.ll = c(z.split.ll, new.ll)
	  z.split = cbind(z.split, new.z)
   
	  pairs = rbind(pairs, c(i, j))
	}  
  }

  select = which.max(z.split.ll) 

  if(select == 1) {
    m = paste0(date(), " ... No additional splitting was performed.") 
  } else {
    m = paste0(date(), " ... Cluster ", pairs[select,1], " was reassigned and cluster ", pairs[select,2], " was split in two.")
  } 

  return(list(z=z.split[,select], message=m))
}



split.each.y = function(counts, y, L, LLFunction, y.prob, min=3, max.clusters.to.try=10, ...) { 

  dot.args = list(...)
  
  ## Identify clusters to split
  y.ta = table(factor(y, levels=1:L))
  y.to.split = which(y.ta >= min)
  y.non.empty = which(y.ta > 0)

  ## Loop through each split-able y and find best split
  clust.split = list()
  for(i in 1:L) {
    if(i %in% y.to.split) {    
      clustLabel = suppressMessages(celda_G(counts[y == i,], L=2, max.iter=3, split.on.iter=-1))
      clust.split = c(clust.split, list(clustLabel$y))
    } else {
      clust.split = c(clust.split, list(NA))
    }  
  }
  
  ## Find second best assignment give current assignments for each cell
  y.prob[cbind(1:nrow(y.prob), y)] = NA
  y.second = apply(y.prob, 1, which.max)
  ll.shuffle = rep(NA, L)
  for(i in y.non.empty) {
    ix = y == i
    new.y = y
    new.y[ix] = y.second[ix]
    params = c(list(counts=counts, y=new.y, L=L), dot.args)
    ll.shuffle[i] = do.call(LLFunction, params)
  }
  y.to.shuffle = head(order(ll.shuffle, decreasing = TRUE, na.last=NA), n = max.clusters.to.try)
  
  if(length(y.to.split) == 0 | length(y.to.shuffle) == 0) {
    m = paste0(date(), " ... Cluster sizes too small. No additional splitting was performed.") 
    return(list(y=y, message=m))
  }
 
  ## Set up variables for holding results with current iteration of y
  params = c(list(counts=counts, y=y, L=L), dot.args)
  y.split.ll = do.call(LLFunction, params)
  y.split = y

  pairs = c(NA, NA)
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
      params = c(list(counts=counts, y=new.y, L=L), dot.args)
      new.ll = do.call(LLFunction, params)

      y.split.ll = c(y.split.ll, new.ll)
      y.split = cbind(y.split, new.y)
      
      pairs = rbind(pairs, c(i, j))
    }
  }

  select = which.max(y.split.ll) 
  
  if(select == 1) {
    m = paste0(date(), " ... No additional splitting was performed.") 
  } else {
    m = paste0(date(), " ... Cluster ", pairs[select,1], " was reassigned and cluster ", pairs[select,2], " was split in two.")
  } 
 
  return(list(y=y.split[,select], message=m))
}





clusterByHC = function(d, min=3, method="ward.D") {
  label = cluster::pam(x = d, k=2)$clustering

  ## If PAM split is too small, perform secondary hclust procedure to split into roughly equal groups
  if(min(table(label)) < min) {
	d.hclust = stats::hclust(d, method=method)

	## Get maximum sample size of each subcluster
	d.hclust.size = sapply(1:length(d.hclust$height), function(i) max(table(stats::cutree(d.hclust, h=d.hclust$height[i]))))
  
	## Find the height of the dendrogram that best splits the samples in "roughly" half
	sample.size = round(length(d.hclust$order)/ 2)
	d.hclust.select = which.min(abs(d.hclust.size - sample.size))
	temp = stats::cutree(d.hclust, h=d.hclust$height[d.hclust.select])
	
	## Set half of the samples as the largest cluster and the other samples to the other cluster
	label.max.cluster = which.max(table(temp))
	label = ifelse(temp == label.max.cluster, 1, 2)
  } 
  
  return(label)
}
