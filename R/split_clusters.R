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
    
    k.dist = cosineDist(counts.norm[,z==k.to.test[i]])
    k.pam = cluster::pam(x=k.dist, k=2)$clustering
    
    ## If PAM split is too small, perform secondary hclust procedure to split into roughly equal groups
    if(min(table(k.pam)) < min.cell) {
      k.hc = hclust(k.dist, method="ward.D")
      
      ## Get maximum sample size of each subcluster
      k.hc.size = sapply(1:length(k.hc$height), function(i) max(table(cutree(k.hc, h=k.hc$height[i]))))
      
      ## Find the height of the dendrogram that best splits the samples in "roughly" half
      sample.size = round(length(k.hc$order)/ 2)
      k.hc.select = which.min(abs(k.hc.size - sample.size))
      k.hc.cut = cutree(k.hc, h=k.hc$height[k.hc.select])
      k.hc.cluster = which.max(table(k.hc.cut))
      
      k.hc.final = ifelse(k.hc.cut == k.hc.cluster, k.to.test[i], empty.K)
      
      ix = (z == k.to.test[i])
      z.split[ix,i] = k.hc.final
      
    } else {
      
      k.pam.final = ifelse(k.pam == 1, k.to.test[i], empty.K)
      ix = (z == k.to.test[i])
      z.split[ix,i] = k.pam.final
      
    }
    params = c(list(counts=counts, z=z.split[,i], K=K), list(...))
    k.split.ll[i] = do.call(LLFunction, params)
  }

  k.to.test.select = sample.ll(k.split.ll)
  
  message(date(), "... Cluster", empty.K, "had", z.ta[empty.K], "cells. Splitting Cluster", k.to.test[k.to.test.select], "into two clusters.\n")
  return(z.split[,k.to.test.select])
}
