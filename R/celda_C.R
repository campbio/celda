calcGibbsProb = function(ix, r, s, z, k, a, b) {
  phi <- c()
  for(j in 1:k) {
  z[ix] <- j
  n <- t(rowsum(t(r), group=z, reorder=TRUE))
  phi <- c(phi, ll.phi.gibbs(n, b))    
  }  
  
  z <- factor(z, levels=1:k)
  m <- table(z[-ix], s[-ix])  
  s.ix <- colnames(m) == s[ix]
  
  final <- log(m[,s.ix] + a) + phi
  return(final)
}


generateCells = function(S=10, C.Range=c(10, 100), N.Range=c(100,5000), 
                         G=5000, k=5, a=1, b=0.1) {
  
  phi <- gtools::rdirichlet(k, rep(b, G))
  theta <- gtools::rdirichlet(S, rep(a, k))
  
  ## Select the number of cells per sample
  nC <- sample(C.Range[1]:C.Range[2], size=S, replace=TRUE)  
  cell.sample <- rep(1:S, nC)
  
  ## Select state of the cells  
  cell.state <- unlist(lapply(1:S, function(i) sample(1:k, size=nC[i], prob=theta[i,], replace=TRUE)))
  
  ## Select number of transcripts per cell
  nN <- sample(N.Range[1]:N.Range[2], size=length(cell.sample), replace=TRUE)
  
  ## Select transcript distribution for each cell
  cell.counts <- sapply(1:length(cell.sample), function(i) rmultinom(1, size=nN[i], prob=phi[cell.state[i],]))
  
  return(list(z=cell.state, counts=cell.counts, sample=cell.sample, k=k, a=a, b=b))
}



celda_C = function(counts, sample, k, a=1, b=0.1, max.iter=25, min.cell=5, 
                   seed=12345, best=TRUE, kick=TRUE, converge=1e-5) {
  
  set.seed(seed)
  require(entropy)
  cat(date(), "... Starting Gibbs sampling\n")
  
  co = counts
  s = sample
  
  z = sample(1:k, ncol(co), replace=TRUE)
  z.all = z
  ll = calcLL(counts=co, s=s, z=z, k=k, alpha=a, beta=b)
  
  z.probs = matrix(NA, nrow=ncol(co), ncol=k)
    
  iter = 1
  continue = TRUE
  while(iter <= max.iter & continue == TRUE) {
    
    ## Determine if any clusters are below the minimum threshold 
    ## and if a kick needs to be performed
    z.ta = table(factor(z, levels=1:k))
    if(min(z.ta) < min.cell & kick==TRUE) {

      all.k.to.kick = which(z.ta < min.cell)
      
      for(j in all.k.to.kick) { 
        all.k.to.test = which(z.ta > 2*min.cell)
        z = kick.z(co, s=s, z=z, k=k, k.to.kick=j, k.to.test=all.k.to.test, 
                   min.cell=min.cell, a=a, b=b)
        z.ta = table(factor(z, levels=1:k))
      }
      
    }
    
    ## Begin process of Gibbs sampling for each cell
    ix = sample(1:ncol(co))
    for(i in ix) {
      probs = calcGibbsProb(i, r=co, s=s, z=z, k=k, a=a, b=b)
      
      z[i] = sample.ll(probs)
      z.probs[i,] = probs
    }

    ## Save Z history
    z.all = cbind(z.all, z)
    
    ## Calculate complete likelihood
    temp.ll = calcLL(counts=co, s=s, z=z, k=k, alpha=a, beta=b)
    ll = c(ll, temp.ll)

    cat(date(), "... Completed iteration:", iter, "| logLik:", temp.ll, "\n")
    
    ## Normalize Z probabilties and test for convergence
    z.probs = exp(sweep(z.probs, 1, apply(z.probs, 1, max), "-"))
    z.probs = sweep(z.probs, 1, rowSums(z.probs), "/")
    f = function(v) sort(v, decreasing=TRUE)[2]
    z.probs.second = max(apply(z.probs, 1, f))
    z.ta = table(z)
    if (z.probs.second < converge & (min(z.ta) >= min.cell | kick==FALSE)) {
      continue = FALSE
      cat("Maximum probability of a cell changing its state is ", z.probs.second, ". Exiting at iteration ", iter, ".", sep="")
    }
    
    iter = iter + 1    
  }
  
  
  if (best == TRUE) {
    ix = which.max(ll)
    z.final = z.all[,ix]
    ll.final = ll[ix]
  } else {
    z.final = z
    ll.final = tail(ll, n=1)
  }
  
  return(list(z=z.final, complete.z=z.all, completeLogLik=ll, 
              finalLogLik=ll.final, z.probability=z.probs))
}



ll.phi.gibbs = function(n, beta) {
  ng = nrow(n)
  nk = ncol(n)
  
  b = sum(lgamma(n+beta))
  d = -sum(lgamma(colSums(n + beta)))
  
  ll = b + d
  return(ll)
}


calcLL = function(counts, s, z, k, alpha, beta) {
  m = table(z, s)
  nk = nrow(m)
  ns = ncol(m)
  
  a = ns*lgamma(nk*alpha)
  b = sum(lgamma(m+alpha))
  c = -ns*nk*lgamma(alpha)
  d = -sum(lgamma(apply(m + alpha, 2, sum)))
  
  theta.ll = a + b + c + d
 
 
  n = sapply(1:k, function(i) apply(counts[,z == i,drop=FALSE], 1, sum))
  ng = nrow(n)
  nk = ncol(n)
  
  a = nk*lgamma(ng*beta)
  b = sum(lgamma(n+beta))
  c = -nk*ng*lgamma(beta)
  d = -sum(lgamma(apply(n + beta, 2, sum)))
  
  phi.ll = a + b + c + d

  final = theta.ll + phi.ll
  return(final)
}




cosineDist <- function(x){
  x = t(x)
  y = as.dist(1 - x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))) 
  return(y)
}

kick.z = function(counts, s, z, k, k.to.kick, k.to.test, min.cell=5, a, b) {
  require(cluster)
  cat(date(), "... Cluster", k.to.kick, "has fewer than", min.cell, "cells. Performing kick by ")
  
  counts.norm = sweep(counts, 2, colSums(counts), "/")
  z.kick = matrix(z, ncol=length(k.to.test), nrow=length(z))
  
  ## Randomly assign clusters to cells with cluster to kick
  z.k.to.kick = sample(1:k, size=sum(z == k.to.kick), replace=TRUE)
  z.kick[z==k.to.kick,] = z.k.to.kick
  
  ## Loop through each cluster, split, and determine logLik
  k.kick.ll = rep(NA, length(k.to.test))
  for(i in 1:length(k.to.test)) {
    k.dist = cosineDist(counts.norm[,z==k.to.test[i]])/2
    k.pam = pam(x=k.dist, k=2)$clustering
    
    ## If PAM split is too small, perform secondary hclust procedure to split into equal groups
    if(min(table(k.pam)) < min.cell) {
      k.hc = hclust(k.dist, method="ward.D")
      
      ## Get maximum sample size of each subcluster
      k.hc.size = sapply(1:length(k.hc$height), function(i) max(table(cutree(k.hc, h=k.hc$height[i]))))
      
      ## Find the height of the dendrogram that best splits the samples in half
      sample.size = round(length(k.hc$order)/ 2)
      k.hc.select = which.min(abs(k.hc.size - sample.size))
      k.hc.cut = cutree(k.hc, h=k.hc$height[k.hc.select])
      k.hc.cluster = which.max(table(k.hc.cut))
      
      k.hc.final = ifelse(k.hc.cut == k.hc.cluster, k.to.test[i], k.to.kick)
      
      ix = (z == k.to.test[i])
      z.kick[ix,i] = k.hc.final
      
    } else {
      
      k.pam.final = ifelse(k.pam == 1, k.to.test[i], k.to.kick)
      ix = (z == k.to.test[i])
      z.kick[ix,i] = k.pam.final
      
    }
    k.kick.ll[i] = calcLL(counts=counts, s=s, z=z.kick[,i], k=k, alpha=a, beta=b)
  }

  k.to.test.select = sample.ll(k.kick.ll)
  
  cat("splitting Cluster", k.to.test[k.to.test.select], "\n")
  return(z.kick[,k.to.test.select])
}


sample.ll = function(ll.probs) {
  probs.sub = exp(ll.probs - max(ll.probs))
  probs.norm = probs.sub / sum(probs.sub)
  probs.select = sample(1:length(ll.probs), size=1, prob=probs.norm)
  return(probs.select)
}

