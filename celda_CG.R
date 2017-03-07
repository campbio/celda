cCG.calcLL = function(counts, s, z, y, K, L, a1, a2, a3, a4) {
  
  ## Calculate for "Theta" component
  m = table(z, s)
  ns = ncol(m)
  
  a = ns*lgamma(K*a1)
  b = sum(lgamma(m+a1))
  c = -ns*K*lgamma(a1)
  d = -sum(lgamma(apply(m + a1, 2, sum)))
  
  theta.ll = a + b + c + d
  
  
  ## Calculate for "Phi" component
  n.phi = rowsum(t(rowsum(counts, group=y, reorder=TRUE)), group=z, reorder=TRUE)
  
  a = K*lgamma(L*a2)
  b = sum(lgamma(n.phi+a2))
  c = -K*L*lgamma(a2)
  d = -sum(lgamma(apply(n.phi + a2, 1, sum)))
  
  phi.ll = a + b + c + d
  
  ## Calculate for "Psi" component
  n.psi = rowSums(counts)
  n.psi.sum = as.numeric(rowsum(n.psi, y))
  
  ny = table(y)
  ny.sum = sum(ny)
  ng = length(n.psi)
  nk = length(n.psi.sum)
  
  a = sum(lgamma(ny * a4))
  b = sum(lgamma(n.psi + a4))
  c = -ny.sum * lgamma(a4)
  d = -sum(lgamma(n.psi.sum + (ny*a4)))
  
  psi.ll = a + b + c + d
  
  
  ## Calculate for "Eta" side
  a = lgamma(L*a3)
  b = sum(lgamma(ny+a3))
  c = -L*lgamma(a3)
  d = -lgamma(sum(ny + a3))
  
  eta.ll = a + b + c + d
  
  final = theta.ll + phi.ll + psi.ll + eta.ll
  return(final)
}

cCG.calcLLliteZ = function(ix, counts, s, z, y, K, L, a1, a2, a3, a4) {
  
  z = factor(z, levels=1:K)
  y = factor(y, levels=1:L)
  
  ## Calculate for "Theta" component
  m = table(z[-ix], s[-ix])
  theta.ll = log(m[z[ix],s[ix]] + a1)
  
  ## Calculate for "Phi" component
  n.phi = rowsum(t(rowsum(counts, group=y, reorder=TRUE)), group=z, reorder=TRUE)
  
  b = sum(lgamma(n.phi+a2))
  d = -sum(lgamma(apply(n.phi + a2, 1, sum)))
  
  phi.ll = b + d
  
  final = theta.ll + phi.ll 
  return(final)
}

cCG.calcLLliteY = function(ix, counts, s, z, y, K, L, a1, a2, a3, a4) {
  
  z = factor(z, levels=1:K)
  y = factor(y, levels=1:L)
  
  ## Calculate for "Phi" component
  n.phi = rowsum(t(rowsum(counts, group=y, reorder=TRUE)), group=z, reorder=TRUE)
  
  b = sum(lgamma(n.phi+a2))
  d = -sum(lgamma(apply(n.phi + a2, 1, sum)))
  
  phi.ll = b + d
  
  ## Calculate for "Psi" component
  n.psi = rowSums(counts)
  n.psi.sum = as.numeric(rowsum(n.psi, y))
  
  ny = table(y)
  ny.sum = sum(ny)
  ng = length(n.psi)
  nk = length(n.psi.sum)
  
  a = sum(lgamma(ny * a4))
  d = -sum(lgamma(n.psi.sum + (ny*a4)))
  
  psi.ll = a + d
  
  
  ## Calculate for "Eta" side
  ng.y.minus = table(y[-ix])
  eta.ll = log(ng.y.minus[y[ix]] + a3)
  
  final = phi.ll + psi.ll + eta.ll
  return(final)
}

cCG.calcGibbsProbZ = function(ix, counts, s, z, y, K, L, a1, a2, a3, a4) {
  
  final = rep(NA, K)
  for(j in 1:K) {
    z[ix] = j
    final[j] = cCG.calcLLliteZ(ix, counts=counts, s=s, z=z, y=y, K=K, L=L, a1=a1, a2=a2, a3=a3, a4=a4)
  }  
  
  return(final)
}

cCG.calcGibbsProbY = function(ix, counts, s, z, y, K, L, a1, a2, a3, a4) {
  
  final = rep(NA, L)
  for(j in 1:L) {
    y[ix] = j
    final[j] = cCG.calcLLliteY(ix, counts=counts, s=s, z=z, y=y, K=K, L=L, a1=a1, a2=a2, a3=a3, a4=a4)
  }  
  
  return(final)
}

cCG.generateCells = function(S=10, C.Range=c(50,100), N.Range=c(500,5000), G=1000, K=3, L=10, a1=1, a2=1, a3=1, a4=1, seed=12345) {
  require(gtools)
  
  set.seed(seed)

  ## Number of cells per sample
  nC = sample(C.Range[1]:C.Range[2], size=S, replace=TRUE)
  nC.sum = sum(nC)
  cell.sample.label = rep(1:S, nC)
  
  ## Select number of transcripts per cell
  nN = sample(N.Range[1]:N.Range[2], size=length(cell.sample.label), replace=TRUE)
  
  ## Generate cell population distribution for each sample
  theta = t(rdirichlet(S, rep(a1, K)))

  ## Assign cells to cellular subpopulations
  z = unlist(lapply(1:S, function(i) sample(1:K, size=nC[i], prob=theta[,i], replace=TRUE)))

  ## Generate transcriptional state distribution for each cell subpopulation
  phi = rdirichlet(K, rep(a2, L))

  ## Assign genes to transcriptional states 
  eta = rdirichlet(1, rep(a3, L))
  y = sample(1:L, size=G, prob=eta, replace=TRUE)
  if(length(table(y)) < L) {
    stop("Some transcriptional states did not receive any genes after sampling. Try increasing G and/or setting gamma > 1.")
  }

  psi = matrix(0, nrow=G, ncol=L)
  for(i in 1:L) {
    ind = y == i
    psi[ind,i] = rdirichlet(1, rep(a4, sum(ind)))
  }
  
  ## Select transcript distribution for each cell
  cell.counts = matrix(0, nrow=G, ncol=nC.sum)
  for(i in 1:nC.sum) {
    #cell.dist = rmultinom(1, size=nN[i], prob=theta[,i])
    
    #transcriptional.state.dist = matrix(0, nrow=K, ncol=L)
    #for(j in 1:K) {
      #transcriptional.state.dist[j,] = rmultinom(1, size=cell.dist[j], prob=phi[i,])
    #}
    transcriptional.state.dist = as.numeric(rmultinom(1, size=nN[i], prob=phi[z[i],]))
    for(j in 1:L) {
      if(transcriptional.state.dist[j] > 0) {
        cell.counts[,i] = cell.counts[,i] + rmultinom(1, size=transcriptional.state.dist[j], prob=psi[,j])
      }
    }  
  }
  
  return(list(z=z, y=y, sample=cell.sample.label, counts=cell.counts, K=K, L=L, C.Range=C.Range, N.Range=N.Range, S=S, a1=a1, a2=a2, a3=a3, a4=a4, theta=theta, phi=phi, psi=psi, eta=eta, seed=seed))
}


celda_CG = function(counts, sample, K, L, a1=1, a2=1, a3=1, a4=1, max.iter=25, min.cell=5, seed=12345, best=TRUE, kick=TRUE, converge=1e-5) {
  set.seed(seed)
  
  cat(date(), "... Starting Gibbs sampling\n")
  
  co = counts
  s = sample
  
  ## Randomly select z and y
  z = sample(1:K, ncol(co), replace=TRUE)
  y = sample(1:L, nrow(co), replace=TRUE)
  z.all = z
  y.all = y
  z.stability = c(0)
  y.stability = c(0)
  
  ll = cCG.calcLL(counts=co, s=s, z=z, y=y, K=K, L=L, a1=a1, a2=a2, a3=a3, a4=a4)
  
  z.probs = matrix(NA, nrow=ncol(co), ncol=K)
  y.probs = matrix(NA, nrow=nrow(co), ncol=L)
  
  iter = 1
  continue = TRUE
  while(iter <= max.iter & continue == TRUE) {
    
    ## Determine if any clusters are below the minimum threshold and if a kick needs to be performed
    z.ta = table(factor(z, levels=1:K))
    if(min(z.ta) < min.cell & kick==TRUE) {
      
      all.k.to.kick = which(z.ta < min.cell)
      
      for(j in all.k.to.kick) { 
        all.k.to.test = which(z.ta > 2*min.cell)
        z = cCG.kick.z(co, s=s, z=z, y=y, K=K, L=L, k.to.kick=j, k.to.test=all.k.to.test, min.cell=min.cell, a1, a2, a3, a4)
        z.ta = table(factor(z, levels=1:K))
      }
    }
    
    ## Begin process of Gibbs sampling for each cell
    ix = sample(1:ncol(co))
    for(i in ix) {
      probs = cCG.calcGibbsProbZ(i, co, s=s, z=z, y=y, K=K, L=L, a1=a1, a2=a2, a3=a3, a4=a4)
      
      z[i] = sample.ll(probs)
      z.probs[i,] = probs
    }
    ## Begin process of Gibbs sampling for each gene
    ix = sample(1:nrow(co))
    for(i in ix) {
      if(sum(y == y[i]) > 1) {
        probs = cCG.calcGibbsProbY(i, co, s=s, z=z, y=y, K=K, L=L, a1=a1, a2=a2, a3=a3, a4=a4)
        y[i] = sample.ll(probs)
        y.probs[i,] = probs
      }  
    }
    
    ## Save Z history
    z.all = cbind(z.all, z)
    y.all = cbind(y.all, y)
    
    ## Calculate complete likelihood
    temp.ll = cCG.calcLL(counts=co, s=s, z=z, y=y, K=K, L=L, a1=a1, a2=a2, a3=a3, a4=a4)
    ll = c(ll, temp.ll)
    
    cat(date(), "... Completed iteration:", iter, "| logLik:", temp.ll, "\n")
    
    ## Normalize Z and Y probabilties and test for convergence
    z.probs = exp(sweep(z.probs, 1, apply(z.probs, 1, max), "-"))
    z.probs = sweep(z.probs, 1, rowSums(z.probs), "/")
    y.probs = exp(sweep(y.probs, 1, apply(y.probs, 1, max), "-"))
    y.probs = sweep(y.probs, 1, rowSums(y.probs), "/")
    
    ## Calculate stability
    z.stab = sum(1-apply(z.probs, 1, entropy) / log(K)) / ncol(co)
    z.stability = c(z.stability, stability(z.probs))
    y.stab = sum(1-apply(y.probs, 1, entropy) / log(L)) / nrow(co)
    y.stability = c(y.stability, stability(y.probs))
    
    iter = iter + 1    
  }
  
  ## Identify which model is the best overall in terms of maximum likelihood
  if(best == TRUE) {
    ix = which.max(ll)
    z.final = z.all[,ix]
    y.final = y.all[,ix]
    ll.final = ll[ix]
  } else {
    z.final = z
    y.final = y
    ll.final = tail(ll, n=1)
  }
  
  return(list(z=z.final, y=y.final, complete.z=z.all, complete.y=y.all, z.stability=z.stability, y.stability=y.stability, completeLogLik=ll, finalLogLik=ll.final, z.prob=z.probs, y.prob=y.probs))
}


cCG.kick.z = function(counts, s, z, y, K, L, k.to.kick, k.to.test, min.cell=5, a1, a2, a3, a4) {
  require(cluster)
  cat(date(), "... Cluster", k.to.kick, "has fewer than", min.cell, "cells. Performing kick by ")
  
  counts.norm = sweep(counts, 2, colSums(counts), "/")
  z.kick = matrix(z, ncol=length(k.to.test), nrow=length(z))
  
  ## Randomly assign clusters to cells with cluster to kick
  z.k.to.kick = sample(1:K, size=sum(z == k.to.kick), replace=TRUE)
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
    k.kick.ll[i] = cCG.calcLL(counts=counts, s=s, z=z.kick[,i], y=y, K=K, L=L, a1=a1, a2=a2, a3=a3, a4=a4)
  }
  
  k.to.test.select = sample.ll(k.kick.ll)
  
  cat("splitting Cluster", k.to.test[k.to.test.select], "\n")
  return(z.kick[,k.to.test.select])
}


