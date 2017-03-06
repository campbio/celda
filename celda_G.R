calcLL = function(counts, z, k, alpha, beta, gamma) {
  n.phi = rowsum(counts, group=z, reorder=TRUE)
  nk = nrow(n.phi)
  nc = ncol(n.phi)
  
  a = nc*lgamma(nk*alpha)
  b = sum(lgamma(n.phi+alpha))
  c = -nc*nk*lgamma(alpha)
  d = -sum(lgamma(apply(n.phi + alpha, 2, sum)))
  
  phi.ll = a + b + c + d

  n.psi = rowSums(counts)
  n.psi.sum = as.numeric(rowsum(n.psi, z))
  
  ng.z = table(z)
  ng = length(n.psi)
  nk = length(n.psi.sum)
  
  a = sum(lgamma(ng.z * beta))
  b = sum(lgamma(n.psi + beta))
  c = -ng * lgamma(beta)
  d = -sum(lgamma(n.psi.sum + (ng.z*beta)))
  
  psi.ll = a + b + c + d
  
  a = lgamma(nk*gamma)
  b = sum(lgamma(ng.z+gamma))
  c = -nk*lgamma(gamma)
  d = -sum(lgamma(sum(ng.z + gamma)))
  
  eta.ll = a + b + c + d

  final = phi.ll + psi.ll + eta.ll
  return(final)
}

calcLLlite = function(ix, counts, z, k, alpha, beta, gamma) {
  
  ng.z.minus = table(z[-ix])
  eta.ll = log(ng.z.minus[z[ix]] + gamma)
  
  n.phi = rowsum(counts, group=z, reorder=TRUE)
  b = sum(lgamma(n.phi+alpha))
  d = -sum(lgamma(apply(n.phi + alpha, 2, sum)))
  phi.ll = b + d
  
  n.psi = rowSums(counts)
  n.psi.sum = as.numeric(rowsum(n.psi, z))
  ng.z = table(z)
  ng = length(n.psi)
  a = sum(lgamma(ng.z*beta))
  d = -sum(lgamma(n.psi.sum + (ng.z*beta)))
  psi.ll = a + d
  
  final = eta.ll + phi.ll + psi.ll
  return(final)
}


calcGibbsProb = function(ix, r, z, k, a, b, g) {
  
  final = rep(NA, k)
  for(j in 1:k) {
    z[ix] = j
    final[j] = calcLLlite(ix, counts=r, z=z, k=k, alpha=a, beta=b, gamma=g)
  }  
  
  return(final)
}





geneCluster = function(counts, k, a=1, b=1, g=1, max.iter=25, min.cell=5, seed=12345, best=TRUE, kick=TRUE, converge=1e-5) {
  set.seed(seed)
  
  cat(date(), "... Starting Gibbs sampling\n")
  
  co = counts

  z = sample(1:k, nrow(co), replace=TRUE)
  z.all = z
  ll = calcLL(counts=co, z=z, k=k, alpha=a, beta=b, gamma=g)
  
  z.probs = matrix(NA, nrow=nrow(co), ncol=k)
  
  iter = 1
  continue = TRUE
  while(iter <= max.iter & continue == TRUE) {
    
    ## Begin process of Gibbs sampling for each cell
    ix = sample(1:nrow(co))
    for(i in ix) {
      probs = calcGibbsProb(i, r=co, z=z, k=k, a=a, b=b, g=g)
      z[i] = sample.ll(probs)
      z.probs[i,] = probs
    }
    
    ## Save Z history
    z.all = cbind(z.all, z)
    
    ## Calculate complete likelihood
    temp.ll = calcLL(counts=co, z=z, k=k, alpha=a, beta=b, gamma=g)
    ll = c(ll, temp.ll)
    
    cat(date(), "... Completed iteration:", iter, "| logLik:", temp.ll, "\n")
    
    ## Normalize Z probabilties and test for convergence
    z.probs = exp(sweep(z.probs, 1, apply(z.probs, 1, max), "-"))
    z.probs = sweep(z.probs, 1, rowSums(z.probs), "/")
    f = function(v) sort(v, decreasing=TRUE)[2]
    z.probs.second = max(apply(z.probs, 1, f))
    z.ta = table(z)
#    if(z.probs.second < converge & (min(z.ta) >= min.cell | kick==FALSE)) {
    if(z.probs.second < converge) {    
      continue = FALSE
      cat("Maximum probability of a cell changing its state is ", z.probs.second, ". Exiting at iteration ", iter, ".", sep="")
    }
    
    iter = iter + 1    
  }
  
  
  if(best == TRUE) {
    ix = which.max(ll)
    z.final = z.all[,ix]
    ll.final = ll[ix]
  } else {
    z.final = z
    ll.final = tail(ll, n=1)
  }
  
  return(list(z=z.final, complete.z=z.all, completeLogLik=ll, finalLogLik=ll.final, z.probability=z.probs))
}









generateCells = function(C=100, N.Range=c(500,5000), G=1000, k=5, a=1, b=1, g=1, seed=12345) {
  require(gtools)
  
  set.seed(seed)
  eta = rdirichlet(1, rep(g, k))
  
  z = sample(1:k, size=G, prob=eta, replace=TRUE)
  if(length(table(z)) < k) {
    stop("Some states did not receive any genes after sampling. Try increasing G and/or setting gamma > 1.")
  }
  
  phi = matrix(0, nrow=G, ncol=k)
  for(i in 1:k) {
    ind = z == i
    phi[ind,i] = rdirichlet(1, rep(b, sum(ind)))
  }
  
  theta = rdirichlet(C, rep(a, k))
  
  ## Select number of transcripts per cell
  nN = sample(N.Range[1]:N.Range[2], size=C, replace=TRUE)
  
  ## Select transcript distribution for each cell
  cell.counts = matrix(0, nrow=G, ncol=C)
  for(i in 1:C) {
    cell.dist = rmultinom(1, size=nN[i], prob=theta[i,])
    for(j in 1:k) {
      cell.counts[,i] = cell.counts[,i] + rmultinom(1, size=cell.dist[j], prob=phi[,j])
    }
  }

  cc.rowsum0 = rowSums(cell.counts) > 0
  cell.counts = cell.counts[cc.rowsum0,]
  z = z[cc.rowsum0]
  
  return(list(z=z, counts=cell.counts, k=k, a=a, b=b, g=g, theta=theta, phi=phi, eta=eta, seed=seed))
}


sample.ll = function(ll.probs) {
  probs.sub = exp(ll.probs - max(ll.probs))
  probs.norm = probs.sub / sum(probs.sub)
  probs.select = sample(1:length(ll.probs), size=1, prob=probs.norm)
  return(probs.select)
}

