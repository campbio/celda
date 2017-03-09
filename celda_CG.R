cCG.calcLL = function(counts, s, m.theta, n.phi, n.psi, z, y, K, L, alpha, beta, gamma, delta) {
  
  ## Calculate for "Theta" component
  #m = table(z, s)
  ns = ncol(m.theta)
  
  a = ns*lgamma(K*alpha)
  b = sum(lgamma(m.theta+alpha))
  c = -ns*K*lgamma(alpha)
  d = -sum(lgamma(apply(m.theta + alpha, 2, sum)))
  
  theta.ll = a + b + c + d
  
  
  ## Calculate for "Phi" component
  #n.phi = rowsum(t(rowsum(counts, group=y, reorder=TRUE)), group=z, reorder=TRUE)
  
  a = K*lgamma(L*beta)
  b = sum(lgamma(n.phi+beta))
  c = -K*L*lgamma(beta)
  d = -sum(lgamma(apply(n.phi + beta, 1, sum)))
  
  phi.ll = a + b + c + d
  
  ## Calculate for "Psi" component
  #n.psi = rowSums(counts)
  n.psi.sum = as.numeric(rowsum(n.psi, y))
  
  ny = table(y)
  ny.sum = sum(ny)
  ng = length(n.psi)
  nk = length(n.psi.sum)
  
  a = sum(lgamma(ny * delta))
  b = sum(lgamma(n.psi + delta))
  c = -ny.sum * lgamma(delta)
  d = -sum(lgamma(n.psi.sum + (ny*delta)))
  
  psi.ll = a + b + c + d
  
  
  ## Calculate for "Eta" side
  a = lgamma(L*gamma)
  b = sum(lgamma(ny+gamma))
  c = -L*lgamma(gamma)
  d = -lgamma(sum(ny + gamma))
  
  eta.ll = a + b + c + d
  
  final = theta.ll + phi.ll + psi.ll + eta.ll
  return(final)
}

cCG.calcGibbsProbZ = function(m.theta, n.phi, alpha, beta) {
  
#  z = factor(z, levels=1:K)
#  y = factor(y, levels=1:L)
  
  ## Calculate for "Theta" component
#  m = table(z[-ix], s[-ix])
  theta.ll = log(m.theta + alpha)
  
  ## Calculate for "Phi" component
#  n.phi = rowsum(t(rowsum(counts, group=y, reorder=TRUE)), group=z, reorder=TRUE)
  
  b = sum(lgamma(n.phi+beta))
  d = -sum(lgamma(apply(n.phi + beta, 1, sum)))
  
  phi.ll = b + d
  
  final = theta.ll + phi.ll 
  return(final)
}

cCG.calcGibbsProbY = function(y, n.phi, n.psi, q.eta, beta, gamma, delta) {
  
  ## Calculate for "Phi" component
  b = sum(lgamma(n.phi+beta))
  d = -sum(lgamma(apply(n.phi + beta, 1, sum)))
  
  phi.ll = b + d
  
  ## Calculate for "Psi" component
  n.psi.y = as.numeric(rowsum(n.psi, y))
  
  ny = table(y)
  ny.sum = sum(ny)
  ng = length(n.psi)
  
  a = sum(lgamma(ny * delta))
  d = -sum(lgamma(n.psi.y + (ny*delta)))
  
  psi.ll = a + d
  
  ## Calculate for "Eta" side
  eta.ll = log(q.eta + gamma)
  
  final = phi.ll + psi.ll + eta.ll
  return(final)
}


cCG.generateCells = function(S=10, C.Range=c(50,100), N.Range=c(500,5000), G=1000, K=3, L=10, alpha=1, beta=1, gamma=1, delta=1, seed=12345) {
  require(gtools)
  
  set.seed(seed)

  ## Number of cells per sample
  nC = sample(C.Range[1]:C.Range[2], size=S, replace=TRUE)
  nC.sum = sum(nC)
  cell.sample.label = rep(1:S, nC)
  
  ## Select number of transcripts per cell
  nN = sample(N.Range[1]:N.Range[2], size=length(cell.sample.label), replace=TRUE)
  
  ## Generate cell population distribution for each sample
  theta = t(rdirichlet(S, rep(alpha, K)))

  ## Assign cells to cellular subpopulations
  z = unlist(lapply(1:S, function(i) sample(1:K, size=nC[i], prob=theta[,i], replace=TRUE)))

  ## Generate transcriptional state distribution for each cell subpopulation
  phi = rdirichlet(K, rep(beta, L))

  ## Assign genes to transcriptional states 
  eta = rdirichlet(1, rep(gamma, L))
  y = sample(1:L, size=G, prob=eta, replace=TRUE)
  if(length(table(y)) < L) {
    stop("Some transcriptional states did not receive any genes after sampling. Try increasing G and/or setting gamma > 1.")
  }

  psi = matrix(0, nrow=G, ncol=L)
  for(i in 1:L) {
    ind = y == i
    psi[ind,i] = rdirichlet(1, rep(delta, sum(ind)))
  }
  
  ## Select transcript distribution for each cell
  cell.counts = matrix(0, nrow=G, ncol=nC.sum)
  for(i in 1:nC.sum) {
    transcriptional.state.dist = as.numeric(rmultinom(1, size=nN[i], prob=phi[z[i],]))
    for(j in 1:L) {
      if(transcriptional.state.dist[j] > 0) {
        cell.counts[,i] = cell.counts[,i] + rmultinom(1, size=transcriptional.state.dist[j], prob=psi[,j])
      }
    }  
  }
  
  return(list(z=z, y=y, sample=cell.sample.label, counts=cell.counts, K=K, L=L, C.Range=C.Range, N.Range=N.Range, S=S, alpha=alpha, beta=beta, gamma=gamma, delta=delta, theta=theta, phi=phi, psi=psi, eta=eta, seed=seed))
}


celda_CG = function(counts, sample, K, L, alpha=1, beta=1, gamma=1, delta=1, max.iter=25, min.cell=5, seed=12345, best=TRUE, kick=TRUE, converge=1e-5) {
  set.seed(seed)
  
  cat(date(), "... Starting Gibbs sampling\n")
  
  co = counts
  s = sample
  
  ## Randomly select z and y
  z = sample(1:K, ncol(co), replace=TRUE)
  y = sample(1:L, nrow(co), replace=TRUE)
  z.all = z
  y.all = y
  z.stability = c(NA)
  y.stability = c(NA)
  
  ## Calculate counts one time up front
  m.theta = table(factor(z, levels=1:K), s)
  n.phi.y = rowsum(co, group=y, reorder=TRUE)
  n.phi.y.z = rowsum(t(n.phi.y), group=z, reorder=TRUE)
  n.psi = rowSums(co)
  q.eta = table(y)
  
  ll = cCG.calcLL(counts=co, s=s, m.theta=m.theta, n.phi=n.phi.y.z, n.psi=n.psi, z=z, y=y, K=K, L=L, alpha=alpha, beta=beta, gamma=gamma, delta=delta)
  
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
        z = cCG.kick.z(co, s=s, z=z, y=y, K=K, L=L, k.to.kick=j, k.to.test=all.k.to.test, min.cell=min.cell, a, beta, gamma, delta)
        z.ta = table(factor(z, levels=1:K))
      }
    }
    
    ## Begin process of Gibbs sampling for each cell
    n.phi.y = rowsum(co, group=y, reorder=TRUE)
    ix = sample(1:ncol(co))
    for(i in ix) {
      # If more than one cell belongs to the current state, then sampling will proceed, otherwise this cell will be skipped
      if(sum(z == z[i]) > 1) {
        
        ## Subtract current cell counts from matrices
        m.theta[z[i],s[i]] = m.theta[z[i],s[i]] - 1
        n.phi.y.z[z[i],] = n.phi.y.z[z[i],] - n.phi.y[,i]
        
        ## Calculate probabilities for each state
        probs = rep(NA, K)
        for(j in 1:K) {
          temp.n.phi.y.z = n.phi.y.z
          temp.n.phi.y.z[j,] = temp.n.phi.y.z[j,] + n.phi.y[,i]
          probs[j] = cCG.calcGibbsProbZ(m.theta=m.theta[j,s[i]], n.phi=temp.n.phi.y.z, alpha=alpha, beta=beta)
        }  
        
        ## Identify next state and add back counts
        z[i] = sample.ll(probs)
        m.theta[z[i],s[i]] = m.theta[z[i],s[i]] + 1
        n.phi.y.z[z[i],] = n.phi.y.z[z[i],] + n.phi.y[,i]
        
      } else {
        probs = rep(0, K)
        probs[z[i]] = 1
      }
      z.probs[i,] = probs
    }
    
    ## Begin process of Gibbs sampling for each gene
    n.phi.z = rowsum(t(co), group=z, reorder=TRUE)
    ix = sample(1:nrow(co))
    for(i in ix) {
      if(sum(y == y[i]) > 1) {
        
        probs = rep(NA, L)
        q.eta[y[i]] = q.eta[y[i]] - 1
        n.phi.y.z[,y[i]] = n.phi.y.z[,y[i]] - n.phi.z[,i]
        
        for(j in 1:L) {
          temp.y = y
          temp.y[i] = j
          temp.n.phi.y.z = n.phi.y.z
          temp.n.phi.y.z[,j] = temp.n.phi.y.z[,j] + n.phi.z[,i]
          probs[j] = cCG.calcGibbsProbY(y=temp.y, n.phi=temp.n.phi.y.z, n.psi=n.psi, q.eta=q.eta[j], beta=beta, gamma=gamma, delta=delta)
        }  
        
        y[i] = sample.ll(probs)
        q.eta[y[i]] = q.eta[y[i]] + 1
        n.phi.y.z[,y[i]] = n.phi.y.z[,y[i]] + n.phi.z[,i]
        
      } else {
        probs = rep(0, L)
        probs[y[i]] = 1
      }
      y.probs[i,] = probs
    }
    
    ## Save Z history
    z.all = cbind(z.all, z)
    y.all = cbind(y.all, y)

    ## Normalize Z and Y marginal probabilties and calculate stability
    z.probs = normalizeLogProbs(z.probs)
    y.probs = normalizeLogProbs(y.probs)
    z.stability = c(z.stability, stability(z.probs))
    y.stability = c(y.stability, stability(y.probs))

    ## Calculate complete likelihood
    temp.ll = cCG.calcLL(counts=co, s=s, m.theta=m.theta, n.phi=n.phi.y.z, n.psi=n.psi, z=z, y=y, K=K, L=L, alpha=alpha, beta=beta, gamma=gamma, delta=delta)
    if((best == TRUE & all(temp.ll > ll)) | iter == 1) {
      z.probs.final = z.probs
      y.probs.final = y.probs
    }
    ll = c(ll, temp.ll)
    
    cat(date(), "... Completed iteration:", iter, "| logLik:", temp.ll, "\n")
    iter = iter + 1    
  }
  
  ## Identify which model is the best overall in terms of maximum likelihood
  if(best == TRUE) {
    ix = which.max(ll)
    z.final = z.all[,ix]
    y.final = y.all[,ix]
    ll.final = ll[ix]
    z.stability.final = z.stability[ix]
    y.stability.final = y.stability[ix]
  } else {
    z.final = z
    y.final = y
    ll.final = tail(ll, n=1)
    z.stability.final = tail(z.stability, n=1)
    y.stability.final = tail(y.stability, n=1)
    z.probs.final = z.probs
    y.probs.final = y.probs
  }
  
  return(list(z=z.final, y=y.final, complete.z=z.all, complete.y=y.all, z.stability=z.stability.final, y.stability=y.stability.final, complete.z.stability=z.stability, complete.y.stability=y.stability, completeLogLik=ll, finalLogLik=ll.final, z.prob=z.probs.final, y.prob=y.probs.final))
}


cCG.kick.z = function(counts, s, z, y, K, L, k.to.kick, k.to.test, min.cell=5, alpha, beta, gamma, delta) {
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
    k.kick.ll[i] = cCG.calcLL(counts=counts, s=s, z=z.kick[,i], y=y, K=K, L=L, alpha=alpha, beta=beta, gamma=gamma, delta=delta)
  }
  
  k.to.test.select = sample.ll(k.kick.ll)
  
  cat("splitting Cluster", k.to.test[k.to.test.select], "\n")
  return(z.kick[,k.to.test.select])
}


