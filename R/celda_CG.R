# -----------------------------------
# Variable description
# -----------------------------------
# C = Cell
# S or s = Sample
# G = Gene
# TS = Transcriptional State
# CP = Cell population
# n = counts of transcripts
# m = counts of cells
# K = Total number of cell populations
# L = Total number of transcriptional states
# nM = Number of cells
# nG = Number of genes
# nS = Number of samples

# -----------------------------------
# Count matrices descriptions
# -----------------------------------

# All n.* variables contain counts of transcripts
# n.CP.by.TS = Number of counts in each Cellular Population per Transcriptional State
# n.TS.by.C = Number of counts in each Transcriptional State per Cell 
# n.CP.by.G = Number of counts in each Cellular Population per Gene
# n.by.G = Number of counts per gene (i.e. rowSums)
# n.by.TS = Number of counts per Transcriptional State

## All m.* variables contain counts of cells
# m.CP.by.S = Number of cells in each Cellular Population per Sample

# nG.by.TS = Number of genes in each Transcriptional State


#' @export
cCG.calcLLFromVariables = function(counts, s, z, y, K, L, alpha, beta, gamma, delta) {
  
  ## Calculate for "Theta" component
  m = table(z, s)
  ns = ncol(m)
  
  a = ns * lgamma(K * alpha)
  b = sum(lgamma(m + alpha))
  c = -ns * K * lgamma(alpha)
  d = -sum(lgamma(colSums(m + alpha)))
  
  theta.ll = a + b + c + d
  
  
  ## Calculate for "Phi" component
  n.CP.by.TS = rowsum(t(rowsum(counts, group=y, reorder=TRUE)), group=z, reorder=TRUE)
  
  a = K * lgamma(L * beta)
  b = sum(lgamma(n.CP.by.TS + beta))
  c = -K * L * lgamma(beta)
  d = -sum(lgamma(rowSums(n.CP.by.TS + beta)))
  
  phi.ll = a + b + c + d
  
  ## Calculate for "Psi" component
  n.by.G = rowSums(counts)
  n.by.TS = as.numeric(rowsum(n.by.G, y))
  
  nG.by.TS = table(y)
  nG = length(y)
  
  a = sum(lgamma(nG.by.TS * delta))
  b = sum(lgamma(n.by.G + delta))
  c = -nG * lgamma(delta)
  d = -sum(lgamma(n.by.TS + (nG.by.TS * delta)))
  
  psi.ll = a + b + c + d
  
  
  ## Calculate for "Eta" side
  a = lgamma(L * gamma)
  b = sum(lgamma(nG.by.TS + gamma))
  c = -L * lgamma(gamma)
  d = -lgamma(sum(nG.by.TS + gamma))
  
  eta.ll = a + b + c + d
  
  final = theta.ll + phi.ll + psi.ll + eta.ll
  return(final)
}

cCG.calcLL = function(K, L, m.CP.by.S, n.CP.by.TS, n.by.G, n.by.TS, nG.by.TS, nS, nG, alpha, beta, gamma, delta) {
  
  ## Calculate for "Theta" component
  a = nS * lgamma(K*alpha)
  b = sum(lgamma(m.CP.by.S+alpha))
  c = -nS * K *lgamma(alpha)
  d = -sum(lgamma(colSums(m.CP.by.S + alpha)))
  
  theta.ll = a + b + c + d
  
  
  ## Calculate for "Phi" component
  a = K * lgamma(L * beta)
  b = sum(lgamma(n.CP.by.TS + beta))
  c = -K * L * lgamma(beta)
  d = -sum(lgamma(rowSums(n.CP.by.TS + beta)))
  
  phi.ll = a + b + c + d
  
  ## Calculate for "Psi" component
  a = sum(lgamma(nG.by.TS * delta))
  b = sum(lgamma(n.by.G + delta))
  c = -nG * lgamma(delta)
  d = -sum(lgamma(n.by.TS + (nG.by.TS * delta)))
  
  psi.ll = a + b + c + d
  
  
  ## Calculate for "Eta" side
  a = lgamma(L*gamma)
  b = sum(lgamma(nG.by.TS + gamma))
  c = -L*lgamma(gamma)
  d = -lgamma(sum(nG.by.TS + gamma))
  
  eta.ll = a + b + c + d
  
  final = theta.ll + phi.ll + psi.ll + eta.ll
  return(final)
}

cCG.calcGibbsProbZ = function(m.CP.by.S, n.CP.by.TS, alpha, beta) {

  ## Calculate for "Theta" component
  theta.ll = log(m.CP.by.S + alpha)
  
  ## Calculate for "Phi" component

  b = sum(lgamma(n.CP.by.TS+beta))
  d = -sum(lgamma(rowSums(n.CP.by.TS + beta)))
  
  phi.ll = b + d
  
  final = theta.ll + phi.ll 
  return(final)
}

cCG.calcGibbsProbY = function(n.CP.by.TS, n.by.TS, nG.by.TS, nG.in.Y, beta, gamma, delta) {
  
  ## Calculate for "Phi" component
  b = sum(lgamma(n.CP.by.TS + beta))
  d = -sum(lgamma(rowSums(n.CP.by.TS + beta)))
  
  phi.ll = b + d
  
  ## Calculate for "Psi" component
  a = sum(lgamma(nG.by.TS * delta))
  d = -sum(lgamma(n.by.TS + (nG.by.TS * delta)))
  
  psi.ll = a + d
  
  ## Calculate for "Eta" side
  eta.ll = log(nG.in.Y + gamma)
  
  final = phi.ll + psi.ll + eta.ll
  return(final)
}


#' @export
simulateCells.celdaCG = function(S=10, C.Range=c(50,100), N.Range=c(500,5000), G=1000, K=3, L=10, alpha=1, beta=1, gamma=1, delta=1, seed=12345) {
  
  set.seed(seed)

  ## Number of cells per sample
  nC = sample(C.Range[1]:C.Range[2], size=S, replace=TRUE)
  nC.sum = sum(nC)
  cell.sample.label = rep(1:S, nC)
  
  ## Select number of transcripts per cell
  nN = sample(N.Range[1]:N.Range[2], size=length(cell.sample.label), replace=TRUE)
  
  ## Generate cell population distribution for each sample
  theta = t(gtools::rdirichlet(S, rep(alpha, K)))

  ## Assign cells to cellular subpopulations
  z = unlist(lapply(1:S, function(i) sample(1:K, size=nC[i], prob=theta[,i], replace=TRUE)))

  ## Generate transcriptional state distribution for each cell subpopulation
  phi = gtools::rdirichlet(K, rep(beta, L))

  ## Assign genes to transcriptional states 
  eta = gtools::rdirichlet(1, rep(gamma, L))
  y = sample(1:L, size=G, prob=eta, replace=TRUE)
  if(length(table(y)) < L) {
    stop("Some transcriptional states did not receive any genes after sampling. Try increasing G and/or setting gamma > 1.")
  }

  psi = matrix(0, nrow=G, ncol=L)
  for(i in 1:L) {
    ind = y == i
    psi[ind,i] = gtools::rdirichlet(1, rep(delta, sum(ind)))
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


#' @export
celda_CG = function(counts, sample.label, K, L, alpha=1, beta=1, gamma=1, delta=1, max.iter=25, min.cell=5, seed=12345, best=TRUE, kick=TRUE, converge=1e-5) {
  set.seed(seed)
  
  message(date(), " ... Starting Gibbs sampling")
  
  if(is.factor(sample.label)) {
    s = as.numeric(sample.label)
  }
  else {
    s = as.numeric(as.factor(sample.label))
  }  
  
  ## Randomly select z and y
  z = sample(1:K, ncol(counts), replace=TRUE)
  y = sample(1:L, nrow(counts), replace=TRUE)
  z.all = z
  y.all = y
  z.stability = c(NA)
  y.stability = c(NA)
  z.probs = matrix(NA, nrow=ncol(counts), ncol=K)
  y.probs = matrix(NA, nrow=nrow(counts), ncol=L)
  
  
  ## Calculate counts one time up front
  m.CP.by.S = table(factor(z, levels=1:K), s)
  n.TS.by.C = rowsum(counts, group=y, reorder=TRUE)
  n.CP.by.TS = rowsum(t(n.TS.by.C), group=z, reorder=TRUE)
  n.by.G = rowSums(counts)
  n.by.TS = as.numeric(rowsum(n.by.G, y))
  nG.by.TS = table(y)
  nS = length(unique(s))
  nG = nrow(counts)
  nM = ncol(counts)
  
  ll = cCG.calcLL(K=K, L=L, m.CP.by.S=m.CP.by.S, n.CP.by.TS=n.CP.by.TS, n.by.G=n.by.G, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, nS=nS, nG=nG, alpha=alpha, beta=beta, gamma=gamma, delta=delta)
  
  iter = 1
  continue = TRUE
  while(iter <= max.iter & continue == TRUE) {
    
    ## Begin process of Gibbs sampling for each cell
    n.TS.by.C = rowsum(counts, group=y, reorder=TRUE)
    ix = sample(1:ncol(counts))
    for(i in ix) {

      ## Subtract current cell counts from matrices
      m.CP.by.S[z[i],s[i]] = m.CP.by.S[z[i],s[i]] - 1
      n.CP.by.TS[z[i],] = n.CP.by.TS[z[i],] - n.TS.by.C[,i]
      
      ## Calculate probabilities for each state
      probs = rep(NA, K)
      for(j in 1:K) {
        temp.n.CP.by.TS = n.CP.by.TS
        temp.n.CP.by.TS[j,] = temp.n.CP.by.TS[j,] + n.TS.by.C[,i]
        probs[j] = cCG.calcGibbsProbZ(m.CP.by.S=m.CP.by.S[j,s[i]], n.CP.by.TS=temp.n.CP.by.TS, alpha=alpha, beta=beta)
      }  
      
      ## Sample next state and add back counts
      previous.z = z
      z[i] = sample.ll(probs)
      m.CP.by.S[z[i],s[i]] = m.CP.by.S[z[i],s[i]] + 1
      n.CP.by.TS[z[i],] = n.CP.by.TS[z[i],] + n.TS.by.C[,i]
      
      ## Perform check for empty clusters; Do not allow on last iteration
      if(sum(z == previous.z[i]) == 0 & iter < max.iter) {
        
        ## Split another cluster into two
        z = split.z(counts=counts, z=z, empty.K=previous.z[i], K=K, LLFunction="cCG.calcLLFromVariables", s=s, y=y, L=L, alpha=alpha, beta=beta, delta=1, gamma=1)
        
        ## Re-calculate variables
        m.CP.by.S = table(factor(z, levels=1:K), s)
        n.CP.by.TS = rowsum(t(n.TS.by.C), group=z, reorder=TRUE)
      }
    
      z.probs[i,] = probs
    }
    
    ## Begin process of Gibbs sampling for each gene
    n.CP.by.G = rowsum(t(counts), group=z, reorder=TRUE)
    ix = sample(1:nrow(counts))
    for(i in ix) {
      if(sum(y == y[i]) > 1) {
        
        ## Subtract current gene counts from matrices
        nG.by.TS[y[i]] = nG.by.TS[y[i]] - 1
        n.CP.by.TS[,y[i]] = n.CP.by.TS[,y[i]] - n.CP.by.G[,i]
        n.by.TS[y[i]] = n.by.TS[y[i]] - n.by.G[i]
        
        probs = rep(NA, L)
        for(j in 1:L) {
          ## Add in counts to each state and determine probability
          temp.n.CP.by.TS = n.CP.by.TS
          temp.n.CP.by.TS[,j] = temp.n.CP.by.TS[,j] + n.CP.by.G[,i]
          temp.n.by.TS = n.by.TS
          temp.n.by.TS[j] = temp.n.by.TS[j] + n.by.G[i]
          temp.nG.by.TS = nG.by.TS
          temp.nG.by.TS[j] = temp.nG.by.TS[j] + 1
          
          probs[j] = cCG.calcGibbsProbY(n.CP.by.TS=temp.n.CP.by.TS, n.by.TS=temp.n.by.TS, nG.by.TS=temp.nG.by.TS, nG.in.Y=nG.by.TS[j], beta=beta, gamma=gamma, delta=delta)
        }  
        
        ## Sample next state and add back counts
        y[i] = sample.ll(probs)
        nG.by.TS[y[i]] = nG.by.TS[y[i]] + 1
        n.CP.by.TS[,y[i]] = n.CP.by.TS[,y[i]] + n.CP.by.G[,i]
        n.by.TS[y[i]] = n.by.TS[y[i]] + n.by.G[i]
        
        
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
    temp.ll = cCG.calcLL(K=K, L=L, m.CP.by.S=m.CP.by.S, n.CP.by.TS=n.CP.by.TS, n.by.G=n.by.G, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, nS=nS, nG=nG, alpha=alpha, beta=beta, gamma=gamma, delta=delta)
    if((best == TRUE & all(temp.ll > ll)) | iter == 1) {
      z.probs.final = z.probs
      y.probs.final = y.probs
    }
    ll = c(ll, temp.ll)
    
    message(date(), " ... Completed iteration: ", iter, " | logLik: ", temp.ll)
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
  
  return(list(z=z.final, y=y.final, complete.z=z.all, complete.y=y.all, 
              z.stability=z.stability.final, y.stability=y.stability.final, 
              complete.z.stability=z.stability, complete.y.stability=y.stability, 
              completeLogLik=ll, finalLogLik=ll.final, z.prob=z.probs.final, y.prob=y.probs.final,
              seed=seed))
}
