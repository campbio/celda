# -----------------------------------
# Variable description
# -----------------------------------
# C = Cell
# S or s = Sample
# G = Gene
# CP = Cell population
# n = counts of transcripts
# m = counts of cells
# K = Total number of cell populations
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
simulateCells.celda_C = function(S=10, C.Range=c(10, 100), N.Range=c(100,5000), 
                         G=500, K=5, alpha=1, beta=1) {
  
  phi <- gtools::rdirichlet(K, rep(beta, G))
  theta <- gtools::rdirichlet(S, rep(alpha, K))
  
  ## Select the number of cells per sample
  nC <- sample(C.Range[1]:C.Range[2], size=S, replace=TRUE)  
  cell.sample <- rep(1:S, nC)
  
  ## Select state of the cells  
  cell.state <- unlist(lapply(1:S, function(i) sample(1:K, size=nC[i], prob=theta[i,], replace=TRUE)))
  
  ## Select number of transcripts per cell
  nN <- sample(N.Range[1]:N.Range[2], size=length(cell.sample), replace=TRUE)
  
  ## Select transcript distribution for each cell
  cell.counts <- sapply(1:length(cell.sample), function(i) rmultinom(1, size=nN[i], prob=phi[cell.state[i],]))
  
  return(list(z=cell.state, counts=cell.counts, sample=cell.sample, K=K, alpha=alpha, beta=beta))
}

#' @export
celda_C = function(counts, sample.label, K, alpha=1, beta=1, max.iter=25, min.cell=5, 
                   seed=12345, best=TRUE, kick=TRUE) {
  
  if(is.factor(sample.label)) {
    s = as.numeric(sample.label)
  }
  else {
    s = as.numeric(as.factor(sample.label))
  }  
  
  set.seed(seed)
  message(date(), " ... Starting Gibbs sampling")
  
  z = sample(1:K, ncol(counts), replace=TRUE)
  z.all = z
  z.stability = c(NA)
  z.probs = matrix(NA, nrow=ncol(counts), ncol=K)
  
  ## Calculate counts one time up front
  m.CP.by.S = table(factor(z, levels=1:K), s)
  n.CP.by.G = rowsum(t(counts), group=z, reorder=TRUE)

  ll = cC.calcLL(m.CP.by.S=m.CP.by.S, n.CP.by.G=n.CP.by.G, s=s, K=K, alpha=alpha, beta=beta)

  iter = 1
  continue = TRUE
  while(iter <= max.iter & continue == TRUE) {
    

    ## Begin process of Gibbs sampling for each cell
    ix = sample(1:ncol(counts))
    for(i in ix) {
      
      ## Subtract current cell counts from matrices
      m.CP.by.S[z[i],s[i]] = m.CP.by.S[z[i],s[i]] - 1
      n.CP.by.G[z[i],] = n.CP.by.G[z[i],] - counts[,i]
        
      ## Calculate probabilities for each state
      probs = rep(NA, K)
      for(j in 1:K) {
        temp.n.CP.by.G = n.CP.by.G
        temp.n.CP.by.G[j,] = temp.n.CP.by.G[j,] + counts[,i]
        probs[j] = cC.calcGibbsProbZ(m.CP.by.S=m.CP.by.S[j,s[i]], n.CP.by.G=temp.n.CP.by.G, alpha=alpha, beta=beta)
      }  

      ## Sample next state and add back counts
      previous.z = z
      z[i] = sample.ll(probs)
      m.CP.by.S[z[i],s[i]] = m.CP.by.S[z[i],s[i]] + 1
      n.CP.by.G[z[i],] = n.CP.by.G[z[i],] + counts[,i]
      
      ## Perform check for empty clusters; Do not allow on last iteration
      if(sum(z == previous.z[i]) == 0 & iter < max.iter) {
        
        ## Split another cluster into two
        z = split.z(counts=counts, z=z, empty.K=previous.z[i], K=K, LLFunction="cC.calcLLFromVariables", s=s, alpha=alpha, beta=beta)
        
        ## Re-calculate variables
        m.CP.by.S = table(factor(z, levels=1:K), s)
        n.CP.by.G = rowsum(t(counts), group=z, reorder=TRUE)
      }
       
      z.probs[i,] = probs
    }  
    

    ## Save history
    z.all = cbind(z.all, z)

    ## Normalize Z and Y marginal probabilties and calculate stability
    z.probs = normalizeLogProbs(z.probs)
    z.stability = c(z.stability, stability(z.probs))

    ## Calculate complete likelihood
    temp.ll = cC.calcLL(m.CP.by.S=m.CP.by.S, n.CP.by.G=n.CP.by.G, s=s, K=K, alpha=alpha, beta=beta)
    if((best == TRUE & all(temp.ll > ll)) | iter == 1) {
      z.probs.final = z.probs
    }
    ll = c(ll, temp.ll)
    
    message(date(), " ... Completed iteration: ", iter, " | logLik: ", temp.ll)
    
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
              finalLogLik=ll.final, z.probability=z.probs, seed=seed))
}


cC.calcGibbsProbZ = function(m.CP.by.S, n.CP.by.G, alpha, beta) {
  
  ## Calculate for "Theta" component
  theta.ll = log(m.CP.by.S + alpha)
  
  ## Calculate for "Phi" component
  b = sum(lgamma(n.CP.by.G + beta))
  d = -sum(lgamma(rowSums(n.CP.by.G + beta)))
  
  phi.ll = b + d
  
  final = theta.ll + phi.ll 
  return(final)
}

#' @export
cC.calcLLFromVariables = function(counts, s, z, K, alpha, beta) {
  
  ## Calculate for "Theta" component
  m.CP.by.S = table(z, s)
  nS = length(unique(s))
  
  a = nS * lgamma(K*alpha)
  b = sum(lgamma(m.CP.by.S + alpha))
  c = -nS * K * lgamma(alpha)
  d = -sum(lgamma(colSums(m.CP.by.S + alpha)))
  
  theta.ll = a + b + c + d
 
  ## Calculate for "Phi" component
  n.CP.by.G = rowsum(t(counts), group=z, reorder=TRUE)
  nG = ncol(n.CP.by.G)
  
  a = K * lgamma(nG * beta)
  b = sum(lgamma(n.CP.by.G + beta))
  c = -K * nG * lgamma(beta)
  d = -sum(lgamma(rowSums(n.CP.by.G + beta)))
  
  phi.ll = a + b + c + d

  final = theta.ll + phi.ll
  return(final)
}

cC.calcLL = function(m.CP.by.S, n.CP.by.G, s, z, K, alpha, beta) {
  
  ## Calculate for "Theta" component
  nS = length(unique(s))
  
  a = nS * lgamma(K * alpha)
  b = sum(lgamma(m.CP.by.S + alpha))
  c = -nS * K * lgamma(alpha)
  d = -sum(lgamma(colSums(m.CP.by.S + alpha)))
  
  theta.ll = a + b + c + d
  
  ## Calculate for "Phi" component
  nG = ncol(n.CP.by.G)
  
  a = K * lgamma(nG * beta)
  b = sum(lgamma(n.CP.by.G + beta))
  c = -K * nG * lgamma(beta)
  d = -sum(lgamma(rowSums(n.CP.by.G + beta)))
  
  phi.ll = a + b + c + d
  
  final = theta.ll + phi.ll
  return(final)
}

