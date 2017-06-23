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

#' Calculate log lileklihood for the celda Cell and Gene clustering model, given a set of cell / gene cluster assignments
#' 
#' @param counts A numeric count matrix
#' @param s A numeric vector of sample labels corresponding to the samples each cell in the count matrix originates from
#' @param z A numeric vector of cell cluster assignments
#' @param y A numeric vector of gene cluster assignments
#' @param K The number of cell clusters
#' @param L The number of gene clusters
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution.
#' @param beta The Dirichlet distribution parameter for Phi; adds a pseudocount to each transcriptional state within each cell.
#' @param delta The Dirichlet distribution parameter for Eta; adds a gene pseudocount to the numbers of genes each state.
#' @param gamma The Dirichlet distribution parameter for Psi; adds a pseudocount to each gene within each transcriptional state.
#' @export
calculate_loglik_from_variables.celda_CG = function(counts, s, z, y, K, L, alpha, beta, delta, gamma) {
  
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

cCG.calcLL = function(K, L, m.CP.by.S, n.CP.by.TS, n.by.G, n.by.TS, nG.by.TS, nS, nG, alpha, beta, delta, gamma) {

  ## Determine if any TS has 0 genes
  ## Need to remove 0 gene states as this will cause the likelihood to fail
  if(sum(nG.by.TS == 0) > 0) {
    ind = which(nG.by.TS > 0)
    L = length(ind)
    n.CP.by.TS = n.CP.by.TS[,ind]
    n.by.TS = n.by.TS[ind]
    nG.by.TS = nG.by.TS[ind]
  }
  
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

cCG.calcGibbsProbY = function(n.CP.by.TS, n.by.TS, nG.by.TS, nG.in.Y, beta, delta, gamma) {
  
  ## Determine if any TS has 0 genes
  ## Need to remove 0 gene states as this will cause the likelihood to fail
  if(sum(nG.by.TS == 0) > 0) {
    ind = which(nG.by.TS > 0)
    n.CP.by.TS = n.CP.by.TS[,ind]
    n.by.TS = n.by.TS[ind]
    nG.by.TS = nG.by.TS[ind]
  }

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

#' Simulate cells from the cell/gene clustering generative model
#' 
#' @param S The number of samples
#' @param C.Range two element vector to specify the lower and upper bound of the counts of cells for each sample
#' @param N.Range two element vector to specify the lower and upper bound of the counts of the transcripts
#' @param G The number of genes
#' @param K The number of cell populations
#' @param L The number of gene clusters being considered
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution
#' @param beta The Dirichlet distribution parameter for Phi; adds a pseudocount to each transcriptional state within each cell
#' @param gamma The Dirichlet distribution parameter for Psi; adds a pseudocount to each gene within each transcriptional state
#' @param delta The Dirichlet distribution parameter for Eta; adds a gene pseudocount to the numbers of genes each state
#' @param seed starting point used for generating simulated data
#' @param ... Additional parameters
#' @export
simulateCells.celda_CG = function(S=10, C.Range=c(50,100), N.Range=c(500,5000), 
                                  G=1000, K=3, L=10, alpha=1, beta=1, gamma=1, 
                                  delta=1, seed=12345, ...) {
  
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
  
  new = reorder.labels.by.size.then.counts(cell.counts, z=z, y=y, K=K, L=L)
  
  ## Ensure that there are no all-0 rows in the counts matrix, which violates a celda modeling
  ## constraint (columns are guarnteed at least one count):
  zero.row.idx = which(rowSums(cell.counts) == 0)
  cell.counts = cell.counts[-zero.row.idx, ]
  new$y = new$y[-zero.row.idx]
  
  return(list(z=new$z, y=new$y, sample.label=cell.sample.label, counts=cell.counts, K=K, L=L, C.Range=C.Range, N.Range=N.Range, S=S, alpha=alpha, beta=beta, gamma=gamma, delta=delta, theta=theta, phi=phi, psi=psi, eta=eta, seed=seed))
}

#' celda Cell and Gene Clustering Model
#' 
#' @param counts A numeric count matrix.
#' @param sample.label A vector indicating the sample for each cell in the count matrix
#' @param K The number of cell populations
#' @param L The number of gene clusters being considered
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution
#' @param beta The Dirichlet distribution parameter for Phi; adds a pseudocount to each transcriptional state within each cell. Default to 1
#' @param delta The Dirichlet distribution parameter for Eta; adds a gene pseudocount to the numbers of genes each state. Default to 1
#' @param gamma The Dirichlet distribution parameter for Psi; adds a pseudocount to each gene within each transcriptional state. Default to 1
#' @param count.checksum An MD5 checksum for the provided counts matrix
#' @param max.iter Maximum iterations of Gibbs sampling to perform. Defaults to 25
#' @param seed Parameter to set.seed() for random number generation
#' @param best Whether to return the cluster assignment with the highest log-likelihood. Defaults to TRUE. Returns last generated cluster assignment when FALSE
#' @param z.split.on.iter On z.split.on.iter-th iterations, a heuristic will be applied using hierarchical clustering to determine if a cell cluster should be merged with another cell cluster and a third cell cluster should be split into two clusters. This helps avoid local optimum during the initialization. Default to be 3
#' @param z.num.splits Maximum number of times to perform the heuristic described in z.split.on.iter
#' @param y.split.on.iter  On every y.split.on.iter iteration, a heuristic will be applied using hierarchical clustering to determine if a gene cluster should be merged with another gene cluster and a third gene cluster should be split into two clusters. This helps avoid local optimum during the initialization. Default to be 3
#' @param y.num.splits Maximum number of times to perform the heuristic described in y.split.on.iter
#' @param thread The thread index, used for logging purposes
#' @param save.history Logical; whether to return the history of cluster assignments. Defaults to FALSE
#' @param save.prob Logical; whether to return the history of cluster assignment probabilities. Defaults to FALSE
#' @param ... Additional parameters
#' @export
celda_CG = function(counts, sample.label=NULL, K, L, alpha=1, beta=1, 
                    delta=1, gamma=1, count.checksum=NULL, max.iter=25,
			              seed=12345, best=TRUE, z.split.on.iter=3, z.num.splits=3,
			              y.split.on.iter=3, y.num.splits=3, thread=1, 
			              save.history=FALSE, save.prob=FALSE, ...) {
  set.seed(seed)
  
  message("Thread ", thread, " ", date(), " ... Starting Gibbs sampling")
  
  if(is.null(sample.label)) {
    s = rep(1, ncol(counts))
    sample.label = s
  } else if(is.factor(sample.label)) {
    s = as.numeric(sample.label)
  } else {
    sample.label = as.factor(sample.label)
    s = as.numeric(sample.label)
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
  nS = length(unique(s))
  m.CP.by.S = matrix(table(factor(z, levels=1:K), s), ncol=nS)
  n.TS.by.C = rowsum(counts, group=y, reorder=TRUE)
  n.CP.by.TS = rowsum(t(n.TS.by.C), group=z, reorder=TRUE)
  n.by.G = rowSums(counts)
  n.by.TS = as.numeric(rowsum(n.by.G, y))
  nG.by.TS = table(y)

  nG = nrow(counts)
  nM = ncol(counts)
  
  ll = cCG.calcLL(K=K, L=L, m.CP.by.S=m.CP.by.S, n.CP.by.TS=n.CP.by.TS, n.by.G=n.by.G, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, nS=nS, nG=nG, alpha=alpha, beta=beta, delta=delta, gamma=gamma)
  
  iter = 1
  continue = TRUE
  z.num.of.splits.occurred = 1
  y.num.of.splits.occurred = 1
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
        z = split.z(counts=counts, z=z, empty.K=previous.z[i], K=K, LLFunction="calculate_loglik_from_variables.celda_CG", s=s, y=y, L=L, alpha=alpha, beta=beta, delta=delta, gamma=gamma)
        
        ## Re-calculate variables
        m.CP.by.S = matrix(table(factor(z, levels=1:K), s), ncol=nS)
        n.CP.by.TS = rowsum(t(n.TS.by.C), group=z, reorder=TRUE)
      }
    
      z.probs[i,] = probs
    }
    
    ## Begin process of Gibbs sampling for each gene
    n.CP.by.G = rowsum(t(counts), group=z, reorder=TRUE)
    ix = sample(1:nrow(counts))
    for(i in ix) {
        
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
	   
		probs[j] = cCG.calcGibbsProbY(n.CP.by.TS=temp.n.CP.by.TS, n.by.TS=temp.n.by.TS, nG.by.TS=temp.nG.by.TS, nG.in.Y=temp.nG.by.TS[j], beta=beta, delta=delta, gamma=gamma)
	  }  

	  ## Sample next state and add back counts
	  previous.y = y
	  y[i] = sample.ll(probs)
	  nG.by.TS[y[i]] = nG.by.TS[y[i]] + 1
	  n.CP.by.TS[,y[i]] = n.CP.by.TS[,y[i]] + n.CP.by.G[,i]
      n.by.TS[y[i]] = n.by.TS[y[i]] + n.by.G[i]
 
      ## Perform check for empty clusters; Do not allow on last iteration
      if(sum(y == previous.y[i]) == 0 & iter < max.iter) {
        
        ## Split another cluster into two
        y = split.y(counts=counts, y=y, empty.L=previous.y[i], L=L, LLFunction="calculate_loglik_from_variables.celda_CG", z=z, s=s, K=K, alpha=alpha, beta=beta, delta=delta, gamma=gamma)
        
        ## Re-calculate variables
        n.TS.by.C = rowsum(counts, group=y, reorder=TRUE)
        n.CP.by.TS = rowsum(t(n.TS.by.C), group=z, reorder=TRUE)
        n.by.TS = as.numeric(rowsum(n.by.G, y))
        nG.by.TS = table(y)
      }
       
      y.probs[i,] = probs
    }
    
    ## Perform split on i-th iteration defined by z.split.on.iter
    if(iter %% z.split.on.iter == 0 & z.num.of.splits.occurred <= z.num.splits & K > 2) {

      message("Thread ", thread, " ", date(), " ... Determining if any cell clusters should be split (", z.num.of.splits.occurred, " of ", z.num.splits, ")")
      z = split.each.z(counts=counts, z=z, y=y, K=K, L=L, alpha=alpha, delta=delta, beta=beta, gamma=gamma, s=s, LLFunction="calculate_loglik_from_variables.celda_CG")
      z.num.of.splits.occurred = z.num.of.splits.occurred + 1

      ## Re-calculate variables
      m.CP.by.S = matrix(table(factor(z, levels=1:K), s), ncol=nS)
      n.TS.by.C = rowsum(counts, group=y, reorder=TRUE)
      n.CP.by.TS = rowsum(t(n.TS.by.C), group=z, reorder=TRUE)
      n.CP.by.G = rowsum(t(counts), group=z, reorder=TRUE)      
    }
    
    ## Perform split if on i-th iteration defined by y.split.on.iter
    if(iter %% y.split.on.iter == 0 & y.num.of.splits.occurred <= y.num.splits & L > 2) {

      message("Thread ", thread, " ", date(), " ... Determining if any gene clusters should be split (", y.num.of.splits.occurred, " of ", y.num.splits, ")")
      y = split.each.y(counts=counts, z=z, y=y, K=K, L=L, alpha=alpha, beta=beta, delta=delta, gamma=gamma, s=s, LLFunction="calculate_loglik_from_variables.celda_CG")
      y.num.of.splits.occurred = y.num.of.splits.occurred + 1

      ## Re-calculate variables
      n.TS.by.C = rowsum(counts, group=y, reorder=TRUE)
      n.CP.by.TS = rowsum(t(n.TS.by.C), group=z, reorder=TRUE)
      n.by.TS = as.numeric(rowsum(n.by.G, y))
      nG.by.TS = table(y)
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
    temp.ll = cCG.calcLL(K=K, L=L, m.CP.by.S=m.CP.by.S, n.CP.by.TS=n.CP.by.TS, n.by.G=n.by.G, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, nS=nS, nG=nG, alpha=alpha, beta=beta, delta=delta, gamma=gamma)
    if((best == TRUE & all(temp.ll > ll)) | iter == 1) {
      z.probs.final = z.probs
      y.probs.final = y.probs
    }
    ll = c(ll, temp.ll)
    
    message("Thread ", thread, " ", date(), " ... Completed iteration: ", iter, " | logLik: ", temp.ll)
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
  
  ## Peform reordering on final Z and Y assigments:
  reordered.labels = reorder.labels.by.size.then.counts(counts, z=z.final, 
                                                        y=y.final, K=K, L=L)
  names = list(row=rownames(counts), column=colnames(counts), 
               sample=levels(sample.label))
  
  
  result = list(z=reordered.labels$z, y=reordered.labels$y, 
                z.stability=z.stability.final, 
                y.stability=y.stability.final,  
                complete.z.stability=z.stability, 
                complete.y.stability=y.stability, completeLogLik=ll, 
                finalLogLik=ll.final, K=K, L=L, alpha=alpha, 
                beta=beta, delta=delta, seed=seed, 
                sample.label=sample.label, names=names,
                count.checksum=count.checksum)
  
  if (save.prob) {
    result$z.prob = z.probs.final
    result$y.prob = y.probs.final
  } 
  
  if (save.history) {
    ## Reorder Z and Y histories based off of final Z/Y reordering above:
    z.all = base::apply(z.all, 2, function(column) reordered.labels$z.map[column])
    y.all = base::apply(y.all, 2, function(column) reordered.labels$y.map[column])
    result$complete.z = z.all
    result$complete.y = y.all
  } 
  
   class(result) = "celda_CG" 
   return(result)
}


#' Generate factorized matrices showing each feature's influence on the celda_CG model clustering 
#' 
#' @param counts A numerix count matrix
#' @param celda.obj object returned from celda_CG function 
#' @param type one of the "counts", "proportion", or "posterior". 
#' @return A list of factorized matrices, of the types requested by the user. NOTE: "population" state matrices are always returned in cell population (rows) x transcriptional states (cols).
#' @export 
factorizeMatrix.celda_CG = function(counts, celda.obj, type=c("counts", "proportion", "posterior")) {

  K = celda.obj$K
  L = celda.obj$L
  z = celda.obj$z
  y = celda.obj$y
  alpha = celda.obj$alpha
  beta = celda.obj$beta
  delta = celda.obj$delta
  sample.label = celda.obj$sample.label
  
  counts.list = c()
  prop.list = c()
  post.list = c()
  res = list()

  counts.list = c()
  prop.list = c()
  post.list = c()
  res = list()
  
  nS = length(unique(sample.label))
  m.CP.by.S = matrix(table(factor(z, levels=1:K), sample.label), ncol=nS)
  n.TS.by.C = rowsum(counts, group=y, reorder=TRUE)
  n.CP.by.TS = rowsum(t(n.TS.by.C), group=z, reorder=TRUE)
  n.by.G = rowSums(counts)
  n.by.TS = as.numeric(rowsum(n.by.G, y))

  n.G.by.TS = matrix(0, nrow=length(y), ncol=L)
  for(i in 1:length(y)) {n.G.by.TS[i,y[i]] = n.by.G[i]}

  L.names = paste0("L", 1:L)
  K.names = paste0("K", 1:K)
  colnames(n.TS.by.C) = celda.obj$names$column
  rownames(n.TS.by.C) = L.names
  colnames(n.G.by.TS) = L.names
  rownames(n.G.by.TS) = celda.obj$names$row
  rownames(m.CP.by.S) = K.names
  colnames(m.CP.by.S) = celda.obj$names$sample
  colnames(n.CP.by.TS) = L.names
  rownames(n.CP.by.TS) = K.names
    
  if(any("counts" %in% type)) {
    counts.list = list(sample.states = m.CP.by.S,
            				   population.states = n.CP.by.TS, 
            				   cell.states = n.TS.by.C,
            				   gene.states = n.G.by.TS)
    res = c(res, list(counts=counts.list))
  }
  if(any("proportion" %in% type)) {
    prop.list = list(sample.states = normalizeCounts(m.CP.by.S, scale.factor=1),
    				   population.states = t(normalizeCounts(t(n.CP.by.TS), scale.factor=1)), 
    				   cell.states = normalizeCounts(n.TS.by.C, scale.factor=1),
    				   gene.states = normalizeCounts(n.G.by.TS, scale.factor=1))
    res = c(res, list(proportions=prop.list))
  }
  if(any("posterior" %in% type)) {
    post.list = list(sample.states = normalizeCounts(m.CP.by.S + alpha, scale.factor=1),
          				   population.states = normalizeCounts(n.CP.by.TS + beta, scale.factor=1), 
          				   gene.states = normalizeCounts(n.G.by.TS + delta, scale.factor=1))
    res = c(res, posterior = list(post.list))						    
  }
  
  return(res)
}


################################################################################
# celda_CG S3 methods                                                          #
################################################################################
#' finalClusterAssignment for the celda Cell and Gene clustering model
#' @param celda.mod A celda model object of "Celda_CG"
#' @export
finalClusterAssignment.celda_CG = function(celda.mod) {
  return(list(z=celda.mod$z, y=celda.mod$y))
}

#' completeClusterHistory for the celda Cell and Gene clustering model
#' @param celda.mod A celda model object of "Celda_CG"
#' @export
completeClusterHistory.celda_CG = function(celda.mod) {
  return(list(complete.z=celda.mod$complete.z, complete.y=celda.mod$complete.y))
}

#' clusterProbabilities for the celda Cell and Gene clustering model
#' @param celda.mod A celda model object of "Celda_CG"
#' @export
clusterProbabilities.celda_CG = function(celda.mod) {
  return(list(z.prob=celda.mod$z.prob, y.prob=celda.mod$y.prob))
}

#' getK for the celda Cell and Gene clustering model
#' @param celda.mod A celda model object of "Celda_CG"
#' @export
getK.celda_CG = function(celda.mod) {
  return(celda.mod$K)
}


#' getL for the celda Cell and Gene clustering model
#' @param celda.mod A celda model object of "Celda_CG"
#' @export
getL.celda_CG = function(celda.mod) {
  return(celda.mod$L)
}


#' celda_heatmap for celda Cell and Gene clustering model
#' @param celda.mod A celda model object of "Celda_CG"
#' @param counts A count matrix
#' @param ... extra parameters passed onto render_celda_heatmap
#' @export
celda_heatmap.celda_CG = function(celda.mod, counts, ...) {
  render_celda_heatmap(counts, z=celda.mod$z, y=celda.mod$y, ...)
}


#' Visualize the performance of celda_CG models grouped by L and K
#' 
#' Plot the performance of a list of celda_CG models returned 
#' from running celda function. For each number of gene clusters 
#' L (cell clusters K), plot the performance of each number of cell
#' clusters K (gene clusters L).
#' @param celda.list A list of celda_CG objects returned from celda function
#' @param method One of "perplexity", "harmonic", or "loglik"
#' @param title Title for the visualize_model_performance
#' @import Rmpfr gridExtra
#' @export
visualize_model_performance.celda_CG = function(celda.list, 
                    method = "perplexity", title="Model Performance (All Chains)") {
  # validate input parameters
  if (class(small.sim.res) != "celda_list") {
    stop("celda.list argument must be of class 'celda_list'")
  } else if (celda.list$content.type != "celda_CG") {
    stop("celda.list must be a 'celda.list' of 'celda_CG' objects")
  } else if (!(method %in% c("perplexity","harmonic","loglik"))) {
    stop("Invalid method, 'method' has to be either 'perplexity', 'harmonic', or 'loglik'")
  } 
  # if every thing went well
  else {
    y.lab = method
    cluster.sizes = unlist(lapply(celda.list$res.list, function(mod) { getK(mod) }))
    log.likelihoods = lapply(celda.list$res.list,
                             function(mod) { completeLogLikelihood(mod) })
    performance.metric = lapply(log.likelihoods, 
                                calculate_performance_metric,
                                method)
    
    # These methods return Rmpfr numbers that are extremely small and can't be 
    # plotted, so log 'em first
    if (method %in% c("perplexity", "harmonic")) {
      performance.metric = lapply(performance.metric, log)
      performance.metric = new("mpfr", unlist(performance.metric))
      performance.metric = as.numeric(performance.metric)
      y.lab = paste0("Log(",method,")")
    } else {
      performance.metric = as.numeric(performance.metric)
    }
    
    plot.df = data.frame(K=cluster.sizes, L=celda.list$run.params$L,
                         metric=performance.metric)
    
    L.list = sort(unique(celda.list$run.params$L))
    K.list = sort(unique(celda.list$run.params$K))
    #chains = sort(unique(celda.list$run.params$chain))
    
    plots = list()
    nc = round(length(L.list)^.5)
    x.lab.K = "K"
    
    for (i in L.list) {
      plots = c(plots, list(ggplot2::ggplot(subset(plot.df, L==i), 
                                            ggplot2::aes(x=K, y=metric, group=K)) + 
                              ggplot2::geom_boxplot(outlier.color=NA, fill=NA) + 
                              ggplot2::geom_point(position=ggplot2::position_jitter(width=0.1, height=0)) +
                              ggplot2::ggtitle(paste0("L = ", i)) + 
                              ggplot2::theme_bw() + 
                              ggplot2::theme(axis.title.x=ggplot2::element_blank(), 
                                             axis.title.y=ggplot2::element_blank())))
    }
    gridExtra::grid.arrange(grobs = plots, ncol = nc, left = grid::textGrob(y.lab, rot = 90), 
                            top = grid::textGrob(title),
                            bottom = grid::textGrob(x.lab.K))
    
    
    plots = list()
    nc = round(length(K.list)^.5)
    x.lab.L = "L"
    
    for (i in K.list) {
      plots = c(plots, list(ggplot2::ggplot(subset(plot.df, K==i), 
                                            ggplot2::aes(x=L, y=metric, group=L)) +
                              ggplot2::geom_boxplot(outlier.color=NA, fill=NA) + 
                              ggplot2::geom_point(position=ggplot2::position_jitter(width=0.1, height=0)) + 
                              ggplot2::ggtitle(paste0("K = ", i)) + 
                              ggplot2::theme_bw() + 
                              ggplot2::theme(axis.title.x=ggplot2::element_blank(), 
                                             axis.title.y=ggplot2::element_blank())))
    }
    gridExtra::grid.arrange(grobs = plots, ncol = nc, left = grid::textGrob(y.lab, rot = 90),
                            top = grid::textGrob(title),
                            bottom = grid::textGrob(x.lab.L))
  }
}
