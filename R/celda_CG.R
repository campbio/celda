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


#' celda Cell and Gene Clustering Model
#' 
#' @param counts A numeric count matrix.
#' @param sample.label A vector indicating the sample for each cell in the count matrix
#' @param K The number of cell populations
#' @param L The number of gene clusters being considered
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution
#' @param beta The Dirichlet distribution parameter for Phi; adds a pseudocount to each transcriptional state within each cell. Default 1.
#' @param delta The Dirichlet distribution parameter for Eta; adds a gene pseudocount to the numbers of genes each state. Default 1.
#' @param gamma The Dirichlet distribution parameter for Psi; adds a pseudocount to each gene within each transcriptional state. Default 1.
#' @param stop.iter Number of iterations without improvement in the log likelihood to stop the Gibbs sampler. Default 10.
#' @param split.on.iter On every 'split.on.iter' iteration, a heuristic will be applied to determine if a gene/cell cluster should be reassigned and another gene/cell cluster should be split into two clusters. Default 10.
#' @param max.iter Maximum iterations of Gibbs sampling to perform regardless of convergence. Default 200.
#' @param seed Parameter to set.seed() for random number generation
#' @param count.checksum An MD5 checksum for the provided counts matrix
#' @param z.init Initial values of z. If NULL, z will be randomly sampled. Default NULL.
#' @param y.init Initial values of y. If NULL, y will be randomly sampled. Default NULL.
#' @param process.counts Whether to cast the counts matrix to integer and round(). Defaults to TRUE.
#' @param logfile The name of the logfile to redirect messages to.
#' @export
celda_CG = function(counts, sample.label=NULL, K, L,
                    alpha=1, beta=1, delta=1, gamma=1, 
                    max.iter=200, stop.iter = 10, split.on.iter=10,
                    seed=12345, count.checksum=NULL,
                    z.init = NULL, y.init = NULL, process.counts=TRUE, logfile=NULL) {
  
  ## Error checking and variable processing
  if (isTRUE(processCounts)) {
    counts = processCounts(counts)
  }
    
  s = processSampleLabels(sample.label, ncol(counts))
  if (is.null(sample.label)) {
    sample.label = s
  }
  
  ## Randomly select z and y or set z/y to supplied initial values
  z = initialize.cluster(K, ncol(counts), initial = z.init, fixed = NULL, seed=seed)
  y = initialize.cluster(L, nrow(counts), initial = y.init, fixed = NULL, seed=seed)
  z.best = z
  y.best = y  
  
  ## Calculate counts one time up front
  nS = length(unique(s))
  m.CP.by.S = matrix(as.integer(table(factor(z, levels=1:K), s)), ncol=nS)
  n.TS.by.C = rowsum.y(counts, y=y, L=L)
  n.CP.by.TS = rowsum.z(n.TS.by.C, z=z, K=K)
  n.CP = as.integer(rowSums(n.CP.by.TS))
  n.by.G = as.integer(rowSums(counts))
  n.by.C = as.integer(colSums(counts))
  n.by.TS = as.integer(rowsum.y(matrix(n.by.G,ncol=1), y=y, L=L))
  nG.by.TS = as.integer(table(factor(y, 1:L)))

  nG = nrow(counts)
  nM = ncol(counts)
  
  ll = cCG.calcLL(K=K, L=L, m.CP.by.S=m.CP.by.S, n.CP.by.TS=n.CP.by.TS, n.by.G=n.by.G, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, nS=nS, nG=nG, alpha=alpha, beta=beta, delta=delta, gamma=gamma)

  set.seed(seed)
  logMessages(date(), "... Starting Gibbs sampling", logfile=logfile, append=FALSE)
  
  iter = 1L
  num.iter.without.improvement = 0L
  do.cell.split = TRUE
  do.gene.split = TRUE  
  while(iter <= max.iter & num.iter.without.improvement <= stop.iter) {
    
    ## Gibbs sampling for each cell
    n.TS.by.C = rowsum.y(counts, y=y, L=L)
    n.TS.by.CP = t(n.CP.by.TS)
	next.z = cC.calcGibbsProbZ(counts=n.TS.by.C, m.CP.by.S=m.CP.by.S, n.G.by.CP=n.TS.by.CP, n.CP=n.CP, n.by.C=n.by.C, z=z, s=s, K=K, nG=L, nM=nM, alpha=alpha, beta=beta)
    m.CP.by.S = next.z$m.CP.by.S
    n.TS.by.CP = next.z$n.G.by.CP
    n.CP = next.z$n.CP
    z = next.z$z

    ## Gibbs sampling for each gene
    n.CP.by.G = rowsum.z(counts, z=z, K=K)
    n.CP.by.TS = t(n.TS.by.CP) 
 	next.y = cG.calcGibbsProbY(counts.t=n.CP.by.G, n.C.by.TS=n.CP.by.TS, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, n.by.G=n.by.G, y=y, L=L, nG=nG, beta=beta, delta=delta, gamma=gamma)
	n.CP.by.TS = next.y$n.C.by.TS
	nG.by.TS = next.y$nG.by.TS
	n.by.TS = next.y$n.by.TS
	y = next.y$y
    
    ## Perform split on i-th iteration defined by split.on.iter
	if(K > 2 & (num.iter.without.improvement == stop.iter | (iter %% split.on.iter == 0 & isTRUE(do.cell.split)))) {
	  logMessages(date(), " ... Determining if any cell clusters should be split.", logfile=logfile, append=TRUE, sep="")
	  res = split.each.z(counts=counts, z=z, y=y, z.prob=t(next.z$probs), K=K, L=L, alpha=alpha, delta=delta, beta=beta, gamma=gamma, s=s, LLFunction="calculateLoglikFromVariables.celda_CG")
	  logMessages(res$message, logfile=logfile, append=TRUE)

	  # Reset convergence counter if a split occured
	  if(!isTRUE(all.equal(z, res$z))) {
		num.iter.without.improvement = 0L
		do.cell.split = TRUE
	  } else {
		do.cell.split = FALSE
	  }

	  ## Re-calculate variables
	  z = res$z      
	  m.CP.by.S = matrix(as.integer(table(factor(z, levels=1:K), s)), ncol=nS)
	  n.TS.by.C = rowsum.y(counts, y=y, L=L)
	  n.CP.by.TS = rowsum.z(n.TS.by.C, z=z, K=K)
	  n.CP = as.integer(rowSums(n.CP.by.TS))
	  n.CP.by.G = rowsum.z(counts, z=z, K=K)     
	}  
	if(L > 2 & (num.iter.without.improvement == stop.iter | (iter %% split.on.iter == 0 & isTRUE(do.gene.split)))) {
	  logMessages(date(), " ... Determining if any gene clusters should be split.", logfile=logfile, append=TRUE, sep="")
	  res = split.each.y(counts=counts, z=z, y=y, y.prob=t(next.y$probs), K=K, L=L, alpha=alpha, beta=beta, delta=delta, gamma=gamma, s=s, LLFunction="calculateLoglikFromVariables.celda_CG")
	  logMessages(res$message, logfile=logfile, append=TRUE)

	  # Reset convergence counter if a split occured	    
	  if(!isTRUE(all.equal(y, res$y))) {
		num.iter.without.improvement = 1L
		do.gene.split = TRUE
	  } else {
		do.gene.split = FALSE
	  }

	  ## Re-calculate variables
	  y = res$y        
	  n.TS.by.C = rowsum.y(counts, y=y, L=L)
	  n.CP.by.TS = rowsum.z(n.TS.by.C, z=z, K=K)
	  n.CP = as.integer(rowSums(n.CP.by.TS))
	  n.by.TS = as.integer(rowsum.y(matrix(n.by.G,ncol=1), y=y, L=L))
	  nG.by.TS = as.integer(table(factor(y, levels=1:L)))
	}      

    ## Calculate complete likelihood
    temp.ll = cCG.calcLL(K=K, L=L, m.CP.by.S=m.CP.by.S, n.CP.by.TS=n.CP.by.TS, n.by.G=n.by.G, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, nS=nS, nG=nG, alpha=alpha, beta=beta, delta=delta, gamma=gamma)
    if((all(temp.ll > ll)) | iter == 1) {
      z.best = z
      y.best = y
      ll.best = temp.ll
      num.iter.without.improvement = 1L
    } else {  
      num.iter.without.improvement = num.iter.without.improvement + 1L   
    }
    ll = c(ll, temp.ll)
  
    logMessages(date(), " ... Completed iteration: ", iter, " | logLik: ", temp.ll, logfile=logfile, append=TRUE, sep="")
    iter = iter + 1L
    
  }
    

  names = list(row=rownames(counts), column=colnames(counts), 
               sample=levels(sample.label))
  
  result = list(z=z.best, y=y.best, completeLogLik=ll, 
                finalLogLik=ll.best, K=K, L=L, alpha=alpha, 
                beta=beta, delta=delta, gamma=gamma, seed=seed, 
                sample.label=sample.label, names=names,
                count.checksum=count.checksum)
  
  class(result) = "celda_CG" 
   
  ## Peform reordering on final Z and Y assigments:
  result = reorder.celda_CG(counts = counts, res = result)
  return(result)
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
#' @param process.counts Whether to cast the counts matrix to integer and round(). Defaults to TRUE.
#' @param ... Unused arguments
#' @param model Dummy parameter for S3 dispatch
#' @export
simulateCells.celda_CG = function(model, S=10, C.Range=c(50,100), N.Range=c(500,5000), 
                                  G=1000, K=3, L=10, alpha=1, beta=1, gamma=1, 
                                  delta=1, seed=12345,
                                  process.counts=TRUE, ...) {
  
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
    transcriptional.state.dist = as.numeric(stats::rmultinom(1, size=nN[i], prob=phi[z[i],]))
    for(j in 1:L) {
      if(transcriptional.state.dist[j] > 0) {
        cell.counts[,i] = cell.counts[,i] + stats::rmultinom(1, size=transcriptional.state.dist[j], prob=psi[,j])
      }
    }  
  }
  
  new = reorder.labels.by.size.then.counts(cell.counts, z=z, y=y, K=K, L=L)
  
  ## Ensure that there are no all-0 rows in the counts matrix, which violates a celda modeling
  ## constraint (columns are guarnteed at least one count):
  zero.row.idx = which(rowSums(cell.counts) == 0)
  if (length(zero.row.idx > 0)) {
    cell.counts = cell.counts[-zero.row.idx, ]
    new$y = new$y[-zero.row.idx]
  } 
 
  ## Assign gene/cell/sample names 
  rownames(cell.counts) = paste0("Gene_", 1:nrow(cell.counts))
  colnames(cell.counts) = paste0("Cell_", 1:ncol(cell.counts))
  cell.sample.label = paste0("Sample_", 1:S)[cell.sample.label]

  return(list(z=new$z, y=new$y, sample.label=cell.sample.label, counts=cell.counts, K=K, L=L, C.Range=C.Range, N.Range=N.Range, S=S, alpha=alpha, beta=beta, gamma=gamma, delta=delta, theta=theta, phi=phi, psi=psi, eta=eta, seed=seed))
}


##' Generate factorized matrices showing each feature's influence on the celda_CG model clustering 
#' 
#' @param counts A numerix count matrix
#' @param celda.mod object returned from celda_CG function 
#' @param type one of the "counts", "proportion", or "posterior". 
#' @return A list of factorized matrices, of the types requested by the user. NOTE: "population" state matrices are always returned in cell population (rows) x transcriptional states (cols).
#' @export 
factorizeMatrix.celda_CG = function(celda.mod, counts, type=c("counts", "proportion", "posterior")) {

  K = celda.mod$K
  L = celda.mod$L
  z = celda.mod$z
  y = celda.mod$y
  alpha = celda.mod$alpha
  beta = celda.mod$beta
  delta = celda.mod$delta
  gamma = celda.mod$gamma
  sample.label = celda.mod$sample.label
  s = as.integer(as.factor(sample.label))
  
  ## Calculate counts one time up front
  nS = length(unique(s))
  m.CP.by.S = matrix(as.integer(table(factor(z, levels=1:K), s)), ncol=nS)
  n.TS.by.C = rowsum.y(counts, y=y, L=L)
  n.CP.by.TS = rowsum.z(n.TS.by.C, z=z, K=K)
  n.by.G = as.integer(rowSums(counts))
  n.by.TS = as.integer(rowsum.y(matrix(n.by.G,ncol=1), y=y, L=L))
  nG = nrow(counts)
  nM = ncol(counts)

  n.G.by.TS = matrix(0, nrow=length(y), ncol=L)
  n.G.by.TS[cbind(1:nG,y)] = n.by.G

  L.names = paste0("L", 1:L)
  K.names = paste0("K", 1:K)
  colnames(n.TS.by.C) = celda.mod$names$column
  rownames(n.TS.by.C) = L.names
  colnames(n.G.by.TS) = L.names
  rownames(n.G.by.TS) = celda.mod$names$row
  rownames(m.CP.by.S) = K.names
  colnames(m.CP.by.S) = celda.mod$names$sample
  colnames(n.CP.by.TS) = L.names
  rownames(n.CP.by.TS) = K.names

  counts.list = c()
  prop.list = c()
  post.list = c()
  res = list()
    
  if(any("counts" %in% type)) {
    counts.list = list(sample.states = m.CP.by.S,
            				   population.states = t(n.CP.by.TS), 
            				   cell.states = n.TS.by.C,
            				   gene.states = n.G.by.TS)
    res = c(res, list(counts=counts.list))
  }
  if(any("proportion" %in% type)) {

    ## Need to avoid normalizing cell/gene states with zero cells/genes
    unique.z = sort(unique(z))
    temp.n.CP.by.TS = t(n.CP.by.TS)
    temp.n.CP.by.TS[,unique.z] = normalizeCounts(temp.n.CP.by.TS[,unique.z], scale.factor=1)

    unique.y = sort(unique(y))
    temp.n.G.by.TS = n.G.by.TS
    temp.n.G.by.TS[,unique.y] = normalizeCounts(temp.n.G.by.TS[,unique.y], scale.factor=1)
    
    prop.list = list(sample.states =  normalizeCounts(m.CP.by.S, scale.factor=1),
    				   population.states = temp.n.CP.by.TS, 
    				   cell.states = normalizeCounts(n.TS.by.C, scale.factor=1),
    				   gene.states = temp.n.G.by.TS)
    res = c(res, list(proportions=prop.list))
  }
  if(any("posterior" %in% type)) {

    gs = n.G.by.TS
    gs[gs > 0] = gs[gs > 0] + delta
    gs = normalizeCounts(gs, scale.factor=1)

    post.list = list(sample.states = normalizeCounts(m.CP.by.S + alpha, scale.factor=1),
          				   population.states = normalizeCounts(t(n.CP.by.TS + beta), scale.factor=1), 
          				   gene.states = gs)
    res = c(res, posterior = list(post.list))						    
  }
  
  return(res)
}

# Calculate the loglikelihood for the celda_CG model
cCG.calcLL = function(K, L, m.CP.by.S, n.CP.by.TS, n.by.G, n.by.TS, nG.by.TS, nS, nG, alpha, beta, delta, gamma) {
  nG.by.TS[nG.by.TS == 0] = 1
  nG = sum(nG.by.TS)
  
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
#' @param ... Additional parameters 
#' @export
calculateLoglikFromVariables.celda_CG = function(counts, s, z, y, K, L, alpha, beta, delta, gamma, ...) {
  
  ## Calculate for "Theta" component
  z = factor(z, 1:K)
  m = table(z, s)
  ns = ncol(m)
  
  a = ns * lgamma(K * alpha)
  b = sum(lgamma(m + alpha))
  c = -ns * K * lgamma(alpha)
  d = -sum(lgamma(colSums(m + alpha)))
  
  theta.ll = a + b + c + d
  
  
  ## Calculate for "Phi" component
  n.TS.by.C = rowsum.y(counts, y=y, L=L)
  n.CP.by.TS = rowsum.z(n.TS.by.C, z=z, K=K)
    
  a = K * lgamma(L * beta)
  b = sum(lgamma(n.CP.by.TS + beta))
  c = -K * L * lgamma(beta)
  d = -sum(lgamma(rowSums(n.CP.by.TS + beta)))
  
  phi.ll = a + b + c + d
  
  ## Calculate for "Psi" component
  n.by.G = rowSums(counts)
  n.by.TS = as.numeric(rowsum.y(matrix(n.by.G,ncol=1), y=y, L=L))

  nG.by.TS = table(factor(y, levels=1:L))
  nG.by.TS[nG.by.TS == 0] = 1
  nG = sum(nG.by.TS)
  
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


#' Calculates the conditional probability of each cell belong to each cluster given all other cluster assignments
#'
#' @param counts The original count matrix used in the model
#' @param celda.mod A model returned from the 'celda_CG' function
#' @param log If FALSE, then the normalized conditional probabilities will be returned. If TRUE, then the unnormalized log probabilities will be returned.  
#' @return A list containging a matrix for the conditional cell cluster probabilities. 
#' @export
clusterProbability.celda_CG = function(counts, celda.mod, log=FALSE) {
  
  s = celda.mod$sample.label
  z = celda.mod$z
  K = celda.mod$K  
  y = celda.mod$y
  L = celda.mod$L
  alpha = celda.mod$alpha
  delta = celda.mod$delta
  beta = celda.mod$beta
  gamma = celda.mod$gamma

  nS = length(unique(s))
  nG = nrow(counts)
  nM = ncol(counts)
  m.CP.by.S = matrix(as.integer(table(factor(z, levels=1:K), s)), ncol=nS)
  n.TS.by.C = rowsum.y(counts, y=y, L=L)
  n.CP.by.TS = rowsum.z(n.TS.by.C, z=z, K=K)
  n.CP = as.integer(rowSums(n.CP.by.TS))
  n.by.G = as.integer(rowSums(counts))
  n.by.C = as.integer(colSums(counts))
  n.by.TS = as.integer(rowsum.y(matrix(n.by.G,ncol=1), y=y, L=L))
  nG.by.TS = as.integer(table(factor(y, 1:L)))
  n.CP.by.G = rowsum.z(counts, z=z, K=K)

  ## Gibbs sampling for each cell
  n.TS.by.CP = t(n.CP.by.TS)
  next.z = cCG.calcGibbsProbZ(m.CP.by.S=m.CP.by.S, n.TS.by.CP=n.TS.by.CP, n.TS.by.C=n.TS.by.C, n.CP=n.CP, n.by.C=n.by.C, z=z, s=s, L=L, K=K, nM=nM, alpha=alpha, beta=beta)
  z.prob = t(next.z$probs)

  ## Gibbs sampling for each gene
  next.y = cCG.calcGibbsProbY(n.CP.by.TS=n.CP.by.TS, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, n.CP.by.G=n.CP.by.G, n.by.G=n.by.G, y=y, nG=nG, L=L, beta=beta, delta=delta, gamma=gamma)
  y.prob = t(next.y$probs)

  if(!isTRUE(log)) {
    z.prob = normalizeLogProbs(z.prob)
    y.prob = normalizeLogProbs(y.prob)
  }
       
  return(list(z.probability=z.prob, y.probability=y.prob))
}


#' @export
calculatePerplexity.celda_CG = function(counts, celda.mod, precision=128) {
  if (!compareCountMatrix(counts, celda.mod)) {
    warning("Provided count matrix was not used to generate the provided celda model.")
  }
  
  factorized = factorizeMatrix(celda.mod, counts, "posterior")
  theta = log(factorized$posterior$sample.states)
  phi   = factorized$posterior$population.states
  psi   = factorized$posterior$gene.states
  sl = celda.mod$sample.label
  
  gene.by.pop.prob = log(psi %*% phi)
  inner.log.prob = (t(gene.by.pop.prob) %*% counts) + theta[, sl]  
  inner.log.prob = Rmpfr::mpfr(inner.log.prob, precision)
  inner.log.prob.exp = exp(inner.log.prob)
  
  log.px = 0
  for(i in 1:ncol(inner.log.prob.exp)) {
    log.px = log.px + Rmpfr::asNumeric(log(sum(inner.log.prob.exp[, i])))
  }
  
  perplexity = exp(-(log.px/sum(counts)))
  return(perplexity)
}


reorder.celda_CG = function(counts,res){
  # Reorder K
  if(res$K > 2) {
    res$z = as.integer(as.factor(res$z))
    fm <- factorizeMatrix(counts = counts, celda.mod = res, type="posterior")
    unique.z = sort(unique(res$z))
    d <- cosineDist(fm$posterior$population.states[,unique.z])
    h <- hclust(d, method = "complete")
    
    res <- recodeClusterZ(res, from = h$order, to = 1:length(h$order))
  }  
  
  # Reorder L
  if(res$L > 2) {
    res$y = as.integer(as.factor(res$y))
    fm <- factorizeMatrix(counts = counts, celda.mod = res, type="posterior")
    unique.y = sort(unique(res$y))
    cs <- prop.table(t(fm$posterior$population.states[unique.y,]), 2)
    d <- cosineDist(cs)
    h <- hclust(d, method = "complete")
    
    res <- recodeClusterY(res, from = h$order, to = 1:length(h$order))
  }
  return(res)
}


#' finalClusterAssignment for the celda Cell and Gene clustering model
#' @param celda.mod A celda model object of "Celda_CG"
#' @export
finalClusterAssignment.celda_CG = function(celda.mod) {
  return(list(z=celda.mod$z, y=celda.mod$y))
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


#' celdaHeatmap for celda Cell and Gene clustering model
#' @param celda.mod A celda model object of "Celda_CG"
#' @param counts A count matrix
#' @param ... extra parameters passed onto renderCeldaHeatmap
#' @export
celdaHeatmap.celda_CG = function(celda.mod, counts, ...) {
  renderCeldaHeatmap(counts, z=celda.mod$z, y=celda.mod$y, ...)
}
