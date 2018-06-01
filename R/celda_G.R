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


#' celda Gene Clustering Model
#'
#' Provides cluster assignments for all genes in a provided single-cell 
#' sequencing count matrix, using the celda Bayesian hierarchical model.
#' 
#' @param counts A numeric count matrix
#' @param L The number of clusters to generate
#' @param beta The Dirichlet distribution parameter for Phi; adds a pseudocount to each transcriptional state within each cell.
#' @param delta The Dirichlet distribution parameter for Eta; adds a gene pseudocount to the numbers of genes each state. Default to 1.
#' @param gamma The Dirichlet distribution parameter for Psi; adds a pseudocount to each gene within each transcriptional state.
#' @param stop.iter Number of iterations without improvement in the log likelihood to stop the Gibbs sampler. Default 10.
#' @param max.iter Maximum iterations of Gibbs sampling to perform regardless of convergence. Default 200.
#' @param split.on.iter On every 'split.on.iter' iteration, a heuristic will be applied to determine if a gene/cell cluster should be reassigned and another gene/cell cluster should be split into two clusters. Default 10.
#' @param split.on.last After the the chain has converged according to 'stop.iter', a heuristic will be applied to determine if a gene/cell cluster should be reassigned and another gene/cell cluster should be split into two clusters. If a split occurs, then 'stop.iter' will be reset. Default TRUE.
#' @param count.checksum An MD5 checksum for the provided counts matrix
#' @param seed Parameter to set.seed() for random number generation.
#' @param y.init Initial values of y. If NULL, y will be randomly sampled. Default NULL.
#' @param logfile The name of the logfile to redirect messages to.
#' @keywords LDA gene clustering gibbs
#' @export
celda_G = function(counts, L, beta=1, delta=1, gamma=1,
					stop.iter=10, max.iter=200, split.on.iter=10, split.on.last=TRUE,
					count.checksum=NULL, seed=12345, 
					y.init=NULL, logfile=NULL) {

  ## Error checking and variable processing
  counts = processCounts(counts)  
  counts.t = t(counts)
  
  if(is.null(count.checksum)) {
    count.checksum = digest::digest(counts, algo="md5")
  }

  ## Randomly select z and y or set z/y to supplied initial values
  y = initialize.cluster(L, nrow(counts), initial = y.init, fixed = NULL, seed=seed)
  y.best = y  

  ## Calculate counts one time up front
  p = cG.decomposeCounts(counts=counts, y=y, L=L)
  n.C.by.TS = p$n.C.by.TS
  n.by.G = p$n.by.G
  n.by.TS = p$n.by.TS
  nG.by.TS = p$nG.by.TS
  nM = p$nM
  nG = p$nG
  rm(p)

  set.seed(seed)
  logMessages(date(), "... Starting celda_G to cluster genes.", logfile=logfile, append=FALSE)
  
  ## Calculate initial log likelihood
  ll <- cG.calcLL(n.C.by.TS=n.C.by.TS, n.by.TS=n.by.TS, n.by.G=n.by.G, nG.by.TS=nG.by.TS, nM=nM, nG=nG, L=L, beta=beta, delta=delta, gamma=gamma)

  iter <- 1L
  num.iter.without.improvement = 0L
  do.gene.split = TRUE
  while(iter <= max.iter & num.iter.without.improvement <= stop.iter) {

	next.y = cG.calcGibbsProbY(counts.t=counts.t, n.C.by.TS=n.C.by.TS, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, n.by.G=n.by.G, y=y, nG=nG, L=L, beta=beta, delta=delta, gamma=gamma)
	n.C.by.TS = next.y$n.C.by.TS
	nG.by.TS = next.y$nG.by.TS
	n.by.TS = next.y$n.by.TS
	y = next.y$y
    
    ## Perform split on i-th iteration of no improvement in log likelihood
    if(L > 2 & (((iter == max.iter | num.iter.without.improvement == stop.iter) & isTRUE(split.on.last)) | (split.on.iter > 0 & iter %% split.on.iter == 0 & isTRUE(do.gene.split)))) {
      logMessages(date(), " ... Determining if any gene clusters should be split.", logfile=logfile, append=TRUE, sep="")
      res = split.each.y(counts=counts, y=y, L=L, y.prob=t(next.y$probs), beta=beta, delta=delta, gamma=gamma, LLFunction="calculateLoglikFromVariables.celda_G")
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
      n.C.by.TS = t(rowSumByGroup(counts, y=y, L=L))
      n.by.TS = as.integer(rowSumByGroup(matrix(n.by.G,ncol=1), y=y, L=L))
      nG.by.TS = as.integer(table(factor(y, 1:L)))
    }
     
    ## Calculate complete likelihood
    temp.ll <- cG.calcLL(n.C.by.TS=n.C.by.TS, n.by.TS=n.by.TS, n.by.G=n.by.G, nG.by.TS=nG.by.TS, nM=nM, nG=nG, L=L, beta=beta, delta=delta, gamma=gamma)
    if((all(temp.ll > ll)) | iter == 1) {
      y.best = y
      ll.best = temp.ll
      num.iter.without.improvement = 1L
    } else {  
      num.iter.without.improvement = num.iter.without.improvement + 1L   
    }
    ll <- c(ll, temp.ll)

    logMessages(date(), " ... Completed iteration: ", iter, " | logLik: ", temp.ll, logfile=logfile, append=TRUE, sep="")

    iter <- iter + 1L
  }
    
  names = list(row=rownames(counts), column=colnames(counts))  

  result = list(y=y.best, completeLogLik=ll, 
                finalLogLik=ll.best, L=L, beta=beta, delta=delta, gamma=gamma,
                count.checksum=count.checksum, seed=seed, names=names)
  class(result) = "celda_G"
  result = reorder.celda_G(counts = counts, res = result)
  
  return(result)
}


# Calculate Log Likelihood For Single Set of Cluster Assignments (Gene Clustering)
#
# This function calculates the log-likelihood of a given set of cluster assigments for the samples
# represented in the provided count matrix.
# 
# 
# @param n.TS.by.C Number of counts in each Transcriptional State per Cell 
# @param n.by.TS Number of counts per Transcriptional State
# @param nG.by.TS Number of genes in each Transcriptional State
# @param nG.in.Y  Number of genes in each of the cell cluster
# @param gamma The Dirichlet distribution parameter for Psi; adds a pseudocount to each gene within each transcriptional state.
# @param delta The Dirichlet distribution parameter for Eta; adds a gene pseudocount to the numbers of genes each state.
# @param beta Vector of non-zero concentration parameters for cluster <-> gene assignment Dirichlet distribution
# @keywords log likelihood
cG.calcGibbsProbY = function(counts.t, n.C.by.TS, n.by.TS, nG.by.TS, n.by.G, y, L, nG, beta, delta, gamma, do.sample=TRUE) {

  ## Set variables up front outside of loop
  probs = matrix(NA, ncol=nG, nrow=L)
  temp.nG.by.TS = nG.by.TS 
  temp.n.by.TS = n.by.TS 
#  temp.n.C.by.TS = n.C.by.TS

  ix <- sample(1:nG)
  for(i in ix) {
	  
	## Subtract current gene counts from matrices
	nG.by.TS[y[i]] = nG.by.TS[y[i]] - 1L
	n.by.TS[y[i]] = n.by.TS[y[i]] - n.by.G[i]
	
	n.C.by.TS_1 = n.C.by.TS
	n.C.by.TS_1[,y[i]] = n.C.by.TS_1[,y[i]] - counts.t[,i]
    n.C.by.TS_1 = .colSums(lgamma(n.C.by.TS_1 + beta), nrow(n.C.by.TS), ncol(n.C.by.TS))
    
	n.C.by.TS_2 = n.C.by.TS
	other.ix = (1:L)[-y[i]]
	n.C.by.TS_2[,other.ix] = n.C.by.TS_2[,other.ix] + counts.t[,i]
    n.C.by.TS_2 = .colSums(lgamma(n.C.by.TS_2 + beta), nrow(n.C.by.TS), ncol(n.C.by.TS))
    	
	## Calculate probabilities for each state
	for(j in 1:L) {
	
	  temp.nG.by.TS = nG.by.TS 
	  temp.n.by.TS = n.by.TS 
#	  temp.n.C.by.TS = n.C.by.TS
	
	  temp.nG.by.TS[j] = temp.nG.by.TS[j] + 1L      
	  temp.n.by.TS[j] = temp.n.by.TS[j] + n.by.G[i]
#	  temp.n.C.by.TS[,j] = temp.n.C.by.TS[,j] + counts.t[,i]


      pseudo.nG.by.TS = temp.nG.by.TS
	  pseudo.nG.by.TS[temp.nG.by.TS == 0L] = 1L
	  pseudo.nG = sum(pseudo.nG.by.TS)

	  other.ix = (1:L)[-j]
      probs[j,i] <- 	sum(lgamma(pseudo.nG.by.TS + gamma)) -
						sum(lgamma(sum(pseudo.nG.by.TS + gamma))) +
						sum(n.C.by.TS_1[other.ix]) + n.C.by.TS_2[j] +
#						sum(lgamma(temp.n.C.by.TS + beta)) +
						sum(lgamma(pseudo.nG.by.TS * delta)) -
						(pseudo.nG * lgamma(delta)) -
						sum(lgamma(temp.n.by.TS + (pseudo.nG.by.TS * delta)))

    }

	## Sample next state and add back counts
	prev.y = y[i]
	if(isTRUE(do.sample)) y[i] = sample.ll(probs[,i])
	
	if(prev.y != y[i]) {
	  n.C.by.TS[,prev.y] = n.C.by.TS[,prev.y] - counts.t[,i]
	  n.C.by.TS[,y[i]] = n.C.by.TS[,y[i]] + counts.t[,i]
	}
	
	nG.by.TS[y[i]] = nG.by.TS[y[i]] + 1L
	n.by.TS[y[i]] = n.by.TS[y[i]] + n.by.G[i]
	#n.C.by.TS[,y[i]] = n.C.by.TS[,y[i]] + counts.t[,i]
  }
  
  return(list(n.C.by.TS=n.C.by.TS, nG.by.TS=nG.by.TS, n.by.TS=n.by.TS, y=y, probs=probs))
}


#' Simulate cells from the gene clustering generative model
#'
#' Generate a simulated count matrix, based off a generative distribution whose 
#' parameters can be provided by the user.
#' 
#' @param model Celda model to use for simulation. One of 'available_models'. 
#' @param C The number of cells
#' @param L The number of transcriptional states
#' @param N.Range Vector of length 2 given the range (min,max) of number of counts for each cell to be randomly generated from the uniform distribution
#' @param G The number of genes for which to simulate counts
#' @param beta The Dirichlet distribution parameter for Phi; adds a pseudocount to each transcriptional state within each cell
#' @param delta The Dirichlet distribution parameter for Psi; adds a pseudocount to each gene within each transcriptional state
#' @param gamma The Dirichlet distribution parameter for Psi; adds a pseudocount to each gene within each transcriptional state.
#' @param seed Parameter to set.seed() for random number generation
#' @param ... Other arguments
#' @export
simulateCells.celda_G = function(model, C=100, N.Range=c(500,5000),  G=1000, 
                                 L=5, beta=1, gamma=1, delta=1, seed=12345, ...) {
  set.seed(seed)
  eta = rdirichlet(1, rep(gamma, L))
  
  y = sample(1:L, size=G, prob=eta, replace=TRUE)
  if(length(table(y)) < L) {
    stop("Some states did not receive any genes after sampling. Try increasing G and/or setting gamma > 1.")
  }
  y = reorder.label.by.size(y, L)$new.labels
  
  psi = matrix(0, nrow=G, ncol=L)
  for(i in 1:L) {
    ind = y == i
    psi[ind,i] = rdirichlet(1, rep(delta, sum(ind)))
  }
  
  phi = rdirichlet(C, rep(beta, L))
  
  ## Select number of transcripts per cell
  nN = sample(N.Range[1]:N.Range[2], size=C, replace=TRUE)
  
  ## Select transcript distribution for each cell
  cell.counts = matrix(0, nrow=G, ncol=C)
  for(i in 1:C) {
    cell.dist = stats::rmultinom(1, size=nN[i], prob=phi[i,])
    for(j in 1:L) {
      cell.counts[,i] = cell.counts[,i] + stats::rmultinom(1, size=cell.dist[j], prob=psi[,j])
    }
  }
  
  ## Ensure that there are no all-0 rows in the counts matrix, which violates a celda modeling
  ## constraint (columns are guarnteed at least one count):
  zero.row.idx = which(rowSums(cell.counts) == 0)
  if (length(zero.row.idx > 0)) {
    cell.counts = cell.counts[-zero.row.idx, ]
    y = y[-zero.row.idx]
  }

  rownames(cell.counts) = paste0("Gene_", 1:nrow(cell.counts))
  colnames(cell.counts) = paste0("Cell_", 1:ncol(cell.counts))
  
  ## Peform reordering on final Z and Y assigments:
  names = list(row=rownames(cell.counts), column=colnames(cell.counts))
  result = list(y=y, completeLogLik=NULL, 
                finalLogLik=NULL, L=L, 
                beta=beta, delta=delta, gamma=gamma, seed=seed, 
                names=names, count.checksum=NULL)
  class(result) = "celda_G" 
  result = reorder.celda_G(counts = cell.counts, res = result)  
  
  return(list(y=result$y, counts=processCounts(cell.counts), L=L, beta=beta, delta=delta, gamma=gamma, phi=phi, psi=psi, eta=eta, seed=seed))
}


#' Generate factorized matrices showing each feature's influence on the celda_G model clustering 
#' 
#' @param counts A numeric count matrix
#' @param celda.mod Object return from celda_C function
#' @param type A character vector containing one or more of "counts", "proportions", or "posterior". "counts" returns the raw number of counts for each entry in each matrix. "proportions" returns the counts matrix where each vector is normalized to a probability distribution. "posterior" returns the posterior estimates which include the addition of the Dirichlet concentration parameter (essentially as a pseudocount).
#' @export
factorizeMatrix.celda_G = function(counts, celda.mod, type=c("counts", "proportion", "posterior")) {

  L = celda.mod$L
  y = celda.mod$y
  beta = celda.mod$beta
  delta = celda.mod$delta
  gamma = celda.mod$gamma
  
  p = cG.decomposeCounts(counts=counts, y=y, L=L)
  n.TS.by.C = t(p$n.C.by.TS)
  n.by.G = p$n.by.G
  n.by.TS = p$n.by.TS
  nG.by.TS = p$nG.by.TS
  nM = p$nM
  nG = p$nG
  rm(p)
  
  nG.by.TS[nG.by.TS == 0] = 1  
  n.G.by.TS = matrix(0, nrow=length(y), ncol=L)
  n.G.by.TS[cbind(1:nG,y)] = n.by.G

  L.names = paste0("L", 1:L)
  colnames(n.TS.by.C) = celda.mod$names$column
  rownames(n.TS.by.C) = L.names
  colnames(n.G.by.TS) = L.names
  rownames(n.G.by.TS) = celda.mod$names$row
  names(nG.by.TS) = L.names
  
  counts.list = c()
  prop.list = c()
  post.list = c()
  res = list()
  
  if(any("counts" %in% type)) {
    counts.list = list(cell.states=n.TS.by.C, gene.states=n.G.by.TS, gene.distribution=nG.by.TS)
    res = c(res, list(counts=counts.list))
  }
  if(any("proportion" %in% type)) {
    ## Need to avoid normalizing cell/gene states with zero cells/genes
    unique.y = sort(unique(y))
    temp.n.G.by.TS = n.G.by.TS
    temp.n.G.by.TS[,unique.y] = normalizeCounts(temp.n.G.by.TS[,unique.y], scale.factor=1)
    temp.nG.by.TS = nG.by.TS/sum(nG.by.TS)
    
    prop.list = list(cell.states = normalizeCounts(n.TS.by.C, scale.factor=1),
    							  gene.states = temp.n.G.by.TS, gene.distribution=temp.nG.by.TS)
    res = c(res, list(proportions=prop.list))
  }
  if(any("posterior" %in% type)) {
  
    gs = n.G.by.TS
    gs[cbind(1:nG,y)] = gs[cbind(1:nG,y)] + delta
    gs = normalizeCounts(gs, scale.factor=1)
    temp.nG.by.TS = (nG.by.TS + gamma)/sum(nG.by.TS + gamma)
    
    post.list = list(cell.states = normalizeCounts(n.TS.by.C + beta, scale.factor=1),
    						    gene.states = gs, gene.distribution=temp.nG.by.TS)
    res = c(res, posterior = list(post.list))						    
  }
  
  return(res)
}


# Calculate log-likelihood of celda_CG model
cG.calcLL = function(n.C.by.TS, n.by.TS, n.by.G, nG.by.TS, nM, nG, L, beta, delta, gamma) {
  
  nG.by.TS[nG.by.TS == 0] = 1
  nG = sum(nG.by.TS)
  
  ## Calculate for "Phi" component
  a <- nM * lgamma(L * beta)
  b <- sum(lgamma(n.C.by.TS + beta))
  c <- -nM * L * lgamma(beta)
  d <- -sum(lgamma(rowSums(n.C.by.TS + beta)))

  phi.ll <- a + b + c + d

  ## Calculate for "Psi" component
  a <- sum(lgamma(nG.by.TS * delta))
  b <- sum(lgamma(n.by.G + delta))
  c <- -nG * lgamma(delta)
  d <- -sum(lgamma(n.by.TS + (nG.by.TS * delta)))
  
  psi.ll <- a + b + c + d

  ## Calculate for "Eta" component
  a <- lgamma(L * gamma)
  b <- sum(lgamma(nG.by.TS + gamma))
  c <- -L * lgamma(gamma)
  d <- -sum(lgamma(sum(nG.by.TS + gamma)))

  eta.ll <- a + b + c + d

  final <- phi.ll + psi.ll + eta.ll
  return(final)
}


#' Calculate the celda_G log likelihood for user-provided cluster assignments
#'
#' This function calculates the log likelihood of each clustering of genes generated
#' over multiple iterations of Gibbs sampling.
#' 
#' @param counts A numeric count matrix
#' @param y A numeric vector of gene cluster assignments
#' @param L The number of clusters being considered
#' @param beta The Dirichlet distribution parameter for Phi; adds a pseudocount to each transcriptional state within each cell
#' @param delta The Dirichlet distribution parameter for Eta; adds a gene pseudocount to the numbers of genes each state
#' @param gamma The Dirichlet distribution parameter for Psi; adds a pseudocount to each gene within each transcriptional state
#' @param ... Additional parameters
#' @keywords log likelihood
#' @return The log likelihood of the provided cluster assignment, as calculated by the celda_G likelihood function
#' @export
calculateLoglikFromVariables.celda_G = function(counts, y, L, beta, delta, gamma) {

  p = cG.decomposeCounts(counts=counts, y=y, L=L)
  final <- cG.calcLL(n.C.by.TS=p$n.C.by.TS, n.by.TS=p$n.by.TS, n.by.G=p$n.by.G, nG.by.TS=p$nG.by.TS, nM=p$nM, nG=p$nG, L=L, beta=beta, delta=delta, gamma=gamma)
  
  return(final)
}


#' Takes raw counts matrix and converts it to a series of matrices needed for log likelihood calculation
#' @param counts A numeric count matrix
#' @param y A numeric vector of gene cluster assignments
#' @param L The number of clusters being considered
cG.decomposeCounts = function(counts, y, L) {

  n.C.by.TS = t(rowSumByGroup(counts, y=y, L=L))
  n.by.G = as.integer(rowSums(counts))
  n.by.TS = as.integer(rowSumByGroup(matrix(n.by.G,ncol=1), y=y, L=L))
  nG.by.TS = as.integer(table(factor(y, 1:L)))
  nM = ncol(counts)
  nG = nrow(counts)
  
  return(list(n.C.by.TS=n.C.by.TS, n.by.G=n.by.G, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, nM=nM, nG=nG))
}


#' Calculates the conditional probability of each cell belong to each cluster given all other cluster assignments
#'
#' @param celda.mod A model returned from the 'celda_G' function
#' @param counts The original count matrix used in the model
#' @param log If FALSE, then the normalized conditional probabilities will be returned. If TRUE, then the unnormalized log probabilities will be returned.  
#' @param ... Other arguments
#' @return A list containging a matrix for the conditional cell cluster probabilities. 
#' @export
clusterProbability.celda_G = function(celda.mod, counts, log=FALSE, ...) {

  y = celda.mod$y
  L = celda.mod$L
  delta = celda.mod$delta
  beta = celda.mod$beta
  gamma = celda.mod$gamma
  
  ## Calculate counts one time up front
  p = cG.decomposeCounts(counts=counts, y=y, L=L)
  next.y = cG.calcGibbsProbY(counts.t=t(counts), n.C.by.TS=p$n.C.by.TS, n.by.TS=p$n.by.TS, nG.by.TS=p$nG.by.TS, n.by.G=p$n.by.G, y=y, nG=p$nG, L=L, beta=beta, delta=delta, gamma=gamma, do.sample=FALSE)  
  y.prob = t(next.y$probs)
  
  if(!isTRUE(log)) {
    y.prob = normalizeLogProbs(y.prob)
  }
  
  return(list(y.probability=y.prob))
}


#' @export
calculatePerplexity.celda_G = function(counts, celda.mod, precision=128) {
 
  factorized = factorizeMatrix(counts = counts, celda.mod = celda.mod, "posterior")
  phi <- factorized$posterior$gene.states
  psi <- factorized$posterior$cell.states
  eta <- factorized$posterior$gene.distribution
  nG.by.TS = factorized$counts$gene.distribution
  
  eta.prob = log(eta) * nG.by.TS
  gene.by.cell.prob = log(phi %*% psi) 
  log.px = sum(eta.prob) + sum(gene.by.cell.prob * counts)
  
  perplexity = exp(-(log.px/sum(counts)))
  return(perplexity)
}


reorder.celda_G = function(counts, res) {
  if(res$L > 2) {
    res$y = as.integer(as.factor(res$y))
    fm <- factorizeMatrix(counts = counts, celda.mod = res)
    unique.y = sort(unique(res$y))
    cs = prop.table(t(fm$posterior$cell.states[unique.y,]), 2)
    d <- cosineDist(cs)
    h <- hclust(d, method = "complete")
    res <- recodeClusterY(res, from = h$order, to = 1:length(h$order))
  }  
  return(res)
}


#' finalClusterAssignment for celda Gene clustering model
#' @param celda.mod A celda model object of class "celda_G"
#' @export
finalClusterAssignment.celda_G = function(celda.mod) {
  return(celda.mod$y)
}


#' getK for celda Gene clustering model
#' @param celda.mod A celda model object of class "celda_G"
#' @export
getK.celda_G = function(celda.mod) { return(NA) }


#' getL for celda Gene clustering model
#' @param celda.mod A celda model object of class "celda_G"
#' @export
getL.celda_G = function(celda.mod) {
  return(celda.mod$L)
}


#' celdaHeatmap for celda Gene clustering model
#' @param celda.mod A celda model object of class "celda_G"
#' @param counts A numeric count matrix
#' @param ... extra parameters passed onto renderCeldaHeatmap
#' @export
celdaHeatmap.celda_G = function(celda.mod, counts, ...) {
  renderCeldaHeatmap(counts, y=celda.mod$y, ...)
}
