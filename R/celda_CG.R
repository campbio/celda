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
#' @param algorithm Use 'EM' or 'Gibbs' sampling for clustering of cells into subpopulations. EM is much faster for larger datasets. Default 'EM'.
#' @param stop.iter Number of iterations without improvement in the log likelihood to stop the Gibbs sampler. Default 10.
#' @param max.iter Maximum iterations of Gibbs sampling to perform regardless of convergence. Default 200.
#' @param split.on.iter On every 'split.on.iter' iteration, a heuristic will be applied to determine if a gene/cell cluster should be reassigned and another gene/cell cluster should be split into two clusters. Default 10.
#' @param split.on.last After the the chain has converged according to 'stop.iter', a heuristic will be applied to determine if a gene/cell cluster should be reassigned and another gene/cell cluster should be split into two clusters. If a split occurs, then 'stop.iter' will be reset. Default TRUE.
#' @param seed Parameter to set.seed() for random number generation. Default 12345. 
#' @param nchains Number of random z/y initializations. Default 3. 
#' @param count.checksum An MD5 checksum for the provided counts matrix. Default NULL.
#' @param z.init Initial values of z. If NULL, z will be randomly sampled. Default NULL.
#' @param y.init Initial values of y. If NULL, y will be randomly sampled. Default NULL.
#' @param logfile The name of the logfile to redirect messages to.
#' @export
celda_CG = function(counts, sample.label=NULL, K, L,
                    alpha=1, beta=1, delta=1, gamma=1, 
                    algorithm = c("EM", "Gibbs"), 
                    stop.iter = 10, max.iter=200, split.on.iter=10, split.on.last=TRUE,
                    seed=12345, nchains=3, count.checksum=NULL,
                    z.init = NULL, y.init = NULL, logfile=NULL) {
 
  if(is.null(count.checksum)) {
    count.checksum = digest::digest(counts, algo="md5")
  }
  counts = processCounts(counts)
    
  s = processSampleLabels(sample.label, ncol(counts))
  if (is.null(sample.label)) {
    sample.label = s
  }

  algorithm <- match.arg(algorithm)
  algorithm.fun <- ifelse(algorithm == "Gibbs", "cC.calcGibbsProbZ", "cC.calcEMProbZ")
  
  if(is.null(count.checksum)) {
    count.checksum = digest::digest(counts, algo="md5")
  }
  
  all.seeds = seed:(seed + nchains - 1)
  
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=FALSE)  
  logMessages("Celda_CG: Clustering cells and genes.", logfile=logfile, append=FALSE)
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=FALSE)  
  
  best.result = NULL  
  for(i in seq_along(all.seeds)) { 

	## Randomly select z and y or set z/y to supplied initial values
    current.seed = all.seeds[i]	
	z = initialize.cluster(K, ncol(counts), initial = z.init, fixed = NULL, seed=current.seed)
	y = initialize.cluster(L, nrow(counts), initial = y.init, fixed = NULL, seed=current.seed)
	z.best = z
	y.best = y  
  
	## Calculate counts one time up front
	p = cCG.decomposeCounts(counts, s, z, y, K, L)
	m.CP.by.S = p$m.CP.by.S
	n.TS.by.C = p$n.TS.by.C
	n.TS.by.CP = p$n.TS.by.CP
	n.CP = p$n.CP
	n.by.G = p$n.by.G
	n.by.C = p$n.by.C
	n.by.TS = p$n.by.TS
	nG.by.TS = p$nG.by.TS
	n.G.by.CP = p$n.G.by.CP
	nM = p$nM
	nG = p$nG
	nS = p$nS
	rm(p)
  
	ll = cCG.calcLL(K=K, L=L, m.CP.by.S=m.CP.by.S, n.TS.by.CP=n.TS.by.CP, n.by.G=n.by.G, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, nS=nS, nG=nG, alpha=alpha, beta=beta, delta=delta, gamma=gamma)

    logMessages(date(), ".. Starting chain", i, "with seed", current.seed, logfile=logfile, append=FALSE)

	set.seed(current.seed)
	iter = 1L
	num.iter.without.improvement = 0L
	do.cell.split = TRUE
	do.gene.split = TRUE  
	while(iter <= max.iter & num.iter.without.improvement <= stop.iter) {

	  ## Gibbs or EM sampling for each cell
	  next.z = do.call(algorithm.fun, list(counts=n.TS.by.C, m.CP.by.S=m.CP.by.S, n.G.by.CP=n.TS.by.CP, n.CP=n.CP, n.by.C=n.by.C, z=z, s=s, K=K, nG=L, nM=nM, alpha=alpha, beta=beta))
	  m.CP.by.S = next.z$m.CP.by.S
	  n.TS.by.CP = next.z$n.G.by.CP
	  n.CP = next.z$n.CP
	  n.G.by.CP = colSumByGroupChange(counts, n.G.by.CP, next.z$z, z, K)
	  z = next.z$z
		
	  ## Gibbs sampling for each gene
	  next.y = cG.calcGibbsProbY(counts=n.G.by.CP, n.TS.by.C=n.TS.by.CP, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, n.by.G=n.by.G, y=y, L=L, nG=nG, beta=beta, delta=delta, gamma=gamma)
	  n.TS.by.CP = next.y$n.TS.by.C
	  nG.by.TS = next.y$nG.by.TS
	  n.by.TS = next.y$n.by.TS
	  n.TS.by.C = rowSumByGroupChange(counts, n.TS.by.C, next.y$y, y, L)
	  y = next.y$y
	
		
	  ## Perform split on i-th iteration defined by split.on.iter
	  if(K > 2 & (((iter == max.iter | num.iter.without.improvement == stop.iter) & isTRUE(split.on.last)) | (split.on.iter > 0 & iter %% split.on.iter == 0 & isTRUE(do.cell.split)))) {
		logMessages(date(), " .... Determining if any cell clusters should be split.", logfile=logfile, append=TRUE, sep="")
		res = cCG.splitZ(counts, m.CP.by.S, n.TS.by.C, n.TS.by.CP, n.by.G, n.by.TS, nG.by.TS, s, z, K, L, nS, nG, alpha, beta, delta, gamma, z.prob=t(next.z$probs), max.clusters.to.try=10, min.cell=3)
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
		m.CP.by.S = res$m.CP.by.S
		n.TS.by.CP = res$n.TS.by.CP
		n.CP = res$n.CP
		n.G.by.CP = colSumByGroup(counts, group=z, K=K)
	  }  
	  if(L > 2 & (((iter == max.iter | num.iter.without.improvement == stop.iter) & isTRUE(split.on.last)) | (split.on.iter > 0 & iter %% split.on.iter == 0 & isTRUE(do.gene.split)))) {
		logMessages(date(), " .... Determining if any gene clusters should be split.", logfile=logfile, append=TRUE, sep="")
		res = cCG.splitY(counts, y, m.CP.by.S, n.G.by.CP, n.TS.by.C, n.TS.by.CP, n.by.G, n.by.TS, nG.by.TS, n.CP, s, z, K, L, nS, nG, alpha, beta, delta, gamma, y.prob=t(next.y$probs), max.clusters.to.try=10, min.cell=3)
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
		n.TS.by.CP = res$n.TS.by.CP
		n.by.TS = res$n.by.TS
		nG.by.TS = res$nG.by.TS
		n.TS.by.C = rowSumByGroup(counts, group=y, L=L)	  
	  }      

	  ## Calculate complete likelihood
	  temp.ll = cCG.calcLL(K=K, L=L, m.CP.by.S=m.CP.by.S, n.TS.by.CP=n.TS.by.CP, n.by.G=n.by.G, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, nS=nS, nG=nG, alpha=alpha, beta=beta, delta=delta, gamma=gamma)
	  if((all(temp.ll > ll)) | iter == 1) {
		z.best = z
		y.best = y
		ll.best = temp.ll
		num.iter.without.improvement = 1L
	  } else {  
		num.iter.without.improvement = num.iter.without.improvement + 1L   
	  }
	  ll = c(ll, temp.ll)
  
	  logMessages(date(), " .... Completed iteration: ", iter, " | logLik: ", temp.ll, logfile=logfile, append=TRUE, sep="")
	  iter = iter + 1L
	}
      
	names = list(row=rownames(counts), column=colnames(counts), 
				 sample=levels(sample.label))
  
	result = list(z=z.best, y=y.best, completeLogLik=ll, 
				  finalLogLik=ll.best, K=K, L=L, alpha=alpha, 
				  beta=beta, delta=delta, gamma=gamma, seed=current.seed, 
				  sample.label=sample.label, names=names,
				  count.checksum=count.checksum)
  
	class(result) = "celda_CG" 

    if(is.null(best.result) || result$finalLogLik > best.result$finalLogLik) {
      best.result = result
    }
    
    logMessages(date(), ".. Finished chain", i, "with seed", current.seed, logfile=logfile, append=FALSE)
  } 
  
  ## Peform reordering on final Z and Y assigments:
  best.result = reorder.celda_CG(counts = counts, res = best.result)
  return(best.result)
}




#' Simulate cells from the cell/gene clustering generative model
#' 
#' @param model Celda model to use for simulation. One of 'available_models'. 
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
#' @param ... Other arguments
#' @export
simulateCells.celda_CG = function(model, S=10, C.Range=c(50,100), N.Range=c(500,5000), 
                                  G=1000, K=3, L=10, alpha=1, beta=1, gamma=5, 
                                  delta=1, seed=12345, ...) {
  
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
    warning("Some transcriptional states did not receive any genes after sampling. Try increasing G and/or making gamma larger.")
    L = length(table(y))
    y = as.integer(as.factor(y))
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
  
  ## Ensure that there are no all-0 rows in the counts matrix, which violates a celda modeling
  ## constraint (columns are guarnteed at least one count):
  zero.row.idx = which(rowSums(cell.counts) == 0)
  if (length(zero.row.idx > 0)) {
    cell.counts = cell.counts[-zero.row.idx, ]
    y = y[-zero.row.idx]
  } 
 
  ## Assign gene/cell/sample names 
  rownames(cell.counts) = paste0("Gene_", 1:nrow(cell.counts))
  colnames(cell.counts) = paste0("Cell_", 1:ncol(cell.counts))
  cell.sample.label = paste0("Sample_", 1:S)[cell.sample.label]

  ## Peform reordering on final Z and Y assigments:
  names = list(row=rownames(cell.counts), column=colnames(cell.counts), 
               sample=unique(cell.sample.label))
  result = list(z=z, y=y, completeLogLik=NULL, 
                finalLogLik=NULL, K=K, L=L, alpha=alpha, 
                beta=beta, delta=delta, gamma=gamma, seed=seed, 
                sample.label=cell.sample.label, names=names,
                count.checksum=NULL)
  class(result) = "celda_CG" 
  cell.counts = processCounts(cell.counts)
  result = reorder.celda_CG(counts = cell.counts, res = result)
  
  return(list(z=result$z, y=result$y, sample.label=cell.sample.label, counts=cell.counts, K=K, L=L, C.Range=C.Range, N.Range=N.Range, S=S, alpha=alpha, beta=beta, gamma=gamma, delta=delta, theta=theta, phi=phi, psi=psi, eta=eta, seed=seed))
}


##' Generate factorized matrices showing each feature's influence on the celda_CG model clustering 
#' 
#' @param counts A numerix count matrix
#' @param celda.mod object returned from celda_CG function 
#' @param type one of the "counts", "proportion", or "posterior". 
#' @param validate.counts Whether to verify that the counts matrix provided was used to generate the results in celda.mod. Defaults to TRUE.
#' @return A list of factorized matrices, of the types requested by the user. NOTE: "population" state matrices are always returned in cell population (rows) x transcriptional states (cols).
#' @export 
factorizeMatrix.celda_CG = function(counts, celda.mod, 
                                    type=c("counts", "proportion", "posterior"),
                                    validate.counts = TRUE) {
  counts = processCounts(counts)
  if (validate.counts) { 
    compareCountMatrix(counts, celda.mod)
  }
  
  K = celda.mod$K
  L = celda.mod$L
  z = celda.mod$z
  y = celda.mod$y
  alpha = celda.mod$alpha
  beta = celda.mod$beta
  delta = celda.mod$delta
  gamma = celda.mod$gamma
  sample.label = celda.mod$sample.label
  s = processSampleLabels(sample.label, ncol(counts))
  
  ## Calculate counts one time up front
  p = cCG.decomposeCounts(counts, s, z, y, K, L)
  nS = p$nS
  nG = p$nG
  nM = p$nM
  m.CP.by.S = p$m.CP.by.S
  n.TS.by.C = p$n.TS.by.C
  n.TS.by.CP = p$n.TS.by.CP
  n.by.G = p$n.by.G
  n.by.TS = p$n.by.TS
  nG.by.TS = p$nG.by.TS
  nG.by.TS[nG.by.TS == 0] = 1


  n.G.by.TS = matrix(0, nrow=length(y), ncol=L)
  n.G.by.TS[cbind(1:nG,y)] = p$n.by.G

  L.names = paste0("L", 1:L)
  K.names = paste0("K", 1:K)
  colnames(n.TS.by.C) = celda.mod$names$column
  rownames(n.TS.by.C) = L.names
  colnames(n.G.by.TS) = L.names
  rownames(n.G.by.TS) = celda.mod$names$row
  rownames(m.CP.by.S) = K.names
  colnames(m.CP.by.S) = celda.mod$names$sample
  colnames(n.TS.by.CP) = K.names
  rownames(n.TS.by.CP) = L.names

  counts.list = c()
  prop.list = c()
  post.list = c()
  res = list()
    
  if(any("counts" %in% type)) {
    counts.list = list(sample.states = m.CP.by.S,
            				   population.states = n.TS.by.CP, 
            				   cell.states = n.TS.by.C,
            				   gene.states = n.G.by.TS,
            				   gene.distribution = nG.by.TS)
    res = c(res, list(counts=counts.list))
  }
  if(any("proportion" %in% type)) {

    ## Need to avoid normalizing cell/gene states with zero cells/genes
    unique.z = sort(unique(z))
    temp.n.TS.by.CP = n.TS.by.CP
    temp.n.TS.by.CP[,unique.z] = normalizeCounts(temp.n.TS.by.CP[,unique.z], scale.factor=1)

    unique.y = sort(unique(y))
    temp.n.G.by.TS = n.G.by.TS
    temp.n.G.by.TS[,unique.y] = normalizeCounts(temp.n.G.by.TS[,unique.y], scale.factor=1)
    temp.nG.by.TS = nG.by.TS/sum(nG.by.TS)
    
    prop.list = list(sample.states =  normalizeCounts(m.CP.by.S, scale.factor=1),
    				   population.states = temp.n.TS.by.CP, 
    				   cell.states = normalizeCounts(n.TS.by.C, scale.factor=1),
    				   gene.states = temp.n.G.by.TS, 
    				   gene.distribution = temp.nG.by.TS)
    res = c(res, list(proportions=prop.list))
  }
  if(any("posterior" %in% type)) {

    gs = n.G.by.TS
    gs[cbind(1:nG,y)] = gs[cbind(1:nG,y)] + delta
    gs = normalizeCounts(gs, scale.factor=1)
	temp.nG.by.TS = (nG.by.TS + gamma)/sum(nG.by.TS + gamma)
	
    post.list = list(sample.states = normalizeCounts(m.CP.by.S + alpha, scale.factor=1),
          				   population.states = normalizeCounts(n.TS.by.CP + beta, scale.factor=1), 
          				   gene.states = gs,
          				   gene.distribution = temp.nG.by.TS)
    res = c(res, posterior = list(post.list))						    
  }
  
  return(res)
}

# Calculate the loglikelihood for the celda_CG model
cCG.calcLL = function(K, L, m.CP.by.S, n.TS.by.CP, n.by.G, n.by.TS, nG.by.TS, nS, nG, alpha, beta, delta, gamma) {
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
  b = sum(lgamma(n.TS.by.CP + beta))
  c = -K * L * lgamma(beta)
  d = -sum(lgamma(colSums(n.TS.by.CP + beta)))
  
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
#' @param sample.label A vector indicating the sample label for each cell (column) in the count matrix
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
calculateLoglikFromVariables.celda_CG = function(counts, sample.label, z, y, K, L, alpha, beta, delta, gamma) {  
  s = processSampleLabels(sample.label, ncol(counts))
  p = cCG.decomposeCounts(counts, s, z, y, K, L)
  final = cCG.calcLL(K=K, L=L, m.CP.by.S=p$m.CP.by.S, n.TS.by.CP=p$n.TS.by.CP, n.by.G=p$n.by.G, n.by.TS=p$n.by.TS, nG.by.TS=p$nG.by.TS, nS=p$nS, nG=p$nG, alpha=alpha, beta=beta, delta=delta, gamma=gamma)
  return(final)
}


#' Takes raw counts matrix and converts it to a series of matrices needed for log likelihood calculation
#' @param counts A numeric count matrix
#' @param s An integer vector indicating the sample label for each cell (column) in the count matrix
#' @param z A numeric vector of cell cluster assignments
#' @param y A numeric vector of gene cluster assignments
#' @param K The number of cell clusters
#' @param L The number of gene clusters
cCG.decomposeCounts = function(counts, s, z, y, K, L) {
  nS = length(unique(s))
  m.CP.by.S = matrix(as.integer(table(factor(z, levels=1:K), s)), ncol=nS)
  n.TS.by.C = rowSumByGroup(counts, group=y, L=L)
  n.TS.by.CP = colSumByGroup(n.TS.by.C, group=z, K=K)
  n.CP = as.integer(colSums(n.TS.by.CP))
  n.by.G = as.integer(rowSums(counts))
  n.by.C = as.integer(colSums(counts))
  n.by.TS = as.integer(rowSumByGroup(matrix(n.by.G,ncol=1), group=y, L=L))
  nG.by.TS = tabulate(y, L)
  n.G.by.CP = colSumByGroup(counts, group=z, K=K)
  
  nG = nrow(counts)
  nM = ncol(counts)
  return(list(m.CP.by.S=m.CP.by.S, n.TS.by.C=n.TS.by.C, n.TS.by.CP=n.TS.by.CP, n.CP=n.CP, n.by.G=n.by.G, n.by.C=n.by.C, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, n.G.by.CP=n.G.by.CP, nM=nM, nG=nG, nS=nS))
}  



#' Calculates the conditional probability of each cell belong to each cluster given all other cluster assignments
#'
#' @param celda.mod A model returned from the 'celda_CG' function
#' @param counts The original count matrix used in the model
#' @param log If FALSE, then the normalized conditional probabilities will be returned. If TRUE, then the unnormalized log probabilities will be returned.  
#' @param ... Other arguments
#' @return A list containging a matrix for the conditional cell cluster probabilities. 
#' @export
clusterProbability.celda_CG = function(celda.mod, counts, log=FALSE, ...) {
  
  s = processSampleLabels(celda.mod$sample.label, ncol(counts))
  z = celda.mod$z
  K = celda.mod$K  
  y = celda.mod$y
  L = celda.mod$L
  alpha = celda.mod$alpha
  delta = celda.mod$delta
  beta = celda.mod$beta
  gamma = celda.mod$gamma

  p = cCG.decomposeCounts(counts, s, z, y, K, L)

  ## Gibbs sampling for each cell
  next.z = cC.calcGibbsProbZ(counts=p$n.TS.by.C, m.CP.by.S=p$m.CP.by.S, n.G.by.CP=p$n.TS.by.CP, n.CP=p$n.CP, n.by.C=p$n.by.C, z=z, s=s, K=K, nG=L, nM=p$nM, alpha=alpha, beta=beta, do.sample=FALSE)
  z.prob = t(next.z$probs)

  ## Gibbs sampling for each gene
  next.y = cG.calcGibbsProbY(counts=p$n.CP.by.G, n.TS.by.C=p$n.TS.by.CP, n.by.TS=p$n.by.TS, nG.by.TS=p$nG.by.TS, n.by.G=p$n.by.G, y=y, L=L, nG=nG, beta=beta, delta=delta, gamma=gamma, do.sample=FALSE)
  y.prob = t(next.y$probs)

  if(!isTRUE(log)) {
    z.prob = normalizeLogProbs(z.prob)
    y.prob = normalizeLogProbs(y.prob)
  }
       
  return(list(z.probability=z.prob, y.probability=y.prob))
}


#' @export
calculatePerplexity.celda_CG = function(counts, celda.mod, validate.counts) {
  
  factorized = factorizeMatrix(counts = counts, celda.mod = celda.mod, 
                               "posterior", validate.counts)
  theta = log(factorized$posterior$sample.states)
  phi   = factorized$posterior$population.states
  psi   = factorized$posterior$gene.states
  sl = celda.mod$sample.label
  eta <- factorized$posterior$gene.distribution
  nG.by.TS = factorized$counts$gene.distribution
  
  eta.prob = log(eta) * nG.by.TS
  gene.by.pop.prob = log(psi %*% phi)
  inner.log.prob = (t(gene.by.pop.prob) %*% counts) + theta[, sl]  
  
  log.px = sum(apply(inner.log.prob, 2, matrixStats::logSumExp)) + sum(eta.prob)
  perplexity = exp(-(log.px/sum(counts)))
  return(perplexity)
}


reorder.celda_CG = function(counts, res){
  # Reorder K
  if(res$K > 2 & isTRUE(length(unique(res$z)) > 1)) {
    res$z = as.integer(as.factor(res$z))
    fm <- factorizeMatrix(counts = counts, celda.mod = res, type="posterior",
                          validate.counts=FALSE)
    unique.z = sort(unique(res$z))
    d <- cosineDist(fm$posterior$population.states[,unique.z])
    h <- hclust(d, method = "complete")
    
    res <- recodeClusterZ(res, from = h$order, to = 1:length(h$order))
  }  
  
  # Reorder L
  if(res$L > 2 & isTRUE(length(unique(res$y)) > 1)) {
    res$y = as.integer(as.factor(res$y))
    fm <- factorizeMatrix(counts = counts, celda.mod = res, type="posterior",
                          validate.counts=FALSE)
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


#' Embeds cells in two dimensions using tSNE based on celda_CG results.
#' 
#' @param counts Counts matrix, should have cell name for column name and gene name for row name.
#' @param celda.mod Celda model to use for tsne. 
#' @param max.cells Integer; Maximum number of cells to plot. Cells will be randomly subsampled if ncol(conts) > max.cells. Larger numbers of cells requires more memory. Default 10000.
#' @param min.cluster.size Integer; Do not subsample cell clusters below this threshold. Default 100. 
#' @param states Numeric vector; determines which gene states to use for tSNE. If NULL, all states will be used. Default NULL.
#' @param perplexity Numeric vector; determines perplexity for tSNE. Default 20.
#' @param max.iter Numeric vector; determines iterations for tsne. Default 1000.
#' @param distance Character vector; determines which distance metric to use for tSNE. One of 'hellinger', 'cosine', 'spearman'.
#' @param seed Seed for random number generation. Default 12345.
#' @param ... Further arguments passed to or from other methods.
#' @export
celdaTsne.celda_CG = function(counts, celda.mod, max.cells=10000, min.cluster.size=100, states=NULL,
								perplexity=20, max.iter=2500, distance="hellinger", seed=12345, ...) {

  fm = factorizeMatrix(counts=counts, celda.mod=celda.mod, type="counts")
    
  states.to.use = 1:nrow(fm$counts$cell.states)
  if (!is.null(states)) {
    if (!all(states %in% states.to.use)) {
      stop("'states' must be a vector of numbers between 1 and ", states.to.use, ".")
    }
    states.to.use = states 
  }
  
  states.to.use = 1:nrow(fm$counts$cell.states)
  if (!is.null(states)) {
	if (!all(states %in% states.to.use)) {
	  stop("'states' must be a vector of numbers between 1 and ", states.to.use, ".")
	}
	states.to.use = states 
  } 
  norm = normalizeCounts(fm$counts$cell.states[states.to.use,], scale.factor=1)

  ## Select a subset of cells to sample if greater than 'max.cells'
  total.cells.to.remove = ncol(norm) - max.cells
  z.include = rep(TRUE, ncol(norm))
  if(total.cells.to.remove > 0) {
	z.ta = tabulate(celda.mod$z, celda.mod$K)
	
	## Number of cells that can be sampled from each cluster without going below the minimum threshold
	cluster.cells.to.sample = z.ta - min.cluster.size          
	cluster.cells.to.sample[cluster.cells.to.sample < 0] = 0
	
	## Number of cells to sample after exluding smaller clusters
	## Rounding can cause number to be off by a few, so ceiling is used with a second round of subtraction
	cluster.n.to.sample = ceiling((cluster.cells.to.sample / sum(cluster.cells.to.sample)) * total.cells.to.remove)
	diff = sum(cluster.n.to.sample) - total.cells.to.remove 
	cluster.n.to.sample[which.max(cluster.n.to.sample)] = cluster.n.to.sample[which.max(cluster.n.to.sample)] - diff

	## Perform sampling for each cluster
	for(i in which(cluster.n.to.sample > 0)) {
	  z.include[sample(which(celda.mod$z == i), cluster.n.to.sample[i])] = FALSE
	}
  }   
  cell.ix = which(z.include)

  res = calculateTsne(norm[,cell.ix], do.pca=FALSE, perplexity=perplexity, max.iter=max.iter, distance=distance, seed=seed)
  final = matrix(NA, nrow=ncol(norm), ncol=2)
  final[cell.ix,] = res
  rownames(final) = colnames(norm)
  colnames(final) = c("tsne_1", "tsne_2")
  return(final)
}
