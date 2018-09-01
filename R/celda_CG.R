#' celda Cell and Gene Clustering Model
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. 
#' @param sample.label Vector or factor. Denotes the sample label for each cell (column) in the count matrix.
#' @param K Integer. Number of cell populations. 
#' @param L Integer. Number of feature modules.  
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount to each cell population in each sample. Default 1. 
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to each feature module in each cell population. Default 1. 
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to each feature in each module. Default 1. 
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to the number of features in each module. Default 1. 
#' @param algorithm String. Algorithm to use for clustering cell subpopulations. One of 'EM' or 'Gibbs'. Default 'EM'.
#' @param stop.iter Integer. Number of iterations without improvement in the log likelihood to stop inference. Default 10.
#' @param max.iter Integer. Maximum number of iterations of Gibbs sampling to perform. Default 200.
#' @param split.on.iter Integer. On every `split.on.iter` iteration, a heuristic will be applied to determine if a cell population or feature module should be reassigned and another cell population or feature module should be split into two clusters. To disable splitting, set to -1. Default 10.
#' @param split.on.last Integer. After the the chain has converged, according to `stop.iter`, a heuristic will be applied to determine if a cell population or feature module should be reassigned and another cell population or feature module should be split into two clusters. If a split occurs, then 'stop.iter' will be reset. Default TRUE.
#' @param seed Integer. Passed to set.seed(). Default 12345.   
#' @param nchains Integer. Number of random cluster initializations. Default 1.  
#' @param count.checksum Character. An MD5 checksum for the `counts` matrix. Default NULL.
#' @param z.init Integer vector. Sets initial starting values of z. If NULL, starting values for each cell will be randomly sampled from 1:K. Default NULL.
#' @param y.init Integer vector. Sets initial starting values of y. If NULL, starting values for each feature will be randomly sampled from 1:L. Default NULL.
#' @param logfile Character. Messages will be redirected to a file named `logfile`. If NULL, messages will be printed to stdout.  Default NULL.
#' @param verbose Logical. Whether to print log messages. Default TRUE. 
#' @return An object of class celda_CG with clustering results and various sampling statistics.
#' @examples
#' celda.sim = simulateCells(model="celda_CG")
#' celda.mod = celda_CG(celda.sim$counts, K=celda.sim$K, L=celda.sim$L,
#'                      sample.label=celda.sim$sample.label,
#'                      max.iter=2, nchains=1)
#' @export
celda_CG = function(counts, sample.label=NULL, K, L,
                    alpha=1, beta=1, delta=1, gamma=1, 
                    algorithm = c("EM", "Gibbs"), 
                    stop.iter = 10, max.iter=200, split.on.iter=10, split.on.last=TRUE,
                    seed=12345, nchains=3, count.checksum=NULL,
                    z.init = NULL, y.init = NULL, logfile=NULL, verbose=TRUE) {

  if(is.null(count.checksum)) {
    count.checksum = digest::digest(counts, algo="md5")
  }
  counts = processCounts(counts)
    
  sample.label = processSampleLabels(sample.label, ncol(counts))
  s = as.integer(sample.label)
  
  algorithm <- match.arg(algorithm)
  algorithm.fun <- ifelse(algorithm == "Gibbs", "cC.calcGibbsProbZ", "cC.calcEMProbZ")
  
  all.seeds = seed:(seed + nchains - 1)
  
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=FALSE, verbose=verbose)  
  logMessages("Starting Celda_CG: Clustering cells and genes.", logfile=logfile, append=TRUE, verbose=verbose)
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose)  
  start.time = Sys.time()
  
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

    logMessages(date(), ".. Starting chain", i, "with seed", current.seed, logfile=logfile, append=FALSE, verbose=verbose)

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
		logMessages(date(), " .... Determining if any cell clusters should be split.", logfile=logfile, append=TRUE, sep="", verbose=verbose)
		res = cCG.splitZ(counts, m.CP.by.S, n.TS.by.C, n.TS.by.CP, n.by.G, n.by.TS, nG.by.TS, n.CP, s, z, K, L, nS, nG, alpha, beta, delta, gamma, z.prob=t(next.z$probs), max.clusters.to.try=K, min.cell=3)
		logMessages(res$message, logfile=logfile, append=TRUE, verbose=verbose)

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
		logMessages(date(), " .... Determining if any gene clusters should be split.", logfile=logfile, append=TRUE, sep="", verbose=verbose)
		res = cCG.splitY(counts, y, m.CP.by.S, n.G.by.CP, n.TS.by.C, n.TS.by.CP, n.by.G, n.by.TS, nG.by.TS, n.CP, s, z, K, L, nS, nG, alpha, beta, delta, gamma, y.prob=t(next.y$probs), max.clusters.to.try=max(L/2, 10), min.cell=3)
		logMessages(res$message, logfile=logfile, append=TRUE, verbose=verbose)

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
  
	  logMessages(date(), " .... Completed iteration: ", iter, " | logLik: ", temp.ll, logfile=logfile, append=TRUE, sep="", verbose=verbose)
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
    
    logMessages(date(), ".. Finished chain", i, "with seed", current.seed, logfile=logfile, append=TRUE, verbose=verbose)
  } 
  
  ## Peform reordering on final Z and Y assigments:
  best.result = reorder.celda_CG(counts = counts, res = best.result)
  
  end.time = Sys.time()
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose)  
  logMessages("Completed Celda_CG. Total time:", format(difftime(end.time, start.time)), logfile=logfile, append=TRUE, verbose=verbose)
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose)  

  return(best.result)
}




#' Simulate cells from the cell/feature bi-clustering generative model
#' 
#' This function generates a list containing a simulated counts matrix, as well as various parameters
#' used in the simulation which can be useful for running celda. 
#' 
#' @param model Character. Options available in `celda::available.models`. 
#' @param S Integer. Number of samples to simulate. 
#' @param C.Range Integer vector. A vector of length 2 that specifies the lower and upper bounds of the number of cells to be generated in each sample. Default c(50, 100). 
#' @param N.Range Integer vector. A vector of length 2 that specifies the lower and upper bounds of the number of counts generated for each cell. Default c(500, 5000). 
#' @param G Numeric. The total number of features to be simulated. 
#' @param K Integer. Number of cell populations. 
#' @param L Integer. Number of feature modules.  
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount to each cell population in each sample. Default 1. 
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to each feature module in each cell population. Default 1. 
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to the number of features in each module. Default 5. 
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to each feature in each module. Default 1. 
#' @param seed Integer. Passed to set.seed(). Default 12345.  
#' @param ... Additional parameters.
#' @return List. Contains the simulated counts matrix, derived cell cluster assignments, the provided parameters, and estimated Dirichlet distribution parameters for the model.
#' @examples
#' celda.cg.sim = simulateCells(model="celda_CG", K=10, L=50)
#' sim.counts = celda.cg.sim$res.listcounts
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
  cell.sample.label = factor(cell.sample.label, levels=paste0("Sample_", 1:S))
  
  ## Peform reordering on final Z and Y assigments:
  cell.counts = processCounts(cell.counts)
  names = list(row=rownames(cell.counts), column=colnames(cell.counts), 
               sample=unique(cell.sample.label))
  result = list(z=z, y=y, completeLogLik=NULL, 
                finalLogLik=NULL, K=K, L=L, alpha=alpha, 
                beta=beta, delta=delta, gamma=gamma, seed=seed, 
                sample.label=cell.sample.label, names=names,
                count.checksum=digest::digest(cell.counts, algo="md5"))
  class(result) = "celda_CG" 
  
  result = reorder.celda_CG(counts = cell.counts, res = result)
  
  return(list(z=result$z, y=result$y, sample.label=cell.sample.label, counts=cell.counts, K=K, L=L, C.Range=C.Range, N.Range=N.Range, S=S, alpha=alpha, beta=beta, gamma=gamma, delta=delta, theta=theta, phi=phi, psi=psi, eta=eta, seed=seed))
}


##' Generate factorized matrices showing each feature's influence on the celda_CG model clustering 
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda model. Options are "celda_C" or "celda_CG". Celda object of class "celda_CG".  
#' @param type Character vector. A vector containing one or more of "counts", "proportion", or "posterior". "counts" returns the raw number of counts for each factorized matrix. "proportions" returns the normalized probabilities for each factorized matrix, which are calculated by dividing the raw counts in each factorized matrix by the total counts in each column. "posterior" returns the posterior estimates. Default `c("counts", "proportion", "posterior")`.  
#' @return A list of factorized matrices, of the types requested by the user. NOTE: "population" state matrices are always returned in cell population (rows) x transcriptional states (cols).
#' @export 
factorizeMatrix.celda_CG = function(counts, celda.mod, 
                                    type=c("counts", "proportion", "posterior")) {                             
  counts = processCounts(counts)
  compareCountMatrix(counts, celda.mod)
  
  K = celda.mod$K
  L = celda.mod$L
  z = celda.mod$z
  y = celda.mod$y
  alpha = celda.mod$alpha
  beta = celda.mod$beta
  delta = celda.mod$delta
  gamma = celda.mod$gamma
  sample.label = celda.mod$sample.label
  s = as.integer(sample.label)
  
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
    temp.n.TS.by.CP[,unique.z] = normalizeCounts(temp.n.TS.by.CP[,unique.z], normalize="proportion")

    unique.y = sort(unique(y))
    temp.n.G.by.TS = n.G.by.TS
    temp.n.G.by.TS[,unique.y] = normalizeCounts(temp.n.G.by.TS[,unique.y], normalize="proportion")
    temp.nG.by.TS = nG.by.TS/sum(nG.by.TS)
    
    prop.list = list(sample.states =  normalizeCounts(m.CP.by.S, normalize="proportion"),
    				   population.states = temp.n.TS.by.CP, 
    				   cell.states = normalizeCounts(n.TS.by.C, normalize="proportion"),
    				   gene.states = temp.n.G.by.TS, 
    				   gene.distribution = temp.nG.by.TS)
    res = c(res, list(proportions=prop.list))
  }
  if(any("posterior" %in% type)) {

    gs = n.G.by.TS
    gs[cbind(1:nG,y)] = gs[cbind(1:nG,y)] + delta
    gs = normalizeCounts(gs, normalize="proportion")
	temp.nG.by.TS = (nG.by.TS + gamma)/sum(nG.by.TS + gamma)
	
    post.list = list(sample.states = normalizeCounts(m.CP.by.S + alpha, normalize="proportion"),
          				   population.states = normalizeCounts(n.TS.by.CP + beta, normalize="proportion"), 
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
#' @param counts Integer matrix. Rows represent features and columns represent cells. 
#' @param sample.label Vector or factor. Denotes the sample label for each cell (column) in the count matrix.
#' @param z Numeric vector. Denotes cell population labels. 
#' @param y Numeric vector. Denotes feature module labels. 
#' @param K Integer. Number of cell populations. 
#' @param L Integer. Number of feature modules.  
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount to each cell population in each sample. Default 1.  
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to each feature module in each cell population. Default 1. 
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to each feature in each module. Default 1. 
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to the number of features in each module. Default 1. 
#' @param ... Additional parameters.
#' @examples
#' celda.sim = simulateCells(model="celda_CG")
#' celda.mod = "celda_CG"
#' loglik = calculateLoglikFromVariables(celda.sim$counts, celda.mod, 
#'                                       sample.label=celda.sim$sample.label,
#'                                       z=celda.sim$z, y=celda.sim$y,
#'                                       K=celda.sim$K, L=celda.sim$L,
#'                                       alpha=celda.sim$alpha, beta=celda.sim$beta,
#'                                       gamma=celda.sim$gamma, delta=celda.sim$delta)
#' @export
calculateLoglikFromVariables.celda_CG = function(counts, sample.label, z, y, K, L, alpha, beta, delta, gamma) {  
  sample.label = processSampleLabels(sample.label, ncol(counts))
  s = as.integer(sample.label)
  p = cCG.decomposeCounts(counts, s, z, y, K, L)
  final = cCG.calcLL(K=K, L=L, m.CP.by.S=p$m.CP.by.S, n.TS.by.CP=p$n.TS.by.CP, n.by.G=p$n.by.G, n.by.TS=p$n.by.TS, nG.by.TS=p$nG.by.TS, nS=p$nS, nG=p$nG, alpha=alpha, beta=beta, delta=delta, gamma=gamma)
  return(final)
}


#' Takes raw counts matrix and converts it to a series of matrices needed for log likelihood calculation
#' @param counts Integer matrix. Rows represent features and columns represent cells. 
#' @param s Integer vector. Contains the sample label for each cell (column) in the count matrix. 
#' @param z Numeric vector. Denotes cell population labels. 
#' @param y Numeric vector. Denotes feature module labels. 
#' @param K Integer. Number of cell populations. 
#' @param L Integer. Number of feature modules.  
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
#' @param celda.mod Celda object of class "celda_CG".
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param log Logical. If FALSE, then the normalized conditional probabilities will be returned. If TRUE, then the unnormalized log probabilities will be returned. Default FALSE.  
#' @param ... Additional parameters.
#' @return A list containging a matrix for the conditional cell and feature cluster probabilities. 
#' @export
clusterProbability.celda_CG = function(counts, celda.mod, log=FALSE, ...) {
  
  s = as.integer(celda.mod$sample.label)
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
  next.y = cG.calcGibbsProbY(counts=p$n.G.by.CP, n.TS.by.C=p$n.TS.by.CP, n.by.TS=p$n.by.TS, nG.by.TS=p$nG.by.TS, n.by.G=p$n.by.G, y=y, L=L, nG=p$nG, beta=beta, delta=delta, gamma=gamma, do.sample=FALSE)
  y.prob = t(next.y$probs)

  if(!isTRUE(log)) {
    z.prob = normalizeLogProbs(z.prob)
    y.prob = normalizeLogProbs(y.prob)
  }
       
  return(list(z.probability=z.prob, y.probability=y.prob))
}


#' Calculate the perplexity from a single celda model
#' 
#' Perplexity can be seen as a measure of how well a provided set of 
#' cluster assignments fit the data being clustered.
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class "celda_C", "celda_G" or "celda_CG".
#' @param new.counts A new counts matrix used to calculate perplexity. If NULL, perplexity will be calculated for the 'counts' matrix. Default NULL.
#' @return Numeric. The perplexity for the provided count data and model.
#' @examples
#' celda.sim = simulateCells(model="celda_CG")
#' celda.mod = celda_CG(celda.sim$counts, K=celda.sim$K, L=celda.sim$L, max.iter=2, nchains=1)
#' perplexity = calculatePerplexity(celda.sim$counts, celda.mod)
#' @export
calculatePerplexity.celda_CG = function(counts, celda.mod, new.counts=NULL) {
  
  counts = processCounts(counts)
  compareCountMatrix(counts, celda.mod)

  if(is.null(new.counts)) {
    new.counts = counts
  } else {
    new.counts = processCounts(new.counts)
  }
  if(nrow(new.counts) != nrow(counts)) {
    stop("new.counts should have the same number of rows as counts.")
  }
  
  factorized = factorizeMatrix(counts = counts, celda.mod = celda.mod, 
                               type=c("posterior", "counts"))
  theta = log(factorized$posterior$sample.states)
  phi   = factorized$posterior$population.states
  psi   = factorized$posterior$gene.states
  s = as.integer(celda.mod$sample.label)
  eta <- factorized$posterior$gene.distribution
  nG.by.TS = factorized$counts$gene.distribution
  
  eta.prob = log(eta) * nG.by.TS
  gene.by.pop.prob = log(psi %*% phi)
  inner.log.prob = (t(gene.by.pop.prob) %*% new.counts) + theta[, s]  
  
  log.px = sum(apply(inner.log.prob, 2, matrixStats::logSumExp)) + sum(eta.prob)
  perplexity = exp(-(log.px/sum(new.counts)))
  return(perplexity)
}


reorder.celda_CG = function(counts, res){
  # Reorder K
  if(res$K > 2 & isTRUE(length(unique(res$z)) > 1)) {
    res$z = as.integer(as.factor(res$z))
    fm <- factorizeMatrix(counts = counts, celda.mod = res, type="posterior")
    unique.z = sort(unique(res$z))
    d <- cosineDist(fm$posterior$population.states[,unique.z])
    h <- hclust(d, method = "complete")
    
    res <- recodeClusterZ(res, from = h$order, to = 1:length(h$order))
  }  
  
  # Reorder L
  if(res$L > 2 & isTRUE(length(unique(res$y)) > 1)) {
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


#' celdaHeatmap for celda Cell and Gene clustering model.
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class "celda_CG". 
#' @param nfeatures Integer. Maximum number of features to select for each module. Default 25.
#' @param ... Additional parameters.
#' @examples 
#' celda.sim = simulateCells("celda_CG")
#' celda.mod = celda_CG(celda.sim$counts, K=celda.sim$K, L=celda.sim$L,
#'                      nchains=1, max.iter=1)
#' celdaHeatmap(celda.sim$counts, celda.mod)
#' @export
celdaHeatmap.celda_CG = function(counts, celda.mod, nfeatures=25, ...) {
  fm = factorizeMatrix(counts, celda.mod, type="proportion")
  top = topRank(fm$proportions$gene.states, n=nfeatures)
  ix = unlist(top$index)
  norm = normalizeCounts(counts, normalize="proportion", transformation.fun=sqrt)
  plotHeatmap(norm[ix,], z=celda.mod$z, y=celda.mod$y[ix], ...)
}


#' Embeds cells in two dimensions using tSNE based on celda_CG results.
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda model. Options available in `celda::available.models`. 
#' @param max.cells Integer. Maximum number of cells to plot. Cells will be randomly subsampled if ncol(conts) > max.cells. Larger numbers of cells requires more memory. Default 10000.
#' @param min.cluster.size Integer. Do not subsample cell clusters below this threshold. Default 100. 
#' @param modules Integer vector. Determines which features modules to use for tSNE. If NULL, all modules will be used. Default NULL.
#' @param perplexity Numeric. Perplexity parameter for tSNE. Default 20.
#' @param max.iter Integer. Maximum number of iterations in tSNE generation. Default 2500.
#' @param seed Integer. Passed to set.seed(). Default 12345.  
#' @param ... Additional parameters.
#' @export
celdaTsne.celda_CG = function(counts, celda.mod, max.cells=10000, min.cluster.size=100, modules=NULL,
								perplexity=20, max.iter=2500, seed=12345, ...) {

  ## Checking if max.cells and min.cluster.size will work
  if(max.cells / min.cluster.size < celda.mod$K) {
    stop(paste0("Cannot distribute ", max.cells, " cells among ", celda.mod$K, " clusters while maintaining a minumum of ", min.cluster.size, " cells per cluster. Try increasing 'max.cells' or decreasing 'min.cluster.size'."))
  }

  fm = factorizeMatrix(counts=counts, celda.mod=celda.mod, type="counts")    
  modules.to.use = 1:nrow(fm$counts$cell.states)
  if (!is.null(modules)) {
    if (!all(modules %in% modules.to.use)) {
      stop("'modules' must be a vector of numbers between 1 and ", states.to.use, ".")
    }
    modules.to.use = modules 
  }
  

  ## Select a subset of cells to sample if greater than 'max.cells'
  total.cells.to.remove = ncol(counts) - max.cells
  z.include = rep(TRUE, ncol(counts))
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

  norm = t(normalizeCounts(fm$counts$cell.states[modules.to.use,cell.ix], normalize="proportion", transformation.fun=sqrt))
  
  res = calculateTsne(norm, do.pca=FALSE, perplexity=perplexity, max.iter=max.iter, seed=seed)
  final = matrix(NA, nrow=ncol(counts), ncol=2)
  final[cell.ix,] = res
  rownames(final) = colnames(counts)
  colnames(final) = c("tsne_1", "tsne_2")
  return(final)
}




#' Renders probability and relative expression heatmaps to visualize the relationship between feature modules and cell populations.
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`. 
#' @param celda.mod Celda object of class "celda_CG".   
#' @param level Character. One of "cell.population" or "sample". "cell.population" will display the absolute probabilities and relative normalized expression of each module in each cell population. "sample" will display the absolute probabilities and relative normalized abundance of each cell population in each sample." Default "cell.population".
#' @param ... Additional parameters.
#' @export 
celdaProbabilityMap.celda_CG <- function(counts, celda.mod, level=c("cell.population", "sample"), ...){
  counts = processCounts(counts)
  compareCountMatrix(counts, celda.mod)
  
  level = match.arg(level)
  factorized <- factorizeMatrix(celda.mod = celda.mod, counts = counts)
  z.include = which(tabulate(celda.mod$z, celda.mod$K) > 0)
  y.include = which(tabulate(celda.mod$y, celda.mod$L) > 0)    

  if(level == "cell.population") {
    pop <- factorized$proportions$population.states[y.include,z.include,drop=FALSE]
    pop.norm = normalizeCounts(pop, normalize="proportion", transformation.fun=sqrt, scale.fun=base::scale)
  
    percentile.9 <- round(quantile(pop,.9), digits = 2) * 100
    col1 <- colorRampPalette(c("#FFFFFF", brewer.pal(n = 9, name = "Blues")))(percentile.9)
    col2 <- colorRampPalette(c("#08306B", c("#006D2C","Yellowgreen","Yellow","Orange","Red")))(100-percentile.9)
    col <- c(col1,col2)
    breaks <-  seq(0, 1, length.out = length(col))     
    
    g1 = plotHeatmap(pop, color.scheme="sequential", scale.row=NULL, cluster.cell=FALSE, cluster.feature=FALSE, show.names.cell=TRUE, show.names.feature=TRUE, breaks = breaks, col=col, main = "Absolute Probability", silent=TRUE)
    g2 = plotHeatmap(pop.norm, color.scheme="divergent", cluster.cell=FALSE, cluster.feature=FALSE, show.names.cell=TRUE, show.names.feature=TRUE, main = "Relative Expression", silent=TRUE) 
    gridExtra::grid.arrange(g1$gtable, g2$gtable, ncol=2)
  } else {
    samp <- factorized$proportions$sample.states
    col <- colorRampPalette(c("white","blue","#08306B","#006D2C","yellowgreen","yellow","orange","red"))(100)
    breaks <-  seq(0, 1, length.out = length(col))     
    g1 = plotHeatmap(samp, color.scheme="sequential", scale.row=NULL, cluster.cell=FALSE, cluster.feature=FALSE, show.names.cell=TRUE, show.names.feature=TRUE, breaks = breaks, col=col, main = "Absolute Probability", silent=TRUE)

    if(ncol(samp) > 1) {
      samp.norm = normalizeCounts(factorized$counts$sample.states, normalize="proportion", transformation.fun=sqrt, scale.fun=base::scale)
      g2 = plotHeatmap(samp.norm, color.scheme="divergent", cluster.cell=FALSE, cluster.feature=FALSE, show.names.cell=TRUE, show.names.feature=TRUE, main = "Relative Abundance", silent=TRUE)   
      gridExtra::grid.arrange(g1$gtable, g2$gtable, ncol=2)
    } else {
      gridExtra::grid.arrange(g1$gtable)
    } 
  }
}


#' Obtain the feature module of a feature of interest
#' 
#' This function will output the gene module of a specific gene(s) from a celda model
#'  
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Model of class "celda_G" or "celda_CG".
#' @param feature Character vector. Identify feature modules for the specified feature names. 
#' @return List. Each entry corresponds to the feature module determined for the provided features
#' @examples
#' celda.mod = celda_CG(celda::pbmc_select, K=10, 
#'                      L=50, max.iter=2, nchains=1)
#' corresponding.module = featureModuleLookup(celda::pbmc_select, celda.mod, c("ENSG00000000938_FGR", "ENSG00000004059_ARF5"))
#' @export
featureModuleLookup.celda_CG = function(counts, celda.mod, feature){
  list <- list()
  for(x in 1:length(feature)){
    if(feature[x] %in% rownames(counts)){
      list[x] <- celda.mod$y[which(rownames(counts) == feature[x])]
    }else{
      list[x] <- paste0("No feature was identified matching '", feature[x], "'.")
    }
  } 
  names(list) <- feature
  return(list)
}
