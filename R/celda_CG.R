#' @title Cell and feature clustering with Celda
#' 
#' @description Clusters the rows and columns of a count matrix containing single-cell data into L modules and K subpopulations, respectively.  
#'
#' @param counts Integer matrix. Rows represent features and columns represent cells. 
#' @param sample.label Vector or factor. Denotes the sample label for each cell (column) in the count matrix.
#' @param K Integer. Number of cell populations. 
#' @param L Integer. Number of feature modules.  
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount to each cell population in each sample. Default 1. 
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to each feature module in each cell population. Default 1. 
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to each feature in each module. Default 1. 
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to the number of features in each module. Default 1. 
#' @param algorithm String. Algorithm to use for clustering cell subpopulations. One of 'EM' or 'Gibbs'. The EM algorithm for cell clustering is faster, especially for larger numbers of cells. However, more chains may be required to ensure a good solution is found. Default 'EM'.
#' @param stop.iter Integer. Number of iterations without improvement in the log likelihood to stop inference. Default 10.
#' @param max.iter Integer. Maximum number of iterations of Gibbs sampling to perform. Default 200.
#' @param split.on.iter Integer. On every `split.on.iter` iteration, a heuristic will be applied to determine if a cell population or feature module should be reassigned and another cell population or feature module should be split into two clusters. To disable splitting, set to -1. Default 10.
#' @param split.on.last Integer. After `stop.iter` iterations have been performed without improvement, a heuristic will be applied to determine if a cell population or feature module should be reassigned and another cell population or feature module should be split into two clusters. If a split occurs, then 'stop.iter' will be reset. Default TRUE.
#' @param seed Integer. Passed to `set.seed()`. Default 12345.   
#' @param nchains Integer. Number of random cluster initializations. Default 3.  
#' @param initialize Chararacter. One of 'random' or 'split'. With 'random', cells and features are randomly assigned to a clusters. With 'split' cell and feature clusters will be recurssively split into two clusters using `celda_C` and `celda_G`, respectively, until the specified K and L is reached. Default 'random'.
#' @param count.checksum Character. An MD5 checksum for the `counts` matrix. Default NULL.
#' @param z.init Integer vector. Sets initial starting values of z. If NULL, starting values for each cell will be randomly sampled from 1:K. 'z.init' can only be used when `initialize' = 'random'`. Default NULL.
#' @param y.init Integer vector. Sets initial starting values of y. If NULL, starting values for each feature will be randomly sampled from 1:L. 'y.init' can only be used when `initialize = 'random'`. Default NULL.
#' @param logfile Character. Messages will be redirected to a file named `logfile`. If NULL, messages will be printed to stdout.  Default NULL.
#' @param verbose Logical. Whether to print log messages. Default TRUE. 
#' @return An object of class `celda_CG` with the cell populations clusters stored in in `z` and feature module clusters stored in `y`.
#' @seealso `celda_G()` for feature clustering and `celda_C()` for clustering cells. `celdaGridSearch()` can be used to run multiple values of K/L and multiple chains in parallel. 
#' @examples
#' celda.mod = celda_CG(celda.CG.sim$counts, K=celda.CG.sim$K, L=celda.CG.sim$L,
#'                      sample.label=celda.CG.sim$sample.label, nchains=1)
#' @export
celda_CG = function(counts, sample.label=NULL, K, L,
                    alpha=1, beta=1, delta=1, gamma=1, 
                    algorithm = c("EM", "Gibbs"), 
                    stop.iter = 10, max.iter=200, split.on.iter=10, split.on.last=TRUE,
                    seed=12345, nchains=3, initialize=c("random", "split"), count.checksum=NULL,
                    z.init = NULL, y.init = NULL, logfile=NULL, verbose=TRUE) {
  validateCounts(counts)
  return(.celda_CG(counts, sample.label, K, L, alpha, beta, delta, gamma,
                   algorithm, stop.iter, max.iter, split.on.iter, split.on.last,
                   seed, nchains, initialize, count.checksum, z.init, y.init,
                   logfile, verbose))
}

.celda_CG = function(counts, sample.label=NULL, K, L,
                    alpha=1, beta=1, delta=1, gamma=1, 
                    algorithm = c("EM", "Gibbs"), 
                    stop.iter = 10, max.iter=200, split.on.iter=10, split.on.last=TRUE,
                    seed=12345, nchains=3, initialize=c("random", "split"), count.checksum=NULL,
                    z.init = NULL, y.init = NULL, logfile=NULL, verbose=TRUE) {

  logMessages("--------------------------------------------------------------------", logfile=logfile, append=FALSE, verbose=verbose)  
  logMessages("Starting Celda_CG: Clustering cells and genes.", logfile=logfile, append=TRUE, verbose=verbose)
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose)  
  start.time = Sys.time()

  counts = processCounts(counts)
  if(is.null(count.checksum)) {
    count.checksum = digest::digest(counts, algo="md5")
  }
    
  sample.label = processSampleLabels(sample.label, ncol(counts))
  s = as.integer(sample.label)
  
  algorithm <- match.arg(algorithm)
  algorithm.fun <- ifelse(algorithm == "Gibbs", "cC.calcGibbsProbZ", "cC.calcEMProbZ")
  initialize = match.arg(initialize)
  
  all.seeds = seed:(seed + nchains - 1)
  
  best.result = NULL  
  for(i in seq_along(all.seeds)) { 
  
	## Initialize cluster labels
    current.seed = all.seeds[i]	
    logMessages(date(), ".. Initializing chain", i, "with", paste0("'",initialize, "' (seed=", current.seed, ")"), logfile=logfile, append=TRUE, verbose=verbose)

    if(initialize == "random") {
  	  z = initialize.cluster(K, ncol(counts), initial = z.init, fixed = NULL, seed=current.seed)
	  y = initialize.cluster(L, nrow(counts), initial = y.init, fixed = NULL, seed=current.seed)
	} else {
	  z = recursive.splitZ(counts, s, K=K, alpha=alpha, beta=beta)
	  y = recursive.splitY(counts, L, beta=beta, delta=delta, gamma=gamma, z=z, 
	                       K=K, K.subclusters=10, min.feature=3, max.cells=100,
	                       seed=seed)
	}  
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
	  temp.ll = cCG.calcLL(K=K, L=L, m.CP.by.S=m.CP.by.S, n.TS.by.CP=n.TS.by.CP, n.by.G=n.by.G, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, nS=nS, nG=nG, alpha=alpha, beta=beta, delta=delta, gamma=gamma)
	  if(((iter == max.iter | (num.iter.without.improvement == stop.iter & all(temp.ll < ll)) & isTRUE(split.on.last)) | (split.on.iter > 0 & iter %% split.on.iter == 0))) {
		if(K > 2 & isTRUE(do.cell.split)) {
		  logMessages(date(), " .... Determining if any cell clusters should be split.", 
		              logfile=logfile, append=TRUE, sep="", verbose=verbose)
		  res = cCG.splitZ(counts, m.CP.by.S, n.TS.by.C, n.TS.by.CP, n.by.G, n.by.TS, 
		                   nG.by.TS, n.CP, s, z, K, L, nS, nG, alpha, beta, delta, 
		                   gamma, z.prob=t(next.z$probs),  max.clusters.to.try=K, 
		                   min.cell=3)
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
		if(L > 2 & isTRUE(do.gene.split)) {
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
  best.result = methods::new("celda_CG", 
                             clusters=list(z=z.best, y=y.best),
                             params=list(K=K, L=L, alpha=alpha, beta=beta, 
                                         delta=delta, gamma=gamma, 
                                         seed=current.seed,
                                         count.checksum=count.checksum),
                             completeLogLik=ll, 
                             finalLogLik=ll.best,
                  				   sample.label=sample.label, 
                  				   names=names)
  best.result = reorder.celda_CG(counts = counts, res = best.result)
  
  end.time = Sys.time()
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose)  
  logMessages("Completed Celda_CG. Total time:", format(difftime(end.time, start.time)), logfile=logfile, append=TRUE, verbose=verbose)
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose)  

  return(best.result)
}




#' @title Simulate cells from the celda_CG model
#' 
#' @description Generates a simulated counts matrix, cell subpopulation clusters, sample labels, and feature module clusters
#' according to the generative process of the celda_CG model. 
#' 
#' @param model Character. Options available in `celda::available.models`. 
#' @param S Integer. Number of samples to simulate. Default 5.  
#' @param C.Range Integer vector. A vector of length 2 that specifies the lower and upper bounds of the number of cells to be generated in each sample. Default c(50, 100). 
#' @param N.Range Integer vector. A vector of length 2 that specifies the lower and upper bounds of the number of counts generated for each cell. Default c(500, 1000). 
#' @param G Integer. The total number of features to be simulated. Default 100.
#' @param K Integer. Number of cell populations. Default 5. 
#' @param L Integer. Number of feature modules. Default 10.
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount to each cell population in each sample. Default 1. 
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to each feature module in each cell population. Default 1. 
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to the number of features in each module. Default 5. 
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to each feature in each module. Default 1. 
#' @param seed Integer. Passed to `set.seed()`. Default 12345.  
#' @param ... Additional parameters.
#' @return List. Contains the simulated matrix `counts`, cell population clusters `z`, feature module clusters `y`, sample assignments `sample.label`, and input parameters.
#' @seealso `celda_C()` for simulating cell subpopulations and `celda_G()` for simulating feature modules. 
#' @examples
#' celda.sim = simulateCells(model="celda_CG")
#' @export
simulateCells.celda_CG = function(model, S=5, C.Range=c(50,100), N.Range=c(500,1000), 
                                  G=100, K=5, L=10, alpha=1, beta=1, gamma=5, 
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

  ## Assign genes to gene modules 
  eta = rdirichlet(1, rep(gamma, L))
  y = sample(1:L, size=G, prob=eta, replace=TRUE)
  if(length(table(y)) < L) {
    warning("Some gene modules did not receive any genes after sampling. Try increasing G and/or making gamma larger.")
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
  result = methods::new("celda_CG", 
                        clusters=list(z=z, y=y),
                        params=list(K=K, L=L, alpha=alpha, beta=beta, delta=delta, 
                                    gamma=gamma, seed=seed,
                                    count.checksum=digest::digest(cell.counts, algo="md5")),
            				    sample.label=cell.sample.label, names=names)
  
  result = reorder.celda_CG(counts = cell.counts, res = result)
  
  return(list(z=result@clusters$z, y=result@clusters$y, sample.label=cell.sample.label, 
              counts=cell.counts, K=K, L=L, C.Range=C.Range, 
              N.Range=N.Range, S=S, alpha=alpha, beta=beta, gamma=gamma, 
              delta=delta, seed=seed))
}


#' @title Matrix factorization for results from celda_CG
#' @description Generates factorized matrices showing the contribution of each feature in each module, each module in each cell and/or cell population, and each cell population in each sample. 
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda model. Options are "celda_C" or "celda_CG". Celda object of class "celda_CG".  
#' @param type Character vector. A vector containing one or more of "counts", "proportion", or "posterior". "counts" returns the raw number of counts for each factorized matrix. "proportions" returns the normalized probabilities for each factorized matrix, which are calculated by dividing the raw counts in each factorized matrix by the total counts in each column. "posterior" returns the posterior estimates. Default `c("counts", "proportion", "posterior")`.  
#' @return A list with elements for `counts`, `proportions`, or `posterior` probabilities. Each element will be a list containing factorized matrices for `module`, `cell.population`, and `sample`. Additionally, the contribution of each module in each individual cell will be included in the `cell` element of `counts` and `proportions` elements. 
#' @seealso `celda_CG()` for clustering features and cells
#' @examples 
#' factorized.matrices = factorizeMatrix(celda.CG.sim$counts, celda.CG.mod, "posterior")
#' @export 
setMethod("factorizeMatrix",
          signature(celda.mod = "celda_CG"),
          function(counts, celda.mod,  
                   type=c("counts", "proportion", "posterior")) {
            counts = processCounts(counts)
            compareCountMatrix(counts, celda.mod)
            
            K = celda.mod@params$K
            L = celda.mod@params$L
            z = celda.mod@clusters$z
            y = celda.mod@clusters$y
            alpha = celda.mod@params$alpha
            beta = celda.mod@params$beta
            delta = celda.mod@params$delta
            gamma = celda.mod@params$gamma
            sample.label = celda.mod@sample.label
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
            colnames(n.TS.by.C) = celda.mod@names$column
            rownames(n.TS.by.C) = L.names
            colnames(n.G.by.TS) = L.names
            rownames(n.G.by.TS) = celda.mod@names$row
            rownames(m.CP.by.S) = K.names
            colnames(m.CP.by.S) = celda.mod@names$sample
            colnames(n.TS.by.CP) = K.names
            rownames(n.TS.by.CP) = L.names
          
            counts.list = c()
            prop.list = c()
            post.list = c()
            res = list()
              
            if(any("counts" %in% type)) {
              counts.list = list(sample = m.CP.by.S,
                      				   cell.population = n.TS.by.CP, 
                      				   cell = n.TS.by.C,
                      				   module = n.G.by.TS,
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
              
              prop.list = list(sample =  normalizeCounts(m.CP.by.S, normalize="proportion"),
              				   cell.population = temp.n.TS.by.CP, 
              				   cell = normalizeCounts(n.TS.by.C, normalize="proportion"),
              				   module = temp.n.G.by.TS, 
              				   gene.distribution = temp.nG.by.TS)
              res = c(res, list(proportions=prop.list))
            }
            if(any("posterior" %in% type)) {
          
              gs = n.G.by.TS
              gs[cbind(1:nG,y)] = gs[cbind(1:nG,y)] + delta
              gs = normalizeCounts(gs, normalize="proportion")
          	temp.nG.by.TS = (nG.by.TS + gamma)/sum(nG.by.TS + gamma)
          	
              post.list = list(sample = normalizeCounts(m.CP.by.S + alpha, normalize="proportion"),
                    				   cell.population = normalizeCounts(n.TS.by.CP + beta, normalize="proportion"), 
                    				   module = gs,
                    				   gene.distribution = temp.nG.by.TS)
              res = c(res, posterior = list(post.list))						    
            }
            
            return(res)
          })

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


#' @title Calculate Celda_CG log likelihood
#' @description Calculates the log likelihood for user-provided cell population and feature module clusters using the `celda_CG()` model.
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
#' @return The log likelihood for the given cluster assignments
#' @seealso `celda_CG()` for clustering features and cells
#' @examples
#' loglik = logLikelihood(celda.CG.sim$counts, model="celda_CG", 
#'                        sample.label=celda.CG.sim$sample.label,
#'                        z=celda.CG.sim$z, y=celda.CG.sim$y,
#'                        K=celda.CG.sim$K, L=celda.CG.sim$L,
#'                        alpha=celda.CG.sim$alpha, beta=celda.CG.sim$beta,
#'                        gamma=celda.CG.sim$gamma, delta=celda.CG.sim$delta)
#' @export
logLikelihood.celda_CG = function(counts, sample.label, z, y, K, L, alpha, beta, delta, gamma) {  
  if (sum(z > K) > 0) stop("An entry in z contains a value greater than the provided K.")
  if (sum(y > L) > 0) stop("An entry in y contains a value greater than the provided L.")
  sample.label = processSampleLabels(sample.label, ncol(counts))
  s = as.integer(sample.label)
  p = cCG.decomposeCounts(counts, s, z, y, K, L)
  final = cCG.calcLL(K=K, L=L, m.CP.by.S=p$m.CP.by.S, n.TS.by.CP=p$n.TS.by.CP, n.by.G=p$n.by.G, n.by.TS=p$n.by.TS, nG.by.TS=p$nG.by.TS, nS=p$nS, nG=p$nG, alpha=alpha, beta=beta, delta=delta, gamma=gamma)
  return(final)
}


# Takes raw counts matrix and converts it to a series of matrices needed for log likelihood calculation
# @param counts Integer matrix. Rows represent features and columns represent cells. 
# @param s Integer vector. Contains the sample label for each cell (column) in the count matrix. 
# @param z Numeric vector. Denotes cell population labels. 
# @param y Numeric vector. Denotes feature module labels. 
# @param K Integer. Number of cell populations. 
# @param L Integer. Number of feature modules.  
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



#' @title Conditional probabilities for cells and features from a Celda_CG model
#' @description Calculates the conditional probability of each cell belonging to each subpopulation given all other cell cluster assignments as well as each feature belonging to each module given all other feature cluster assignments in a `celda_CG()` result. 
#'
#' @param celda.mod Celda object of class `celda_CG`.
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param log Logical. If FALSE, then the normalized conditional probabilities will be returned. If TRUE, then the unnormalized log probabilities will be returned. Default FALSE.  
#' @param ... Additional parameters.
#' @return A list containging a matrix for the conditional cell and feature cluster probabilities. 
#' @seealso `celda_CG()` for clustering features and cells
#' @examples
#' cluster.prob = clusterProbability(celda.CG.sim$counts, celda.CG.mod)
#' @export
setMethod("clusterProbability",
          signature(celda.mod = "celda_CG"),
          function(counts, celda.mod, log=FALSE, ...) {
  
            s = as.integer(celda.mod@sample.label)
            z = celda.mod@clusters$z
            K = celda.mod@params$K  
            y = celda.mod@clusters$y
            L = celda.mod@params$L
            alpha = celda.mod@params$alpha
            delta = celda.mod@params$delta
            beta = celda.mod@params$beta
            gamma = celda.mod@params$gamma
          
            p = cCG.decomposeCounts(counts, s, z, y, K, L)
          
            ## Gibbs sampling for each cell
            next.z = cC.calcGibbsProbZ(counts=p$n.TS.by.C, m.CP.by.S=p$m.CP.by.S, 
                                       n.G.by.CP=p$n.TS.by.CP, n.CP=p$n.CP, 
                                       n.by.C=p$n.by.C, z=z, s=s, K=K, nG=L, 
                                       nM=p$nM, alpha=alpha, beta=beta, 
                                       do.sample=FALSE)
            z.prob = t(next.z$probs)
          
            ## Gibbs sampling for each gene
            next.y = cG.calcGibbsProbY(counts=p$n.G.by.CP, n.TS.by.C=p$n.TS.by.CP, 
                                       n.by.TS=p$n.by.TS, nG.by.TS=p$nG.by.TS, 
                                       n.by.G=p$n.by.G, y=y, L=L, nG=p$nG, 
                                       beta=beta, delta=delta, gamma=gamma, 
                                       do.sample=FALSE)
            y.prob = t(next.y$probs)
          
            if(!isTRUE(log)) {
              z.prob = normalizeLogProbs(z.prob)
              y.prob = normalizeLogProbs(y.prob)
            }
                 
            return(list(z.probability=z.prob, y.probability=y.prob))
         })


#' @title Calculate the perplexity on new data with a celda_CG model
#' @description Perplexity is a statistical measure of how well a probability model can predict new data. Lower perplexity indicates a better model. 
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class "celda_C", "celda_G" or "celda_CG".
#' @param new.counts A new counts matrix used to calculate perplexity. If NULL, perplexity will be calculated for the 'counts' matrix. Default NULL.
#' @return Numeric. The perplexity for the provided count data and model.
#' @seealso `celda_CG()` for clustering features and cells
#' @examples
#' perplexity = perplexity(celda.CG.sim$counts, celda.CG.mod)
#' @export
setMethod("perplexity",
          signature(celda.mod = "celda_CG"),
          function(counts, celda.mod, new.counts=NULL) {
            if (!("celda_CG" %in% class(celda.mod))) stop("The celda.mod provided was not of class celda_CG.")
            
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
            theta = log(factorized$posterior$sample)
            phi   = factorized$posterior$cell.population
            psi   = factorized$posterior$module
            s = as.integer(celda.mod@sample.label)
            eta <- factorized$posterior$gene.distribution
            nG.by.TS = factorized$counts$gene.distribution
            
            eta.prob = log(eta) * nG.by.TS
            gene.by.pop.prob = log(psi %*% phi)
            inner.log.prob = (t(gene.by.pop.prob) %*% new.counts) + theta[, s]  
            
            log.px = sum(apply(inner.log.prob, 2, matrixStats::logSumExp)) + sum(eta.prob)
            perplexity = exp(-(log.px/sum(new.counts)))
            return(perplexity)
          })


reorder.celda_CG = function(counts, res){
  # Reorder K
  if(res@params$K > 2 & isTRUE(length(unique(res@clusters$z)) > 1)) {
    res@clusters$z = as.integer(as.factor(res@clusters$z))
    fm <- factorizeMatrix(counts = counts, celda.mod = res, type="posterior")
    unique.z = sort(unique(res@clusters$z))
    d <- cosineDist(fm$posterior$cell.population[,unique.z])
    h <- stats::hclust(d, method = "complete")
    
    res <- recodeClusterZ(res, from = h$order, to = 1:length(h$order))
  }  
  
  # Reorder L
  if(res@params$L > 2 & isTRUE(length(unique(res@clusters$y)) > 1)) {
    res@clusters$y = as.integer(as.factor(res@clusters$y))
    fm <- factorizeMatrix(counts = counts, celda.mod = res, type="posterior")
    unique.y = sort(unique(res@clusters$y))
    cs <- prop.table(t(fm$posterior$cell.population[unique.y,]), 2)
    d <- cosineDist(cs)
    h <- stats::hclust(d, method = "complete")
    
    res <- recodeClusterY(res, from = h$order, to = 1:length(h$order))
  }
  return(res)
}


#' @title Heatmap for celda_CG
#' @description Renders an expression heatmap to visualize `celda_CG()` results. The top `nfeatures` for each module will be included in the heatmap. 
#'
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class `celda_CG`. 
#' @param nfeatures Integer. Maximum number of features to select for each module. Default 25.
#' @param ... Additional parameters.
#' @seealso `celda_CG()` for clustering features and cells and `celdaTsne()` for generating 2-dimensional coordinates
#' @examples 
#' celdaHeatmap(celda.CG.sim$counts, celda.CG.mod)
#' @return A list containing dendrograms and the heatmap grob
#' @export
setMethod("celdaHeatmap",
          signature(celda.mod = "celda_CG"),
          function(counts, celda.mod, nfeatures=25, ...) {
            fm = factorizeMatrix(counts, celda.mod, type="proportion")
            top = celda::topRank(fm$proportions$module, n=nfeatures)
            ix = unlist(top$index)
            norm = normalizeCounts(counts, normalize="proportion", transformation.fun=sqrt)
            plotHeatmap(norm[ix,], z=celda.mod@clusters$z, y=celda.mod@clusters$y[ix], ...)
          })


#' @title tSNE for celda_CG
#' @description Embeds cells in two dimensions using tSNE based on a `celda_CG` model. tSNE is run on module probabilities to reduce the number of features instead of using PCA. Module probabilities square-root trasformed before applying tSNE. 
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class `celda_CG`. 
#' @param max.cells Integer. Maximum number of cells to plot. Cells will be randomly subsampled if ncol(counts) > max.cells. Larger numbers of cells requires more memory. Default 25000.
#' @param min.cluster.size Integer. Do not subsample cell clusters below this threshold. Default 100. 
#' @param modules Integer vector. Determines which features modules to use for tSNE. If NULL, all modules will be used. Default NULL.
#' @param perplexity Numeric. Perplexity parameter for tSNE. Default 20.
#' @param max.iter Integer. Maximum number of iterations in tSNE generation. Default 2500.
#' @param seed Integer. Passed to `set.seed()`. Default 12345.  
#' @param ... Additional parameters.
#' @seealso `celda_CG()` for clustering features and cells  and `celdaHeatmap()` for displaying expression
#' @examples
#' tsne.res = celdaTsne(celda.CG.sim$counts, celda.CG.mod)
#' @return A two column matrix of t-SNE coordinates
#' @export
setMethod("celdaTsne",
          signature(celda.mod = "celda_CG"),
          function(counts, celda.mod, max.cells=25000, min.cluster.size=100,
                     initial.dims=20, modules=NULL, perplexity=20, max.iter=2500, 
                     seed=12345, ...) {

            ## Checking if max.cells and min.cluster.size will work
            if((max.cells < ncol(counts)) & (max.cells / min.cluster.size < celda.mod@params$K)) {
              stop(paste0("Cannot distribute ", max.cells, " cells among ",
                          celda.mod@params$K, " clusters while maintaining a minumum of ", 
                          min.cluster.size, " cells per cluster. Try increasing 'max.cells' or decreasing 'min.cluster.size'."))
            }
            
            fm = factorizeMatrix(counts=counts, celda.mod=celda.mod, 
                                 type="counts")    
            modules.to.use = 1:nrow(fm$counts$cell)
            if (!is.null(modules)) {
              if (!all(modules %in% modules.to.use)) {
                stop("'modules' must be a vector of numbers between 1 and ", 
                     modules.to.use, ".")
              }
              modules.to.use = modules 
            }
            
          
            ## Select a subset of cells to sample if greater than 'max.cells'
            total.cells.to.remove = ncol(counts) - max.cells
            z.include = rep(TRUE, ncol(counts))
            if(total.cells.to.remove > 0) {
            	z.ta = tabulate(celda.mod@clusters$z, celda.mod@params$K)
            	
            	## Number of cells that can be sampled from each cluster without 
            	## going below the minimum threshold
            	cluster.cells.to.sample = z.ta - min.cluster.size          
            	cluster.cells.to.sample[cluster.cells.to.sample < 0] = 0
            	
            	## Number of cells to sample after exluding smaller clusters
            	## Rounding can cause number to be off by a few, so ceiling is used 
            	## with a second round of subtraction
            	cluster.n.to.sample = ceiling((cluster.cells.to.sample / sum(cluster.cells.to.sample)) * total.cells.to.remove)
            	diff = sum(cluster.n.to.sample) - total.cells.to.remove 
            	cluster.n.to.sample[which.max(cluster.n.to.sample)] = cluster.n.to.sample[which.max(cluster.n.to.sample)] - diff
            
            	## Perform sampling for each cluster
            	for(i in which(cluster.n.to.sample > 0)) {
            	  z.include[sample(which(celda.mod@clusters$z == i), cluster.n.to.sample[i])] = FALSE
            	}
            }   
            cell.ix = which(z.include)
          
            norm = t(normalizeCounts(fm$counts$cell[modules.to.use,cell.ix], normalize="proportion", transformation.fun=sqrt))
            
            res = calculateTsne(norm, do.pca=FALSE, perplexity=perplexity, max.iter=max.iter, seed=seed)
            final = matrix(NA, nrow=ncol(counts), ncol=2)
            final[cell.ix,] = res
            rownames(final) = colnames(counts)
            colnames(final) = c("tsne_1", "tsne_2")
            return(final)
          })



#' @title Probability map for a celda_CG model
#' @description Renders probability and relative expression heatmaps to visualize the relationship between features and cell populations or cell populations and samples.
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`. 
#' @param celda.mod Celda object of class `celda_CG`.   
#' @param level Character. One of 'cell.population' or 'sample'. 'cell.population' will display the absolute probabilities and relative normalized expression of each module in each cell population. 'sample' will display the absolute probabilities and relative normalized abundance of each cell population in each sample. Default 'cell.population'.
#' @param ... Additional parameters.
#' @examples
#' celdaProbabilityMap(celda.CG.sim$counts, celda.CG.mod)
#' @return A grob containing the specified plots
#' @seealso `celda_CG()` for clustering features and cells
#' @export 
setMethod("celdaProbabilityMap",
          signature(celda.mod = "celda_CG"),
          function(counts, celda.mod, level=c("cell.population", "sample"), 
                   ...){
            counts = processCounts(counts)
            compareCountMatrix(counts, celda.mod)
            
            level = match.arg(level)
            factorized <- factorizeMatrix(celda.mod = celda.mod, counts = counts)
            z.include = which(tabulate(celda.mod@clusters$z, celda.mod@params$K) > 0)
            y.include = which(tabulate(celda.mod@clusters$y, celda.mod@params$L) > 0)    
          
            if(level == "cell.population") {
              pop <- factorized$proportions$cell.population[y.include,z.include,drop=FALSE]
              pop.norm = normalizeCounts(pop, normalize="proportion", transformation.fun=sqrt, scale.fun=base::scale)
            
              percentile.9 <- round(stats::quantile(pop,.9), digits = 2) * 100
              col1 <- colorRampPalette(c("#FFFFFF", brewer.pal(n = 9, name = "Blues")))(percentile.9)
              col2 <- colorRampPalette(c("#08306B", c("#006D2C","Yellowgreen","Yellow","Orange","Red")))(100-percentile.9)
              col <- c(col1,col2)
              breaks <-  seq(0, 1, length.out = length(col))     
              
              g1 = plotHeatmap(pop, color.scheme="sequential", scale.row=NULL, cluster.cell=FALSE, cluster.feature=FALSE, show.names.cell=TRUE, show.names.feature=TRUE, breaks = breaks, col=col, main = "Absolute Probability", silent=TRUE)
              g2 = plotHeatmap(pop.norm, color.scheme="divergent", cluster.cell=FALSE, cluster.feature=FALSE, show.names.cell=TRUE, show.names.feature=TRUE, main = "Relative Expression", silent=TRUE) 
              gridExtra::grid.arrange(g1$gtable, g2$gtable, ncol=2)
            } else {
              samp <- factorized$proportions$sample
              col <- colorRampPalette(c("white","blue","#08306B","#006D2C","yellowgreen","yellow","orange","red"))(100)
              breaks <-  seq(0, 1, length.out = length(col))     
              g1 = plotHeatmap(samp, color.scheme="sequential", scale.row=NULL, cluster.cell=FALSE, cluster.feature=FALSE, show.names.cell=TRUE, show.names.feature=TRUE, breaks = breaks, col=col, main = "Absolute Probability", silent=TRUE)
          
              if(ncol(samp) > 1) {
                samp.norm = normalizeCounts(factorized$counts$sample, normalize="proportion", transformation.fun=sqrt, scale.fun=base::scale)
                g2 = plotHeatmap(samp.norm, color.scheme="divergent", cluster.cell=FALSE, cluster.feature=FALSE, show.names.cell=TRUE, show.names.feature=TRUE, main = "Relative Abundance", silent=TRUE)   
                gridExtra::grid.arrange(g1$gtable, g2$gtable, ncol=2)
              } else {
                gridExtra::grid.arrange(g1$gtable)
              } 
            }
          })


#' @title Lookup the module of a feature
#' @description Finds the module assignments of given features in a `celda_G()` model
#'  
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Model of class `celda_CG`.
#' @param feature Character vector. The module assignemnts will be found for feature names in this vector. 
#' @param exact.match Logical. Whether an exact match or a partial match using `grep()` is used to look up the feature in the rownames of the counts matrix. Default TRUE. 
#' @return List. Each element contains the module of the provided feature.
#' @seealso `celda_CG()` for clustering features and cells
#' @examples
#' module = featureModuleLookup(celda.CG.sim$counts, celda.CG.mod, 
#'                              c("Gene_1", "Gene_XXX"))
#' @export
setMethod("featureModuleLookup",
          signature(celda.mod = "celda_CG"),
          function(counts, celda.mod, feature, exact.match = TRUE){
            list <- list()
            if(!isTRUE(exact.match)){
              feature.grep <- c()
              for(x in 1:length(feature)){
                feature.grep <- c(feature.grep, rownames(counts)[grep(feature[x],
                                                                      rownames(counts))]) 
              }
              feature <- feature.grep
            }
            for(x in 1:length(feature)){
              if(feature[x] %in% rownames(counts)){
                list[x] <- celda.mod@clusters$y[which(rownames(counts) == feature[x])]
              }else{
                list[x] <- paste0("No feature was identified matching '", feature[x], "'.")
              }
            } 
            names(list) <- feature
            return(list)
         })
