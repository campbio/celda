#' @title Feature clustering with Celda
#'
#' @description Clusters the rows of a count matrix containing single-cell data into L modules. 
#'
#' @param counts Integer matrix. Rows represent features and columns represent cells. 
#' @param L Integer. Number of feature modules.  
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to each feature module in each cell. Default 1. 
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to each feature in each module. Default 1. 
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to the number of features in each module. Default 1. 
#' @param stop.iter Integer. Number of iterations without improvement in the log likelihood to stop inference. Default 10.
#' @param max.iter Integer. Maximum number of iterations of Gibbs sampling to perform. Default 200.
#' @param split.on.iter Integer. On every `split.on.iter` iteration, a heuristic will be applied to determine if a feature module should be reassigned and another feature module should be split into two clusters. To disable splitting, set to -1. Default 10.
#' @param split.on.last Integer. After `stop.iter` iterations have been performed without improvement, a heuristic will be applied to determine if a cell population should be reassigned and another cell population should be split into two clusters. If a split occurs, then `stop.iter` will be reset. Default TRUE.
#' @param seed Integer. Passed to `set.seed()`. Default 12345.  
#' @param nchains Integer. Number of random cluster initializations. Default 3.  
#' @param initialize Chararacter. One of 'random' or 'split'. With 'random', features are randomly assigned to a clusters. With 'split' cell and feature clusters will be recurssively split into two clusters using `celda_G()` until the specified L is reached. Default 'random'.
#' @param count.checksum Character. An MD5 checksum for the `counts` matrix. Default NULL.
#' @param y.init Integer vector. Sets initial starting values of y. If NULL, starting values for each feature will be randomly sampled from `1:L`. `y.init` can only be used when `initialize = 'random'`. Default NULL.
#' @param logfile Character. Messages will be redirected to a file named `logfile`. If NULL, messages will be printed to stdout.  Default NULL.
#' @param verbose Logical. Whether to print log messages. Default TRUE. 
#' @return An object of class `celda_G` with the feature module clusters stored in `y`.
#' @seealso `celda_C()` for cell clustering and `celda_CG()` for simultaneous clustering of features and cells. `celdaGridSearch()` can be used to run multiple values of L and multiple chains in parallel. 
#' @examples
#' celda.mod = celda_G(celda.G.sim$counts, L=celda.G.sim$L)
#' @export
celda_G = function(counts, L, beta=1, delta=1, gamma=1,
					stop.iter=10, max.iter=200, split.on.iter=10, split.on.last=TRUE,
					seed=12345, nchains=3, initialize=c("random", "split"), count.checksum=NULL, 
					y.init=NULL, logfile=NULL, verbose=TRUE) {
  
  validateCounts(counts)
  return(.celda_G(counts, L, beta, delta, gamma, stop.iter, max.iter, split.on.iter,
                  split.on.last, seed, nchains, initialize, count.checksum,
                  y.init, logfile, verbose))
}

.celda_G = function(counts, L, beta=1, delta=1, gamma=1,
					stop.iter=10, max.iter=200, split.on.iter=10, split.on.last=TRUE,
					seed=12345, nchains=3, initialize=c("random", "split"), count.checksum=NULL, 
					y.init=NULL, logfile=NULL, verbose=TRUE) {

  logMessages("--------------------------------------------------------------------", logfile=logfile, append=FALSE, verbose=verbose)  
  logMessages("Starting Celda_G: Clustering genes.", logfile=logfile, append=TRUE, verbose=verbose)
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose)  
  start.time = Sys.time()

  ## Error checking and variable processing
  counts = processCounts(counts)
  if(is.null(count.checksum)) {
    count.checksum = digest::digest(counts, algo="md5")
  }
  initialize = match.arg(initialize)
   
  all.seeds = seed:(seed + nchains - 1)
    
  best.result = NULL  
  for(i in seq_along(all.seeds)) {   

	## Randomly select y or y to supplied initial values
	## Initialize cluster labels
    current.seed = all.seeds[i]	
    logMessages(date(), ".. Initializing chain", i, "with", paste0("'",initialize, "' (seed=", current.seed, ")"), logfile=logfile, append=TRUE, verbose=verbose)

    if(initialize == "random") {
	  y = initialize.cluster(L, nrow(counts), initial = y.init, fixed = NULL, seed=current.seed)
	} else {
	  y = recursive.splitY(counts, L, beta=beta, delta=delta, gamma=gamma, z=NULL, K=NULL, K.subclusters=NULL, min.feature=3, max.cells=100, seed=seed)
	}  
	y.best = y  

	## Calculate counts one time up front
	p = cG.decomposeCounts(counts=counts, y=y, L=L)
	n.TS.by.C = p$n.TS.by.C
	n.by.G = p$n.by.G
	n.by.TS = p$n.by.TS
	nG.by.TS = p$nG.by.TS
	nM = p$nM
	nG = p$nG
	rm(p)

	set.seed(seed)
  
	## Calculate initial log likelihood
	ll <- cG.calcLL(n.TS.by.C=n.TS.by.C, n.by.TS=n.by.TS, n.by.G=n.by.G, nG.by.TS=nG.by.TS, nM=nM, nG=nG, L=L, beta=beta, delta=delta, gamma=gamma)

	iter <- 1L
	num.iter.without.improvement = 0L
	do.gene.split = TRUE
	while(iter <= max.iter & num.iter.without.improvement <= stop.iter) {

	  next.y = cG.calcGibbsProbY(counts=counts, n.TS.by.C=n.TS.by.C, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, n.by.G=n.by.G, y=y, nG=nG, L=L, beta=beta, delta=delta, gamma=gamma)
	  n.TS.by.C = next.y$n.TS.by.C
	  nG.by.TS = next.y$nG.by.TS
	  n.by.TS = next.y$n.by.TS
	  y = next.y$y
	
	  ## Perform split on i-th iteration of no improvement in log likelihood
	  temp.ll <- cG.calcLL(n.TS.by.C=n.TS.by.C, n.by.TS=n.by.TS, n.by.G=n.by.G, nG.by.TS=nG.by.TS, nM=nM, nG=nG, L=L, beta=beta, delta=delta, gamma=gamma)
	  if(L > 2 & (((iter == max.iter | (num.iter.without.improvement == stop.iter & all(temp.ll < ll))) & isTRUE(split.on.last)) | (split.on.iter > 0 & iter %% split.on.iter == 0 & isTRUE(do.gene.split)))) {
		logMessages(date(), " .... Determining if any gene clusters should be split.", logfile=logfile, append=TRUE, sep="", verbose=verbose)
		res = cG.splitY(counts, y, n.TS.by.C, n.by.TS, n.by.G, nG.by.TS, nM, nG, L, 
		                beta, delta, gamma, y.prob=t(next.y$probs), 
		                min.feature=3, max.clusters.to.try=max(L/2, 10))
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
		n.TS.by.C = res$n.TS.by.C
		n.by.TS = res$n.by.TS
		nG.by.TS = res$nG.by.TS
	  }
	 
	  ## Calculate complete likelihood
	  temp.ll <- cG.calcLL(n.TS.by.C=n.TS.by.C, n.by.TS=n.by.TS, n.by.G=n.by.G, nG.by.TS=nG.by.TS, nM=nM, nG=nG, L=L, beta=beta, delta=delta, gamma=gamma)
	  if((all(temp.ll > ll)) | iter == 1) {
		y.best = y
		ll.best = temp.ll
		num.iter.without.improvement = 1L
	  } else {  
		num.iter.without.improvement = num.iter.without.improvement + 1L   
	  }
	  ll <- c(ll, temp.ll)

	  logMessages(date(), ".... Completed iteration:", iter, "| logLik:", temp.ll, logfile=logfile, append=TRUE, verbose=verbose)
	  iter = iter + 1    	  
    }
    
	names = list(row=rownames(counts), column=colnames(counts))  

	result = list(y=y.best, completeLogLik=ll, 
				  finalLogLik=ll.best, L=L, beta=beta, delta=delta, gamma=gamma,
				  count.checksum=count.checksum, seed=current.seed, names=names)
	
	if (is.null(best.result) || result$finalLogLik > best.result$finalLogLik) {
      best.result = result
    }
    
    logMessages(date(), ".. Finished chain", i, "with seed", current.seed, logfile=logfile, append=TRUE, verbose=verbose)
  } 
  
  
  best.result = methods::new("celda_G", 
                             clusters=list(y=y.best),
                             params=list(L=L, beta=beta, delta=delta, gamma=gamma,
                                         count.checksum=count.checksum, seed=current.seed),
                             completeLogLik=ll, finalLogLik=ll.best, 
                             names=names)
  best.result = reorder.celda_G(counts = counts, res = best.result) 
  
  end.time = Sys.time()
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose)  
  logMessages("Completed Celda_G. Total time:", format(difftime(end.time, start.time)), logfile=logfile, append=TRUE, verbose=verbose)
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose)  

  return(best.result)
}


# Calculate Log Likelihood For Single Set of Cluster Assignments (Gene Clustering)
#
# This function calculates the log-likelihood of a given set of cluster assigments for the samples
# represented in the provided count matrix.
# 
# 
# @param n.TS.by.C Number of counts in each Transcriptional State per Cell.
# @param n.by.TS Number of counts per Transcriptional State.
# @param nG.by.TS Number of genes in each Transcriptional State.
# @param nG.in.Y  Number of genes in each of the cell cluster.
# @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to the number of features in each module. Default 1. 
# @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to each feature in each module. Default 1. 
# @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to each feature module in each cell. Default 1. 
# @keywords log likelihood
cG.calcGibbsProbY = function(counts, n.TS.by.C, n.by.TS, nG.by.TS, n.by.G, y, L, nG, beta, delta, gamma, do.sample=TRUE) {

  ## Set variables up front outside of loop
  probs = matrix(NA, ncol=nG, nrow=L)
  temp.nG.by.TS = nG.by.TS 
  temp.n.by.TS = n.by.TS 
#  temp.n.TS.by.C = n.TS.by.C

  ix <- sample(1:nG)
  for(i in ix) {
	  
	## Subtract current gene counts from matrices
	nG.by.TS[y[i]] = nG.by.TS[y[i]] - 1L
	n.by.TS[y[i]] = n.by.TS[y[i]] - n.by.G[i]
	
	n.TS.by.C_1 = n.TS.by.C
	n.TS.by.C_1[y[i],] = n.TS.by.C_1[y[i],] - counts[i,]
    n.TS.by.C_1 = .rowSums(lgamma(n.TS.by.C_1 + beta), nrow(n.TS.by.C), ncol(n.TS.by.C))
    
	n.TS.by.C_2 = n.TS.by.C
	other.ix = (1:L)[-y[i]]
	n.TS.by.C_2[other.ix,] = n.TS.by.C_2[other.ix,] + counts[rep(i, length(other.ix)),]
    n.TS.by.C_2 = .rowSums(lgamma(n.TS.by.C_2 + beta), nrow(n.TS.by.C), ncol(n.TS.by.C))
    
	## Calculate probabilities for each state
	for(j in 1:L) {
	
	  temp.nG.by.TS = nG.by.TS 
	  temp.n.by.TS = n.by.TS 
#	  temp.n.TS.by.C = n.TS.by.C
	
	  temp.nG.by.TS[j] = temp.nG.by.TS[j] + 1L      
	  temp.n.by.TS[j] = temp.n.by.TS[j] + n.by.G[i]
#	  temp.n.TS.by.C[j,] = temp.n.TS.by.C[j,] + counts[i,]


      pseudo.nG.by.TS = temp.nG.by.TS
	  pseudo.nG.by.TS[temp.nG.by.TS == 0L] = 1L
	  pseudo.nG = sum(pseudo.nG.by.TS)

	  other.ix = (1:L)[-j]
      probs[j,i] <- 	sum(lgamma(pseudo.nG.by.TS + gamma)) -
						sum(lgamma(sum(pseudo.nG.by.TS + gamma))) +
						sum(n.TS.by.C_1[other.ix]) + n.TS.by.C_2[j] +
#						sum(lgamma(temp.n.TS.by.C + beta)) +
						sum(lgamma(pseudo.nG.by.TS * delta)) -
						(pseudo.nG * lgamma(delta)) -
						sum(lgamma(temp.n.by.TS + (pseudo.nG.by.TS * delta)))

    }

	## Sample next state and add back counts
	prev.y = y[i]
	if(isTRUE(do.sample)) y[i] = sample.ll(probs[,i])
	
	if(prev.y != y[i]) {
	  n.TS.by.C[prev.y,] = n.TS.by.C[prev.y,] - counts[i,]
	  n.TS.by.C[y[i],] = n.TS.by.C[y[i],] + counts[i,]
	}
	
	nG.by.TS[y[i]] = nG.by.TS[y[i]] + 1L
	n.by.TS[y[i]] = n.by.TS[y[i]] + n.by.G[i]
	#n.TS.by.C[y[i],] = n.TS.by.C[y[i],] + counts[i,]
  }
  
  return(list(n.TS.by.C=n.TS.by.C, nG.by.TS=nG.by.TS, n.by.TS=n.by.TS, y=y, probs=probs))
}

#' @title Simulate cells from the celda_G model
#' 
#' @description Generates a simulated counts matrix and feature module clusters
#' according to the generative process of the celda_G model.  
#' 
#' @param model Character. Options available in `celda::available.models`. 
#' @param C Integer. Number of cells to simulate. Default 100. 
#' @param L Integer. Number of feature modules. Default 10.
#' @param N.Range Integer vector. A vector of length 2 that specifies the lower and upper bounds of the number of counts generated for each cell. Default c(500, 5000). 
#' @param G Integer. The total number of features to be simulated. Default 100. 
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to each feature module in each cell. Default 1. 
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to each feature in each module. Default 1. 
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to the number of features in each module. Default 5. 
#' @param seed Integer. Passed to `set.seed()`. Default 12345.  
#' @param ... Additional parameters.
#' @return List. Contains the simulated matrix `counts`, feature module clusters `y`, and input parameters.
#' @seealso `celda_C()` for simulating cell subpopulations and `celda_CG()` for simulating feature modules and cell populations. 
#' @examples
#' celda.g.sim = simulateCells(model="celda_G")
#' @export
simulateCells.celda_G = function(model, C=100, N.Range=c(500,1000), G=100, 
                                 L=10, beta=1, gamma=5, delta=1, seed=12345, ...) {
  set.seed(seed)
  eta = rdirichlet(1, rep(gamma, L))
  
  y = sample(1:L, size=G, prob=eta, replace=TRUE)
  if(length(table(y)) < L) {
    stop("Some states did not receive any genes after sampling. Try increasing G and/or setting gamma > 1.")
  }
  
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
  cell.counts = processCounts(cell.counts)
  names = list(row=rownames(cell.counts), column=colnames(cell.counts))
  result = methods::new("celda_G", clusters=list(y=y), 
                        params=list(L=L, beta=beta, delta=delta, gamma=gamma,
                                    seed=seed,
                                    count.checksum=digest::digest(cell.counts, algo="md5")),
                        names=names)
  result = reorder.celda_G(counts = cell.counts, res = result)  
  
  return(list(y=result@clusters$y, counts=processCounts(cell.counts), L=L, 
              beta=beta, delta=delta, gamma=gamma, seed=seed))
}


#' @title Matrix factorization for results from celda_G
#' @description Generates factorized matrices showing the contribution of each feature in each module and each module in each cell.
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class "celda_G". 
#' @param type Character vector. A vector containing one or more of "counts", "proportion", or "posterior". "counts" returns the raw number of counts for each factorized matrix. "proportions" returns the normalized probabilities for each factorized matrix, which are calculated by dividing the raw counts in each factorized matrix by the total counts in each column. "posterior" returns the posterior estimates. Default `c("counts", "proportion", "posterior")`. 
#' @return A list with elements for `counts`, `proportions`, or `posterior` probabilities. Each element will be a list containing factorized matrices for `module` and `cell`.
#' @seealso `celda_G()` for clustering features
#' @examples 
#' factorized.matrices = factorizeMatrix(celda.G.sim$counts, celda.G.mod, "posterior")
#' @export
setMethod("factorizeMatrix",
          signature(celda.mod = "celda_G"),
           function(counts, celda.mod, 
                    type=c("counts", "proportion", "posterior")) {
            counts = processCounts(counts)
            compareCountMatrix(counts, celda.mod)
            
            L = celda.mod@params$L
            y = celda.mod@clusters$y
            beta = celda.mod@params$beta
            delta = celda.mod@params$delta
            gamma = celda.mod@params$gamma
            
            p = cG.decomposeCounts(counts=counts, y=y, L=L)
            n.TS.by.C = p$n.TS.by.C
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
            colnames(n.TS.by.C) = celda.mod@names$column
            rownames(n.TS.by.C) = L.names
            colnames(n.G.by.TS) = L.names
            rownames(n.G.by.TS) = celda.mod@names$row
            names(nG.by.TS) = L.names
            
            counts.list = c()
            prop.list = c()
            post.list = c()
            res = list()
            
            if(any("counts" %in% type)) {
              counts.list = list(cell=n.TS.by.C, module=n.G.by.TS, 
                                 gene.distribution=nG.by.TS)
              res = c(res, list(counts=counts.list))
            }
            if(any("proportion" %in% type)) {
              ## Need to avoid normalizing cell/gene states with zero cells/genes
              unique.y = sort(unique(y))
              temp.n.G.by.TS = n.G.by.TS
              temp.n.G.by.TS[,unique.y] = normalizeCounts(temp.n.G.by.TS[,unique.y], 
                                                          normalize="proportion")
              temp.nG.by.TS = nG.by.TS/sum(nG.by.TS)
              
              prop.list = list(cell = normalizeCounts(n.TS.by.C, normalize="proportion"),
              							  module = temp.n.G.by.TS, gene.distribution=temp.nG.by.TS)
              res = c(res, list(proportions=prop.list))
            }
            if(any("posterior" %in% type)) {
            
              gs = n.G.by.TS
              gs[cbind(1:nG,y)] = gs[cbind(1:nG,y)] + delta
              gs = normalizeCounts(gs, normalize="proportion")
              temp.nG.by.TS = (nG.by.TS + gamma)/sum(nG.by.TS + gamma)
              
              post.list = list(cell = normalizeCounts(n.TS.by.C + beta, normalize="proportion"),
              						     module = gs, gene.distribution=temp.nG.by.TS)
              res = c(res, posterior = list(post.list))						    
            }
            
            return(res)
          })


# Calculate log-likelihood of celda_CG model
cG.calcLL = function(n.TS.by.C, n.by.TS, n.by.G, nG.by.TS, nM, nG, L, beta, delta, gamma) {
  
  nG.by.TS[nG.by.TS == 0] = 1
  nG = sum(nG.by.TS)
  
  ## Calculate for "Phi" component
  a <- nM * lgamma(L * beta)
  b <- sum(lgamma(n.TS.by.C + beta))
  c <- -nM * L * lgamma(beta)
  d <- -sum(lgamma(colSums(n.TS.by.C + beta)))

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


#' @title Calculate Celda_G log likelihood
#' @description Calculates the log likelihood for user-provided feature module clusters using the `celda_G()` model.
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. 
#' @param y Numeric vector. Denotes feature module labels. 
#' @param L Integer. Number of feature modules.  
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to each feature module in each cell. Default 1. 
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to each feature in each module. Default 1. 
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to the number of features in each module. Default 1. 
#' @param ... Additional parameters.
#' @keywords log likelihood
#' @return The log-likelihood for the given cluster assignments
#' @seealso `celda_G()` for clustering features
#' @examples
#' loglik = logLikelihood(celda.G.sim$counts, model="celda_G", 
#'                        y=celda.G.sim$y, L=celda.G.sim$L,
#'                        beta=celda.G.sim$beta, delta=celda.G.sim$delta,
#'                        gamma=celda.G.sim$gamma)
#' @export
logLikelihood.celda_G = function(counts, y, L, beta, delta, gamma) {
  if (sum(y > L) > 0) stop("An entry in y contains a value greater than the provided L.")
  p = cG.decomposeCounts(counts=counts, y=y, L=L)
  final <- cG.calcLL(n.TS.by.C=p$n.TS.by.C, n.by.TS=p$n.by.TS, n.by.G=p$n.by.G, nG.by.TS=p$nG.by.TS, nM=p$nM, nG=p$nG, L=L, beta=beta, delta=delta, gamma=gamma)
  
  return(final)
}


# Takes raw counts matrix and converts it to a series of matrices needed for log likelihood calculation
# @param counts Integer matrix. Rows represent features and columns represent cells. 
# @param y Numeric vector. Denotes feature module labels. 
# @param L Integer. Number of feature modules.  
cG.decomposeCounts = function(counts, y, L) {

  n.TS.by.C = rowSumByGroup(counts, group=y, L=L)
  n.by.G = as.integer(.rowSums(counts, nrow(counts), ncol(counts)))
  n.by.TS = as.integer(rowSumByGroup(matrix(n.by.G,ncol=1), group=y, L=L))
  nG.by.TS = tabulate(y, L)
  nM = ncol(counts)
  nG = nrow(counts)
  
  return(list(n.TS.by.C=n.TS.by.C, n.by.G=n.by.G, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, nM=nM, nG=nG))
}


cG.reDecomposeCounts = function(counts, y, previous.y, n.TS.by.C, n.by.G, L) {
  ## Recalculate counts based on new label
  n.TS.by.C = rowSumByGroupChange(counts, n.TS.by.C, y, previous.y, L)
  n.by.TS = as.integer(rowSumByGroup(matrix(n.by.G,ncol=1), group=y, L=L))
  nG.by.TS = tabulate(y, L)

  return(list(n.TS.by.C=n.TS.by.C, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS))  
}


#' @title Conditional probabilities for features in modules from a Celda_G model
#' @description Calculates the conditional probability of each feature belonging to each module given all other feature cluster assignments in a `celda_G()` result. 
#'
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class `celda_G`. 
#' @param log Logical. If FALSE, then the normalized conditional probabilities will be returned. If TRUE, then the unnormalized log probabilities will be returned. Default FALSE.  
#' @param ... Additional parameters.
#' @return A list containging a matrix for the conditional cell cluster probabilities. 
#' @seealso `celda_G()` for clustering features
#' @examples
#' cluster.prob = clusterProbability(celda.G.sim$counts, celda.G.mod)
#' @export
setMethod("clusterProbability",
           signature(celda.mod = "celda_G"),
           function(counts, celda.mod, log=FALSE, ...) {
              y = celda.mod@clusters$y
              L = celda.mod@params$L
              delta = celda.mod@params$delta
              beta = celda.mod@params$beta
              gamma = celda.mod@params$gamma
              
              ## Calculate counts one time up front
              p = cG.decomposeCounts(counts=counts, y=y, L=L)
              next.y = cG.calcGibbsProbY(counts=counts, n.TS.by.C=p$n.TS.by.C, 
                                         n.by.TS=p$n.by.TS, nG.by.TS=p$nG.by.TS, 
                                         n.by.G=p$n.by.G, y=y, nG=p$nG, L=L, 
                                         beta=beta, delta=delta, gamma=gamma, 
                                         do.sample=FALSE)  
              y.prob = t(next.y$probs)
              
              if(!isTRUE(log)) {
                y.prob = normalizeLogProbs(y.prob)
              }
              
              return(list(y.probability=y.prob))
          })


#' @title Calculate the perplexity on new data with a celda_G model
#' @description Perplexity is a statistical measure of how well a probability model can predict new data. Lower perplexity indicates a better model. 
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class "celda_C"
#' @param new.counts A new counts matrix used to calculate perplexity. If NULL, perplexity will be calculated for the 'counts' matrix. Default NULL.
#' @return Numeric. The perplexity for the provided count data and model.
#' @seealso `celda_G()` for clustering features
#' @examples
#' perplexity = perplexity(celda.G.sim$counts, celda.G.mod)
#' @export
setMethod("perplexity",
          signature(celda.mod = "celda_G"),
          function(counts, celda.mod, new.counts=NULL) {
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
            phi <- factorized$posterior$module
            psi <- factorized$posterior$cell
            eta <- factorized$posterior$gene.distribution
            nG.by.TS = factorized$counts$gene.distribution
            
            eta.prob = log(eta) * nG.by.TS
            gene.by.cell.prob = log(phi %*% psi) 
            log.px = sum(eta.prob) + sum(gene.by.cell.prob * new.counts)
            
            perplexity = exp(-(log.px/sum(new.counts)))
            return(perplexity)
          })


reorder.celda_G = function(counts, res) {
  if(res@params$L > 2 & isTRUE(length(unique(res@clusters$y)) > 1)) {
    res@clusters$y = as.integer(as.factor(res@clusters$y))
    fm <- factorizeMatrix(counts = counts, celda.mod = res)
    unique.y = sort(unique(res@clusters$y))
    cs = prop.table(t(fm$posterior$cell[unique.y,]), 2)
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
#' @param celda.mod Celda object of class `celda_G`. 
#' @param nfeatures Integer. Maximum number of features to select for each module. Default 25.
#' @param ... Additional parameters.
#' @seealso `celda_G()` for clustering features and `celdaTsne()` for generating 2-dimensional coordinates
#' @examples
#' celdaHeatmap(celda.G.sim$counts, celda.G.mod)
#' @return list A list containing the dendrograms and the heatmap grob
#' @export
setMethod("celdaHeatmap",
          signature(celda.mod = "celda_G"),
          function(counts, celda.mod, nfeatures=25, ...) {
            fm = factorizeMatrix(counts, celda.mod, type="proportion")
            top = topRank(fm$proportions$module, n=nfeatures)
            ix = unlist(top$index)
            norm = normalizeCounts(counts, normalize="proportion", transformation.fun=sqrt)
            plotHeatmap(norm[ix,], y=celda.mod@clusters$y[ix], ...)
          })

#' @title tSNE for celda_G
#' @description Embeds cells in two dimensions using tSNE based on a `celda_G` model. tSNE is run on module probabilities to reduce the number of features instead of using PCA. Module probabilities square-root trasformed before applying tSNE. 
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class `celda_G`.  
#' @param max.cells Integer. Maximum number of cells to plot. Cells will be randomly subsampled if ncol(conts) > max.cells. Larger numbers of cells requires more memory. Default 10000.
#' @param modules Integer vector. Determines which feature modules to use for tSNE. If NULL, all modules will be used. Default NULL.
#' @param perplexity Numeric. Perplexity parameter for tSNE. Default 20.
#' @param max.iter Integer. Maximum number of iterations in tSNE generation. Default 2500.
#' @param seed Integer. Passed to `set.seed()`. Default 12345.  
#' @param ... Additional parameters.
#' @seealso `celda_G()` for clustering features and `celdaHeatmap()` for displaying expression
#' @examples
#' tsne.res = celdaTsne(celda.G.sim$counts, celda.G.mod)
#' @return A two column matrix of t-SNE coordinates
#' @export
setMethod("celdaTsne",
          signature(celda.mod = "celda_G"),
          function(counts, celda.mod, max.cells=25000, min.cluster.size=100,
                    initial.dims=20, modules=NULL, perplexity=20, max.iter=2500, 
                    seed=12345, ...) {
  
            if(max.cells > ncol(counts)) {
              max.cells = ncol(counts)
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
             
            cell.ix = sample(1:ncol(counts), max.cells)
            norm = t(normalizeCounts(fm$counts$cell[modules.to.use,cell.ix], 
                                     normalize="proportion", 
                                     transformation.fun=sqrt))
          
            res = calculateTsne(norm, do.pca=FALSE, perplexity=perplexity, 
                                max.iter=max.iter, seed=seed)
            rownames(res) = colnames(counts)
            colnames(res) = c("tsne_1", "tsne_2")
            return(res)
          })


#' @title Lookup the module of a feature
#' @description Finds the module assignments of given features in a `celda_G()` model
#'  
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Model of class `celda_G`.
#' @param feature Character vector. The module assignemnts will be found for feature names in this vector. 
#' @param exact.match Logical. Whether an exact match or a partial match using `grep()` is used to look up the feature in the rownames of the counts matrix. Default TRUE. 
#' @return List. Each element contains the module of the provided feature.
#' @seealso `celda_G()` for clustering features
#' @examples
#' module = featureModuleLookup(celda.G.sim$counts, celda.G.mod, 
#'                              c("Gene_1", "Gene_XXX"))
#' @export
setMethod("featureModuleLookup",
          signature(celda.mod = "celda_G"),
          function(counts, celda.mod, feature, exact.match = TRUE){
            list <- list()
            if(!isTRUE(exact.match)){
              feature.grep <- c()
              for(x in 1:length(feature)){
                feature.grep <- c(feature.grep, 
                                  rownames(counts)[grep(feature[x],rownames(counts))]) 
              }
              feature <- feature.grep
            }
            for(x in 1:length(feature)){
              if(feature[x] %in% rownames(counts)){
                list[x] <- celda.mod@clusters$y[which(rownames(counts) == feature[x])]
              }else{
                list[x] <- paste0("No feature was identified matching '", 
                                  feature[x], "'.")
              }
            } 
            names(list) <- feature
            return(list)
         })
