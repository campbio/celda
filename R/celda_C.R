#' @title Cell clustering with Celda
#' @description Clusters the columns of a count matrix containing single-cell data into K subpopulations. 
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. 
#' @param sample.label Vector or factor. Denotes the sample label for each cell (column) in the count matrix.
#' @param K Integer. Number of cell populations. 
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount to each cell population in each sample. Default 1. 
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to each feature in each cell population. Default 1. 
#' @param algorithm String. Algorithm to use for clustering cell subpopulations. One of 'EM' or 'Gibbs'. The EM algorithm is faster, especially for larger numbers of cells. However, more chains may be required to ensure a good solution is found. Default 'EM'.
#' @param stop.iter Integer. Number of iterations without improvement in the log likelihood to stop inference. Default 10.
#' @param max.iter Integer. Maximum number of iterations of Gibbs sampling or EM to perform. Default 200.
#' @param split.on.iter Integer. On every `split.on.iter` iteration, a heuristic will be applied to determine if a cell population should be reassigned and another cell population should be split into two clusters. To disable splitting, set to -1. Default 10.
#' @param split.on.last Integer. After `stop.iter` iterations have been performed without improvement, a heuristic will be applied to determine if a cell population should be reassigned and another cell population should be split into two clusters. If a split occurs, then `stop.iter` will be reset. Default TRUE.
#' @param seed Integer. Passed to `set.seed()`. Default 12345.  
#' @param nchains Integer. Number of random cluster initializations. Default 3.  
#' @param initialize Chararacter. One of 'random' or 'split'. With 'random', cells are randomly assigned to a clusters. With 'split' cell clusters will be recurssively split into two clusters using `celda_C` until the specified K is reached. Default 'random'.
#' @param count.checksum "Character. An MD5 checksum for the `counts` matrix. Default NULL.
#' @param z.init Integer vector. Sets initial starting values of z. If NULL, starting values for each cell will be randomly sampled from `1:K`. 'z.init' can only be used when `initialize = 'random'`. Default NULL.
#' @param logfile Character. Messages will be redirected to a file named `logfile`. If NULL, messages will be printed to stdout.  Default NULL.
#' @param verbose Logical. Whether to print log messages. Default TRUE. 
#' @return An object of class `celda_C` with the cell population clusters stored in in `z`.
#' @seealso `celda_G()` for feature clustering and `celda_CG()` for simultaneous clustering of features and cells. `celdaGridSearch()` can be used to run multiple values of K and multiple chains in parallel. 
#' @examples
#' celda.mod = celda_C(celda.C.sim$counts, K=celda.C.sim$K, 
#'                     sample.label=celda.C.sim$sample.label)
#' @export
celda_C = function(counts, sample.label=NULL, K, alpha=1, beta=1,
  					        algorithm = c("EM", "Gibbs"), 
                   	stop.iter = 10, max.iter=200, split.on.iter=10, split.on.last=TRUE,
                   	seed=12345, nchains=3, initialize=c("random", "split"), count.checksum=NULL, 
                   	z.init = NULL, logfile=NULL, verbose=TRUE) {
  validateCounts(counts)
  return(.celda_C(counts, sample.label, K, alpha, beta, algorithm, stop.iter,
                  max.iter, split.on.iter, split.on.last, seed, nchains,
                  initialize, count.checksum, z.init, logfile, verbose))
}

.celda_C = function(counts, sample.label=NULL, K, alpha=1, beta=1,
  					        algorithm = c("EM", "Gibbs"), 
                   	stop.iter = 10, max.iter=200, split.on.iter=10, split.on.last=TRUE,
                   	seed=12345, nchains=3, initialize=c("random", "split"), count.checksum=NULL, 
                   	z.init = NULL, logfile=NULL, verbose=TRUE) {
  
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=FALSE, verbose=verbose)  
  logMessages("Starting Celda_C: Clustering cells.", logfile=logfile, append=TRUE, verbose=verbose)
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose)  
  start.time = Sys.time()
                 	 
  ## Error checking and variable processing
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
	logMessages(date(), ".. Initializing chain", i, "with", paste0("'", initialize, "' (seed=", current.seed, ")"), logfile=logfile, append=TRUE, verbose=verbose)
	
    if(initialize == "random") {
  	  z = initialize.cluster(K, ncol(counts), initial = z.init, fixed = NULL, seed=current.seed)
	} else {
	  z = recursive.splitZ(counts, s, K=K, alpha=alpha, beta=beta)
	}  
	z.best = z
  
	## Calculate counts one time up front
	p = cC.decomposeCounts(counts, s, z, K)
	nS = p$nS
	nG = p$nG
	nM = p$nM
	m.CP.by.S = p$m.CP.by.S
	n.G.by.CP = p$n.G.by.CP
	n.CP = p$n.CP
	n.by.C = p$n.by.C
  
	ll = cC.calcLL(m.CP.by.S=m.CP.by.S, n.G.by.CP=n.G.by.CP, s=s, K=K, nS=nS, nG=nG, alpha=alpha, beta=beta)

    set.seed(seed)
	iter = 1L
	num.iter.without.improvement = 0L
	do.cell.split = TRUE
	while(iter <= max.iter & num.iter.without.improvement <= stop.iter) {
	
	  next.z = do.call(algorithm.fun, list(counts=counts, m.CP.by.S=m.CP.by.S, n.G.by.CP=n.G.by.CP, n.by.C=n.by.C, n.CP=n.CP, z=z, s=s, K=K, nG=nG, nM=nM, alpha=alpha, beta=beta))

	  m.CP.by.S = next.z$m.CP.by.S
	  n.G.by.CP = next.z$n.G.by.CP
	  n.CP = next.z$n.CP
	  z = next.z$z

	  ## Perform split on i-th iteration of no improvement in log likelihood
	  temp.ll = cC.calcLL(m.CP.by.S=m.CP.by.S, n.G.by.CP=n.G.by.CP, s=s, K=K, nS=nS, nG=nG, alpha=alpha, beta=beta)
	  if(K > 2 & (((iter == max.iter | (num.iter.without.improvement == stop.iter & all(temp.ll < ll))) & isTRUE(split.on.last)) | (split.on.iter > 0 & iter %% split.on.iter == 0 & isTRUE(do.cell.split)))) {

		logMessages(date(), " .... Determining if any cell clusters should be split.", logfile=logfile, append=TRUE, sep="", verbose=verbose)
		res = cC.splitZ(counts, m.CP.by.S, n.G.by.CP, n.CP, s, z, K, nS, nG, alpha, beta, z.prob=t(next.z$probs), max.clusters.to.try=K, min.cell=3)
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
		n.G.by.CP = res$n.G.by.CP
		n.CP = res$n.CP
	  }


	  ## Calculate complete likelihood
	  temp.ll = cC.calcLL(m.CP.by.S=m.CP.by.S, n.G.by.CP=n.G.by.CP, s=s, K=K, nS=nS, nG=nG, alpha=alpha, beta=beta)
	  if((all(temp.ll > ll)) | iter == 1) {
		z.best = z
		ll.best = temp.ll
		num.iter.without.improvement = 1L
	  } else {  
		num.iter.without.improvement = num.iter.without.improvement + 1L   
	  }
	  ll = c(ll, temp.ll)
	
	  logMessages(date(), ".... Completed iteration:", iter, "| logLik:", temp.ll, logfile=logfile, append=TRUE, verbose=verbose)
	  iter = iter + 1    
	}
	
	names = list(row=rownames(counts), column=colnames(counts), sample=levels(sample.label))

	result = list(z=z.best, completeLogLik=ll,  
      				  finalLogLik=ll.best, seed=current.seed, K=K, 
      				  sample.label=sample.label, alpha=alpha, 
      				  beta=beta, count.checksum=count.checksum,
      				  names=names)
	
    if(is.null(best.result) || result$finalLogLik > best.result$finalLogLik) {
      best.result = result
    }
    
    logMessages(date(), ".. Finished chain", i, "with seed", current.seed, logfile=logfile, append=TRUE, verbose=verbose)
  }  
  best.result = methods::new("celda_C",
                             clusters=list(z=best.result$z),
                             params=list(K=best.result$K, alpha=best.result$alpha, beta=best.result$beta,
                                         count.checksum=best.result$count.checksum,
                                         seed=best.result$seed),
                             sample.label=best.result$sample.label, 
                             completeLogLik=best.result$completeLogLik,
                             finalLogLik=best.result$finalLogLik,
                             names=best.result$names)
  best.result = reorder.celda_C(counts = counts, res = best.result)
  
  end.time = Sys.time()
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose)  
  logMessages("Completed Celda_C. Total time:", format(difftime(end.time, start.time)), logfile=logfile, append=TRUE, verbose=verbose)
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose)  

  return(best.result)
}


# Gibbs sampling for the celda_C Model
cC.calcGibbsProbZ = function(counts, m.CP.by.S, n.G.by.CP, n.by.C, n.CP, z, s, K, nG, nM, alpha, beta, do.sample=TRUE) {

  ## Set variables up front outside of loop  
  probs = matrix(NA, ncol=nM, nrow=K)

  ix = sample(1:nM)
  for(i in ix) {

	## Subtract cell counts from current population assignment
	n.G.by.CP_1 = n.G.by.CP
	n.G.by.CP_1[,z[i]] = n.G.by.CP[,z[i]] - counts[,i]
	n.G.by.CP_1 = .colSums(lgamma(n.G.by.CP_1 + beta), nrow(n.G.by.CP), ncol(n.G.by.CP))

	n.CP_1 = n.CP
	n.CP_1[z[i]] = n.CP_1[z[i]] - n.by.C[i]
	n.CP_1 = lgamma(n.CP_1 + (nG * beta))
	
	## Add cell counts to all other populations
	n.G.by.CP_2 = n.G.by.CP
	other.ix = (1:K)[-z[i]]
	n.G.by.CP_2[,other.ix] = n.G.by.CP_2[,other.ix] + counts[,i]
	n.G.by.CP_2 = .colSums(lgamma(n.G.by.CP_2 + beta), nrow(n.G.by.CP), ncol(n.G.by.CP))

	n.CP_2 = n.CP
	n.CP_2[other.ix] = n.CP_2[other.ix] + n.by.C[i]
	n.CP_2 = lgamma(n.CP_2 + (nG * beta))


	m.CP.by.S[z[i],s[i]] = m.CP.by.S[z[i],s[i]] - 1L	

	## Calculate probabilities for each state
	for(j in 1:K) {
      other.ix = (1:K)[-j]
	  probs[j,i] = 	log(m.CP.by.S[j,s[i]] + alpha) +		## Theta simplified
				  sum(n.G.by.CP_1[other.ix]) +				## Phi Numerator (other cells)
				  n.G.by.CP_2[j] -							## Phi Numerator (current cell)
				  sum(n.CP_1[other.ix]) -					## Phi Denominator (other cells)
  				  n.CP_2[j]									## Phi Denominator (current cell)
	}  

	## Sample next state and add back counts
	prev.z = z[i]
	if(isTRUE(do.sample)) z[i] = sample.ll(probs[,i])
	
	if(prev.z != z[i]) { 
	    n.G.by.CP[,prev.z] = n.G.by.CP[,prev.z] - counts[,i]
	  	n.G.by.CP[,z[i]] = n.G.by.CP[,z[i]] + counts[,i]
	  	
	  	n.CP[prev.z] = n.CP[prev.z] - n.by.C[i]
	  	n.CP[z[i]] = n.CP[z[i]] + n.by.C[i]
	}
	m.CP.by.S[z[i],s[i]] = m.CP.by.S[z[i],s[i]] + 1L
  }  

  return(list(m.CP.by.S=m.CP.by.S, n.G.by.CP=n.G.by.CP, n.CP=n.CP, z=z, probs=probs))
}

cC.calcEMProbZ = function(counts, m.CP.by.S, n.G.by.CP, n.by.C, n.CP, z, s, K, nG, nM, alpha, beta, do.sample=TRUE) {

  ## Expectation given current cell population labels
  theta = fastNormPropLog(m.CP.by.S, alpha)
  phi = fastNormPropLog(n.G.by.CP, beta)
  #theta = log(normalizeCounts(m.CP.by.S + alpha, normalize="proportion"))
  #phi = log(normalizeCounts(n.G.by.CP + beta, normalize="proportion"))
  
  ## Maximization to find best label for each cell
  probs = eigenMatMultInt(phi, counts) + theta[, s]  
  #probs = (t(phi) %*% counts) + theta[, s]  
  
  z.previous = z
  z = apply(probs, 2, which.max)

  ## Recalculate counts based on new label
  #p = cC.decomposeCounts(counts, s, z, K)
  p = cC.reDecomposeCounts(counts, s, z, z.previous, n.G.by.CP, K)
  m.CP.by.S = p$m.CP.by.S
  n.G.by.CP = p$n.G.by.CP
  n.CP = p$n.CP

  return(list(m.CP.by.S=m.CP.by.S, n.G.by.CP=n.G.by.CP, n.CP=n.CP, z=z, probs=probs))
}

#' @title Simulate cells from the celda_C model
#' 
#' @description Generates a simulated counts matrix, cell subpopulation clusters, and sample labels
#' according to the generative process of the celda_C model. 
#' 
#' @param model Character. Options available in `celda::available.models`. 
#' @param S Integer. Number of samples to simulate. Default 5.
#' @param C.Range Vector of length 2 given the range (min,max) of number of cells for each sample to be randomly generated from the uniform distribution. Default c(50, 100).
#' @param N.Range Integer vector. A vector of length 2 that specifies the lower and upper bounds of the number of counts generated for each cell. Default c(500, 1000). 
#' @param G Integer. The total number of features to be simulated. Default 100.
#' @param K Integer. Number of cell populations. Default 5.
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount to each cell population in each sample. Default 1. 
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to each feature in each cell population. Default 1. 
#' @param seed Integer. Passed to `set.seed()`. Default 12345.  
#' @param ... Additional parameters.
#' @return List. Contains the simulated matrix `counts`, cell population clusters `z`, sample assignments `sample.label`, and input parameters.
#' @seealso `celda_G()` for simulating feature modules and `celda_CG()` for simulating feature modules and cell populations. 
#' @examples
#' celda.c.sim = simulateCells(model="celda_C", K=10)
#' sim.counts = celda.c.sim$counts
#' @export
simulateCells.celda_C = function(model, S=5, C.Range=c(50, 100), N.Range=c(500,1000), 
                         G=100, K=5, alpha=1, beta=1, seed=12345, ...) {
 
  set.seed(seed) 
    
  phi <- rdirichlet(K, rep(beta, G))
  theta <- rdirichlet(S, rep(alpha, K))
  
  ## Select the number of cells per sample
  nC <- sample(C.Range[1]:C.Range[2], size=S, replace=TRUE)  
  cell.sample.label <- rep(1:S, nC)
  
  ## Select state of the cells  
  z <- unlist(lapply(1:S, function(i) sample(1:K, size=nC[i], prob=theta[i,], replace=TRUE)))
    
  ## Select number of transcripts per cell
  nN <- sample(N.Range[1]:N.Range[2], size=length(cell.sample.label), replace=TRUE)
  
  ## Select transcript distribution for each cell
  cell.counts <- sapply(1:length(cell.sample.label), function(i) stats::rmultinom(1, size=nN[i], prob=phi[z[i],]))
  
  rownames(cell.counts) = paste0("Gene_", 1:nrow(cell.counts))
  colnames(cell.counts) = paste0("Cell_", 1:ncol(cell.counts)) 
  cell.sample.label = paste0("Sample_", 1:S)[cell.sample.label]
  cell.sample.label = factor(cell.sample.label, levels=paste0("Sample_", 1:S))
  
  ## Peform reordering on final Z and Y assigments:
  cell.counts = processCounts(cell.counts) 
  names = list(row=rownames(cell.counts), column=colnames(cell.counts), 
               sample=unique(cell.sample.label))
  result = methods::new("celda_C", clusters=list(z=z),
                        params=list(alpha=alpha, beta=beta, seed=seed,
                                    count.checksum=digest::digest(cell.counts, algo="md5"),
                                    K=K),
                        sample.label=cell.sample.label, 
                        names=names)
  class(result) = "celda_C" 
  result = reorder.celda_C(counts = cell.counts, res = result)
  
  return(list(z=result@clusters$z, counts=processCounts(cell.counts), 
              sample.label=cell.sample.label, K=K, alpha=alpha, 
              beta=beta, C.Range=C.Range, N.Range=N.Range, S=S))
}


#' @title Matrix factorization for results from celda_C()
#' @description Generates factorized matrices showing the contribution of each feature in each cell population or each cell population in each sample. 
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class "celda_C".
#' @param type Character vector. A vector containing one or more of "counts", "proportion", or "posterior". "counts" returns the raw number of counts for each factorized matrix. "proportions" returns the normalized probabilities for each factorized matrix, which are calculated by dividing the raw counts in each factorized matrix by the total counts in each column. "posterior" returns the posterior estimates. Default `c("counts", "proportion", "posterior")`. 
#' @examples 
#' factorized.matrices = factorizeMatrix(celda.C.sim$counts, celda.C.mod, "posterior")
#' @return A list with elements for `counts`, `proportions`, or `posterior` probabilities. Each element will be a list containing factorized matrices for `module` and `sample`.
#' @seealso `celda_C()` for clustering cells
#' @export
setMethod("factorizeMatrix", 
          signature(celda.mod = "celda_C"),
          function(counts, celda.mod, 
                   type=c("counts", "proportion", "posterior")) {
            counts = processCounts(counts) 
            compareCountMatrix(counts, celda.mod)
            
            K = celda.mod@params$K
            z = celda.mod@clusters$z
            alpha = celda.mod@params$alpha
            beta = celda.mod@params$beta
            sample.label = celda.mod@sample.label
            s = as.integer(sample.label)
                  
            p = cC.decomposeCounts(counts, s, z, K)
            m.CP.by.S = p$m.CP.by.S
            n.G.by.CP = p$n.G.by.CP
              
            K.names = paste0("K", 1:K)
            rownames(n.G.by.CP) = celda.mod@names$row
            colnames(n.G.by.CP) = K.names
            rownames(m.CP.by.S) = K.names
            colnames(m.CP.by.S) = celda.mod@names$sample
          
            counts.list = c()
            prop.list = c()
            post.list = c()
            res = list()
                          
            if(any("counts" %in% type)) {
              counts.list = list(sample=m.CP.by.S, module=n.G.by.CP)
              res = c(res, list(counts=counts.list))
            }
            if(any("proportion" %in% type)) {
              ## Need to avoid normalizing cell/gene states with zero cells/genes
              unique.z = sort(unique(z))
              temp.n.G.by.CP = n.G.by.CP
              temp.n.G.by.CP[,unique.z] = normalizeCounts(temp.n.G.by.CP[,unique.z], 
                                                          normalize="proportion")
          
              prop.list = list(sample = normalizeCounts(m.CP.by.S, 
                                                        normalize="proportion"),
                               module = temp.n.G.by.CP)
              res = c(res, list(proportions=prop.list))
            }
            if(any("posterior" %in% type)) {
              post.list = list(sample = normalizeCounts(m.CP.by.S + alpha, 
                                                        normalize="proportion"),
                               module = normalizeCounts(n.G.by.CP + beta, 
                                                        normalize="proportion"))
              res = c(res, posterior = list(post.list))                           
            }
          
            return(res)
})


# Calculate log-likelihood for celda_C model
cC.calcLL = function(m.CP.by.S, n.G.by.CP, s, z, K, nS, nG, alpha, beta) {
  
  ## Calculate for "Theta" component
  a = nS * lgamma(K * alpha)
  b = sum(lgamma(m.CP.by.S + alpha))
  c = -nS * K * lgamma(alpha)
  d = -sum(lgamma(colSums(m.CP.by.S + alpha)))
  
  theta.ll = a + b + c + d
  
  ## Calculate for "Phi" component
  a = K * lgamma(nG * beta)
  b = sum(lgamma(n.G.by.CP + beta))
  c = -K * nG * lgamma(beta)
  d = -sum(lgamma(colSums(n.G.by.CP + beta)))
  
  phi.ll = a + b + c + d
  
  final = theta.ll + phi.ll
  return(final)
}


#' @title Calculate Celda_C log likelihood
#' @description Calculates the log likelihood for user-provided cell population clusters using the `celda_C()` model.
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. 
#' @param sample.label Vector or factor. Denotes the sample label for each cell (column) in the count matrix.
#' @param z Numeric vector. Denotes cell population labels. 
#' @param K Integer. Number of cell populations. 
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount to each cell population in each sample. Default 1. 
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to each feature in each cell population. Default 1. 
#' @param ... Additional parameters.
#' @return Numeric. The log likelihood for the given cluster assignments
#' @seealso `celda_C()` for clustering cells
#' @examples
#' loglik = logLikelihood(celda.C.sim$counts, model="celda_C", 
#'                        sample.label=celda.C.sim$sample.label,
#'                        z=celda.C.sim$z, K=celda.C.sim$K,
#'                        alpha=celda.C.sim$alpha, beta=celda.C.sim$beta)
#' @export
logLikelihood.celda_C = function(counts, model, sample.label, z, K, 
                                 alpha, beta) {
  if (sum(z > K) > 0) stop("An entry in z contains a value greater than the provided K.")
  sample.label = processSampleLabels(sample.label, ncol(counts))
  s = as.integer(sample.label)
  p = cC.decomposeCounts(counts, s, z, K)  
  final = cC.calcLL(m.CP.by.S=p$m.CP.by.S, n.G.by.CP=p$n.G.by.CP, 
                    s=s, z=z, K=K, nS=p$nS, nG=p$nG, alpha=alpha,
                    beta=beta)
  return(final)
}


# Takes raw counts matrix and converts it to a series of matrices needed for log likelihood calculation
# @param counts Integer matrix. Rows represent features and columns represent cells. 
# @param s Integer vector. Contains the sample label for each cell (column) in the count matrix. 
# @param z Numeric vector. Denotes cell population labels. 
# @param K Integer. Number of cell populations. 
cC.decomposeCounts = function(counts, s, z, K) {
  nS = length(unique(s))
  nG = nrow(counts)
  nM = ncol(counts)

  m.CP.by.S = matrix(as.integer(table(factor(z, levels=1:K), s)), ncol=nS)
  n.G.by.CP = colSumByGroup(counts, group=z, K=K)
  n.CP = as.integer(colSums(n.G.by.CP))
  n.by.C = as.integer(colSums(counts))

  return(list(m.CP.by.S=m.CP.by.S, n.G.by.CP=n.G.by.CP, n.CP=n.CP, n.by.C=n.by.C, nS=nS, nG=nG, nM=nM))
}

cC.reDecomposeCounts = function(counts, s, z, previous.z, n.G.by.CP, K) {

  ## Recalculate counts based on new label
  n.G.by.CP = colSumByGroupChange(counts, n.G.by.CP, z, previous.z, K)
  nS = length(unique(s))
  m.CP.by.S = matrix(as.integer(table(factor(z, levels=1:K), s)), ncol=nS)
  n.CP = as.integer(colSums(n.G.by.CP))

  return(list(m.CP.by.S=m.CP.by.S, n.G.by.CP=n.G.by.CP, n.CP=n.CP))  
}


#' @title Conditional probabilities for cells in subpopulations from a Celda_C model
#' @description Calculates the conditional probability of each cell belonging to each subpopulation given all other cell cluster assignments in a `celda_C()` result. 
#'
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class `celda_C`.
#' @param log Logical. If FALSE, then the normalized conditional probabilities will be returned. If TRUE, then the unnormalized log probabilities will be returned. Default FALSE.  
#' @param ... Additional parameters.
#' @return A list containging a matrix for the conditional cell subpopulation cluster probabilities. 
#' @seealso `celda_C()` for clustering cells
#' @examples
#' cluster.prob = clusterProbability(celda.C.sim$counts, celda.C.mod)
#' @export
setMethod("clusterProbability", 
          signature(celda.mod = "celda_C"),
          function(counts, celda.mod, log=FALSE, ...) {
            z = celda.mod@clusters$z
            sample.label = celda.mod@sample.label
            s = as.integer(sample.label)
            
            K = celda.mod@params$K
            alpha = celda.mod@params$alpha
            beta = celda.mod@params$beta
            
            p = cC.decomposeCounts(counts, s, z, K)  
            
            next.z = cC.calcGibbsProbZ(counts=counts, m.CP.by.S=p$m.CP.by.S, 
                                       n.G.by.CP=p$n.G.by.CP, n.by.C=p$n.by.C,
                                       n.CP=p$n.CP, z=z, s=s, K=K, nG=p$nG, 
                                       nM=p$nM, alpha=alpha, beta=beta, 
                                       do.sample=FALSE)  
            z.prob = t(next.z$probs)
            
            if(!isTRUE(log)) {
              z.prob = normalizeLogProbs(z.prob)
            }
             
            return(list(z.probability=z.prob))
         })


#' @title Calculate the perplexity on new data with a celda_C model
#' @description Perplexity is a statistical measure of how well a probability model can predict new data. Lower perplexity indicates a better model. 
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class "celda_C"
#' @param new.counts A new counts matrix used to calculate perplexity. If NULL, perplexity will be calculated for the 'counts' matrix. Default NULL.
#' @return Numeric. The perplexity for the provided count data and model.
#' @seealso `celda_C()` for clustering cells
#' @examples
#' perplexity = perplexity(celda.C.sim$counts, celda.C.mod)
#' @export
setMethod("perplexity",
          signature(celda.mod = "celda_C"),
          function(counts, celda.mod, new.counts=NULL) {
            compareCountMatrix(counts, celda.mod)
            if (!("celda_C" %in% class(celda.mod))) { 
              stop("The celda.mod provided was not of class celda_C.")
            }
            
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
                                         type="posterior")
            theta = log(factorized$posterior$sample)
            phi = log(factorized$posterior$module)
            s = as.integer(celda.mod@sample.label)
            
            inner.log.prob = (t(phi) %*% new.counts) + theta[, s]  
            log.px = sum(apply(inner.log.prob, 2, matrixStats::logSumExp))
            
            perplexity = exp(-(log.px/sum(new.counts)))
            return(perplexity)
         })  


reorder.celda_C = function(counts, res){
  if(res@params$K > 2 & isTRUE(length(unique(res@clusters$z)) > 1)) {
    res@clusters$z = as.integer(as.factor(res@clusters$z))
    fm <- factorizeMatrix(counts = counts, celda.mod = res)
    unique.z = sort(unique(res@clusters$z))
    d <- cosineDist(fm$posterior$module[,unique.z])
    h <- stats::hclust(d, method = "complete")
    res <- recodeClusterZ(res, from = h$order, to = 1:length(h$order))
  }  
  return(res)
}


#' @title Heatmap for celda_C
#' @description Renders an expression heatmap to visualize `celda_C()` results. Features to include in the heatmap must be supplied. 
#'  
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class `celda_C`.
#' @param feature.ix Integer vector. Indices of features to plot, such the top features from a differential expression analysis. 
#' @param ... Additional parameters.
#' @seealso `celda_C()` for clustering cells and `celdaTsne()` for generating 2-dimensional coordinates
#' @examples 
#' celdaHeatmap(celda.C.sim$counts, celda.C.mod)
#' @return list A list containing dendrograms and the heatmap grob
#' @export
setMethod("celdaHeatmap", 
          signature(celda.mod = "celda_C"),
          function(counts, celda.mod, feature.ix, ...) {
            norm = normalizeCounts(counts, normalize="proportion", transformation.fun=sqrt)
            plotHeatmap(norm[feature.ix,], z=celda.mod@clusters$z, ...)
          })


#' @title tSNE for celda_C
#' @description Embeds cells in two dimensions using tSNE based on a `celda_C` model. PCA on the normalized counts is used to reduce the number of features before applying tSNE. 
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class `celda_C`. 
#' @param max.cells Integer. Maximum number of cells to plot. Cells will be randomly subsampled if ncol(counts) > max.cells. Larger numbers of cells requires more memory. Default 25000.
#' @param min.cluster.size Integer. Do not subsample cell clusters below this threshold. Default 100. 
#' @param initial.dims Integer. PCA will be used to reduce the dimentionality of the dataset. The top 'initial.dims' principal components will be used for tSNE. Default 20.
#' @param perplexity Numeric. Perplexity parameter for tSNE. Default 20.
#' @param max.iter Integer. Maximum number of iterations in tSNE generation. Default 2500.
#' @param seed Integer. Passed to `set.seed()`. Default 12345.  
#' @param ... Additional parameters.
#' @seealso `celda_C()` for clustering cells and `celdaHeatmap()` for displaying expression
#' @examples
#' tsne.res = celdaTsne(celda.C.sim$counts, celda.C.mod)
#' @return A two column matrix of t-SNE coordinates
#' @export
setMethod("celdaTsne",
          signature(celda.mod = "celda_C"),
            function(counts, celda.mod, max.cells=25000, min.cluster.size=100,
                     initial.dims=20, modules=NULL, perplexity=20, max.iter=2500, 
                     seed=12345, ...) {

            counts = processCounts(counts)
            compareCountMatrix(counts, celda.mod)
            
            ## Checking if max.cells and min.cluster.size will work
            if((max.cells < ncol(counts)) & (max.cells / min.cluster.size < celda.mod@params$K)) {
              stop(paste0("Cannot distribute ", max.cells, " cells among ", 
                          celda.mod@params$K, " clusters while maintaining a minumum of ", 
                          min.cluster.size, 
                          " cells per cluster. Try increasing 'max.cells' or decreasing 'min.cluster.size'."))
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
          	## Rounding can cause number to be off by a few, so ceiling is 
          	## used with a second round of subtraction
          	cluster.n.to.sample = ceiling((cluster.cells.to.sample / sum(cluster.cells.to.sample)) * total.cells.to.remove)
          	diff = sum(cluster.n.to.sample) - total.cells.to.remove 
          	cluster.n.to.sample[which.max(cluster.n.to.sample)] = cluster.n.to.sample[which.max(cluster.n.to.sample)] - diff
          
          	## Perform sampling for each cluster
          	for(i in which(cluster.n.to.sample > 0)) {
          	  z.include[sample(which(celda.mod@clusters$z == i), cluster.n.to.sample[i])] = FALSE
          	}
            }   
            cell.ix = which(z.include)
          
            norm = t(normalizeCounts(counts[,cell.ix], normalize="proportion", 
                                     transformation.fun=sqrt))
            res = calculateTsne(norm, perplexity=perplexity, max.iter=max.iter, 
                                seed=seed, do.pca=TRUE, 
                                initial.dims = initial.dims)
            final = matrix(NA, nrow=ncol(counts), ncol=2)
            final[cell.ix,] = res
            rownames(final) = colnames(counts)
            colnames(final) = c("tsne_1", "tsne_2")
            return(final)
          })



#' @title Probability map for a celda_C model
#' @description Renders probability and relative expression heatmaps to visualize the relationship between cell populations and samples.
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`. 
#' @param celda.mod Celda object of class `celda_C`.   
#' @param level Character. 'sample' will display the absolute probabilities and relative normalized abundance of each cell population in each sample. Default 'sample'.
#' @param ... Additional parameters.
#' @seealso `celda_C()` for clustering cells
#' @examples
#' celdaProbabilityMap(celda.C.sim$counts, celda.C.mod)
#' @return A grob containing the specified plots
#' @export 
setMethod("celdaProbabilityMap",
          signature(celda.mod = "celda_C"),
          function(counts, celda.mod, level=c("sample"), ...){
            counts = processCounts(counts)
            compareCountMatrix(counts, celda.mod)
            
            z.include = which(tabulate(celda.mod@clusters$z, celda.mod@params$K) > 0)
            
            level = match.arg(level)
            factorized <- factorizeMatrix(celda.mod = celda.mod, 
                                          counts = counts)
            
            samp <- factorized$proportions$sample[z.include,,drop=FALSE]
            col <- colorRampPalette(c("white","blue","#08306B","#006D2C",
                                      "yellowgreen","yellow","orange",
                                      "red"))(100)
            breaks <-  seq(0, 1, length.out = length(col))     
            g1 = plotHeatmap(samp, color.scheme="sequential", scale.row=NULL, 
                             cluster.cell=FALSE, cluster.feature=FALSE, 
                             show.names.cell=TRUE, show.names.feature=TRUE, 
                             breaks = breaks, col=col, 
                             main = "Absolute Probability", silent=TRUE)
          
            if(ncol(samp) > 1) {
            	samp.norm = normalizeCounts(samp, normalize="proportion", 
            	                            transformation.fun=sqrt, 
            	                            scale.fun=base::scale)
            	g2 = plotHeatmap(samp.norm, color.scheme="divergent", 
            	                 cluster.cell=FALSE, cluster.feature=FALSE, 
            	                 show.names.cell=TRUE, show.names.feature=TRUE, 
            	                 main = "Relative Abundance", silent=TRUE)   
            	gridExtra::grid.arrange(g1$gtable, g2$gtable, ncol=2)
            } else {
          	  gridExtra::grid.arrange(g1$gtable)
            } 
        })


#' @title Lookup the module of a feature
#' @description Finds the module assignments of given features in a `celda_C()` model
#'  
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Model of class `celda_C`.
#' @param feature Character vector. The module assignemnts will be found for feature names in this vector. 
#' @param exact.match Logical. Whether an exact match or a partial match using `grep()` is required to look up the feature in the rownames of the counts matrix. Default TRUE. 
#' @return List. Each element contains the module of the provided feature.
#' @export
setMethod("featureModuleLookup",
          signature(celda.mod = "celda_C"),
          function(counts, celda.mod, feature, exact.match){
            stop("Celda_C models do not contain feature modules.")
          })
