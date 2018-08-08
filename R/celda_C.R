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


#' celda Cell Clustering Model
#' 
#' @param counts A numeric count matrix
#' @param sample.label A vector indicating the sample for each cell (column) in the count matrix
#' @param K An integer or range of integers indicating the desired number of cell clusters (for celda_C / celda_CG models)
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution
#' @param beta Non-zero concentration parameter for gene Dirichlet distribution
#' @param algorithm Use 'EM' or 'Gibbs' sampling for clustering of cells into subpopulations. EM is much faster for larger datasets. Default 'EM'.
#' @param stop.iter Number of iterations without improvement in the log likelihood to stop the inference algorithm. Default 10.
#' @param max.iter Maximum iterations of inference algorithm to perform regardless of convergence. Default 200.
#' @param split.on.iter On every 'split.on.iter' iteration, a heuristic will be applied to determine if a gene/cell cluster should be reassigned and another gene/cell cluster should be split into two clusters. Default 10.
#' @param split.on.last After the the chain has converged according to 'stop.iter', a heuristic will be applied to determine if a gene/cell cluster should be reassigned and another gene/cell cluster should be split into two clusters. If a split occurs, then 'stop.iter' will be reset. Default TRUE.
#' @param seed Parameter to set.seed() for random number generation. Default 12345.
#' @param nchains Number of random z initializations. Default 3. 
#' @param count.checksum An MD5 checksum for the provided counts matrix. Default NULL.
#' @param z.init Initial values of z. If NULL, z will be randomly sampled. Default NULL.
#' @param logfile If NULL, messages will be displayed as normal. If set to a file name, messages will be redirected messages to the file. Default NULL.
#' @return An object of class celda_C with clustering results and Gibbs sampling statistics
#' @export
celda_C = function(counts, sample.label=NULL, K, alpha=1, beta=1,
					 algorithm = c("EM", "Gibbs"), 
                 	 stop.iter = 10, max.iter=200, split.on.iter=10, split.on.last=TRUE,
                 	 seed=12345, nchains=3, count.checksum=NULL, 
                 	 z.init = NULL, logfile=NULL) {
  
  ## Error checking and variable processing
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

  all.seeds = seed:(seed + nchains - 1)
  
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=FALSE)  
  logMessages("Celda_C: Clustering cells.", logfile=logfile, append=FALSE)
  logMessages("--------------------------------------------------------------------", logfile=logfile, append=FALSE)  

  best.result = NULL  
  for(i in seq_along(all.seeds)) { 
  
	## Randomly select z or set z to supplied initial values
	current.seed = all.seeds[i]	
	z = initialize.cluster(K, ncol(counts), initial = z.init, fixed = NULL, seed=current.seed)
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

    logMessages(date(), ".. Starting chain", i, "with seed", current.seed, logfile=logfile, append=FALSE)

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
	  if(K > 2 & (((iter == max.iter | num.iter.without.improvement == stop.iter) & isTRUE(split.on.last)) | (split.on.iter > 0 & iter %% split.on.iter == 0 & isTRUE(do.cell.split)))) {

		logMessages(date(), " .... Determining if any cell clusters should be split.", logfile=logfile, append=TRUE, sep="")
		res = cC.splitZ(counts, m.CP.by.S, n.G.by.CP, s, z, K, nS, nG, alpha, beta, z.prob=t(next.z$probs), max.clusters.to.try=10, min.cell=3)
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
	
	  logMessages(date(), ".... Completed iteration:", iter, "| logLik:", temp.ll, logfile=logfile, append=TRUE)
	  iter = iter + 1    
	}
	
	names = list(row=rownames(counts), column=colnames(counts), sample=levels(sample.label))

	result = list(z=z.best, completeLogLik=ll,  
				  finalLogLik=ll.best, seed=current.seed, K=K, 
				  sample.label=sample.label, alpha=alpha, 
				  beta=beta, count.checksum=count.checksum, 
				  names=names)
  
	class(result) = "celda_C"
	
    if(is.null(best.result) || result$finalLogLik > best.result$finalLogLik) {
      best.result = result
    }
    
    logMessages(date(), ".. Finished chain", i, "with seed", current.seed, logfile=logfile, append=FALSE)
  }  
  
  best.result = reorder.celda_C(counts = counts, res = best.result)
  return(best.result)
}


# Gibbs sampling for the celda_C Model
cC.calcGibbsProbZ = function(counts, m.CP.by.S, n.G.by.CP, n.by.C, n.CP, z, s, K, nG, nM, alpha, beta, do.sample=TRUE) {

  ## Set variables up front outside of loop  
  probs = matrix(NA, ncol=nM, nrow=K)
#  temp.n.G.by.CP = n.G.by.CP
#  temp.n.CP = n.CP
#  n.G.by.CP_1 = n.G.by.CP
#  n.G.by.CP_2 = n.G.by.CP

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
  theta = log(normalizeCounts(m.CP.by.S + alpha, scale.factor=1))
  phi = log(normalizeCounts(n.G.by.CP + beta, scale.factor=1))
  
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

#' Simulate cells from the cell clustering generative model
#' 
#' @param model Celda model to use for simulation. One of 'available_models'. 
#' @param S Total number of samples
#' @param C.Range Vector of length 2 given the range (min,max) of number of cells for each sample to be randomly generated from the uniform distribution
#' @param N.Range Vector of length 2 given the range (min,max) of number of counts for each cell to be randomly generated from the uniform distribution
#' @param G Total number of Genes to be simulated
#' @param K An integer or range of integers indicating the desired number of cell clusters (for celda_C / celda_CG models)
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution
#' @param beta Non-zero concentration parameter for gene Dirichlet distribution
#' @param seed starting point used for generating simulated data
#' @param ... Other arguments
#' @export
simulateCells.celda_C = function(model, S=10, C.Range=c(10, 100), N.Range=c(100,5000), 
                         G=500, K=5, alpha=1, beta=1, seed=12345, ...) {
 
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

  ## Peform reordering on final Z and Y assigments:
  names = list(row=rownames(cell.counts), column=colnames(cell.counts), 
               sample=unique(cell.sample.label))
  result = list(z=z, completeLogLik=NULL, 
                finalLogLik=NULL, K=K, 
                alpha=alpha, beta=beta, seed=seed, 
                sample.label=cell.sample.label, names=names,
                count.checksum=NULL)
  class(result) = "celda_C" 
  result = reorder.celda_C(counts = cell.counts, res = result)
  
  return(list(z=result$z, counts=processCounts(cell.counts), sample.label=cell.sample.label, K=K, alpha=alpha, beta=beta, C.Range=C.Range, N.Range=N.Range, S=S))
}


#' Generate factorized matrices showing each feature's influence on the celda_C model clustering 
#' 
#' @param counts A numeric count matrix
#' @param celda.mod Object return from celda_C function
#' @param type A character vector containing one or more of "counts", "proportions", or "posterior". "counts" returns the raw number of counts for each entry in each matrix. "proportions" returns the counts matrix where each vector is normalized to a probability distribution. "posterior" returns the posterior estimates which include the addition of the Dirichlet concentration parameter (essentially as a pseudocount).
#' @param validate.counts Whether to verify that the counts matrix provided was used to generate the results in celda.mod. Defaults to TRUE.
#' @export
factorizeMatrix.celda_C = function(counts, celda.mod, type=c("counts", "proportion", "posterior")) {

  counts = processCounts(counts)
  K = celda.mod$K
  z = celda.mod$z
  alpha = celda.mod$alpha
  beta = celda.mod$beta
  sample.label = celda.mod$sample.label
  s = processSampleLabels(sample.label, ncol(counts))
        
  p = cC.decomposeCounts(counts, s, z, K)
  m.CP.by.S = p$m.CP.by.S
  n.G.by.CP = p$n.G.by.CP
    
  K.names = paste0("K", 1:K)
  rownames(n.G.by.CP) = celda.mod$names$row
  colnames(n.G.by.CP) = K.names
  rownames(m.CP.by.S) = K.names
  colnames(m.CP.by.S) = celda.mod$names$sample

  counts.list = c()
  prop.list = c()
  post.list = c()
  res = list()
                
  if(any("counts" %in% type)) {
    counts.list = list(sample.states=m.CP.by.S, gene.states=n.G.by.CP)
    res = c(res, list(counts=counts.list))
  }
  if(any("proportion" %in% type)) {
    ## Need to avoid normalizing cell/gene states with zero cells/genes
    unique.z = sort(unique(z))
    temp.n.G.by.CP = n.G.by.CP
    temp.n.G.by.CP[,unique.z] = normalizeCounts(temp.n.G.by.CP[,unique.z], scale.factor=1)

    prop.list = list(sample.states = normalizeCounts(m.CP.by.S, scale.factor=1),
                     gene.states = temp.n.G.by.CP)
    res = c(res, list(proportions=prop.list))
  }
  if(any("posterior" %in% type)) {
    post.list = list(sample.states = normalizeCounts(m.CP.by.S + alpha, scale.factor=1),
                     gene.states = normalizeCounts(n.G.by.CP + beta, scale.factor=1))
    res = c(res, posterior = list(post.list))                           
  }

  return(res)
}


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


#' Calculate the celda_C log likelihood for user-provided cluster assignments
#' 
#' @param counts A numeric count matrix
#' @param sample.label A vector indicating the sample label for each cell (column) in the count matrix
#' @param z A numeric vector of cluster assignments
#' @param K The total number of clusters in z
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution
#' @param beta Non-zero concentration parameter for gene Dirichlet distribution
#' @param ... Additional parameters
#' @export
calculateLoglikFromVariables.celda_C = function(counts, sample.label, z, K, alpha, beta) {
  s = processSampleLabels(sample.label, ncol(counts))
  p = cC.decomposeCounts(counts, s, z, K)  
  final = cC.calcLL(m.CP.by.S=p$m.CP.by.S, n.G.by.CP=p$n.G.by.CP, s=s, z=z, K=K, nS=p$nS, nG=p$nG, alpha=alpha, beta=beta)
  return(final)
}


#' Takes raw counts matrix and converts it to a series of matrices needed for log likelihood calculation
#' @param counts A numeric count matrix
#' @param s An integer vector indicating the sample label for each cell (column) in the count matrix
#' @param z A numeric vector of cluster assignments
#' @param K The total number of clusters in z
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


#' Calculates the conditional probability of each cell belong to each cluster given all other cluster assignments
#'
#' @param celda.mod A model returned from the 'celda_C' function
#' @param counts The original count matrix used in the model
#' @param log If FALSE, then the normalized conditional probabilities will be returned. If TRUE, then the unnormalized log probabilities will be returned.  
#' @param ... Other arguments
#' @return A list containging a matrix for the conditional cell cluster probabilities. 
#' @export
clusterProbability.celda_C = function(celda.mod, counts, log=FALSE, ...) {

  z = celda.mod$z
  sample.label = celda.mod$sample.label
  s = processSampleLabels(sample.label, ncol(counts))
  
  K = celda.mod$K
  alpha = celda.mod$alpha
  beta = celda.mod$beta
  
  p = cC.decomposeCounts(counts, s, z, K)  
  
  next.z = cC.calcGibbsProbZ(counts=counts, m.CP.by.S=p$m.CP.by.S, n.G.by.CP=p$n.G.by.CP, n.by.C=p$n.by.C, n.CP=p$n.CP, z=z, s=s, K=K, nG=p$nG, nM=p$nM, alpha=alpha, beta=beta, do.sample=FALSE)  
  z.prob = t(next.z$probs)
  
  if(!isTRUE(log)) {
    z.prob = normalizeLogProbs(z.prob)
  }
   
  return(list(z.probability=z.prob))
}


#' @export
calculatePerplexity.celda_C = function(counts, celda.mod, validate.counts) {
  
  factorized = factorizeMatrix(counts = counts, celda.mod = celda.mod, 
                               "posterior", validate.counts)
  theta = log(factorized$posterior$sample.states)
  phi = log(factorized$posterior$gene.states)
  sl = celda.mod$sample.label
  
  inner.log.prob = (t(phi) %*% counts) + theta[, sl]  
  log.px = sum(apply(inner.log.prob, 2, matrixStats::logSumExp))
  
  perplexity = exp(-(log.px/sum(counts)))
  return(perplexity)
}  


reorder.celda_C = function(counts, res){
  if(res$K > 2 & isTRUE(length(unique(res$z)) > 1)) {
    res$z = as.integer(as.factor(res$z))
    fm <- factorizeMatrix(counts = counts, celda.mod = res,
                          validate.counts = FALSE)
    unique.z = sort(unique(res$z))
    d <- cosineDist(fm$posterior$gene.states[,unique.z])
    h <- hclust(d, method = "complete")
    res <- recodeClusterZ(res, from = h$order, to = 1:length(h$order))
  }  
  return(res)
}


#' finalClusterAssignment for celda Cell clustering funciton 
#' @param celda.mod A celda model object of class "celda_C"
#' @export
finalClusterAssignment.celda_C = function(celda.mod) {
  return(celda.mod$z)
}


#' getK for celda Cell clustering function 
#' @param celda.mod A celda model object of class "celda_C"
#' @export
getK.celda_C = function(celda.mod) {
  return(celda.mod$K)
}


#' getL for celda Cell clustering function 
#' @param celda.mod A celda model object of class "celda_C"
#' @export
getL.celda_C = function(celda.mod) { return(NA) }


#' celdaHeatmap for celda Cell clustering function 
#' @param celda.mod A celda model object of class "celda_C"
#' @param counts A numeric count matrix
#' @param ... extra parameters passed onto the renderCeldaHeatmap
#' @export
celdaHeatmap.celda_C = function(celda.mod, counts, ...) {
  renderCeldaHeatmap(counts, z=celda.mod$z, ...)
}


#' Embeds cells in two dimensions using tSNE based on celda_C results.
#' 
#' @param counts Counts matrix, should have cell name for column name and gene name for row name.
#' @param celda.mod Celda model to use for tsne. 
#' @param max.cells Integer; Maximum number of cells to plot. Cells will be randomly subsampled if ncol(conts) > max.cells. Larger numbers of cells requires more memory. Default 10000.
#' @param min.cluster.size Integer; Do not subsample cell clusters below this threshold. Default 100. 
#' @param initial.dims PCA will be used to reduce the dimentionality of the dataset. The top 'initial.dims' principal components will be used for tSNE.
#' @param perplexity Numeric vector; determines perplexity for tSNE. Default 20.
#' @param max.iter Numeric vector; determines iterations for tsne. Default 1000.
#' @param seed Seed for random number generation. Default 12345.
#' @param ... Further arguments passed to or from other methods.
#' @export
celdaTsne.celda_C = function(counts, celda.mod,  
							 max.cells=10000, min.cluster.size=100, initial.dims=20,
							 perplexity=20, max.iter=2500, seed=12345, ...) {

  norm = sqrt(normalizeCounts(counts, scale.factor=1))

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

  res = calculateTsne(norm[,cell.ix], perplexity=perplexity, max.iter=max.iter, distance="euclidean", seed=seed, do.pca=TRUE, initial.dims = initial.dims)
  final = matrix(NA, nrow=ncol(norm), ncol=2)
  final[cell.ix,] = res
  rownames(final) = colnames(norm)
  colnames(final) = c("tsne_1", "tsne_2")
  return(final)
}
