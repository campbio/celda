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
#' @param stop.iter Number of iterations without improvement in the log likelihood to stop the Gibbs sampler. Default 10.
#' @param max.iter Maximum iterations of Gibbs sampling to perform regardless of convergence. Default 200.
#' @param split.on.iter On every 'split.on.iter' iteration, a heuristic will be applied to determine if a gene/cell cluster should be reassigned and another gene/cell cluster should be split into two clusters. Default 10.
#' @param split.on.last After the the chain has converged according to 'stop.iter', a heuristic will be applied to determine if a gene/cell cluster should be reassigned and another gene/cell cluster should be split into two clusters. If a split occurs, then 'stop.iter' will be reset. Default TRUE.
#' @param count.checksum An MD5 checksum for the provided counts matrix
#' @param seed Parameter to set.seed() for random number generation
#' @param z.init Initial values of z. If NULL, z will be randomly sampled. Default NULL.
#' @param process.counts Whether to cast the counts matrix to integer and round(). Defaults to TRUE.
#' @param logfile If NULL, messages will be displayed as normal. If set to a file name, messages will be redirected messages to the file. Default NULL.
#' @return An object of class celda_C with clustering results and Gibbs sampling statistics
#' @export
celda_C = function(counts, sample.label=NULL, K, alpha=1, beta=1, 
                 	 stop.iter = 10, max.iter=200, split.on.iter=10, split.on.last=TRUE,
                 	 count.checksum=NULL, seed=12345,
                 	 z.init = NULL, process.counts=TRUE, logfile=NULL) {
  
  ## Error checking and variable processing
  if (isTRUE(process.counts)) {
    counts = processCounts(counts)  
  }
    
  s = processSampleLabels(sample.label, ncol(counts))
  if (is.null(sample.label)) {
    sample.label = s
  }
  
  ## Randomly select z and y or set z/y to supplied initial values
  z = initialize.cluster(K, ncol(counts), initial = z.init, fixed = NULL, seed=seed)
  z.best = z
  
  ## Calculate counts one time up front
  nS = length(unique(s))
  nG = nrow(counts)
  nM = ncol(counts)

  m.CP.by.S = matrix(as.integer(table(factor(z, levels=1:K), s)), ncol=nS)
  n.G.by.CP = t(rowsum.z(counts, z=z, K=K))
  n.CP = as.integer(colSums(n.G.by.CP))
  n.by.C = as.integer(colSums(counts))
  
  ll = cC.calcLL(m.CP.by.S=m.CP.by.S, n.G.by.CP=n.G.by.CP, s=s, K=K, nS=nS, nG=nG, alpha=alpha, beta=beta)

  set.seed(seed)
  logMessages(date(), "... Starting Gibbs sampling", logfile=logfile, append=FALSE)
  
  iter = 1L
  num.iter.without.improvement = 0L
  do.cell.split = TRUE
  while(iter <= max.iter & num.iter.without.improvement <= stop.iter) {
    
    next.z = cC.calcGibbsProbZ(counts=counts, m.CP.by.S=m.CP.by.S, n.G.by.CP=n.G.by.CP, n.by.C=n.by.C, n.CP=n.CP, z=z, s=s, K=K, nG=nG, nM=nM, alpha=alpha, beta=beta)
    m.CP.by.S = next.z$m.CP.by.S
    n.G.by.CP = next.z$n.G.by.CP
    n.CP = next.z$n.CP
    z = next.z$z
    
    ## Perform split on i-th iteration of no improvement in log likelihood
    if(K > 2 & (((iter == max.iter | num.iter.without.improvement == stop.iter) & isTRUE(split.on.last)) | (iter %% split.on.iter == 0 & isTRUE(do.cell.split)))) {

      logMessages(date(), " ... Determining if any cell clusters should be split.", logfile=logfile, append=TRUE, sep="")
      res = split.each.z(counts=counts, z=z, K=K, z.prob=t(next.z$probs), alpha=alpha, beta=beta, s=s, LLFunction="calculateLoglikFromVariables.celda_C")
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
      n.G.by.CP = t(rowsum.z(counts, z=z, K=K))
      n.CP = as.integer(colSums(n.G.by.CP))
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
    
    logMessages(date(), "... Completed iteration:", iter, "| logLik:", temp.ll, logfile=logfile, append=TRUE)
    iter = iter + 1    
  }
    
  names = list(row=rownames(counts), column=colnames(counts), sample=levels(sample.label))

  result = list(z=z.best, completeLogLik=ll,  
                finalLogLik=ll.best, seed=seed, K=K, 
                sample.label=sample.label, alpha=alpha, 
                beta=beta, count.checksum=count.checksum, 
                names=names)
  
  class(result) = "celda_C"
  result = reorder.celda_C(counts = counts, res = result)
  
  return(result)
}


# Gibbs sampling for the celda_C Model
cC.calcGibbsProbZ = function(counts, m.CP.by.S, n.G.by.CP, n.by.C, n.CP, z, s, K, nG, nM, alpha, beta, do.sample=TRUE) {

  ## Set variables up front outside of loop  
  probs = matrix(NA, ncol=nM, nrow=K)
  temp.n.G.by.CP = n.G.by.CP
  temp.n.CP = n.CP
  
  ix = sample(1:nM)
  for(i in ix) {
	
	## Subtract current cell counts from matrices
	m.CP.by.S[z[i],s[i]] = m.CP.by.S[z[i],s[i]] - 1L
	n.G.by.CP[,z[i]] = n.G.by.CP[,z[i]] - counts[,i]
	n.CP[z[i]] = n.CP[z[i]] - n.by.C[i]
  
	## Calculate probabilities for each state
	for(j in 1:K) {
	  temp.n.G.by.CP = n.G.by.CP
	  temp.n.G.by.CP[,j] = temp.n.G.by.CP[,j] + counts[,i]
	  temp.n.CP = n.CP
	  temp.n.CP[j] = temp.n.CP[j] + n.by.C[i]

	  probs[j,i] = 	log(m.CP.by.S[j,s[i]] + alpha) +		## Theta simplified
				  sum(lgamma(temp.n.G.by.CP + beta)) -	## Phi Numerator
				  sum(lgamma(temp.n.CP + (nG * beta)))	## Phi Denominator
	}  

	## Sample next state and add back counts
	if(isTRUE(do.sample)) z[i] = sample.ll(probs[,i])
	
	m.CP.by.S[z[i],s[i]] = m.CP.by.S[z[i],s[i]] + 1L
	n.G.by.CP[,z[i]] = n.G.by.CP[,z[i]] + counts[,i]
	n.CP[z[i]] = n.CP[z[i]] + n.by.C[i]
  }  

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
#' @export
simulateCells.celda_C = function(model, S=10, C.Range=c(10, 100), N.Range=c(100,5000), 
                         G=500, K=5, alpha=1, beta=1, seed=12345) {
 
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
  
  return(list(z=result$z, counts=cell.counts, sample.label=cell.sample.label, K=K, alpha=alpha, beta=beta, C.Range=C.Range, N.Range=N.Range, S=S))
}


#' Generate factorized matrices showing each feature's influence on the celda_C model clustering 
#' 
#' @param counts A numeric count matrix
#' @param celda.mod Object return from celda_C function
#' @param type A character vector containing one or more of "counts", "proportions", or "posterior". "counts" returns the raw number of counts for each entry in each matrix. "proportions" returns the counts matrix where each vector is normalized to a probability distribution. "posterior" returns the posterior estimates which include the addition of the Dirichlet concentration parameter (essentially as a pseudocount).
#' @export
factorizeMatrix.celda_C = function(celda.mod, counts, type=c("counts", "proportion", "posterior")) {

  K = celda.mod$K
  z = celda.mod$z
  alpha = celda.mod$alpha
  beta = celda.mod$beta
  sample.label = celda.mod$sample.label
  s = as.integer(as.factor(sample.label))
        
  nS = length(unique(s))
  nG = nrow(counts)
  nM = ncol(counts)

  m.CP.by.S = matrix(as.integer(table(factor(z, levels=1:K), s)), ncol=nS)
  n.G.by.CP = t(rowsum.z(counts, z=z, K=K))
  
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
#' @param s A vector indicating the sample for each cell (column) in the count matrix
#' @param z A numeric vector of cluster assignments
#' @param K The total number of clusters in z
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution
#' @param beta Non-zero concentration parameter for gene Dirichlet distribution
#' @param ... Additional parameters
#' @export
calculateLoglikFromVariables.celda_C = function(counts, s, z, K, alpha, beta) {
  
  ## Calculate for "Theta" component
  z = factor(z, 1:K)
  m.CP.by.S = table(z, s)
  nS = length(unique(s))
  n.CP.by.G = rowsum(t(counts), group=z, reorder=TRUE)
  nG = ncol(n.CP.by.G)
  
  final = cC.calcLL(m.CP.by.S=m.CP.by.S, n.G.by.CP=t(n.CP.by.G), s=s, z=z, K=K, nS=nS, nG=nG, alpha=alpha, beta=beta)
  return(final)
}


#' Calculates the conditional probability of each cell belong to each cluster given all other cluster assignments
#'
#' @param counts The original count matrix used in the model
#' @param celda.mod A model returned from the 'celda_C' function
#' @param log If FALSE, then the normalized conditional probabilities will be returned. If TRUE, then the unnormalized log probabilities will be returned.  
#' @return A list containging a matrix for the conditional cell cluster probabilities. 
#' @export
clusterProbability.celda_C = function(counts, celda.mod, log=FALSE) {

  z = celda.mod$z
  s = celda.mod$sample.label
  K = celda.mod$K
  alpha = celda.mod$alpha
  beta = celda.mod$beta
  
  nS = length(unique(s))
  nG = nrow(counts)
  nM = ncol(counts)
  m.CP.by.S = matrix(as.integer(table(factor(z, levels=1:K), s)), nrow=K, ncol=nS)
  n.G.by.CP = t(rowsum.z(counts, z=z, K=K))
  n.CP = as.integer(colSums(n.G.by.CP))
  n.by.C = as.integer(colSums(counts))
  
  next.z = cC.calcGibbsProbZ(counts=counts, m.CP.by.S=m.CP.by.S, n.G.by.CP=n.G.by.CP, n.by.C=n.by.C, n.CP=n.CP, z=z, s=s, K=K, nG=nG, nM=nM, alpha=alpha, beta=beta, do.sample=FALSE)  
  z.prob = t(next.z$probs)
  
  if(!isTRUE(log)) {
    z.prob = normalizeLogProbs(z.prob)
  }
   
  return(list(z.probability=z.prob))
}


#' @export
calculatePerplexity.celda_C = function(counts, celda.mod, precision=128) {
  if (!compareCountMatrix(counts, celda.mod)) {
    warning("Provided count matrix was not used to generate the provided celda model.")
  }
  
  # TODO Can try to turn into a single giant matrix multiplication by duplicating
  #     phi / theta / sl
  # TODO Cast to sparse matrices?
  factorized = factorizeMatrix(celda.mod, counts, "posterior")
  theta = log(factorized$posterior$sample.states)
  phi = log(factorized$posterior$gene.states)
  sl = celda.mod$sample.label
  
  inner.log.prob = (t(phi) %*% counts) + theta[, sl]  
  inner.log.prob = Rmpfr::mpfr(inner.log.prob, precision)
  inner.log.prob.exp = exp(inner.log.prob)
  
  log.px = 0
  for(i in 1:ncol(inner.log.prob.exp)) {
    log.px = log.px + Rmpfr::asNumeric(log(sum(inner.log.prob.exp[, i])))
  }
  
  perplexity = exp(-(log.px/sum(counts)))
  return(perplexity)
}  


reorder.celda_C = function(counts,res){
  if(res$K > 2) {
    res$z = as.integer(as.factor(res$z))
    fm <- factorizeMatrix(counts = counts, celda.mod = res)
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


