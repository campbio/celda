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


#' Simulate cells from the cell clustering generative model
#' 
#' @param S Total number of samples
#' @param C.Range Vector of length 2 given the range (min,max) of number of cells for each sample to be randomly generated from the uniform distribution
#' @param N.Range Vector of length 2 given the range (min,max) of number of counts for each cell to be randomly generated from the uniform distribution
#' @param G Total number of Genes to be simulated
#' @param K An integer or range of integers indicating the desired number of cell clusters (for celda_C / celda_CG models)
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution
#' @param beta Non-zero concentration parameter for gene Dirichlet distribution
#' @param model Dummy parameter for S3 dispatch
#' @param ... Unused arguments
#' @export
simulateCells.celda_C = function(model, S=10, C.Range=c(10, 100), N.Range=c(100,5000), 
                         G=500, K=5, alpha=1, beta=1, ...) {
  
  phi <- gtools::rdirichlet(K, rep(beta, G))
  theta <- gtools::rdirichlet(S, rep(alpha, K))
  
  ## Select the number of cells per sample
  nC <- sample(C.Range[1]:C.Range[2], size=S, replace=TRUE)  
  cell.sample <- rep(1:S, nC)
  
  ## Select state of the cells  
  cell.state <- unlist(lapply(1:S, function(i) sample(1:K, size=nC[i], prob=theta[i,], replace=TRUE)))
  cell.state = reorder.label.by.size(cell.state, K)$new.labels
    
  ## Select number of transcripts per cell
  nN <- sample(N.Range[1]:N.Range[2], size=length(cell.sample), replace=TRUE)
  
  ## Select transcript distribution for each cell
  cell.counts <- sapply(1:length(cell.sample), function(i) stats::rmultinom(1, size=nN[i], prob=phi[cell.state[i],]))
  
  rownames(cell.counts) = paste0("Gene_", 1:nrow(cell.counts))
  colnames(cell.counts) = paste0("Cell_", 1:ncol(cell.counts)) 
  cell.sample = paste0("Sample_", 1:S)[cell.sample]

  return(list(z=cell.state, counts=cell.counts, sample.label=cell.sample, K=K, alpha=alpha, beta=beta))
}


#' celda Cell Clustering Model
#' 
#' @param counts A numeric count matrix
#' @param sample.label A vector indicating the sample for each cell (column) in the count matrix
#' @param K An integer or range of integers indicating the desired number of cell clusters (for celda_C / celda_CG models)
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution
#' @param beta Non-zero concentration parameter for gene Dirichlet distribution
#' @param stop.iter Number of iterations without improvement in the log likelihood to stop the Gibbs sampler. Default 10.
#' @param split.on.iter On every 'split.on.iter' iteration, a heuristic will be applied to determine if a gene/cell cluster should be reassigned and another gene/cell cluster should be split into two clusters. Default 10.
#' @param max.iter Maximum iterations of Gibbs sampling to perform regardless of convergence. Default 200.
#' @param count.checksum An MD5 checksum for the provided counts matrix
#' @param seed Parameter to set.seed() for random number generation
#' @param z.init Initial values of z. If NULL, z will be randomly sampled. Default NULL.
#' @param process.counts Whether to cast the counts matrix to integer and round(). Defaults to TRUE.
#' @param logfile If NULL, messages will be displayed as normal. If set to a file name, messages will be redirected messages to the file. Default NULL.
#' @param ... additonal parameters
#' @return An object of class celda_C with clustering results and Gibbs sampling statistics
#' @export
celda_C = function(counts, sample.label=NULL, K,
					alpha=1, beta=1, 
                   	stop.iter = 10, split.on.iter=10, max.iter=200, 
                   	count.checksum=NULL, seed=12345,
                   	z.init = NULL, process.counts=TRUE, logfile=NULL, ...) {

  ## Error checking and variable processing
  if (process.counts) {
    counts = processCounts(counts)  
  }
    
  if(is.null(sample.label)) {
    s = rep(1, ncol(counts))
    sample.label = s 
  } else if(is.factor(sample.label)) {
    s = as.numeric(sample.label)
  } else {
    sample.label = as.factor(sample.label)
    s = as.integer(sample.label)
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
    if(K > 2 & (num.iter.without.improvement == stop.iter | (iter %% split.on.iter == 0 & isTRUE(do.cell.split)))) {

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







#' Calculate the celda_C log likelihood for user-provided cluster assignments
#' 
#' @param counts A numeric count matrix
#' @param s A vector indicating the sample for each cell (column) in the count matrix
#' @param z A numeric vector of cluster assignments
#' @param K The total number of clusters in z
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution
#' @param beta Non-zero concentration parameter for gene Dirichlet distribution
#' @param ... Additional parameters
calculateLoglikFromVariables.celda_C = function(counts, s, z, K, alpha, beta, ...) {
  
  ## Calculate for "Theta" component
  z = factor(z, 1:K)
  m.CP.by.S = table(z, s)
  nS = length(unique(s))
  
  a = nS * lgamma(K*alpha)
  b = sum(lgamma(m.CP.by.S + alpha))
  c = -nS * K * lgamma(alpha)
  d = -sum(lgamma(colSums(m.CP.by.S + alpha)))
  
  theta.ll = a + b + c + d
 
  ## Calculate for "Phi" component
  n.CP.by.G = rowsum(t(counts), group=z, reorder=TRUE)
  nG = ncol(n.CP.by.G)
  
  a = K * lgamma(nG * beta)
  b = sum(lgamma(n.CP.by.G + beta))
  c = -K * nG * lgamma(beta)
  d = -sum(lgamma(rowSums(n.CP.by.G + beta)))
  
  phi.ll = a + b + c + d

  final = theta.ll + phi.ll
  return(final)
}

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




################################################################################
# celda_C S3 methods                                                           #
################################################################################
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


#' visualizeModelPerformance for celda Cell clustering function
#' @param celda.list A celda_list object returned from celda()
#' @param method One of "perplexity", "loglik"
#' @param title Title for the plot
#' @param log Currently not working for celda_C objects
#' @import Rmpfr
#' @export
visualizeModelPerformance.celda_C = function(celda.list, method="perplexity", 
                                               title="Model Performance (All Chains)",
                                               log = F) {
  
  cluster.sizes = unlist(lapply(celda.list$res.list, function(mod) { getK(mod) }))
  log.likelihoods = lapply(celda.list$res.list,
                           function(mod) { completeLogLikelihood(mod) })
  performance.metric = lapply(log.likelihoods, 
                              calculatePerformanceMetric,
                              method)
  
  # These methods return Rmpfr numbers that are extremely small and can't be 
  # plotted, so log 'em first
  if (method %in% c("perplexity")) {
    performance.metric = lapply(performance.metric, log)
    performance.metric = methods::new("mpfr", unlist(performance.metric))
    performance.metric = as.numeric(performance.metric)
  }
  
  plot.df = data.frame(size=cluster.sizes,
                       metric=performance.metric)
  return(renderModelPerformancePlot(plot.df, "K", method, title))
}
