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
calculateLoglikFromVariables.celda_G = function(counts, y, L, beta, delta, gamma, ...) {
  n.TS.by.C <- rowsum.y(counts, y=y, L=L)
  
  nM <- ncol(n.TS.by.C)
  
  a <- nM * lgamma(L * beta)
  b <- sum(lgamma(n.TS.by.C + beta))
  c <- -nM * L * lgamma(beta)
  d <- -sum(lgamma(colSums(n.TS.by.C + beta)))
  
  phi.ll <- a + b + c + d

  n.by.G <- rowSums(counts)
  n.by.TS = as.numeric(rowsum.y(matrix(n.by.G,ncol=1), y=y, L=L))
  
  nG.by.TS = table(factor(y, 1:L))
  nG.by.TS[nG.by.TS == 0] = 1
  nG <- nrow(counts)

  a <- sum(lgamma(nG.by.TS * delta))
  b <- sum(lgamma(n.by.G + delta))
  c <- -nG * lgamma(delta)
  d <- -sum(lgamma(n.by.TS + (nG.by.TS * delta)))
  
  psi.ll <- a + b + c + d

  a <- lgamma(L * gamma)
  b <- sum(lgamma(nG.by.TS + gamma))
  c <- -L * lgamma(gamma)
  d <- -sum(lgamma(sum(nG.by.TS + gamma)))

  eta.ll <- a + b + c + d

  final <- phi.ll + psi.ll + eta.ll
  
  return(final)
}


cG.calcLL = function(n.C.by.TS, n.by.TS, n.by.G, nG.by.TS, nM, nG, L, beta, delta, gamma) {
  
  nG.by.TS[nG.by.TS == 0] = 1

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

reorder.celdaG = function(counts,res){
  #Reorder L
  fm <- factorizeMatrix(counts = counts, celda.mod = res)
  fm.norm <- t(normalizeCounts(fm$proportions$gene.states,scale.factor = 1))
  d <- dist((fm.norm),diag = TRUE, upper = TRUE)
  h <- hclust(d, method = "complete")
  res <- recodeClusterY(res,from = h$order,
                        to = c(1:nrow(fm$counts$cell.states)))
  return(res)
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
cG.calcGibbsProbY = function(n.C.by.TS, n.by.TS, nG.by.TS, nG.in.Y, beta, delta, gamma) {
  nG.by.TS[nG.by.TS == 0] = 1

  ## Calculate for "Eta" component
  b <- sum(lgamma(nG.by.TS + gamma))
  d <- -sum(lgamma(sum(nG.by.TS + gamma)))
  eta.ll <- b + d
    
  ## Calculate for "Phi" component
  phi.ll <- sum(lgamma(n.C.by.TS + beta))
  
  ## Calculate for "Psi" component
  a <- sum(lgamma(nG.by.TS * delta))
  d <- -sum(lgamma(n.by.TS + (nG.by.TS * delta)))
  psi.ll <- a + d
  
  final <- eta.ll + phi.ll + psi.ll
  return(final)
}



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
#' @param max.iter Maximum iterations of Gibbs sampling to perform. Defaults to 25.
#' @param count.checksum An MD5 checksum for the provided counts matrix
#' @param split.on.iter  On every split.on.iter iteration, a heuristic will be applied to determine if a gene cluster should be reassigned and another gene cluster should be split into two clusters. Default to be 5. 
#' @param num.splits Maximum number of times to perform the heuristic described in split.on.iter. Default 5.
#' @param seed Parameter to set.seed() for random number generation.
#' @param y.init Initial values of y. If NULL, y will be randomly sampled. Default NULL.
#' @param logfile The name of the logfile to redirect messages to.
#' @param ...  Additional parameters
#' @keywords LDA gene clustering gibbs
#' @export
celda_G = function(counts, L, beta=1, delta=1, gamma=1, max.iter=50,
                   count.checksum=NULL, seed=12345, 
                   split.on.iter=5,  num.splits=5, 
                   y.init=NULL, logfile=NULL, ...) {

  ## Error checking and variable processing
  counts = processCounts(counts)  
  counts.t = t(counts)
  
  ## Randomly select z and y or set z/y to supplied initial values
  y = initialize.cluster(L, nrow(counts), initial = y.init, fixed = NULL, seed=seed)
  y.best = y  

  ## Calculate counts one time up front
  n.C.by.TS = t(rowsum.y(counts, y=y, L=L))
  n.by.G = as.integer(rowSums(counts))
  n.by.TS = as.integer(rowsum.y(matrix(n.by.G,ncol=1), y=y, L=L))
  nG.by.TS = as.integer(table(factor(y, 1:L)))
  nM = ncol(counts)
  nG = nrow(counts)

  set.seed(seed)
  logMessages(date(), "... Starting Gibbs sampling", logfile=logfile, append=FALSE)
  
  ## Calculate initial log likelihood
  ll <- cG.calcLL(n.C.by.TS=n.C.by.TS, n.by.TS=n.by.TS, n.by.G=n.by.G, nG.by.TS=nG.by.TS, nM=nM, nG=nG, L=L, beta=beta, delta=delta, gamma=gamma)

  iter <- 1L
  continue = TRUE
  num.of.splits.occurred = 1L
  while(iter <= max.iter & continue == TRUE) {
    
    ## Begin process of Gibbs sampling for each cell
    ix <- sample(1:nrow(counts))
    for(i in ix) {
        
      ## Subtract current gene counts from matrices
      nG.by.TS[y[i]] = nG.by.TS[y[i]] - 1L
      n.by.TS[y[i]] = n.by.TS[y[i]] - n.by.G[i]
      n.C.by.TS[,y[i]] = n.C.by.TS[,y[i]] - counts.t[,i]

      ## Calculate probabilities for each state
      probs = rep(NA, L)
      for(j in 1:L) {
      
        temp.nG.by.TS = nG.by.TS 
        temp.n.by.TS = n.by.TS 
        temp.n.C.by.TS = n.C.by.TS
	  
        temp.nG.by.TS[j] = temp.nG.by.TS[j] + 1L
        temp.n.by.TS[j] = temp.n.by.TS[j] + n.by.G[i]
        temp.n.C.by.TS[,j] = temp.n.C.by.TS[,j] + counts.t[,i]

        #probs[j] <- cG.calcGibbsProbY(n.C.by.TS=temp.n.C.by.TS,
        #              n.by.TS=temp.n.by.TS, 
        #              nG.by.TS=temp.nG.by.TS, 
        #              nG.in.Y=temp.nG.by.TS[j], 
        #              beta=beta, delta=delta, gamma=gamma)
		pseudo.nG.by.TS = temp.nG.by.TS
		pseudo.nG.by.TS[temp.nG.by.TS == 0L] = 1L
		
		probs[j] = 	sum(lgamma(pseudo.nG.by.TS + gamma)) -					## Eta Numerator
					sum(lgamma(sum(pseudo.nG.by.TS + gamma))) +				## Eta Denominator
					sum(lgamma(temp.n.C.by.TS + beta)) +					## Phi Numerator
					sum(lgamma(pseudo.nG.by.TS * delta)) -					## Psi Numerator
					sum(lgamma(temp.n.by.TS + (pseudo.nG.by.TS * delta)))   ## Psi Denominator
      }
	
      ## Sample next state and add back counts
      y[i] <- sample.ll(probs)
      nG.by.TS[y[i]] = nG.by.TS[y[i]] + 1L
      n.by.TS[y[i]] = n.by.TS[y[i]] + n.by.G[i]
      n.C.by.TS[,y[i]] = n.C.by.TS[,y[i]] + counts.t[,i]
    }

    ## Perform split if on i-th iteration defined by y.split.on.iter
    if(iter %% split.on.iter == 0 & num.of.splits.occurred <= num.splits & L > 2) {
      logMessages(date(), " ... Determining if any gene clusters should be split (", num.of.splits.occurred, " of ", num.splits, ")", logfile=logfile, append=TRUE, sep="")
      res = split.each.y(counts=counts, y=y, L=L, beta=beta, delta=delta, gamma=gamma, LLFunction="calculateLoglikFromVariables.celda_G")
      logMessages(res$message, logfile=logfile, append=TRUE)
      
      y = res$y
      num.of.splits.occurred = num.of.splits.occurred + 1L

      ## Re-calculate variables
      n.C.by.TS = t(rowsum.y(counts, y=y, L=L))
      n.by.TS = as.integer(rowsum.y(matrix(n.by.G,ncol=1), y=y, L=L))
      nG.by.TS = as.integer(table(factor(y, 1:L)))
    }
     
    ## Calculate complete likelihood
    temp.ll <- cG.calcLL(n.C.by.TS=n.C.by.TS, n.by.TS=n.by.TS, n.by.G=n.by.G, nG.by.TS=nG.by.TS, nM=nM, nG=nG, L=L, beta=beta, delta=delta, gamma=gamma)
    if((all(temp.ll > ll)) | iter == 1) {
      y.best = y
      ll.best = temp.ll
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
  
  #result = reorder.celdaG(counts = counts, res = result)
  
  return(result)
}


#' Simulate cells from the gene clustering generative model
#'
#' Generate a simulated count matrix, based off a generative distribution whose 
#' parameters can be provided by the user.
#' 
#' @param C The number of cells
#' @param L The number of transcriptional states
#' @param N.Range Vector of length 2 given the range (min,max) of number of counts for each cell to be randomly generated from the uniform distribution
#' @param G The number of genes for which to simulate counts
#' @param beta The Dirichlet distribution parameter for Phi; adds a pseudocount to each transcriptional state within each cell
#' @param delta The Dirichlet distribution parameter for Psi; adds a pseudocount to each gene within each transcriptional state
#' @param gamma The Dirichlet distribution parameter for Psi; adds a pseudocount to each gene within each transcriptional state.
#' @param seed Parameter to set.seed() for random number generation
#' @param model Dummy parameter for S3 dispatch
#' @param ... Unused arguments
#' @export
simulateCells.celda_G = function(model, C=100, N.Range=c(500,5000),  G=1000, 
                                 L=5, beta=1, gamma=1, delta=1, seed=12345, ...) {
  set.seed(seed)
  eta = gtools::rdirichlet(1, rep(gamma, L))
  
  y = sample(1:L, size=G, prob=eta, replace=TRUE)
  if(length(table(y)) < L) {
    stop("Some states did not receive any genes after sampling. Try increasing G and/or setting gamma > 1.")
  }
  y = reorder.label.by.size(y, L)$new.labels
  
  psi = matrix(0, nrow=G, ncol=L)
  for(i in 1:L) {
    ind = y == i
    psi[ind,i] = gtools::rdirichlet(1, rep(delta, sum(ind)))
  }
  
  phi = gtools::rdirichlet(C, rep(beta, L))
  
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
  
  return(list(y=y, counts=cell.counts, L=L, beta=beta, delta=delta, gamma=gamma, phi=phi, psi=psi, eta=eta, seed=seed))
}



#' Calculates the conditional probability of each cell belong to each cluster given all other cluster assignments
#'
#' @param counts The original count matrix used in the model
#' @param celda.mod A model returned from the 'celda_G' function
#' @param log If FALSE, then the normalized conditional probabilities will be returned. If TRUE, then the unnormalized log probabilities will be returned.  
#' @return A list containging a matrix for the conditional cell cluster probabilities. 
#' @export
clusterProbability.celda_G = function(counts, celda.mod, log=FALSE) {

  y = celda.mod$y
  L = celda.mod$L
  delta = celda.mod$delta
  beta = celda.mod$beta
  gamma = celda.mod$gamma
  
  ## Calculate counts one time up front
  n.TS.by.C = rowsum.y(counts, y=y, L=L)
  n.by.G = rowSums(counts)
  n.by.TS = as.numeric(rowsum.y(matrix(n.by.G,ncol=1), y=y, L=L))
  nG.by.TS = table(factor(y, 1:L))
  nM = ncol(counts)
  nG = nrow(counts)

  y.prob = matrix(NA, ncol=L, nrow=nrow(counts))
  for(i in 1:nrow(counts)) {
	## Subtract current gene counts from matrices
	nG.by.TS[y[i]] = nG.by.TS[y[i]] - 1
	n.by.TS[y[i]] = n.by.TS[y[i]] - n.by.G[i]
	n.TS.by.C[y[i],] = n.TS.by.C[y[i],] - counts[i,]

	## Calculate probabilities for each state
	for(j in 1:L) {
	  temp.nG.by.TS = nG.by.TS 
	  temp.n.by.TS = n.by.TS 
	  temp.n.TS.by.C = n.TS.by.C 
	
	  temp.nG.by.TS[j] = temp.nG.by.TS[j] + 1
	  temp.n.by.TS[j] = temp.n.by.TS[j] + n.by.G[i]
	  temp.n.TS.by.C[j,] = temp.n.TS.by.C[j,] + counts[i,]
	
	  y.prob[i,j] <- cG.calcGibbsProbY(n.TS.by.C=temp.n.TS.by.C,
					n.by.TS=temp.n.by.TS, 
					nG.by.TS=temp.nG.by.TS, 
					nG.in.Y=temp.nG.by.TS[j], 
					beta=beta, delta=delta, gamma=gamma)
	}
	
	nG.by.TS[y[i]] = nG.by.TS[y[i]] + 1
    n.by.TS[y[i]] = n.by.TS[y[i]] + n.by.G[i]
    n.TS.by.C[y[i],] = n.TS.by.C[y[i],] + counts[i,]

  }
  
  if(!isTRUE(log)) {
    y.prob = normalizeLogProbs(y.prob)
  }
  
  return(list(y.probability=y.prob))
}



#' Generate factorized matrices showing each feature's influence on the celda_G model clustering 
#' 
#' @param counts A numeric count matrix
#' @param celda.mod Object return from celda_C function
#' @param type A character vector containing one or more of "counts", "proportions", or "posterior". "counts" returns the raw number of counts for each entry in each matrix. "proportions" returns the counts matrix where each vector is normalized to a probability distribution. "posterior" returns the posterior estimates which include the addition of the Dirichlet concentration parameter (essentially as a pseudocount).
#' @export
factorizeMatrix.celda_G = function(celda.mod, counts, type=c("counts", "proportion", "posterior")) {

  L = celda.mod$L
  y = celda.mod$y
  beta = celda.mod$beta
  delta = celda.mod$delta
  
  counts.list = c()
  prop.list = c()
  post.list = c()
  res = list()
  
  n.TS.by.C = rowsum(counts, group=y, reorder=TRUE)
  n.by.G = rowSums(counts)
  n.by.TS = as.numeric(rowsum(n.by.G, y))

  n.G.by.TS = matrix(0, nrow=length(y), ncol=L)
  for(i in 1:length(y)) {n.G.by.TS[i,y[i]] = n.by.G[i]}

  L.names = paste0("L", 1:L)
  colnames(n.TS.by.C) = celda.mod$names$column
  rownames(n.TS.by.C) = L.names
  colnames(n.G.by.TS) = L.names
  rownames(n.G.by.TS) = celda.mod$names$row
  
  if(any("counts" %in% type)) {
    counts.list = list(cell.states=n.TS.by.C, gene.states=n.G.by.TS)
    res = c(res, list(counts=counts.list))
  }
  if(any("proportion" %in% type)) {
    prop.list = list(cell.states = normalizeCounts(n.TS.by.C, scale.factor=1),
    							  gene.states = normalizeCounts(n.G.by.TS, scale.factor=1))
    res = c(res, list(proportions=prop.list))
  }
  if(any("posterior" %in% type)) {
  
    gs = n.G.by.TS
    gs[gs > 0] = gs[gs > 0] + delta
    gs = normalizeCounts(gs, scale.factor=1)

    post.list = list(cell.states = normalizeCounts(n.TS.by.C + beta, scale.factor=1),
    						    gene.states = gs)
    res = c(res, posterior = list(post.list))						    
  }
  
  return(res)
}



################################################################################
# celda_G S3 methods                                                           #
################################################################################
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


# TODO DRYer implementation in concert with celda_C
#' visualizeModelPerformance for the celda Gene function
#' @param celda.list A celda_list object returned from celda()
#' @param method One of "perplexity" or "loglik"
#' @param title Title for the plot
#' @param log Currently not working for celda.G objects
#' @import Rmpfr
#' @export
visualizeModelPerformance.celda_G = function(celda.list, method="perplexity",
                                               title="Model Performance (All Chains)",
                                               log = F) {
  
  cluster.sizes = unlist(lapply(celda.list$res.list, function(mod) { getL(mod) }))
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
  return(renderModelPerformancePlot(plot.df, "L", method, title))
}
