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


#' Calculate Log Likelihood For A Set of Gene Clusterings (Gene Clustering)
#'
#' This function calculates the log likelihood of each clustering of genes generated
#' over multiple iterations of Gibbs sampling.
#' 
#' @param counts A numeric count matrix
#' @param y A numeric vector of gene cluster assignments
#' @param L The number of clusters being considered
#' @param beta The Dirichlet distribution parameter for Phi; adds a pseudocount to each transcriptional state within each cell
#' @param gamma The Dirichlet distribution parameter for Psi; adds a pseudocount to each gene within each transcriptional state
#' @param delta The Dirichlet distribution parameter for Eta; adds a gene pseudocount to the numbers of genes each state
#' @keywords log likelihood
#' @examples TODO
#' @export
cG.calcLLFromVariables = function(counts, y, L, beta, delta) {
  n.TS.by.C <- rowsum(counts, group=y, reorder=TRUE)
  
  nM <- ncol(n.TS.by.C)
  
  a <- nM * lgamma(L * beta)
  b <- sum(lgamma(n.TS.by.C + beta))
  c <- -nM * L * lgamma(beta)
  d <- -sum(lgamma(colSums(n.TS.by.C + beta)))
  
  phi.ll <- a + b + c + d

  n.by.G <- rowSums(counts)
  n.by.TS <- as.numeric(rowsum(n.by.G, y))
  
  nG.by.TS <- table(y)
  nG <- nrow(counts)

  a <- sum(lgamma(nG.by.TS * delta))
  b <- sum(lgamma(n.by.G + delta))
  c <- -nG * lgamma(delta)
  d <- -sum(lgamma(n.by.TS + (nG.by.TS * delta)))
  
  psi.ll <- a + b + c + d

  final <- phi.ll + psi.ll
  return(final)
}


cG.calcLL = function(n.TS.by.C, n.by.TS, n.by.G, nG.by.TS, nM, nG, L, beta, delta) {
  
  ## Determine if any TS has 0 genes
  ## Need to remove 0 gene states as this will cause the likelihood to fail
  if(sum(nG.by.TS == 0) > 0) {
    ind = which(nG.by.TS > 0)
    L = length(ind)
    n.TS.by.C = n.TS.by.C[ind,]
    n.by.TS = n.by.TS[ind]
    nG.by.TS = nG.by.TS[ind]
  }

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

  final <- phi.ll + psi.ll  
  return(final)
}




#' Calculate Log Likelihood For Single Set of Cluster Assignments (Gene Clustering)
#'
#' This function calculates the log-likelihood of a given set of cluster assigments for the samples
#' represented in the provided count matrix.
#' 
#' @param ix The index of the cell being assigned a cluster during the current iteration of Gibbs sampling
#' @param counts A numeric count matrix
#' @param z A numeric vector of cluster assignments
#' @param k The number of clusters being considered
#' @param alpha Vector of non-zero concentration parameters for sample <-> cluster assignment Dirichlet distribution
#' @param beta Vector of non-zero concentration parameters for cluster <-> gene assignment Dirichlet distribution
#' @keywords log likelihood
#' @examples TODO
cG.calcGibbsProbY = function(n.TS.by.C, n.by.TS, nG.by.TS, beta, delta) {
 
  ## Determine if any TS has 0 genes
  ## Need to remove 0 gene states as this will cause the likelihood to fail
  if(sum(nG.by.TS == 0) > 0) {
    ind = which(nG.by.TS > 0)
    n.TS.by.C = n.TS.by.C[ind,]
    n.by.TS = n.by.TS[ind]
    nG.by.TS = nG.by.TS[ind]
  }
 
  ## Calculate for "Phi" component
  b <- sum(lgamma(n.TS.by.C + beta))
  d <- -sum(lgamma(colSums(n.TS.by.C + beta)))
  phi.ll <- b + d
  
  ## Calculate for "Psi" component
  a <- sum(lgamma(nG.by.TS * delta))
  d <- -sum(lgamma(n.by.TS + (nG.by.TS * delta)))
  psi.ll <- a + d
  
  final <- phi.ll + psi.ll
  return(final)
}



#' Cluster Genes from Single Cell Sequencing Data
#'
#' geneCluster provides cluster assignments for all genes in a provided single-cell 
#' sequencing count matrix, using the celda Bayesian hierarchical model.
#' 
#' @param counts A numeric count matrix
#' @param k The number of clusters to generate
#' @param a Vector of non-zero concentration parameters for sample <-> cluster assignment Dirichlet distribution
#' @param b Vector of non-zero concentration parameters for cluster <-> gene assignment Dirichlet distribution
#' @param g Number of cell states ("topics")
#' @param max.iter Maximum iterations of Gibbs sampling to perform. Defaults to 25.
#' @param min.cell Desired minimum number of cells per cluster
#' @param seed Parameter to set.seed() for random number generation
#' @param best Whether to return the cluster assignment with the highest log-likelihood. Defaults to TRUE. Returns last generated cluster assignment when FALSE.
#' @param kick Whether to randomize cluster assignments when a cluster has fewer than min.cell cells assigned to it during Gibbs sampling. (TODO param currently unused?)
#' @param converge Threshold at which to consider the Markov chain converged
#' @keywords LDA gene clustering gibbs
#' @examples TODO
#' @export
celda_G = function(counts, L, beta=1, delta=1, max.iter=25,
                       seed=12345, best=TRUE, y.split.on.iter=3, y.num.splits=3) {
  
  set.seed(seed)
  message(date(), " ... Starting Gibbs sampling")

  y <- sample(1:L, nrow(counts), replace=TRUE)
  y.all <- y

  ## Calculate counts one time up front
  n.TS.by.C = rowsum(counts, group=y, reorder=TRUE)
  n.by.G = rowSums(counts)
  n.by.TS = as.numeric(rowsum(n.by.G, y))
  nG.by.TS = table(y)
  nM = ncol(counts)
  nG = nrow(counts)
  
  ## Calculate initial log likelihood
  ll <- cG.calcLL(n.TS.by.C=n.TS.by.C, n.by.TS=n.by.TS, n.by.G=n.by.G, nG.by.TS=nG.by.TS, nM=nM, nG=nG, L=L, beta=beta, delta=delta)
  
  y.probs <- matrix(NA, nrow=nrow(counts), ncol=L)
  iter <- 1
  continue = TRUE
  y.num.of.splits.occurred = 1
  while(iter <= max.iter & continue == TRUE) {
    
    ## Begin process of Gibbs sampling for each cell
    ix <- sample(1:nrow(counts))
    for(i in ix) {
        
	  ## Subtract current gene counts from matrices
	  nG.by.TS[y[i]] = nG.by.TS[y[i]] - 1
	  n.by.TS[y[i]] = n.by.TS[y[i]] - n.by.G[i]
	  n.TS.by.C[y[i],] = n.TS.by.C[y[i],] - counts[i,]

	  ## Calculate probabilities for each state
	  probs = rep(NA, L)
	  for(j in 1:L) {
		temp.nG.by.TS = nG.by.TS
		temp.n.by.TS = n.by.TS
		temp.n.TS.by.C = n.TS.by.C
	  
		temp.nG.by.TS[j] = temp.nG.by.TS[j] + 1
		temp.n.by.TS[j] = temp.n.by.TS[j] + n.by.G[i]
		temp.n.TS.by.C[j,] = temp.n.TS.by.C[j,] + counts[i,]
	  
		probs[j] <- cG.calcGibbsProbY(n.TS.by.C=temp.n.TS.by.C, n.by.TS=temp.n.by.TS, nG.by.TS=temp.nG.by.TS, beta=beta, delta=delta)
	  }
	
	  ## Sample next state and add back counts
	  previous.y = y
	  y[i] <- sample.ll(probs)
	  nG.by.TS[y[i]] = nG.by.TS[y[i]] + 1
	  n.by.TS[y[i]] = n.by.TS[y[i]] + n.by.G[i]
	  n.TS.by.C[y[i],] = n.TS.by.C[y[i],] + counts[i,]

      ## Perform check for empty clusters; Do not allow on last iteration
      if(sum(y == previous.y[i]) == 0 & iter < max.iter & L > 2) {
        
        ## Split another cluster into two
        y = split.y(counts=counts, y=y, empty.L=previous.y[i], L=L, LLFunction="cG.calcLLFromVariables", beta=beta, delta=delta)
        
        ## Re-calculate variables
        n.TS.by.C = rowsum(counts, group=y, reorder=TRUE)
        n.by.TS = as.numeric(rowsum(n.by.G, y))
        nG.by.TS = table(y)
      }

      y.probs[i,] <- probs
    }

    ## Perform split if on i-th iteration defined by y.split.on.iter
    if(iter %% y.split.on.iter == 0 & y.num.of.splits.occurred <= y.num.splits & L > 2) {

      message(date(), " ... Determining if any gene clusters should be split (", y.num.of.splits.occurred, " of ", y.num.splits, ")")
      y = split.each.y(counts=counts, y=y, L=L, beta=beta, delta=delta, LLFunction="cG.calcLLFromVariables")
      y.num.of.splits.occurred = y.num.of.splits.occurred + 1

      ## Re-calculate variables
      n.TS.by.C = rowsum(counts, group=y, reorder=TRUE)
      n.by.TS = as.numeric(rowsum(n.by.G, y))
      nG.by.TS = table(y)
   }
    
    ## Save history
    y.all <- cbind(y.all, y)
    
    ## Calculate complete likelihood
    temp.ll <- cG.calcLL(n.TS.by.C=n.TS.by.C, n.by.TS=n.by.TS, n.by.G=n.by.G, nG.by.TS=nG.by.TS, nM=nM, nG=nG, L=L, beta=beta, delta=delta)
    if((best == TRUE & all(temp.ll > ll)) | iter == 1) {
      y.probs.final = y.probs
    }
    ll <- c(ll, temp.ll)

    message(date(), " ... Completed iteration: ", iter, " | logLik: ", temp.ll)

    iter <- iter + 1    
  }
  
  
  if(best == TRUE) {
    ix <- which.max(ll)
    y.final <- y.all[,ix]
    ll.final <- ll[ix]
  } else {
    y.final <- y
    ll.final <- tail(ll, n=1)
  }
  
  y.final.order = reorder.label.by.size(y.final, L)
  
  return(list(y=y.final.order, complete.y=y.all, completeLogLik=ll, 
              finalLogLik=ll.final, y.probability=y.probs, L=L, beta=beta, delta=delta,
              seed=seed))
}


#' Generate Count Data
#'
#' Generate a simulated count matrix, based off a generative distribution whose 
#' parameters can be provided by the user.
#' 
#' @param C The number of cells
#' @param L The number of transcriptional states
#' @param N.Range The range of counts each gene should have
#' @param G The number of genes for which to simulate counts
#' @param beta The Dirichlet distribution parameter for Phi; adds a pseudocount to each transcriptional state within each cell
#' @param delta The Dirichlet distribution parameter for Psi; adds a pseudocount to each gene within each transcriptional state
#' @param seed Parameter to set.seed() for random number generation
#' @examples TODO
#' @export
simulateCells.celda_G = function(C=100, N.Range=c(500,5000),  G=1000, 
                                         L=5, beta=1, gamma=1, delta=1, seed=12345) {
  set.seed(seed)
  eta = gtools::rdirichlet(1, rep(gamma, L))
  
  y = sample(1:L, size=G, prob=eta, replace=TRUE)
  if(length(table(y)) < L) {
    stop("Some states did not receive any genes after sampling. Try increasing G and/or setting gamma > 1.")
  }
  y = reorder.label.by.size(y, L)
  
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
    cell.dist = rmultinom(1, size=nN[i], prob=phi[i,])
    for(j in 1:L) {
      cell.counts[,i] = cell.counts[,i] + rmultinom(1, size=cell.dist[j], prob=psi[,j])
    }
  }

  return(list(y=y, counts=cell.counts, L=L, beta=beta, delta=delta, gamma=gamma, phi=phi, psi=psi, eta=eta, seed=seed))
}

