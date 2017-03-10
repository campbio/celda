#' Calculate Log Likelihood For A Set of Gene Clusterings (Gene Clustering)
#'
#' This function calculates the log likelihood of each clustering of genes generated
#' over multiple iterations of Gibbs sampling.
#' 
#' @param counts A numeric count matrix
#' @param z A numeric vector of cluster assignments
#' @param k The number of clusters being considered
#' @param alpha Vector of non-zero concentration parameters for sample <-> cluster assignment Dirichlet distribution
#' @param beta Vector of non-zero concentration parameters for cluster <-> gene assignment Dirichlet distribution
#' @param gamma The number of cell states ("topics")
#' @keywords log likelihood
#' @examples TODO
calcLL_gene_clustering = function(counts, z, k, alpha, beta, gamma) {
  n.phi <- rowsum(counts, group=z, reorder=TRUE)
  nk <- nrow(n.phi)
  nc <- ncol(n.phi)
  
  a <- nc*lgamma(nk*alpha)
  b <- sum(lgamma(n.phi+alpha))
  c <- -nc*nk*lgamma(alpha)
  d <- -sum(lgamma(apply(n.phi + alpha, 2, sum)))
  
  phi.ll <- a + b + c + d

  n.psi <- rowSums(counts)
  n.psi.sum <- as.numeric(rowsum(n.psi, z))
  
  ng.z <- table(z)
  ng <- length(n.psi)
  nk <- length(n.psi.sum)
  
  a <- sum(lgamma(ng.z * beta))
  b <- sum(lgamma(n.psi + beta))
  c <- -ng * lgamma(beta)
  d <- -sum(lgamma(n.psi.sum + (ng.z*beta)))
  
  psi.ll <- a + b + c + d
  
  a <- lgamma(nk*gamma)
  b <- sum(lgamma(ng.z+gamma))
  c <- -nk*lgamma(gamma)
  d <- -sum(lgamma(sum(ng.z + gamma)))
  
  eta.ll <- a + b + c + d

  final <- phi.ll + psi.ll + eta.ll
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
#' @param gamma The number of cell states ("topics")
#' @keywords log likelihood
#' @examples TODO
calcLLlite_gene_clustering = function(ix, counts, z, k, alpha, beta, gamma) {
  
  ng.z.minus <- table(z[-ix])
  eta.ll <- log(ng.z.minus[z[ix]] + gamma)
  
  n.phi <- rowsum(counts, group=z, reorder=TRUE)
  b <- sum(lgamma(n.phi+alpha))
  d <- -sum(lgamma(apply(n.phi + alpha, 2, sum)))
  phi.ll <- b + d
  
  n.psi <- rowSums(counts)
  n.psi.sum <- as.numeric(rowsum(n.psi, z))
  ng.z <- table(z)
  ng <- length(n.psi)
  a <- sum(lgamma(ng.z*beta))
  d <- -sum(lgamma(n.psi.sum + (ng.z*beta)))
  psi.ll <- a + d
  
  final <- eta.ll + phi.ll + psi.ll
  return(final)
}


#' Calculate Log Likelihood For Single Set of Cluster Assignments (Gene Clustering)
#'
#' This function calculates the log-likelihood of a cell's membership to each possible clusters,
#' given the cluster assignment for all other cells.
#' 
#' @param ix The index in z corresponding to the cell currently being considered during Gibbs sampling
#' @param r A numeric count matrix
#' @param z A numeric vector of cluster assignments
#' @param k The number of clusters being considered
#' @param a Vector of non-zero concentration parameters for sample <-> cluster assignment Dirichlet distribution
#' @param b Vector of non-zero concentration parameters for cluster <-> gene assignment Dirichlet distribution
#' @param g The number of cell states ("topics")
#' @examples TODO 
calcGibbsProb = function(ix, r, z, k, a, b, g) {
  final <- rep(NA, k)
  for(j in 1:k) {
    z[ix] <- j
    final[j] <- calcLLlite_gene_clustering(ix, counts=r, z=z, k=k, 
                                           alpha=a, beta=b, gamma=g)
  }  
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
geneCluster = function(counts, k, a=1, b=1, g=1, max.iter=25,  min.cell=5, 
                       seed=12345, best=TRUE, kick=TRUE, converge=1e-5) {
  
  set.seed(seed)
  cat(date(), "... Starting Gibbs sampling\n")
  
  co <- counts

  z <- sample(1:k, nrow(co), replace=TRUE)
  z.all <- z
  ll <- calcLL_gene_clustering(counts=co, z=z, k=k, alpha=a, beta=b, gamma=g)
  
  z.probs <- matrix(NA, nrow=nrow(co), ncol=k)
  
  iter <- 1
  continue <- TRUE
  while(iter <= max.iter & continue == TRUE) {
    
    ## Begin process of Gibbs sampling for each cell
    ix <- sample(1:nrow(co))
    for(i in ix) {
      probs <- calcGibbsProb(i, r=co, z=z, k=k, a=a, b=b, g=g)
      z[i] <- sample.ll(probs)
      z.probs[i,] <- probs
    }
    
    ## Save Z history
    z.all <- cbind(z.all, z)
    
    ## Calculate complete likelihood
    temp.ll <- calcLL_gene_clustering(counts=co, z=z, k=k, alpha=a, beta=b, gamma=g)
    ll <- c(ll, temp.ll)
    
    cat(date(), "... Completed iteration:", iter, "| logLik:", temp.ll, "\n")
    
    ## Normalize Z probabilties and test for convergence
    z.probs <- exp(sweep(z.probs, 1, apply(z.probs, 1, max), "-"))
    z.probs <- sweep(z.probs, 1, rowSums(z.probs), "/")
    f <- function(v) sort(v, decreasing=TRUE)[2]
    z.probs.second <- max(apply(z.probs, 1, f))
    z.ta <- table(z)
#    if(z.probs.second < converge & (min(z.ta) >= min.cell | kick==FALSE)) {
    if(z.probs.second < converge) {    
      continue <- FALSE
      cat("Maximum probability of a cell changing its state is ",  z.probs.second, 
          ". Exiting at iteration ", iter, ".", sep="")
    }
    
    iter <- iter + 1    
  }
  
  
  if(best == TRUE) {
    ix <- which.max(ll)
    z.final <- z.all[,ix]
    ll.final <- ll[ix]
  } else {
    z.final <- z
    ll.final <- tail(ll, n=1)
  }
  
  return(list(z=z.final, complete.z=z.all, completeLogLik=ll, 
              finalLogLik=ll.final, z.probability=z.probs))
}


#' Generate Count Data
#'
#' Generate a simulated count matrix, based off a generative distribution whose 
#' parameters can be provided by the user.
#' 
#' @param C The number of cells
#' @param N.Range The range of counts each gene should have
#' @param G The number of genes for which to simulate counts
##' @param a Vector of non-zero concentration parameters for sample <-> cluster assignment Dirichlet distribution
#' @param b Vector of non-zero concentration parameters for cluster <-> gene assignment Dirichlet distribution
#' @param g The number of cell states ("topics")' @param k The number of gene clusters to simulate from
#' @param seed Parameter to set.seed() for random number generation
#' @examples TODO
#' @export
generateCells_gene_clustering = function(C=100, N.Range=c(500,5000),  G=1000, 
                                         k=5, a=1, b=1, g=1, seed=12345) {
  set.seed(seed)
  eta = gtools::rdirichlet(1, rep(g, k))
  
  z = sample(1:k, size=G, prob=eta, replace=TRUE)
  if(length(table(z)) < k) {
    stop("Some states did not receive any genes after sampling. Try increasing G and/or setting gamma > 1.")
  }
  
  phi = matrix(0, nrow=G, ncol=k)
  for(i in 1:k) {
    ind = z == i
    phi[ind,i] = gtools::rdirichlet(1, rep(b, sum(ind)))
  }
  
  theta = gtools::rdirichlet(C, rep(a, k))
  
  ## Select number of transcripts per cell
  nN = sample(N.Range[1]:N.Range[2], size=C, replace=TRUE)
  
  ## Select transcript distribution for each cell
  cell.counts = matrix(0, nrow=G, ncol=C)
  for(i in 1:C) {
    cell.dist = rmultinom(1, size=nN[i], prob=theta[i,])
    for(j in 1:k) {
      cell.counts[,i] = cell.counts[,i] + rmultinom(1, size=cell.dist[j], prob=phi[,j])
    }
  }

  cc.rowsum0 = rowSums(cell.counts) > 0
  cell.counts = cell.counts[cc.rowsum0,]
  z = z[cc.rowsum0]
  
  return(list(z=z, counts=cell.counts, k=k, a=a, b=b, g=g, theta=theta, phi=phi, eta=eta, seed=seed))
}


#' Sample log-likelihood probabilities
#'
#' Given a set of log-likelihoods for cluster membership, return a single cluster assignment.
#' 
#' @param counts A numeric count matrix
#' @keywords log likelihood sample
#' @return A single integer in 1:k corresponding to a cluster assignment
#' @examples TODO
sample.ll = function(ll.probs) {
  probs.sub <- exp(ll.probs - max(ll.probs))
  probs.norm <- probs.sub / sum(probs.sub)
  probs.select <- sample(1:length(ll.probs), size=1, prob=probs.norm)
  return(probs.select)
}

