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


#' simulateCells for celda Cell clustering function 
#' @param S Total number of samples
#' @param C.Range Vector of length 2 given the range (min,max) of number of cells for each sample to be randomly generated from the uniform distribution
#' @param N.Range Vector of length 2 given the range (min,max) of number of counts for each cell to be randomly generated from the uniform distribution
#' @param G Total number of Genes to be simulated
#' @param K An integer or range of integers indicating the desired number of cell clusters (for celda_C / celda_CG models)
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution
#' @param beta Non-zero concentration parameter for gene Dirichlet distribution
#' @export
simulateCells.celda_C = function(S=10, C.Range=c(10, 100), N.Range=c(100,5000), 
                         G=500, K=5, alpha=1, beta=1) {
  
  phi <- gtools::rdirichlet(K, rep(beta, G))
  theta <- gtools::rdirichlet(S, rep(alpha, K))
  
  ## Select the number of cells per sample
  nC <- sample(C.Range[1]:C.Range[2], size=S, replace=TRUE)  
  cell.sample <- rep(1:S, nC)
  
  ## Select state of the cells  
  cell.state <- unlist(lapply(1:S, function(i) sample(1:K, size=nC[i], prob=theta[i,], replace=TRUE)))
  cell.state = reorder.label.by.size(cell.state, K)
    
  ## Select number of transcripts per cell
  nN <- sample(N.Range[1]:N.Range[2], size=length(cell.sample), replace=TRUE)
  
  ## Select transcript distribution for each cell
  cell.counts <- sapply(1:length(cell.sample), function(i) rmultinom(1, size=nN[i], prob=phi[cell.state[i],]))
  
  ## Ensure that there are no all-0 rows in the counts matrix, which violates a celda modeling
  ## constraint (columns are guarnteed at least one count):
  zero.row.idx = which(rowSums(cell.counts) == 0)
  cell.counts = cell.counts[-zero.row.idx, ]
  
  return(list(z=cell.state, counts=cell.counts, sample=cell.sample, K=K, alpha=alpha, beta=beta))
}


#' celda Cell clustering function 
#' @param counts A numeric count matrix
#' @param sample.label A vector indicating the sample for each cell (column) in the count matrix
#' @param K An integer or range of integers indicating the desired number of cell clusters (for celda_C / celda_CG models)
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution
#' @param beta Non-zero concentration parameter for gene Dirichlet distribution
#' @param max.iter Maximum iterations of Gibbs sampling to perform. Defaults to 25 
#' @param seed Parameter to set.seed() for random number generation
#' @param best Whether to return the cluster assignment with the highest log-likelihood. Defaults to TRUE. Returns last generated cluster assignment when FALSE.
#' @param z.split.on.iter On every "z.split.on.iter" iteration, a heuristic will be applied using hierarchical clustering to determine if a cell cluster should be merged with another cell cluster and a third cell cluster should be split into two clusters. This helps avoid local optimum during the initialization.
#' @param z.num.splits Maximum number of times to perform the heuristic described in “z.split.on.iter”.
#' @param thread The thread index, used for logging purposes
#' @param save.history Logical; whether to return the history of cluster assignments. Defaults to FALSE
#' @param save.prob Logical; whether to return the history of cluster assignment probabilities. Defaults to FALSE
#' @export
celda_C = function(counts, sample.label=NULL, K, alpha=1, beta=1, max.iter=25, 
                   seed=12345, best=TRUE, z.split.on.iter=3, z.num.splits=3, 
                   thread=1, save.history=FALSE, save.prob=FALSE, ...) {
  
  if(is.null(sample.label)) {
    s = rep(1, ncol(counts))
  } else if(is.factor(sample.label)) {
    s = as.numeric(sample.label)
  } else {
    sample.label = as.factor(sample.label)
    s = as.numeric(sample.label)
  }  
  
  set.seed(seed)
  message("Thread ", thread, " ", date(), " ... Starting Gibbs sampling")
  
  z = sample(1:K, ncol(counts), replace=TRUE)
  z.all = z
  z.stability = c(NA)
  z.probs = matrix(NA, nrow=ncol(counts), ncol=K)
  
  ## Calculate counts one time up front
  nS = length(unique(s))
  m.CP.by.S = matrix(table(factor(z, levels=1:K), s), ncol=nS)
  n.CP.by.G = rowsum(t(counts), group=z, reorder=TRUE)

  
  ll = cC.calcLL(m.CP.by.S=m.CP.by.S, n.CP.by.G=n.CP.by.G, s=s, K=K, nS=nS, alpha=alpha, beta=beta)

  iter = 1
  continue = TRUE
  z.num.of.splits.occurred = 1
  while(iter <= max.iter & continue == TRUE) {
    

    ## Begin process of Gibbs sampling for each cell
    ix = sample(1:ncol(counts))
    for(i in ix) {
      
      ## Subtract current cell counts from matrices
      m.CP.by.S[z[i],s[i]] = m.CP.by.S[z[i],s[i]] - 1
      n.CP.by.G[z[i],] = n.CP.by.G[z[i],] - counts[,i]
        
      ## Calculate probabilities for each state
      probs = rep(NA, K)
      for(j in 1:K) {
        temp.n.CP.by.G = n.CP.by.G
        temp.n.CP.by.G[j,] = temp.n.CP.by.G[j,] + counts[,i]
        probs[j] = cC.calcGibbsProbZ(m.CP.by.S=m.CP.by.S[j,s[i]], n.CP.by.G=temp.n.CP.by.G, alpha=alpha, beta=beta)
      }  

      ## Sample next state and add back counts
      previous.z = z
      z[i] = sample.ll(probs)
      m.CP.by.S[z[i],s[i]] = m.CP.by.S[z[i],s[i]] + 1
      n.CP.by.G[z[i],] = n.CP.by.G[z[i],] + counts[,i]
      
      ## Perform check for empty clusters; Do not allow on last iteration
      if(sum(z == previous.z[i]) == 0 & iter < max.iter & K > 2) {
        
        ## Split another cluster into two
        z = split.z(counts=counts, z=z, empty.K=previous.z[i], K=K, LLFunction="cC.calcLLFromVariables", s=s, alpha=alpha, beta=beta)
        
        ## Re-calculate variables
        m.CP.by.S = matrix(table(factor(z, levels=1:K), s), ncol=nS)
        n.CP.by.G = rowsum(t(counts), group=z, reorder=TRUE)
      }
       
      z.probs[i,] = probs
    }  
    
    ## Perform split if on i-th iteration defined by split.on.iter
    if(iter %% z.split.on.iter == 0 & z.num.of.splits.occurred <= z.num.splits & K > 2) {

      message("Thread ", thread, " ", date(), " ... Determining if any cell clusters should be split (", z.num.of.splits.occurred, " of ", z.num.splits, ")")
      z = split.each.z(counts=counts, z=z, K=K, alpha=alpha, beta=beta, s=s, LLFunction="cC.calcLLFromVariables")
      z.num.of.splits.occurred = z.num.of.splits.occurred + 1

      ## Re-calculate variables
      m.CP.by.S = matrix(table(factor(z, levels=1:K), s), ncol=nS)
      n.CP.by.G = rowsum(t(counts), group=z, reorder=TRUE)
    }

    ## Save history
    z.all = cbind(z.all, z)

    ## Normalize Z and Y marginal probabilties and calculate stability
    z.probs = normalizeLogProbs(z.probs)
    z.stability = c(z.stability, stability(z.probs))

    ## Calculate complete likelihood
    temp.ll = cC.calcLL(m.CP.by.S=m.CP.by.S, n.CP.by.G=n.CP.by.G, s=s, K=K, nS=nS, alpha=alpha, beta=beta)
    if((best == TRUE & all(temp.ll > ll)) | iter == 1) {
      z.probs.final = z.probs
    }
    ll = c(ll, temp.ll)
    
    message("Thread ", thread, " ",date(), " ... Completed iteration: ", iter, " | logLik: ", temp.ll)
    
    iter = iter + 1    
  }
  
  if (best == TRUE) {
    ix = which.max(ll)
    z.final = z.all[,ix]
    ll.final = ll[ix]
  } else {
    z.final = z
    ll.final = tail(ll, n=1)
  }
  
  z.final.reorder = reorder.label.by.size(z.final, K)
  names = list(row=rownames(counts), column=colnames(counts), sample=levels(sample.label))

  result = list(z=z.final.reorder, completeLogLik=ll,  finalLogLik=ll.final, seed=seed, K=K, sample.label=sample.label, alpha=alpha, beta=beta, names=names)
  if (save.prob) result$z.probability = z.probs else result$z.probability = NA
  if (save.history) result$complete.z = z.all else result$complete.z = NA
  
  class(result) = "celda_C"
  return(result)
}


cC.calcGibbsProbZ = function(m.CP.by.S, n.CP.by.G, alpha, beta) {
  
  ## Calculate for "Theta" component
  theta.ll = log(m.CP.by.S + alpha)
  
  ## Calculate for "Phi" component
  b = sum(lgamma(n.CP.by.G + beta))
  d = -sum(lgamma(rowSums(n.CP.by.G + beta)))
  
  phi.ll = b + d
  
  final = theta.ll + phi.ll 
  return(final)
}


#' calcLLFromVariables for celda Cell clustering function 
#' @param counts A numeric count matrix
#' @param s A vector indicating the sample for each cell (column) in the count matrix
#' @param z A numeric vector of cluster assignments
#' @param K An integer or range of integers indicating the desired number of cell clusters (for celda_C / celda_CG models)
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution
#' @param beta Non-zero concentration parameter for gene Dirichlet distribution
#' @export
cC.calcLLFromVariables = function(counts, s, z, K, alpha, beta) {
  
  ## Calculate for "Theta" component
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

cC.calcLL = function(m.CP.by.S, n.CP.by.G, s, z, K, nS, alpha, beta) {
  
  ## Calculate for "Theta" component
  a = nS * lgamma(K * alpha)
  b = sum(lgamma(m.CP.by.S + alpha))
  c = -nS * K * lgamma(alpha)
  d = -sum(lgamma(colSums(m.CP.by.S + alpha)))
  
  theta.ll = a + b + c + d
  
  ## Calculate for "Phi" component
  nG = ncol(n.CP.by.G)
  
  a = K * lgamma(nG * beta)
  b = sum(lgamma(n.CP.by.G + beta))
  c = -K * nG * lgamma(beta)
  d = -sum(lgamma(rowSums(n.CP.by.G + beta)))
  
  phi.ll = a + b + c + d
  
  final = theta.ll + phi.ll
  return(final)
}


#' factorizeMatrix for celda Cell clustering funciton 
#' @param counts A numeric count matrix
#' @param celda.obj Object return from celda_C function
#' @param type A character vector containing one or more of "counts", "proportions", or "posterior". "counts" returns the raw number of counts for each entry in each matrix. "proportions" returns the counts matrix where each vector is normalized to a probability distribution. "posterior" returns the posterior estimates which include the addition of the Dirichlet concentration parameter (essentially as a pseudocount).
#' @export
factorizeMatrix.celda_C = function(counts, celda.obj, type=c("counts", "proportion", "posterior")) {

  K = celda.obj$K
  z = celda.obj$z
  alpha = celda.obj$alpha
  beta = celda.obj$beta
  sample.label = celda.obj$sample.label

  counts.list = c()
  prop.list = c()
  post.list = c()
  res = list()
        
  nS = length(unique(sample.label))
  m.CP.by.S = matrix(table(factor(z, levels=1:K), sample.label), ncol=nS)
  n.G.by.CP = t(rowsum(t(counts), group=z, reorder=TRUE))

  K.names = paste0("K", 1:K)
  rownames(n.G.by.CP) = celda.obj$names$row
  colnames(n.G.by.CP) = K.names
  rownames(m.CP.by.S) = K.names
  colnames(m.CP.by.S) = celda.obj$names$sample
              
  if(any("counts" %in% type)) {
    counts.list = list(sample.states=m.CP.by.S, gene.states=n.G.by.CP)
    res = c(res, list(counts=counts.list))
  }
  if(any("proportion" %in% type)) {
    prop.list = list(sample.states = normalizeCounts(m.CP.by.S, scale.factor=1),
                     gene.states = normalizeCounts(n.G.by.CP, scale.factor=1))
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

#' completeClusterHistory for celda Cell clustering funciton 
#' @param celda.mod A celda model object of class "celda_C"
#' @export
completeClusterHistory.celda_C = function(celda.mod) {
  return(celda.mod$complete.z)
}

#' clusterProbabilities for celda Cell clustering function 
#' @param celda.mod A celda model object of class "celda_C"
#' @export
clusterProbabilities.celda_C = function(celda.mod) {
  return(celda.mod$z.probability)
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

#' celda_heatmap for celda Cell clustering function 
#' @param celda.mod A celda model object of class "celda_C"
#' @param counts A numeric count matrix
#' @export
celda_heatmap.celda_C = function(celda.mod, counts, ...) {
  render_celda_heatmap(counts, z=celda.mod$z, ...)
}


#' visualize_model_performance for celda Cell clustering function
#' @param celda.list A celda_list object returned from celda()
#' @param method One of “perplexity”, “harmonic”, or “loglik”
#' @param title Title for the plot
#' @export
#' @import Rmpfr
visualize_model_performance.celda_C = function(celda.list, method="perplexity", 
                                               title="Model Performance (All Chains)") {
  
  cluster.sizes = unlist(lapply(celda.list$res.list, function(mod) { getK(mod) }))
  log.likelihoods = lapply(celda.list$res.list,
                           function(mod) { completeLogLikelihood(mod) })
  performance.metric = lapply(log.likelihoods, 
                              calculate_performance_metric,
                              method)
  
  # These methods return Rmpfr numbers that are extremely small and can't be 
  # plotted, so log 'em first
  if (method %in% c("perplexity", "harmonic")) {
    performance.metric = lapply(performance.metric, log)
    performance.metric = new("mpfr", unlist(performance.metric))
    performance.metric = as.numeric(performance.metric)
  }
  
  plot.df = data.frame(size=cluster.sizes,
                       metric=performance.metric)
  return(celda::render_model_performance_plot(plot.df, "K", method, title))
}
