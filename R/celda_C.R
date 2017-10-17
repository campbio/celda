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
simulateCells.celda_C = function(model, S=10, C.Range=c(10, 100), N.Range=c(100,5000), 
                         G=500, K=5, alpha=1, beta=1) {
  
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
#' @param count.checksum An MD5 checksum for the provided counts matrix
#' @param max.iter Maximum iterations of Gibbs sampling to perform. Defaults to 25 
#' @param seed Parameter to set.seed() for random number generation
#' @param z.split.on.iter On every "z.split.on.iter" iteration, a heuristic will be applied using hierarchical clustering to determine if a cell cluster should be merged with another cell cluster and a third cell cluster should be split into two clusters. This helps avoid local optimum during the initialization.
#' @param z.num.splits Maximum number of times to perform the heuristic described in z.split.on.iter
#' @param z.init Initial values of z. If NULL, z will be randomly sampled. Default NULL.
#' @param logfile If NULL, messages will be displayed as normal. If set to a file name, messages will be redirected messages to the file. Default NULL.
#' @param ... additonal parameters
#' @return An object of class celda_C with clustering results and Gibbs sampling statistics
#' @export
celda_C = function(counts, sample.label=NULL, K, alpha=1, beta=1, 
                   count.checksum=NULL, max.iter=25, seed=12345,
                   z.split.on.iter=3, z.num.splits=3, 
                   z.init = NULL, logfile=NULL, ...) {
    
  if(is.null(sample.label)) {
    s = rep(1, ncol(counts))
    sample.label = s 
  } else if(is.factor(sample.label)) {
    s = as.numeric(sample.label)
  } else {
    sample.label = as.factor(sample.label)
    s = as.numeric(sample.label)
  }  
  
  set.seed(seed)
  log_messages(date(), "... Starting Gibbs sampling", logfile=logfile, append=FALSE)
  
  ## Randomly select z or set z to supplied initial values
  if(is.null(z.init)) {
    z = sample(1:K, ncol(counts), replace=TRUE)
  } else {
    if(length(unique(z.init)) != K || length(z.init) != ncol(counts)) {
      stop("'z.init' needs to be a combination of K unique values that is the same length as the number of columns in 'counts' matrix.")
    }
    z = as.numeric(as.factor(z.init))
  }  

  z.best = z
  
  ## Calculate counts one time up front
  nS = length(unique(s))
  nG = nrow(counts)
  nM = ncol(counts)

  m.CP.by.S = matrix(table(factor(z, levels=1:K), s), ncol=nS)
  n.CP.by.G = rowsum(t(counts), group=z, reorder=TRUE)
  n.CP = rowSums(n.CP.by.G)  

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
      n.CP[z[i]] = n.CP[z[i]] - sum(counts[,i])
    
      ## Calculate probabilities for each state
      probs = rep(NA, K)
      for(j in 1:K) {
        temp.n.CP.by.G = n.CP.by.G
        temp.n.CP.by.G[j,] = temp.n.CP.by.G[j,] + counts[,i]
        temp.n.CP = n.CP
        temp.n.CP[j] = temp.n.CP[j] + sum(counts[,i])

        probs[j] = cC.calcGibbsProbZ(m.CP.by.S=m.CP.by.S[j,s[i]], n.CP.by.G=temp.n.CP.by.G, n.CP=temp.n.CP, nG=nG, alpha=alpha, beta=beta)
      }  

      ## Sample next state and add back counts
      previous.z = z
      z[i] = sample.ll(probs)
      m.CP.by.S[z[i],s[i]] = m.CP.by.S[z[i],s[i]] + 1
      n.CP.by.G[z[i],] = n.CP.by.G[z[i],] + counts[,i]
      n.CP[z[i]] = n.CP[z[i]] + sum(counts[,i])

      ## Perform check for empty clusters; Do not allow on last iteration
      if(sum(z == previous.z[i]) == 0 & K > 2) {
        
        ## Split another cluster into two
        res = split.z(counts=counts, z=z, empty.K=previous.z[i], K=K, LLFunction="calculate_loglik_from_variables.celda_C", s=s, alpha=alpha, beta=beta)
        log_messages(res$message, logfile=logfile, append=TRUE)
        z = res$z
        
        ## Re-calculate variables
        m.CP.by.S = matrix(table(factor(z, levels=1:K), s), ncol=nS)
        n.CP.by.G = rowsum(t(counts), group=z, reorder=TRUE)
        n.CP = rowSums(n.CP.by.G)
      }
    }  
    
    ## Perform split if on i-th iteration defined by split.on.iter
    if(iter %% z.split.on.iter == 0 & z.num.of.splits.occurred <= z.num.splits & K > 2) {

      log_messages(date(), " ... Determining if any cell clusters should be split (", z.num.of.splits.occurred, " of ", z.num.splits, ")", logfile=logfile, append=TRUE, sep="")
      res = split.each.z(counts=counts, z=z, K=K, alpha=alpha, beta=beta, s=s, LLFunction="calculate_loglik_from_variables.celda_C")
      log_messages(res$message, logfile=logfile, append=TRUE)
      
      z = res$z
      z.num.of.splits.occurred = z.num.of.splits.occurred + 1

      ## Re-calculate variables
      m.CP.by.S = matrix(table(factor(z, levels=1:K), s), ncol=nS)
      n.CP.by.G = rowsum(t(counts), group=z, reorder=TRUE)
      n.CP = rowSums(n.CP.by.G)
    }


    ## Calculate complete likelihood
    temp.ll = cC.calcLL(m.CP.by.S=m.CP.by.S, n.CP.by.G=n.CP.by.G, s=s, K=K, nS=nS, alpha=alpha, beta=beta)
    if((all(temp.ll > ll)) | iter == 1) {
      z.best = z
      ll.best = temp.ll
    }
    ll = c(ll, temp.ll)
    
    log_messages(date(), "... Completed iteration:", iter, "| logLik:", temp.ll, logfile=logfile, append=TRUE)
    
    iter = iter + 1    
  }
    
  reordered.labels = reorder.label.by.size(z.best, K)
  z.final.reorder = reordered.labels$new.labels
  names = list(row=rownames(counts), column=colnames(counts), sample=levels(sample.label))

  result = list(z=z.final.reorder, completeLogLik=ll,  
                finalLogLik=ll.best, seed=seed, K=K, 
                sample.label=sample.label, alpha=alpha, 
                beta=beta, count.checksum=count.checksum, 
                names=names)
  
  class(result) = "celda_C"
  
  return(result)
}


cC.calcGibbsProbZ = function(m.CP.by.S, n.CP.by.G, n.CP, nG, alpha, beta) {
  
  ## Calculate for "Theta" component
  theta.ll = log(m.CP.by.S + alpha)
  
  ## Calculate for "Phi" component
  b = sum(lgamma(n.CP.by.G + beta))
  d = -sum(lgamma(n.CP + (nG*beta)))
  
  phi.ll = b + d
  
  final = theta.ll + phi.ll 
  return(final)
}


#' Calculates the conditional probability of each cell belong to each cluster given all other cluster assignments
#'
#' @param counts The original count matrix used in the model
#' @param celda.mod A model returned from the 'celda_C' function
#' @return A list containging a matrix for the conditional cell cluster probabilities. 
#' @export
cluster_probability.celda_C = function(counts, celda.mod) {

  z = celda.mod$z
  s = celda.mod$sample.label
  K = celda.mod$K
  alpha = celda.mod$alpha
  beta = celda.mod$beta
  
  nS = length(unique(s))
  nG = nrow(counts)
  nM = ncol(counts)
  m.CP.by.S = matrix(table(factor(z, levels=1:K), s), ncol=nS)
  n.CP.by.G = rowsum(t(counts), group=z, reorder=TRUE)
  n.CP = rowSums(n.CP.by.G)  

  z.prob = matrix(NA, ncol=K, nrow=ncol(counts))
  for(i in 1:ncol(counts)) {
  
	## Subtract current cell counts from matrices
	m.CP.by.S[z[i],s[i]] = m.CP.by.S[z[i],s[i]] - 1
	n.CP.by.G[z[i],] = n.CP.by.G[z[i],] - counts[,i]
	n.CP[z[i]] = n.CP[z[i]] - sum(counts[,i])

	## Calculate probabilities for each state
	for(j in 1:K) {
	  temp.n.CP.by.G = n.CP.by.G
	  temp.n.CP.by.G[j,] = temp.n.CP.by.G[j,] + counts[,i]
	  temp.n.CP = n.CP
	  temp.n.CP[j] = temp.n.CP[j] + sum(counts[,i])

	  z.prob[i,j] = cC.calcGibbsProbZ(m.CP.by.S=m.CP.by.S[j,s[i]], n.CP.by.G=temp.n.CP.by.G, n.CP=temp.n.CP, nG=nG, alpha=alpha, beta=beta)
	}  

	## Add back counts
	m.CP.by.S[z[i],s[i]] = m.CP.by.S[z[i],s[i]] + 1
	n.CP.by.G[z[i],] = n.CP.by.G[z[i],] + counts[,i]
	n.CP[z[i]] = n.CP[z[i]] + sum(counts[,i])
  }
  
  return(list(z.probability=normalizeLogProbs(z.prob)))
}







#' Calculate the celda_C log likelihood for user-provided cluster assignments
#' 
#' @param counts A numeric count matrix
#' @param s A vector indicating the sample for each cell (column) in the count matrix
#' @param z A numeric vector of cluster assignments
#' @param K The total number of clusters in z
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution
#' @param beta Non-zero concentration parameter for gene Dirichlet distribution
calculateLoglikFromVariables.celda_C = function(counts, s, z, K, alpha, beta, ...) {
  
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


#' Generate factorized matrices showing each feature's influence on the celda_C model clustering 
#' 
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
#' @param ... extra parameters passed onto the render_celda_heatmap
#' @export
celda_heatmap.celda_C = function(celda.mod, counts, ...) {
  render_celda_heatmap(counts, z=celda.mod$z, ...)
}


#' visualize_model_performance for celda Cell clustering function
#' @param celda.list A celda_list object returned from celda()
#' @param method One of "perplexity", "loglik"
#' @param title Title for the plot
#' @param log Currently not working for celda_C objects
#' @import Rmpfr
#' @export
visualize_model_performance.celda_C = function(celda.list, method="perplexity", 
                                               title="Model Performance (All Chains)",
                                               log = F) {
  
  cluster.sizes = unlist(lapply(celda.list$res.list, function(mod) { getK(mod) }))
  log.likelihoods = lapply(celda.list$res.list,
                           function(mod) { completeLogLikelihood(mod) })
  performance.metric = lapply(log.likelihoods, 
                              calculate_performance_metric,
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
  return(render_model_performance_plot(plot.df, "K", method, title))
}
