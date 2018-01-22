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

#' Calculate log lileklihood for the celda Cell and Gene clustering model, given a set of cell / gene cluster assignments
#' 
#' @param counts A numeric count matrix
#' @param s A numeric vector of sample labels corresponding to the samples each cell in the count matrix originates from
#' @param z A numeric vector of cell cluster assignments
#' @param y A numeric vector of gene cluster assignments
#' @param K The number of cell clusters
#' @param L The number of gene clusters
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution.
#' @param beta The Dirichlet distribution parameter for Phi; adds a pseudocount to each transcriptional state within each cell.
#' @param delta The Dirichlet distribution parameter for Eta; adds a gene pseudocount to the numbers of genes each state.
#' @param gamma The Dirichlet distribution parameter for Psi; adds a pseudocount to each gene within each transcriptional state.
#' @param ... Additional parameters 
calculateLoglikFromVariables.celda_CG = function(counts, s, z, y, K, L, alpha, beta, delta, gamma, ...) {
  
  ## Calculate for "Theta" component
  z = factor(z, 1:K)
  m = table(z, s)
  ns = ncol(m)
  
  a = ns * lgamma(K * alpha)
  b = sum(lgamma(m + alpha))
  c = -ns * K * lgamma(alpha)
  d = -sum(lgamma(colSums(m + alpha)))
  
  theta.ll = a + b + c + d
  
  
  ## Calculate for "Phi" component
  n.TS.by.C = rowsum.y(counts, y=y, L=L)
  n.CP.by.TS = rowsum.z(n.TS.by.C, z=z, K=K)
    
  a = K * lgamma(L * beta)
  b = sum(lgamma(n.CP.by.TS + beta))
  c = -K * L * lgamma(beta)
  d = -sum(lgamma(rowSums(n.CP.by.TS + beta)))
  
  phi.ll = a + b + c + d
  
  ## Calculate for "Psi" component
  n.by.G = rowSums(counts)
  n.by.TS = as.numeric(rowsum.y(matrix(n.by.G,ncol=1), y=y, L=L))

  nG.by.TS = table(factor(y, levels=1:L))
  nG.by.TS[nG.by.TS == 0] = 1
  nG = length(y)
  
  a = sum(lgamma(nG.by.TS * delta))
  b = sum(lgamma(n.by.G + delta))
  c = -nG * lgamma(delta)
  d = -sum(lgamma(n.by.TS + (nG.by.TS * delta)))
  
  psi.ll = a + b + c + d
    
  ## Calculate for "Eta" side
  a = lgamma(L * gamma)
  b = sum(lgamma(nG.by.TS + gamma))
  c = -L * lgamma(gamma)
  d = -lgamma(sum(nG.by.TS + gamma))

  eta.ll = a + b + c + d

  final = theta.ll + phi.ll + psi.ll + eta.ll
  return(final)
}

.calcLL = function(K, L, m.CP.by.S, n.CP.by.TS, n.by.G, n.by.TS, nG.by.TS, nS, nG, alpha, beta, delta, gamma) {
  nG.by.TS[nG.by.TS == 0] = 1
  
  ## Calculate for "Theta" component
  a = nS * lgamma(K*alpha)
  b = sum(lgamma(m.CP.by.S+alpha))
  c = -nS * K *lgamma(alpha)
  d = -sum(lgamma(colSums(m.CP.by.S + alpha)))
  
  theta.ll = a + b + c + d
  
  
  ## Calculate for "Phi" component
  a = K * lgamma(L * beta)
  b = sum(lgamma(n.CP.by.TS + beta))
  c = -K * L * lgamma(beta)
  d = -sum(lgamma(rowSums(n.CP.by.TS + beta)))
  
  phi.ll = a + b + c + d
  
  ## Calculate for "Psi" component
  a = sum(lgamma(nG.by.TS * delta))
  b = sum(lgamma(n.by.G + delta))
  c = -nG * lgamma(delta)
  d = -sum(lgamma(n.by.TS + (nG.by.TS * delta)))
  
  psi.ll = a + b + c + d
    
  ## Calculate for "Eta" side
  a = lgamma(L*gamma)
  b = sum(lgamma(nG.by.TS + gamma))
  c = -L*lgamma(gamma)
  d = -lgamma(sum(nG.by.TS + gamma))

  eta.ll = a + b + c + d

  final = theta.ll + phi.ll + psi.ll + eta.ll
  return(final)
}

cCG.calcGibbsProbZ = function(m.CP.by.S, n.TS.by.CP, n.TS.by.C, n.CP, n.by.C, z, s, L, K, nM, alpha, beta, do.sample=TRUE) {

  probs = matrix(NA, ncol=nM, nrow=K)
  ix = sample(1:nM)
  for(i in ix) {

	## Subtract current cell counts from matrices
	m.CP.by.S[z[i],s[i]] = m.CP.by.S[z[i],s[i]] - 1L
	n.TS.by.CP[,z[i]] = n.TS.by.CP[,z[i]] - n.TS.by.C[,i]
	n.CP[z[i]] = n.CP[z[i]] - n.by.C[i]

	## Calculate probabilities for each state
	for(j in 1:K) {
	  temp.n.TS.by.CP = n.TS.by.CP
	  temp.n.TS.by.CP[,j] = temp.n.TS.by.CP[,j] + n.TS.by.C[,i]
	  temp.n.CP = n.CP
	  temp.n.CP[j] = temp.n.CP[j] + n.by.C[i]

	  probs[j,i] = 	log(m.CP.by.S[j,s[i]] + alpha) +		## Theta simplified
				  sum(lgamma(temp.n.TS.by.CP + beta)) -		## Phi Numerator
				  sum(lgamma(temp.n.CP + (L * beta)))		## Phi Denominator

	}  

	## Sample next state and add back counts
	if(isTRUE(do.sample)) z[i] = sample.ll(probs[,i])
	
	m.CP.by.S[z[i],s[i]] = m.CP.by.S[z[i],s[i]] + 1L
	n.TS.by.CP[,z[i]] = n.TS.by.CP[,z[i]] + n.TS.by.C[,i]
	n.CP[z[i]] = n.CP[z[i]] + n.by.C[i]
  }
  
  return(list(m.CP.by.S=m.CP.by.S, n.TS.by.CP=n.TS.by.CP, n.CP=n.CP, z=z, probs=probs))
}

cCG.calcGibbsProbY = function(n.CP.by.TS, n.by.TS, nG.by.TS, n.CP.by.G, n.by.G, y, nG, L, beta, delta, gamma, do.sample=TRUE) {

  probs = matrix(NA, ncol=nG, nrow=L)
  ix = sample(1:nG)
  for(i in ix) {
	  
	## Subtract current gene counts from matrices
	nG.by.TS[y[i]] = nG.by.TS[y[i]] - 1L
	n.CP.by.TS[,y[i]] = n.CP.by.TS[,y[i]] - n.CP.by.G[,i]
	n.by.TS[y[i]] = n.by.TS[y[i]] - n.by.G[i]
   
	for(j in 1:L) {
	  ## Add in counts to each state and determine probability
	  temp.n.CP.by.TS = n.CP.by.TS 
	  temp.n.CP.by.TS[,j] = temp.n.CP.by.TS[,j] + n.CP.by.G[,i]
	  temp.n.by.TS = n.by.TS 
	  temp.n.by.TS[j] = temp.n.by.TS[j] + n.by.G[i]
	  temp.nG.by.TS = nG.by.TS 
	  temp.nG.by.TS[j] = temp.nG.by.TS[j] + 1L

	  pseudo.nG.by.TS = temp.nG.by.TS
	  pseudo.nG.by.TS[temp.nG.by.TS == 0L] = 1L

	  probs[j,i] <-	sum(lgamma(pseudo.nG.by.TS + gamma)) -					## Eta Numerator
				  sum(lgamma(sum(pseudo.nG.by.TS + gamma))) +				## Eta Denominator
				  sum(lgamma(temp.n.CP.by.TS + beta)) +						## Phi Numerator
				  sum(lgamma(pseudo.nG.by.TS * delta)) -					## Psi Numerator
				  sum(lgamma(temp.n.by.TS + (pseudo.nG.by.TS * delta)))		## Psi Denominator
	}  

	## Sample next state and add back counts
	if(isTRUE(do.sample)) y[i] = sample.ll(probs[,i])
	
	nG.by.TS[y[i]] = nG.by.TS[y[i]] + 1L
	n.CP.by.TS[,y[i]] = n.CP.by.TS[,y[i]] + n.CP.by.G[,i]
	n.by.TS[y[i]] = n.by.TS[y[i]] + n.by.G[i]
  }

  return(list(n.CP.by.TS=n.CP.by.TS, nG.by.TS=nG.by.TS, n.by.TS=n.by.TS, y=y, probs=probs))
}

#' Simulate cells from the cell/gene clustering generative model
#' 
#' @param S The number of samples
#' @param C.Range two element vector to specify the lower and upper bound of the counts of cells for each sample
#' @param N.Range two element vector to specify the lower and upper bound of the counts of the transcripts
#' @param G The number of genes
#' @param K The number of cell populations
#' @param L The number of gene clusters being considered
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution
#' @param beta The Dirichlet distribution parameter for Phi; adds a pseudocount to each transcriptional state within each cell
#' @param gamma The Dirichlet distribution parameter for Psi; adds a pseudocount to each gene within each transcriptional state
#' @param delta The Dirichlet distribution parameter for Eta; adds a gene pseudocount to the numbers of genes each state
#' @param seed starting point used for generating simulated data
#' @param ... Unused arguments
#' @param model Dummy parameter for S3 dispatch
#' @export
simulateCells.celda_CG = function(model, S=10, C.Range=c(50,100), N.Range=c(500,5000), 
                                  G=1000, K=3, L=10, alpha=1, beta=1, gamma=1, 
                                  delta=1, seed=12345, ...) {
  
  set.seed(seed)

  ## Number of cells per sample
  nC = sample(C.Range[1]:C.Range[2], size=S, replace=TRUE)
  nC.sum = sum(nC)
  cell.sample.label = rep(1:S, nC)
  
  ## Select number of transcripts per cell
  nN = sample(N.Range[1]:N.Range[2], size=length(cell.sample.label), replace=TRUE)
  
  ## Generate cell population distribution for each sample
  theta = t(gtools::rdirichlet(S, rep(alpha, K)))

  ## Assign cells to cellular subpopulations
  z = unlist(lapply(1:S, function(i) sample(1:K, size=nC[i], prob=theta[,i], replace=TRUE)))

  ## Generate transcriptional state distribution for each cell subpopulation
  phi = gtools::rdirichlet(K, rep(beta, L))

  ## Assign genes to transcriptional states 
  eta = gtools::rdirichlet(1, rep(gamma, L))
  y = sample(1:L, size=G, prob=eta, replace=TRUE)
  if(length(table(y)) < L) {
    stop("Some transcriptional states did not receive any genes after sampling. Try increasing G and/or setting gamma > 1.")
  }

  psi = matrix(0, nrow=G, ncol=L)
  for(i in 1:L) {
    ind = y == i
    psi[ind,i] = gtools::rdirichlet(1, rep(delta, sum(ind)))
  }
  
  ## Select transcript distribution for each cell
  cell.counts = matrix(0, nrow=G, ncol=nC.sum)
  for(i in 1:nC.sum) {
    transcriptional.state.dist = as.numeric(stats::rmultinom(1, size=nN[i], prob=phi[z[i],]))
    for(j in 1:L) {
      if(transcriptional.state.dist[j] > 0) {
        cell.counts[,i] = cell.counts[,i] + stats::rmultinom(1, size=transcriptional.state.dist[j], prob=psi[,j])
      }
    }  
  }
  
  new = reorder.labels.by.size.then.counts(cell.counts, z=z, y=y, K=K, L=L)
  
  ## Ensure that there are no all-0 rows in the counts matrix, which violates a celda modeling
  ## constraint (columns are guarnteed at least one count):
  zero.row.idx = which(rowSums(cell.counts) == 0)
  if (length(zero.row.idx > 0)) {
    cell.counts = cell.counts[-zero.row.idx, ]
    new$y = new$y[-zero.row.idx]
  } 
 
  ## Assign gene/cell/sample names 
  rownames(cell.counts) = paste0("Gene_", 1:nrow(cell.counts))
  colnames(cell.counts) = paste0("Cell_", 1:ncol(cell.counts))
  cell.sample.label = paste0("Sample_", 1:S)[cell.sample.label]

  return(list(z=new$z, y=new$y, sample.label=cell.sample.label, counts=cell.counts, K=K, L=L, C.Range=C.Range, N.Range=N.Range, S=S, alpha=alpha, beta=beta, gamma=gamma, delta=delta, theta=theta, phi=phi, psi=psi, eta=eta, seed=seed))
}

#' celda Cell and Gene Clustering Model
#' 
#' @param counts A numeric count matrix.
#' @param sample.label A vector indicating the sample for each cell in the count matrix
#' @param K The number of cell populations
#' @param L The number of gene clusters being considered
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution
#' @param beta The Dirichlet distribution parameter for Phi; adds a pseudocount to each transcriptional state within each cell. Default to 1
#' @param delta The Dirichlet distribution parameter for Eta; adds a gene pseudocount to the numbers of genes each state. Default to 1
#' @param gamma The Dirichlet distribution parameter for Psi; adds a pseudocount to each gene within each transcriptional state. Default to 1
#' @param count.checksum An MD5 checksum for the provided counts matrix
#' @param max.iter Maximum iterations of Gibbs sampling to perform. Defaults to 25
#' @param seed Parameter to set.seed() for random number generation
#' @param split.on.iter  On every split.on.iter iteration, a heuristic will be applied to determine if a gene or cell cluster should be reassigned and another gene or cell cluster should be split into two clusters. Default to be 10. 
#' @param num.splits Maximum number of times to perform the heuristic described in split.on.iter. Default 3.
#' @param z.init Initial values of z. If NULL, z will be randomly sampled. Default NULL.
#' @param y.init Initial values of y. If NULL, y will be randomly sampled. Default NULL.
#' @param logfile The name of the logfile to redirect messages to.
#' @param ... Additional parameters
#' @export
celda_CG = function(counts, sample.label=NULL, K, L, alpha=1, beta=1, 
                    delta=1, gamma=1, count.checksum=NULL, max.iter=50,
			              seed=12345, split.on.iter=10, num.splits=3,
			              z.init = NULL, y.init = NULL, logfile=NULL, ...) {
  
  ## Error checking and variable processing
  counts = processCounts(counts)
    
  if(is.null(sample.label)) {
    s = rep(1, ncol(counts))
    sample.label = s
  } else if(is.factor(sample.label)) {
    s = as.numeric(sample.label)
  } else {
    sample.label = as.factor(sample.label)
    s = as.numeric(sample.label)
  }  
  
  ## Randomly select z and y or set z/y to supplied initial values
  z = initialize.cluster(K, ncol(counts), initial = z.init, fixed = NULL, seed=seed)
  y = initialize.cluster(L, nrow(counts), initial = y.init, fixed = NULL, seed=seed)
  z.best = z
  y.best = y  
  
  ## Calculate counts one time up front
  nS = length(unique(s))
  m.CP.by.S = matrix(as.integer(table(factor(z, levels=1:K), s)), ncol=nS)
  n.TS.by.C = rowsum.y(counts, y=y, L=L)
  n.CP.by.TS = rowsum.z(n.TS.by.C, z=z, K=K)
  n.CP = as.integer(rowSums(n.CP.by.TS))
  n.by.G = as.integer(rowSums(counts))
  n.by.C = as.integer(colSums(counts))
  n.by.TS = as.integer(rowsum.y(matrix(n.by.G,ncol=1), y=y, L=L))
  nG.by.TS = as.integer(table(factor(y, 1:L)))

  nG = nrow(counts)
  nM = ncol(counts)
  
  ll = .calcLL(K=K, L=L, m.CP.by.S=m.CP.by.S, n.CP.by.TS=n.CP.by.TS, n.by.G=n.by.G, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, nS=nS, nG=nG, alpha=alpha, beta=beta, delta=delta, gamma=gamma)

  set.seed(seed)
  logMessages(date(), "... Starting Gibbs sampling", logfile=logfile, append=FALSE)
  
  iter = 1L
  continue = TRUE
  num.of.splits.occurred = 1L
  while(iter <= max.iter & continue == TRUE) {
    
    ## Gibbs sampling for each cell
    n.TS.by.C = rowsum.y(counts, y=y, L=L)
    n.TS.by.CP = t(n.CP.by.TS)
	next.z = cCG.calcGibbsProbZ(m.CP.by.S=m.CP.by.S, n.TS.by.CP=n.TS.by.CP, n.TS.by.C=n.TS.by.C, n.CP=n.CP, n.by.C=n.by.C, z=z, s=s, L=L, K=K, nM=nM, alpha=alpha, beta=beta)
    m.CP.by.S = next.z$m.CP.by.S
    n.TS.by.CP = next.z$n.TS.by.CP
    n.CP = next.z$n.CP
    z = next.z$z

    ## Gibbs sampling for each gene
    n.CP.by.G = rowsum.z(counts, z=z, K=K)
    n.CP.by.TS = t(n.TS.by.CP) 
	next.y = cCG.calcGibbsProbY(n.CP.by.TS=n.CP.by.TS, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, n.CP.by.G=n.CP.by.G, n.by.G=n.by.G, y=y, nG=nG, L=L, beta=beta, delta=delta, gamma=gamma)
	n.CP.by.TS = next.y$n.CP.by.TS
	nG.by.TS = next.y$nG.by.TS
	n.by.TS = next.y$n.by.TS
	y = next.y$y
    
    ## Perform split on i-th iteration defined by split.on.iter
    if(iter %% split.on.iter == 0 & num.of.splits.occurred <= num.splits) {
      if(K > 2) {
		logMessages(date(), " ... Determining if any cell clusters should be split (", num.of.splits.occurred, " of ", num.splits, ")", logfile=logfile, append=TRUE, sep="")
		res = split.each.z(counts=counts, z=z, y=y, z.prob=t(next.z$probs), K=K, L=L, alpha=alpha, delta=delta, beta=beta, gamma=gamma, s=s, LLFunction="calculateLoglikFromVariables.celda_CG")
		logMessages(res$message, logfile=logfile, append=TRUE)
		z = res$z      

		## Re-calculate variables
		m.CP.by.S = matrix(as.integer(table(factor(z, levels=1:K), s)), ncol=nS)
		n.TS.by.C = rowsum.y(counts, y=y, L=L)
		n.CP.by.TS = rowsum.z(n.TS.by.C, z=z, K=K)
		n.CP = as.integer(rowSums(n.CP.by.TS))
		n.CP.by.G = rowsum.z(counts, z=z, K=K)      
	  }  
      if(L > 2) {
        logMessages(date(), " ... Determining if any gene clusters should be split (", num.of.splits.occurred, " of ", num.splits, ")", logfile=logfile, append=TRUE, sep="")
        res = split.each.y(counts=counts, z=z, y=y, y.prob=t(next.y$probs), K=K, L=L, alpha=alpha, beta=beta, delta=delta, gamma=gamma, s=s, LLFunction="calculateLoglikFromVariables.celda_CG")
	    logMessages(res$message, logfile=logfile, append=TRUE)
        y = res$y

        ## Re-calculate variables
        n.TS.by.C = rowsum.y(counts, y=y, L=L)
        n.CP.by.TS = rowsum.z(n.TS.by.C, z=z, K=K)
        n.CP = as.integer(rowSums(n.CP.by.TS))
        n.by.TS = as.integer(rowsum.y(matrix(n.by.G,ncol=1), y=y, L=L))
        nG.by.TS = as.integer(table(factor(y, levels=1:L)))
      }
      num.of.splits.occurred = num.of.splits.occurred + 1L
    }

    ## Calculate complete likelihood
    temp.ll = .calcLL(K=K, L=L, m.CP.by.S=m.CP.by.S, n.CP.by.TS=n.CP.by.TS, n.by.G=n.by.G, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, nS=nS, nG=nG, alpha=alpha, beta=beta, delta=delta, gamma=gamma)
    if((all(temp.ll > ll)) | iter == 1) {
      z.best = z
      y.best = y
      ll.best = temp.ll
    }
    ll = c(ll, temp.ll)
    
    logMessages(date(), " ... Completed iteration: ", iter, " | logLik: ", temp.ll, logfile=logfile, append=TRUE, sep="")
    iter = iter + 1    
  }
    

  names = list(row=rownames(counts), column=colnames(counts), 
               sample=levels(sample.label))
  
  result = list(z=z.best, y=y.best, completeLogLik=ll, 
                finalLogLik=ll.best, K=K, L=L, alpha=alpha, 
                beta=beta, delta=delta, gamma=gamma, seed=seed, 
                sample.label=sample.label, names=names,
                count.checksum=count.checksum)
  
  class(result) = "celda_CG" 
   
  ## Peform reordering on final Z and Y assigments:
  result = reorder.celda_CG(counts = counts, res = result)
  return(result)
}


#' Generate factorized matrices showing each feature's influence on the celda_CG model clustering 
#' 
#' @param counts A numerix count matrix
#' @param celda.mod object returned from celda_CG function 
#' @param type one of the "counts", "proportion", or "posterior". 
#' @return A list of factorized matrices, of the types requested by the user. NOTE: "population" state matrices are always returned in cell population (rows) x transcriptional states (cols).
#' @export 
factorizeMatrix.celda_CG = function(celda.mod, counts, type=c("counts", "proportion", "posterior")) {

  K = celda.mod$K
  L = celda.mod$L
  z = celda.mod$z
  y = celda.mod$y
  alpha = celda.mod$alpha
  beta = celda.mod$beta
  delta = celda.mod$delta
  gamma = celda.mod$gamma
  sample.label = celda.mod$sample.label
  
  counts.list = c()
  prop.list = c()
  post.list = c()
  res = list()

  counts.list = c()
  prop.list = c()
  post.list = c()
  res = list()
  
  nS = length(unique(sample.label))
  m.CP.by.S = matrix(table(factor(z, levels=1:K), sample.label), ncol=nS)
  n.TS.by.C = rowsum(counts, group=y, reorder=TRUE)
  n.CP.by.TS = rowsum(t(n.TS.by.C), group=z, reorder=TRUE)
  n.by.G = rowSums(counts)
  n.by.TS = as.numeric(rowsum(n.by.G, y))

  n.G.by.TS = matrix(0, nrow=length(y), ncol=L)
  for(i in 1:length(y)) {n.G.by.TS[i,y[i]] = n.by.G[i]}

  L.names = paste0("L", 1:L)
  K.names = paste0("K", 1:K)
  colnames(n.TS.by.C) = celda.mod$names$column
  rownames(n.TS.by.C) = L.names
  colnames(n.G.by.TS) = L.names
  rownames(n.G.by.TS) = celda.mod$names$row
  rownames(m.CP.by.S) = K.names
  colnames(m.CP.by.S) = celda.mod$names$sample
  colnames(n.CP.by.TS) = L.names
  rownames(n.CP.by.TS) = K.names
    
  if(any("counts" %in% type)) {
    counts.list = list(sample.states = m.CP.by.S,
            				   population.states = t(n.CP.by.TS), 
            				   cell.states = n.TS.by.C,
            				   gene.states = n.G.by.TS)
    res = c(res, list(counts=counts.list))
  }
  if(any("proportion" %in% type)) {
    prop.list = list(sample.states = normalizeCounts(m.CP.by.S, scale.factor=1),
    				   population.states = normalizeCounts(t(n.CP.by.TS), scale.factor=1), 
    				   cell.states = normalizeCounts(n.TS.by.C, scale.factor=1),
    				   gene.states = normalizeCounts(n.G.by.TS, scale.factor=1))
    res = c(res, list(proportions=prop.list))
  }
  if(any("posterior" %in% type)) {

    gs = n.G.by.TS
    gs[gs > 0] = gs[gs > 0] + delta
    gs = normalizeCounts(gs, scale.factor=1)

    post.list = list(sample.states = normalizeCounts(m.CP.by.S + alpha, scale.factor=1),
          				   population.states = normalizeCounts(t(n.CP.by.TS + beta), scale.factor=1), 
          				   gene.states = gs)
    res = c(res, posterior = list(post.list))						    
  }
  
  return(res)
}


#' Calculates the conditional probability of each cell belong to each cluster given all other cluster assignments
#'
#' @param counts The original count matrix used in the model
#' @param celda.mod A model returned from the 'celda_CG' function
#' @param log If FALSE, then the normalized conditional probabilities will be returned. If TRUE, then the unnormalized log probabilities will be returned.  
#' @return A list containging a matrix for the conditional cell cluster probabilities. 
#' @export
clusterProbability.celda_CG = function(counts, celda.mod, log=FALSE) {
  
  s = celda.mod$sample.label
  z = celda.mod$z
  K = celda.mod$K  
  y = celda.mod$y
  L = celda.mod$L
  alpha = celda.mod$alpha
  delta = celda.mod$delta
  beta = celda.mod$beta
  gamma = celda.mod$gamma

  nS = length(unique(s))
  nG = nrow(counts)
  nM = ncol(counts)
  m.CP.by.S = matrix(as.integer(table(factor(z, levels=1:K), s)), ncol=nS)
  n.TS.by.C = rowsum.y(counts, y=y, L=L)
  n.CP.by.TS = rowsum.z(n.TS.by.C, z=z, K=K)
  n.CP = as.integer(rowSums(n.CP.by.TS))
  n.by.G = as.integer(rowSums(counts))
  n.by.C = as.integer(colSums(counts))
  n.by.TS = as.integer(rowsum.y(matrix(n.by.G,ncol=1), y=y, L=L))
  nG.by.TS = as.integer(table(factor(y, 1:L)))
  n.CP.by.G = rowsum.z(counts, z=z, K=K)

  ## Gibbs sampling for each cell
  n.TS.by.CP = t(n.CP.by.TS)
  next.z = cCG.calcGibbsProbZ(m.CP.by.S=m.CP.by.S, n.TS.by.CP=n.TS.by.CP, n.TS.by.C=n.TS.by.C, n.CP=n.CP, n.by.C=n.by.C, z=z, s=s, L=L, K=K, nM=nM, alpha=alpha, beta=beta)
  z.prob = t(next.z$probs)

  ## Gibbs sampling for each gene
  next.y = cCG.calcGibbsProbY(n.CP.by.TS=n.CP.by.TS, n.by.TS=n.by.TS, nG.by.TS=nG.by.TS, n.CP.by.G=n.CP.by.G, n.by.G=n.by.G, y=y, nG=nG, L=L, beta=beta, delta=delta, gamma=gamma)
  y.prob = t(next.y$probs)

  if(!isTRUE(log)) {
    z.prob = normalizeLogProbs(z.prob)
    y.prob = normalizeLogProbs(y.prob)
  }
       
  return(list(z.probability=z.prob, y.probability=y.prob))
}

reorder.celda_CG = function(counts,res){
  # Reorder K
  if(res$K > 2) {
    fm <- factorizeMatrix(counts = counts, celda.mod = res)
    d <- cosineDist(fm$proportions$population.states)
    h <- hclust(d, method = "complete")
    res <- recodeClusterZ(res, from = h$order, to = 1:res$K)
  }  
  
  # Reorder L
  if(res$L > 2) {
    fm <- factorizeMatrix(counts = counts, celda.mod = res)
    cs <- prop.table(t(fm$proportions$population.states), 2)
    d <- cosineDist(cs)
    h <- hclust(d, method = "complete")
    res <- recodeClusterY(res, from = h$order, to = 1:res$L)
  }
  return(res)
}



################################################################################
# celda_CG S3 methods                                                          #
################################################################################
#' finalClusterAssignment for the celda Cell and Gene clustering model
#' @param celda.mod A celda model object of "Celda_CG"
#' @export
finalClusterAssignment.celda_CG = function(celda.mod) {
  return(list(z=celda.mod$z, y=celda.mod$y))
}


#' getK for the celda Cell and Gene clustering model
#' @param celda.mod A celda model object of "Celda_CG"
#' @export
getK.celda_CG = function(celda.mod) {
  return(celda.mod$K)
}


#' getL for the celda Cell and Gene clustering model
#' @param celda.mod A celda model object of "Celda_CG"
#' @export
getL.celda_CG = function(celda.mod) {
  return(celda.mod$L)
}


#' celdaHeatmap for celda Cell and Gene clustering model
#' @param celda.mod A celda model object of "Celda_CG"
#' @param counts A count matrix
#' @param ... extra parameters passed onto renderCeldaHeatmap
#' @export
celdaHeatmap.celda_CG = function(celda.mod, counts, ...) {
  renderCeldaHeatmap(counts, z=celda.mod$z, y=celda.mod$y, ...)
}


#' Visualize the performance of celda_CG models grouped by L and K
#' 
#' Plot the performance of a list of celda_CG models returned 
#' from running celda function. For each number of gene clusters 
#' L (cell clusters K), plot the performance of each number of cell
#' clusters K (gene clusters L).
#' @param celda.list A list of celda_CG objects returned from celda function
#' @param method One of "perplexity" or "loglik"
#' @param title Title for the visualizeModelPerformance
#' @param log Set log to TRUE to visualize the log(perplexity) of Celda_CG objects.
#' @import Rmpfr
#' @export
visualizeModelPerformance.celda_CG = function(celda.list, method="perplexity",
                                                title="Model Performance (All Chains)",
                                                log = FALSE) {
 
  validateKLPlotParameters(celda.list, method)
 
  y.lab = method
  cluster.sizes = unlist(lapply(celda.list$res.list, function(mod) { getK(mod) }))
  log.likelihoods = lapply(celda.list$res.list,
                           function(mod) { completeLogLikelihood(mod) })
  performance.metric = lapply(log.likelihoods, 
                              calculatePerformanceMetric,
                              method, log)
  
  # These methods return Rmpfr numbers that are extremely small and can't be 
  # plotted, so log 'em first
  if (method == "perplexity") {
    if (!(log)) {
      performance.metric = lapply(performance.metric, log)
      performance.metric = methods::new("mpfr", unlist(performance.metric))
    }
    y.lab = paste0("Log(",method,")")
  } 
  
  performance.metric = as.numeric(performance.metric)

  plot.df = data.frame(K=cluster.sizes, L=celda.list$run.params$L,
                       metric=performance.metric)
  
  L.list = sort(unique(celda.list$run.params$L))
  K.list = sort(unique(celda.list$run.params$K))
  #chains = sort(unique(celda.list$run.params$chain))
  
  plots = list()
  nc = round(length(L.list)^.5)
  x.lab.K = "K"
  
  for (i in L.list) {
    plots = c(plots, list(ggplot2::ggplot(subset(plot.df, L==i), 
                          ggplot2::aes(x=K, y=metric, group=K)) + 
                            ggplot2::geom_boxplot(outlier.color=NA, fill=NA) + 
                            ggplot2::geom_point(position=ggplot2::position_jitter(width=0.1, height=0)) +
                            ggplot2::ggtitle(paste0("L = ", i)) + 
                            ggplot2::theme_bw() + 
                            ggplot2::theme(axis.title.x=ggplot2::element_blank(), 
                                           axis.title.y=ggplot2::element_blank()) +
                            ggplot2::scale_y_continuous(limits = c(min(plot.df$metric),
                                                                   max(plot.df$metric)))))
  }
  gridExtra::grid.arrange(grobs = plots, ncol = nc, left = grid::textGrob(y.lab, rot = 90), 
                          top = grid::textGrob(title),
                          bottom = grid::textGrob(x.lab.K))
  
  
  plots = list()
  nc = round(length(K.list)^.5)
  x.lab.L = "L"
  
  for (i in K.list) {
    plots = c(plots, list(ggplot2::ggplot(subset(plot.df, K==i), 
                                          ggplot2::aes(x=L, y=metric, group=L)) +
                            ggplot2::geom_boxplot(outlier.color=NA, fill=NA) + 
                            ggplot2::geom_point(position=ggplot2::position_jitter(width=0.1, height=0)) + 
                            ggplot2::ggtitle(paste0("K = ", i)) + 
                            ggplot2::theme_bw() + 
                            ggplot2::theme(axis.title.x=ggplot2::element_blank(), 
                                           axis.title.y=ggplot2::element_blank()) +
                            ggplot2::scale_y_continuous(limits = c(min(plot.df$metric),
                                                                   max(plot.df$metric)))))
  }
  gridExtra::grid.arrange(grobs = plots, ncol = nc, left = grid::textGrob(y.lab, rot = 90),
                          top = grid::textGrob(title),
                          bottom = grid::textGrob(x.lab.L))
}


#' Render an interactive figure demonstrating the clustering performance of different K/L parameters
#' 
#' Plot the performance of a list of celda_CG models returned 
#' from running the celda function. This is visualized by rendering a boxplot for each chain
#' run for each combination of K/L (cell/gene).
#' 
#' @param celda.list A list of celda_CG objects returned from celda function
#' @param method One of "perplexity" or "loglik", passed through to calculatePerformanceMetric()
#' @param title The plot title
#' @import Rmpfr 
#' @export
renderInteractiveKLPlot = function(celda.list,  method="perplexity", 
                                   title="Model Performance (All Chains)") {
  
  validateKLPlotParameters(celda.list, method)
  
  chain.ks = unlist(lapply(celda.list$res.list, function(mod) { getK(mod) })) # TODO celda_list getter
  chain.ls = unlist(lapply(celda.list$res.list, function(mod) { getL(mod) })) # TODO celda_list getter
  
  
  log.likelihoods = lapply(celda.list$res.list,
                           function(mod) { completeLogLikelihood(mod) })
  performance.metric = lapply(log.likelihoods, 
                              calculatePerformanceMetric,
                              method)
  
  # The performance metric methods return Rmpfr numbers that are extremely small and can't be 
  # plotted via ggplot2, so log 'em first. 
  # TODO: celda_list getter that calculates these metrics.
  if (method %in% c("perplexity")) {
    performance.metric = lapply(performance.metric, log)
    performance.metric = methods::new("mpfr", unlist(performance.metric))
    performance.metric = as.numeric(performance.metric)
  } else {
    performance.metric = as.numeric(performance.metric)
  }
  
  figure.df = data.frame(K=chain.ks, L=chain.ls, metric=performance.metric)
  figure.df$key = paste(figure.df$K, figure.df$L, sep=",")
  
  # Aggregate all the K/L combinations and order them by median, 
  # in order to arrange the boxplots for each K/L combination by descending median:
  grouped.figure.df    = dplyr::group_by(figure.df, key)
  sorted.key.by.median = dplyr::arrange(dplyr::summarize(grouped.figure.df, median(metric)), 
                                        dplyr::desc(`median(metric)`))
  figure.df$key = factor(figure.df$key, levels=sorted.key.by.median$key)
  
  # TODO: Return plot or nah?
  method.label = method 
  if (method %in% c("perplexity")) {
     method.label = paste("log(", method, ")", sep="")
  }
  k.l.plot = ggplot2::ggplot(figure.df, ggplot2::aes(x=key, y=metric, label=key)) +
               ggplot2::geom_boxplot(outlier.color=NA, fill=NA) + 
               ggplot2::geom_point(position=ggplot2::position_jitter(width=0.1, height=0)) +
               ggplot2::xlab("K,L Value") + ggplot2::ylab(method.label) +
               ggplot2::ggtitle(title) +  ggplot2::theme_bw() +
               ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=1)) +
               ggplot2::theme(axis.text.y=ggplot2::element_text(hjust=1))
  plotly::ggplotly(k.l.plot)
}


# Sanity checks for parameters to the model performance plotting functions for
# celda_CG models
validateKLPlotParameters = function(celda.list, method) {
 if (class(celda.list) != "celda_list") {
    stop("celda.list argument must be of class 'celda_list'")
 } else if (celda.list$content.type != "celda_CG") {
    stop("celda.list must be a 'celda.list' of 'celda_CG' objects")
 } else if (!(method %in% c("perplexity", "loglik"))) {
    stop("Invalid method, 'method' has to be either 'perplexity' or 'loglik'")
 } 
}
