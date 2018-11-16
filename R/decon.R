#' This function generates a list containing two count matrices -- one for real expression, the other one for contamination, as well as other parameters 
#' used in the simulation which can be useful for running decontamination
#'
#' @param C Integer. Number of cells to be simulated. Default to be 300
#' @param G Integer. Number of genes to be simulated. Default to be 100
#' @param K Integer. Number of cell populations to be simulated. Default to be 3
#' @param N.Range Integer vector. A vector of length 2 that specifies the lower and upper bounds of the number of counts generated for each cell. Default to be c(500, 1000)
#' @param beta Numeric. Concentration parameter for Phi. Default to be 0.5
#' @param delta Numeric / Numeric vector. Concentration parameter for Theta. If input as a single numeric value, symmetric values for beta distribution are specified; if input as a vector of lenght 2, the two values will be the shape1 and shape2 paramters of the beta distribution respectively
#' @param seed Integer. Passed to set.seed(). Default to be 12345
#' @export
simulateObservedMatrix = function(C=300, G=100, K=3, N.Range=c(500,1000), beta = 0.5, delta=c(1,2),  seed=12345) {
  
    set.seed(seed) 

    if(length(delta)==1) { 
	    cp.byC = rbeta(n=C, shape1=delta, shape2=delta) 
    } else {
            cp.byC = rbeta(n=C, shape1=delta[1], shape2=delta[2] ) 
    } 
    
    z = sample(1:K, size=C, replace=TRUE) 
    N.byC = sample(min(N.Range):max(N.Range), size=C , replace=TRUE  ) 
    cN.byC = sapply(1:C, function(i)  rbinom(n=1, size=N.byC[i], p=cp.byC[i] )  ) 
    rN.byC = N.byC - cN.byC 
    
    phi = rdirichlet(K, rep(beta, G)) 
    
    
    # sample real expressed count matrix 
    cell.rmat = sapply(1:C, function(i) stats::rmultinom(1, size=rN.byC[i], prob=phi[z[i],] )  )

    rownames(cell.rmat) = paste0("Gene_", 1:G) 
    colnames(cell.rmat) = paste0("Cell_", 1:C) 
    
    
    # sample contamination count matrix 
    n.G.by.K = rowSums(cell.rmat) - colSumByGroup(cell.rmat, group=z, K=K) 
    eta = normalizeCounts(counts=n.G.by.K, normalize="proportion") 

    cell.cmat = sapply(1:C, function(i) stats::rmultinom(1, size=cN.byC[i], prob=eta[,z[i]] )  )
    rownames(cell.cmat) = paste0("Gene_", 1:G) 
    colnames(cell.cmat) = paste0("Cell_", 1:C)

    return(list("rmat"=cell.rmat, "cmat"=cell.cmat, "N.by.C"=N.byC, "z"=z, "eta"=eta , "phi"=t(phi) ) ) 

}




