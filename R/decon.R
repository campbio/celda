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


#' This function calculates the log-likelihood
#' 
#' @param omat Numeric/Integer matrix. Observed count matrix, rows represent features and columns represent cells
#' @param z Integer vector. Cell population labels
#' @param phi Numeric matrix. Rows represent features and columns represent cell populations
#' @param eta Numeric matrix. Rows represent features and columns represent cell populations 
#' @param theta Numeric vector. Proportion of truely expressed transcripts
decon.calcLL = function(omat, z,  phi, eta, theta){
  #ll = sum( t(omat) * log( (1-conP )*geneDist[z,] + conP * conDist[z, ] + 1e-20 ) )  # when dist_mat are K x G matrices
  ll = sum( t(omat) * log( theta * phi[z,] + (1-theta) * eta[z,]  + 1e-20  ))   
  return(ll)
}


#' This function updates decontamination 
#' 
#' @param omat Numeric/Integer matrix. Observed count matrix, rows represent features and columns represent cells
#' @param phi Numeric matrix. Rows represent features and columns represent cell populations
#' @param eta Numeric matrix. Rows represent features and columns represent cell populations
#' @param theta Numeric vector. Proportion of truely expressed transctripts
#' @param z Integer vector. Cell population labels
#' @param K Integer. Number of cell populations
#' @param beta Numeric. Concentration parameter for Phi
#' @param delta Numeric. Concentration parameter for Theta
cD.calcEMDecontamination = function(omat, phi, eta, theta,  z, K, beta, delta ) {

   log.Pr = log(  t(phi)[z,] + 1e-20) + log(theta + 1e-20 )
   log.Pc = log(  t(eta)[z,] + 1e-20) + log(1-theta + 1e-20)  

   Pr = exp(log.Pr) / (  exp(log.Pr) + exp(log.Pc) )   
   
   est.rmat = t(Pr) * omat 
   rn.G.by.K = colSumByGroup.numeric(est.rmat, z, K)  
   cn.G.by.K = rowSums(rn.G.by.K) - rn.G.by.K 

   #update parameters 
   theta = (colSums(est.rmat) + delta) / (colSums(omat) + 2*delta)  
   phi  = normalizeCounts(rn.G.by.K, pseudocount.normalize =beta ) 
   eta  = normalizeCounts(cn.G.by.K, pseudocount.normalize = beta) 

   return( list("est.rmat"=est.rmat,  "theta"=theta, "phi"=phi, "eta"=eta ) ) 
}


#' This function updates decontamination 
#' 
#' @param omat Numeric/Integer matrix. Observed count matrix, rows represent features and columns represent cells. 
#' @param z Integer vector. Cell population labels. 
#' @param max.iter Integer. Maximum iterations of EM algorithm. Default to be 200 
#' @param beta Numeric. Concentration parameter for Phi. Default to be 1e-6
#' @param delta Numeric / Numeric vector. Concentration parameter for Theta. Default to be 10 
#' @export 
DeconX = function(omat, z, max.iter=200, beta=1e-6, delta=10 ) {

  nG = nrow(omat)
  nC = ncol(omat)
  K =  length(unique(z)) 

  # initialization
  theta  = runif(nC, min = 0.1, max = 0.9)  
  est.rmat = t (t(omat) * theta )       
  phi =   colSumByGroup.numeric(est.rmat, z, K)
  eta =   rowSums(phi) - phi 
  phi = normalizeCounts(phi, pseudocount.normalize =beta )
  eta = normalizeCounts(eta, pseudocount.normalize = beta)
  ll = c()

  # EM updates
  for (iteration in 1:max.iter) {


    next.decon = cD.calcEMDecontamination(omat=omat, phi=phi, eta=eta, theta=theta, z=z, K=K,  beta=beta, delta=delta) 

    theta = next.decon$theta 
    phi = next.decon$phi
    eta = next.decon$eta 

    # Calculate log-likelihood
    ll.temp = decon.calcLL(omat=omat, z=z, phi = t(phi), eta = t(eta), theta=theta )
    ll  = c(ll, ll.temp)
  }

  res.conp  = 1- colSums(next.decon$est.rmat) / colSums(omat)

  run.params = list("beta" = beta, "delta" = delta, "iteration"=max.iter)
  res.list = list("logLikelihood" = ll, "est.rmat"=next.decon$est.rmat , "est.conp"= res.conp, "theta"=theta ,  "est.GeneDist"=phi,  "est.ConDist"=eta )
  return(list("run.params"=run.params, "res.list"=res.list ))
}
