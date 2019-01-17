#' This function generates a list containing two count matrices -- one for real expression, the other one for contamination, as well as other parameters 
#' used in the simulation which can be useful for running decontamination
#'
#' @param C Integer. Number of cells to be simulated. Default to be 300
#' @param G Integer. Number of genes to be simulated. Default to be 100
#' @param K Integer. Number of cell populations to be simulated. Default to be 3
#' @param N.Range Integer vector. A vector of length 2 that specifies the lower and upper bounds of the number of counts generated for each cell. Default to be c(500, 1000)
#' @param beta Numeric. Concentration parameter for Phi. Default to be 0.5
#' @param delta Numeric / Numeric vector. Concentration parameter for Theta. If input as a single numeric value, symmetric values for beta distribution are specified; if input as a vector of lenght 2, the two values will be the shape1 and shape2 paramters of the beta distribution respectively
#' @param seed Integer. Passed to set.seed(). Default to be 12345. If NULL, no calls to `set.seed()` are made.
#' @examples 
#' contamination.sim =  simulateObservedMatrix(  K=3,  delta=c(1,9)) 
#' contamination.sim =  simulateObservedMatrix(  K=3,  delta = 1) 
#' @export
simulateObservedMatrix = function(C=300, G=100, K=3, N.Range=c(500,1000), beta = 0.5, delta=c(1,2),  seed=12345) {
  
    if (!is.null(seed)) {
      set.seed(seed) 
    }

    if(length(delta)==1) { 
	    cp.byC = rbeta(n=C, shape1=delta, shape2=delta) 
    } else {
            cp.byC = rbeta(n=C, shape1=delta[1], shape2=delta[2] ) 
    } 
 
    z = sample(1:K, size=C, replace=TRUE) 
    if( length(unique(z)) < K) {
        cat("Only", length(unique(z)), "clusters are simulated. Try to increase numebr of cells 'C' if more clusters are needed", "\n")  
        K = length(unique(z)) 
        z = plyr::mapvalues(z, unique(z), 1:length(unique(z)) )
    }


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

    return(list("rmat"=cell.rmat, "cmat"=cell.cmat, "N.by.C"=N.byC, "z"=z, "eta"=eta , "phi"=t(phi)  ) ) 

}


# This function calculates the log-likelihood
# 
# omat Numeric/Integer matrix. Observed count matrix, rows represent features and columns represent cells
# z Integer vector. Cell population labels
# phi Numeric matrix. Rows represent features and columns represent cell populations
# eta Numeric matrix. Rows represent features and columns represent cell populations 
# theta Numeric vector. Proportion of truely expressed transcripts
decon.calcLL = function(omat, z,  phi, eta, theta){
  #ll = sum( t(omat) * log( (1-conP )*geneDist[z,] + conP * conDist[z, ] + 1e-20 ) )  # when dist_mat are K x G matrices
  ll = sum( t(omat) * log( theta * t(phi)[z,] + (1-theta) * t(eta)[z,]  + 1e-20  ))   
  return(ll)
}

# This function calculates the log-likelihood of background distribution decontamination
#
# bgDist Numeric matrix. Rows represent feature and columns are the times that the background-distribution has been replicated. 
bg.calcLL = function(omat, cellDist, bgDist, theta){
  ll = sum( t(omat) * log( theta * t(cellDist) + (1-theta) * t(bgDist) + 1e-20) ) 
  return(ll) 
}


# This function updates decontamination 
# 
#  phi Numeric matrix. Rows represent features and columns represent cell populations
#  eta Numeric matrix. Rows represent features and columns represent cell populations
#  theta Numeric vector. Proportion of truely expressed transctripts
cD.calcEMDecontamination = function(omat, phi, eta, theta,  z, K, beta, delta ) {

   log.Pr = log(  t(phi)[z,] + 1e-20) + log(theta + 1e-20 )
   log.Pc = log(  t(eta)[z,] + 1e-20) + log(1-theta + 1e-20)  

   Pr = exp(log.Pr) / (  exp(log.Pr) + exp(log.Pc) )   
   
   est.rmat = t(Pr) * omat 
   rn.G.by.K = colSumByGroup.numeric(est.rmat, z, K)  
   cn.G.by.K = rowSums(rn.G.by.K) - rn.G.by.K 

   #update parameters 
   theta = (colSums(est.rmat) + delta) / (colSums(omat) + 2*delta)  
   phi  = normalizeCounts(rn.G.by.K, normalize="proportion",  pseudocount.normalize =beta ) 
   eta  = normalizeCounts(cn.G.by.K, normalize="proportion",  pseudocount.normalize = beta) 

   return( list("est.rmat"=est.rmat,  "theta"=theta, "phi"=phi, "eta"=eta ) ) 
}


# This function updates decontamination using background distribution 
# 
cD.calcEMbgDecontamination = function(omat, cellDist, bgDist, theta, beta, delta){

  meanN.by.C = apply(omat, 2, mean) 
  log.Pr = log( t(cellDist) + 1e-20) + log( theta + 1e-20) # + log( t(omat) / meanN.by.C )   # better when without panelty 
  log.Pc = log( t(bgDist)  + 1e-20) + log( 1-theta + 2e-20) 

  Pr = exp( log.Pr) / (exp(log.Pr) + exp( log.Pc) ) 

  est.rmat = t(Pr) * omat

  # Update paramters 
  theta = (colSums(est.rmat) + delta) / (colSums(omat) + 2*delta) 
  cellDist = normalizeCounts(est.rmat, normalize="proportion", pseudocount.normalize =beta) 

  return( list("est.rmat"=est.rmat, "theta"=theta, "cellDist"=cellDist) ) 
}




#' This function updates decontamination 
#' 
#' @param omat Numeric/Integer matrix. Observed count matrix, rows represent features and columns represent cells. 
#' @param z Integer vector. Cell population labels. Default NULL 
#' @param max.iter Integer. Maximum iterations of EM algorithm. Default to be 200 
#' @param beta Numeric. Concentration parameter for Phi. Default to be 1e-6
#' @param delta Numeric. Symmetric concentration parameter for Theta. Default to be 10 
#' @param logfile Character. Messages will be redirected to a file named `logfile`. If NULL, messages will be printed to stdout.  Default NULL
#' @param verbose Logical. Whether to print log messages. Default TRUE
#' @param seed Integer. Passed to set.seed(). Default to be 1234567. If NULL, no calls to `set.seed()` are made.
#' @examples 
#' decon.c = DecontX( omat = contamination.sim$rmat + contamination.sim$cmat, z=contamination.sim$z, max.iter=3)
#' decon.bg = DecontX( omat=contamination.sim$rmat + contamination.sim$cmat, max.iter=3 ) 
#' @export 
DecontX = function(omat, z=NULL, max.iter=200, beta=1e-6, delta=10, logfile=NULL, verbose=TRUE, seed=1234567 ) {

  checkCounts.decon(omat)
  checkParameters.decon(proportion.prior=delta, distribution.prior=beta) 

  nG = nrow(omat)
  nC = ncol(omat)
  K =  length(unique(z)) 
 
  if ( is.null(z) ) {
      decon.method = "background"
  } else {
      decon.method = "clustering"
      z = processCellLabels(z, num.cells= nC) 
  }

  logMessages("----------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose) 
  logMessages("Start DecontX. Decontamination", logfile=logfile, append=TRUE, verbose=verbose) 
  logMessages("----------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose) 
  start.time = Sys.time()

  if( decon.method == "clustering") {

  # initialization
  if (!is.null(seed)) {
    set.seed(seed)
  }
  theta  = runif(nC, min = 0.1, max = 0.5)  
  est.rmat = t (t(omat) * theta )       
  phi =   colSumByGroup.numeric(est.rmat, z, K)
  eta =   rowSums(phi) - phi 
  phi = normalizeCounts(phi, normalize="proportion", pseudocount.normalize =beta )
  eta = normalizeCounts(eta, normalize="proportion", pseudocount.normalize = beta)
  ll = c()

  # EM updates
  for (iteration in 1:max.iter) {


    next.decon = cD.calcEMDecontamination(omat=omat, phi=phi, eta=eta, theta=theta, z=z, K=K,  beta=beta, delta=delta) 

    theta = next.decon$theta 
    phi = next.decon$phi
    eta = next.decon$eta 

    # Calculate log-likelihood
    ll.temp = decon.calcLL(omat=omat, z=z, phi = phi, eta = eta, theta=theta )
    ll  = c(ll, ll.temp)
  }

  }

  if ( decon.method == "background") {

  # Initialization
  if (!is.null(seed)) {
    set.seed(seed) 
  }
  theta = runif( nC, min =0.1, max=0.5) 
  est.rmat = t( t(omat) *theta) 
  bgDist = rowSums( omat ) / sum( omat) 
  bgDist = matrix( rep( bgDist, nC), ncol=nC) 
  cellDist = normalizeCounts( est.rmat, normalize="proportion", pseudocount.normalize=beta) 
  ll =c()

  # EM updates
  for (iteration in 1:max.iter) {


    next.decon = cD.calcEMbgDecontamination(omat=omat, cellDist=cellDist, bgDist=bgDist, theta=theta, beta=beta, delta=delta)  

    theta = next.decon$theta
    cellDist = next.decon$cellDist

    # Calculate log-likelihood
    ll.temp = bg.calcLL(omat=omat, cellDist=cellDist, bgDist=bgDist, theta=theta) 
    ll  = c(ll, ll.temp) 
  }

  }

  res.conp  = 1- colSums(next.decon$est.rmat) / colSums(omat)


  end.time = Sys.time() 
  logMessages("----------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose) 
  logMessages("Completed DecontX. Total time:", format(difftime(end.time, start.time)), logfile=logfile, append=TRUE, verbose=verbose) 
  logMessages("----------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose) 

  run.params = list("beta" = beta, "delta" = delta, "iteration"=max.iter)

  res.list = list("logLikelihood" = ll, "est.rmat"=next.decon$est.rmat , "est.conp"= res.conp, "theta"=theta   )
  if( decon.method=="clustering" ) {
    posterior.params = list( "est.GeneDist"=phi,  "est.ConDist"=eta  ) 
    res.list = append( res.list , posterior.params ) 
  }
  
  return(list("run.params"=run.params, "res.list"=res.list, "method"=decon.method  ))
}


### checking 
# Make sure provided parameters are the right type and value range 
checkParameters.decon = function(proportion.prior, distribution.prior) {
  if( length(proportion.prior) > 1 | proportion.prior <= 0   )  {
    stop("'delta' should be a single positive value.")
  }
  if( length( distribution.prior) > 1 | distribution.prior <=0  ) {
    stop("'beta' should be a single positive value.") 
  }
}

# Make sure provided rmat is the right type
checkCounts.decon = function(omat) {
  if ( sum(is.na(omat)) >0   ) {
    stop("Missing value in 'omat' matrix.") 
  }
}

# Make sure provided cell labels are the right type
processCellLabels = function(z,  num.cells) {
  if ( length(z) != num.cells )  {
    stop("'z' must be of the same length as the number of cells in the 'counts' matrix.") 
  }
  if( length(unique(z)) <2 ) {
    stop("'z' must have at least 2 different values.") # Even though everything runs smoothly when length(unique(z)) == 1, result is not trustful  
  }
  if( !is.factor(z) ) {
    z = plyr::mapvalues(z, unique(z), 1:length(unique(z)) )
    z = as.factor(z) 
  }
  return(z)
}
