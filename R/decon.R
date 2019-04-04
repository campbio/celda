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
#' @return A list object containing the real expression matrix and contamination
#' expression matrix as well as other parameters used in the simulation.
#' @examples 
#' contamination.sim =  simulateContaminatedMatrix(  K=3,  delta=c(1,9)) 
#' contamination.sim =  simulateContaminatedMatrix(  K=3,  delta = 1) 
#' @export
simulateContaminatedMatrix = function(C=300, G=100, K=3, N.Range=c(500,1000), beta = 0.5, delta=c(1,2),  seed=12345) {
  
    setSeed(seed)

    if(length(delta) == 1 ) { 
        cp.byC = rbeta(n=C, shape1=delta, shape2=delta) 
    } else {
        cp.byC = rbeta(n=C, shape1=delta[1], shape2=delta[2] ) 
    } 
 
    z = sample(1:K, size=C, replace=TRUE) 
    if( length(unique(z)) < K) {
        warning("Only ", length(unique(z)), " clusters are simulated. Try to increase numebr of cells 'C' if more clusters are needed", "\n")  
        K = length(unique(z)) 
        z = plyr::mapvalues(z, unique(z), 1:length(unique(z)) )
    }


    N.byC = sample(min(N.Range):max(N.Range), size=C, replace=TRUE  ) 
    cN.byC = sapply(1:C, function(i)  rbinom(n=1, size=N.byC[i], p=cp.byC[i] )  ) 
    rN.byC = N.byC - cN.byC 
    
    phi = rdirichlet(K, rep(beta, G)) 
 
 
    ## sample real expressed count matrix 
    cell.rmat = sapply(1:C, function(i) stats::rmultinom(1, size=rN.byC[i], prob=phi[z[i],] )  )

    rownames(cell.rmat) = paste0("Gene_", 1:G) 
    colnames(cell.rmat) = paste0("Cell_", 1:C) 
    
    
    ## sample contamination count matrix 
    n.G.by.K = rowSums(cell.rmat) - colSumByGroup(cell.rmat, group=z, K=K) 
    eta = normalizeCounts(counts=n.G.by.K, normalize="proportion") 

    cell.cmat = sapply(1:C, function(i) stats::rmultinom(1, size=cN.byC[i], prob=eta[,z[i]] )  )
    cell.omat = cell.rmat + cell.cmat 

    rownames(cell.omat) = paste0("Gene_", 1:G) 
    colnames(cell.omat) = paste0("Cell_", 1:C)

    return(list("nativeCounts"=cell.rmat, "observedCounts"=cell.omat, "N.by.C"=N.byC, "z"=z, "eta"=eta , "phi"=t(phi)  ) ) 

}


# This function calculates the log-likelihood
# 
# counts Numeric/Integer matrix. Observed count matrix, rows represent features and columns represent cells
# z Integer vector. Cell population labels
# phi Numeric matrix. Rows represent features and columns represent cell populations
# eta Numeric matrix. Rows represent features and columns represent cell populations 
# theta Numeric vector. Proportion of truely expressed transcripts
decon.calcLL = function(counts, z,  phi, eta, theta){
    #ll = sum( t(counts) * log( (1-conP )*geneDist[z,] + conP * conDist[z, ] + 1e-20 ) )  # when dist_mat are K x G matrices
    ll = sum( t(counts) * log( theta * t(phi)[z,] + (1-theta) * t(eta)[z,]  + 1e-20  ))   
    return(ll)
}

# This function calculates the log-likelihood of background distribution decontamination
#
# bgDist Numeric matrix. Rows represent feature and columns are the times that the background-distribution has been replicated. 
bg.calcLL = function(counts, cellDist, bgDist, theta){
    ll = sum( t(counts) * log( theta * t(cellDist) + (1-theta) * t(bgDist) + 1e-20) ) 
    return(ll) 
}


# This function updates decontamination 
# 
#  phi Numeric matrix. Rows represent features and columns represent cell populations
#  eta Numeric matrix. Rows represent features and columns represent cell populations
#  theta Numeric vector. Proportion of truely expressed transctripts
cD.calcEMDecontamination = function(counts, phi, eta, theta,  z, K, beta, delta ) {

    ## Notes: use fix-point iteration to update prior for theta, no need to feed delta anymore
    log.Pr = log(  t(phi)[z,] + 1e-20) + log(theta + 1e-20 )
    log.Pc = log(  t(eta)[z,] + 1e-20) + log(1-theta + 1e-20)  

    Pr = exp(log.Pr) / (  exp(log.Pr) + exp(log.Pc) )   
    Pc = 1 - Pr 
    delta.v2 = MCMCprecision::fit_dirichlet( matrix( c( Pr, Pc) , ncol = 2 ) )$alpha 

    est.rmat = t(Pr) * counts 
    rn.G.by.K = colSumByGroup.numeric(est.rmat, z, K)  
    cn.G.by.K = rowSums(rn.G.by.K) - rn.G.by.K 

    ## Update parameters 
    theta = (colSums(est.rmat) + delta.v2[1]) / (colSums(counts) + sum(delta.v2)  )  
    phi  = normalizeCounts(rn.G.by.K, normalize="proportion",  pseudocount.normalize =beta ) 
    eta  = normalizeCounts(cn.G.by.K, normalize="proportion",  pseudocount.normalize = beta) 

    return( list("est.rmat"=est.rmat,  "theta"=theta, "phi"=phi, "eta"=eta, "delta"=delta.v2 ) ) 
}


# This function updates decontamination using background distribution 
# 
cD.calcEMbgDecontamination = function(counts, cellDist, bgDist, theta, beta){

    meanN.by.C = apply(counts, 2, mean) 
    log.Pr = log( t(cellDist) + 1e-20) + log( theta + 1e-20) # + log( t(counts) / meanN.by.C )   # better when without panelty 
    log.Pc = log( t(bgDist)  + 1e-20) + log( 1-theta + 2e-20) 

    Pr = exp( log.Pr) / (exp(log.Pr) + exp( log.Pc) ) 
    Pc = 1- Pr 
    delta.v2 = MCMCprecision::fit_dirichlet( matrix( c( Pr, Pc) , ncol = 2 ) )$alpha

    est.rmat = t(Pr) * counts

    ## Update paramters 
    theta = (colSums(est.rmat) + delta.v2[1] ) / (colSums(counts) + sum(delta.v2)  ) 
    cellDist = normalizeCounts(est.rmat, normalize="proportion", pseudocount.normalize =beta) 

    return( list("est.rmat"=est.rmat, "theta"=theta, "cellDist"=cellDist, "delta"=delta.v2 ) ) 
}


#' This function updates decontamination on dataset with multiple batches 
#' @param counts Numeric/Integer matrix. Observed count matrix, rows represent features and columns represent cells. 
#' @param z Integer vector. Cell population labels. Default NULL 
#' @param batch Integer vector. Cell batch labels. Default NULL
#' @param max.iter Integer. Maximum iterations of EM algorithm. Default to be 200 
#' @param beta Numeric. Concentration parameter for Phi. Default to be 1e-6
#' @param delta Numeric. Symmetric concentration parameter for Theta. Default to be 10 
#' @param logfile Character. Messages will be redirected to a file named `logfile`. If NULL, messages will be printed to stdout.  Default NULL
#' @param verbose Logical. Whether to print log messages. Default TRUE
#' @param seed Integer. Passed to set.seed(). Default to be 1234567. If NULL, no calls to `set.seed()` are made.
#' @return A list object which contains the decontaminated count matrix and related parameters
#' @examples 
#' decon.c = DecontX( counts = contamination.sim$rmat + contamination.sim$cmat, z=contamination.sim$z, max.iter=3)
#' decon.bg = DecontX( counts=contamination.sim$rmat + contamination.sim$cmat, max.iter=3 ) 
#' @export 
DecontX = function( counts, z=NULL, batch=NULL,  max.iter=200, beta=1e-6, delta=10, logfile=NULL, verbose=TRUE, seed=1234567 ) {
  

    if ( !is.null(batch) ) { 

        ## Set result lists upfront for all cells from different batches 
        logLikelihood = c() 
        est.rmat = matrix(NA, ncol = ncol(counts), nrow = nrow(counts ), dimnames = list( rownames( counts), colnames( counts ) )      ) 
        theta = rep(NA, ncol(counts ) ) 
        est.conp = rep(NA, ncol(counts) ) 

        batch.index = unique(batch) 

        for ( BATCH in batch.index ) {  
            counts.BATCH = counts[, batch == BATCH  ]
            if( !is.null(z) ) { z.BATCH = z[ batch == BATCH ]  }  else { z.BATCH = z } 
            res.BATCH = DecontXoneBatch( counts = counts.BATCH, z = z.BATCH, batch = BATCH  ,max.iter = max.iter, beta = beta, delta = delta, logfile=logfile, verbose=verbose, seed=seed ) 

            est.rmat[, batch == BATCH ] = res.BATCH$res.list$est.nativeCounts
            est.conp[ batch == BATCH ] = res.BATCH$res.list$est.conp 
            theta[  batch == BATCH ] = res.BATCH$res.list$theta 

            if( is.null(logLikelihood)  )   { 
                logLikelihood = res.BATCH$res.list$logLikelihood  
            } else { 
                logLikelihood = addLogLikelihood( logLikelihood, res.BATCH$res.list$logLikelihood ) 
            }
        } 
      
        run.params = res.BATCH$run.params 
        method = res.BATCH$method 
        res.list = list( "logLikelihood"=logLikelihood, "est.nativeCounts"=est.rmat,"est.conp"= est.conp,  "theta"=theta ) 

        return( list("run.params"=run.params, "res.list"=res.list, "method"=method )  ) 
    } 

    return( DecontXoneBatch( counts = counts, z = z, max.iter = max.iter, beta = beta, delta = delta, logfile=logfile, verbose=verbose, seed=seed )  ) 
}




# This function updates decontamination for one batch  
#
DecontXoneBatch = function(counts, z=NULL, batch=NULL, max.iter=200, beta=1e-6, delta=10, logfile=NULL, verbose=TRUE, seed=1234567 ) {

    checkCounts.decon(counts)
    checkParameters.decon(proportionPrior=delta, distributionPrior=beta) 

    nG = nrow(counts)
    nC = ncol(counts)
    K =  length(unique(z)) 
 
    if ( is.null(z) ) {
        decon.method = "background"
    } else {
        decon.method = "clustering"
        z = processCellLabels(z, num.cells= nC) 
    }

    ## Set seed for initialization 
    setSeed(seed)


    iter = 1L 
    num.iter.without.improvement = 0L
    stop.iter = 3L 


    logMessages("----------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose) 
    logMessages("Start DecontX. Decontamination", logfile=logfile, append=TRUE, verbose=verbose) 
    if ( !is.null(batch) ) {  logMessages("batch: ",  batch, logfile=logfile, append=TRUE, verbose=verbose)    }
    logMessages("----------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose) 
    start.time = Sys.time()

    if( decon.method == "clustering") {

        ## Initialization
        delta.init = delta 
        #theta  = stats::runif(nC, min = 0.1, max = 0.5)  
        theta = rbeta( n = nC, shape1 = delta.init, shape2 = delta.init ) 
        est.rmat = t (t(counts) * theta )       
        phi =   colSumByGroup.numeric(est.rmat, z, K)
        eta =   rowSums(phi) - phi 
        phi = normalizeCounts(phi, normalize="proportion", pseudocount.normalize =beta )
        eta = normalizeCounts(eta, normalize="proportion", pseudocount.normalize = beta)
        ll = c()
  

        ll.round =  decon.calcLL(counts=counts, z=z, phi = phi, eta = eta, theta=theta ) 


        ## EM updates
        while (iter <=  max.iter &  num.iter.without.improvement <= stop.iter  ) {

            next.decon = cD.calcEMDecontamination(counts=counts, phi=phi, eta=eta, theta=theta, z=z, K=K,  beta=beta, delta=delta) 

            theta = next.decon$theta 
            phi = next.decon$phi
            eta = next.decon$eta 
            delta = next.decon$delta

            ## Calculate log-likelihood
            ll.temp = decon.calcLL(counts=counts, z=z, phi = phi, eta = eta, theta=theta )
            ll  = c(ll, ll.temp)
            ll.round = c( ll.round, round( ll.temp, 2) )  

            if ( round( ll.temp, 2) > ll.round[ iter ]  |  iter == 1   )    {
                num.iter.without.improvement = 1L 
            }   else { 
                num.iter.without.improvement = num.iter.without.improvement + 1L 
            } 
            iter = iter + 1L 

        }

    }

    if ( decon.method == "background") {

        ## Initialization
        delta.init = delta 
        theta = rbeta( n = nC, shape1 = delta.init,  shape2 = delta.init ) 
        est.rmat = t( t(counts) *theta) 
        bgDist = rowSums( counts ) / sum( counts) 
        bgDist = matrix( rep( bgDist, nC), ncol=nC) 
        cellDist = normalizeCounts( est.rmat, normalize="proportion", pseudocount.normalize=beta) 
        ll =c()


        ll.round = bg.calcLL(counts=counts, cellDist=cellDist, bgDist=bgDist, theta=theta)  

        ## EM updates
        while (iter <=  max.iter &  num.iter.without.improvement <= stop.iter  ) { 


            next.decon = cD.calcEMbgDecontamination(counts=counts, cellDist=cellDist, bgDist=bgDist, theta=theta, beta=beta)  

            theta = next.decon$theta
            cellDist = next.decon$cellDist
            delta = next.decon$delta 


            ## Calculate log-likelihood
            ll.temp = bg.calcLL(counts=counts, cellDist=cellDist, bgDist=bgDist, theta=theta) 
            ll  = c(ll, ll.temp) 
            ll.round = c( ll.round, round( ll.temp, 2) )

            if ( round( ll.temp, 2) > ll.round[ iter ]  |  iter == 1   )    {
                num.iter.without.improvement = 1L
            }   else {
                num.iter.without.improvement = num.iter.without.improvement + 1L
            }
            iter = iter + 1L
        }

    }

    res.conp  = 1- colSums(next.decon$est.rmat) / colSums(counts)

    end.time = Sys.time() 
    logMessages("----------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose) 
    logMessages("Completed DecontX. Total time:", format(difftime(end.time, start.time)), logfile=logfile, append=TRUE, verbose=verbose) 
    if ( !is.null(batch) ) {  logMessages("batch: ",  batch, logfile=logfile, append=TRUE, verbose=verbose)    }
    logMessages("----------------------------------------------------------------------", logfile=logfile, append=TRUE, verbose=verbose) 

    run.params = list("beta.init"=beta, "delta.init"=delta.init, "iteration"=iter-1L, "seed"=seed)

    res.list = list("logLikelihood" = ll, "est.nativeCounts"=next.decon$est.rmat , "est.conp"= res.conp, "theta"=theta , "delta"=delta)
    #if( decon.method=="clustering" ) {
    #    posterior.params = list( "est.GeneDist"=phi,  "est.ConDist"=eta  ) 
    #    res.list = append( res.list , posterior.params ) 
    #}
  
    return(list("run.params"=run.params, "res.list"=res.list, "method"=decon.method  ))
}



## Make sure provided parameters are the right type and value range 
checkParameters.decon = function(proportionPrior, distributionPrior) {
    if( length(proportionPrior) > 1 | any(proportionPrior <= 0)   )  {
        stop("'delta' should be a single positive value.")
    }
    if( length( distributionPrior) > 1 | any(distributionPrior <=0)  ) {
        stop("'beta' should be a single positive value.") 
    }
}


## Make sure provided count matrix is the right type
checkCounts.decon = function(counts) {
    if ( sum(is.na(counts)) >0   ) {
        stop("Missing value in 'counts' matrix.") 
    }
}


## Make sure provided cell labels are the right type
processCellLabels = function(z, num.cells) {
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


## Add two (veried-length) vectors of logLikelihood 
addLogLikelihood = function( ll.a, ll.b   ) { 
    length.a = length(ll.a ) 
    length.b = length(ll.b) 

    if( length.a >= length.b  )  { 
        ll.b = c( ll.b, rep( ll.b[length.b] , length.a - length.b )  ) 
        ll = ll.a + ll.b 
    } else { 
        ll.a = c( ll.a, rep( ll.a[length.a], length.b - length.a ) ) 
        ll = ll.a + ll.b
    } 
  
    return( ll ) 
}



 
