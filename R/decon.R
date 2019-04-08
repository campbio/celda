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
        cat("Only", length(unique(z)), "clusters are simulated. Try to increase numebr of cells 'C' if more clusters are needed", "\n")  
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
bg.calcLL = function(counts, global.z, cb.z, phi, eta,theta){
    #ll = sum( t(counts) * log( theta * t(cellDist) + (1-theta) * t(bgDist) + 1e-20) ) 
    ll = sum( t(counts) * log( theta * t(phi)[cb.z,] + (1-theta) * t(eta)[global.z,] +1e-20 )) 
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
cD.calcEMbgDecontamination = function(counts, global.z, cb.z, tr.z, phi, eta, theta, beta, vctr.idf){

    log.Pr = log( t(phi)[ cb.z,] + 1e-20) + log( theta + 1e-20 )  
    log.Pc = log( t(eta)[ global.z,] + 1e-20) + log(1-theta + 1e-20) 

    Pr = exp( log.Pr) / (exp(log.Pr) + exp( log.Pc) ) 
    Pc = 1- Pr 
    delta.v2 = MCMCprecision::fit_dirichlet( matrix( c( Pr, Pc) , ncol = 2 ) )$alpha

    est.rmat = t(Pr) * counts
    phi.unnormalized = colSumByGroup.numeric( est.rmat, cb.z, max(cb.z) ) 
    eta.unnormalized = rowSums(phi.unnormalized) - colSumByGroup.numeric( phi.unnormalized, tr.z, max(tr.z) ) 

    ## Update paramters 
    theta = (colSums(est.rmat) + delta.v2[1] ) / (colSums(counts) + sum(delta.v2)  ) 
    phi = normalizeCounts(phi.unnormalized, normalize="proportion",  pseudocount.normalize =beta ) 
    eta = normalizeCounts(eta.unnormalized, normalize="proportion",  pseudocount.normalize =beta ) 

    # tf-idf 
    K = length(unique(global.z)) 
    for(k in 1:K) { 
        eta[,k] = eta[,k] %*% diag(vctr.idf[,k]) 
        eta[,k] = eta[,k] / sum(eta[,k]) 
    }

    return( list("est.rmat"=est.rmat, "theta"=theta, "phi"=phi, "eta"=eta, "delta"=delta.v2 ) ) 
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

        ## Initialize cell label
        initial.label = decontx.initializeZ(counts=counts) 
        global.z = initial.label$global.z
        cb.z = initial.label$cb.z
        tr.z = initial.label$tr.z


        ## Initialization
        delta.init = delta 
        theta = runif( n = nC, min = 0.1, max = 0.5 ) 
        #theta = rbeta( n = nC, shape1 = delta.init,  shape2 = delta.init ) 
        est.rmat = t( t(counts) *theta) 

        phi = colSumByGroup.numeric( est.rmat, cb.z, max(cb.z)  ) 
        eta = rowSums(phi) - colSumByGroup.numeric( phi,  tr.z,  max(tr.z)  ) 
        phi = normalizeCounts( phi, normalize="proportion", pseudocount.normalize =beta ) 
        eta = normalizeCounts( eta, normalize="proportion", pseudocount.normalize = beta) 
        ll =c()

        #tf-idf
        K = length(unique(global.z)) 
        vctr.idf = matrix(NA, ncol=K, nrow=nG)
        for (k in 1:K) {
            # what if there is only one cell in cluster k, might be a bug here
            vctr.idf[,k] = log10( (sum( global.z==k ) + 1) / apply(counts[,global.z==k], 1, FUN=function(v) {return(sum(v>0)+1)}) ) + 1
        }
        # then during the iteration: matrix-rize each vector from this matrix and right multiply with the eta, L1 normalize the result as the new weighted eta. The calculation of log-likelihood should also use tf-idf weighted eta.
        for( k in 1:K) {
            # tf.idf weighted eta
            eta[,k] = eta[,k] %*% diag(vctr.idf[,k])  # tf.idf weight to adjust eta
            eta[,k] = eta[,k] / sum(eta[,k])  # L1 normalization
        }


        ll.round = bg.calcLL(counts = counts, global.z = global.z, cb.z = cb.z, phi = phi, eta = eta, theta = theta ) 

        ## EM updates
        while (iter <=  max.iter &  num.iter.without.improvement <= stop.iter  ) { 


            next.decon = cD.calcEMbgDecontamination(counts=counts, global.z=global.z, cb.z=cb.z, tr.z=tr.z,  phi=phi, eta=eta, theta=theta, beta=beta, vctr.idf=vctr.idf)  

            theta = next.decon$theta
            phi = next.decon$phi
            eta = next.decon$eta 
            delta = next.decon$delta


            ## Calculate log-likelihood
            ll.temp = bg.calcLL(counts = counts, global.z = global.z, cb.z = cb.z, phi = phi, eta = eta, theta = theta) 
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
    if( length(proportionPrior) > 1 | proportionPrior <= 0   )  {
        stop("'delta' should be a single positive value.")
    }
    if( length( distributionPrior) > 1 | distributionPrior <=0  ) {
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




## Initialization of cell labels for DecontX when they are not given
decontx.initializeZ = function( counts, K=10, min.cell=3, seed=1111) {
    nC = ncol(counts) 
    if( nC<100 )  K = ceiling(sqrt( nC )) 

    global.z = initialize.splitZ(counts, K=K, K.subcluster=NULL, alpha=1, beta=1, min.cell = 3, seed=seed) 
    global.K = max( global.z) 

    local.z = rep(NA, nC) 
    for( k in 1:global.K) { 
        local.counts = counts[, global.z == k ] 
        local.K = min( K, ceiling( sqrt(ncol(local.counts))   )   )
        local.z[ global.z ==k ] = initialize.splitZ( local.counts, K=local.K, K.subcluster=NULL, alpha=1, beta=1, min.cell = 3, seed=seed)    
    }

    cb.z  = interaction( global.z, local.z, lex.order=TRUE, drop=TRUE)   # combined z label  
    tr.z = as.integer( sub("\\..*", "", levels(cb.z), perl=TRUE) )  # transitional z label
    cb.z = as.integer( plyr::mapvalues( cb.z, from=levels(cb.z), to=1:length(levels(cb.z)))  )


    return( list( "global.z"=global.z, "local.z"=local.z, "tr.z"=tr.z, "cb.z"=cb.z  )  )
} 
