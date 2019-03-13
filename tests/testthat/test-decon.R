  # DecontXoneBatch 
  library(celda) 
  context("Testing Deconx")


  Decon.sim = simulateContaminatedMatrix(K=10, delta=c(1,5), seed = 123)
  model_DecontXoneBatch = DecontXoneBatch( Decon.sim$rmat + Decon.sim$cmat, z= Decon.sim$z, max.iter=2, seed=1234567) 
  model_DecontXoneBatch.iter1 = DecontXoneBatch( Decon.sim$rmat + Decon.sim$cmat, z= Decon.sim$z, max.iter=1, seed=1234567)

  model_DecontXoneBatchbg = DecontX( Decon.sim$rmat+Decon.sim$cmat, max.iter=2, seed=1234567) 

  Decon.sim2 = simulateContaminatedMatrix(K=10, delta=c(1,5), seed = 74) 
  batch_DecontX = DecontX( cbind( Decon.sim$rmat+Decon.sim$cmat, Decon.sim2$rmat+Decon.sim2$cmat ) ,z=c( Decon.sim$z, Decon.sim2$z) , batch = rep( 1:2, each = ncol(Decon.sim$rmat) )  , max.iter=2, seed=1234567)  
  batch_DecontX.bg = DecontX( cbind( Decon.sim$rmat+Decon.sim$cmat, Decon.sim2$rmat+Decon.sim2$cmat ) , batch = rep( 1:2, each = ncol(Decon.sim$rmat) )  , max.iter=2, seed=1234567)  


  # simulateContaminatedMatrix
  test_that( desc = "Testing simulateContaminatedMatrix", { 
    expect_equivalent( object= colSums(Decon.sim$rmat) + colSums(Decon.sim$cmat),  expected=Decon.sim$N.by.C)
    expect_equal( object=dim(Decon.sim$phi), expected=dim(Decon.sim$eta) ) 
    expect_equal( typeof(Decon.sim$rmat), "integer")
    expect_equal( typeof(Decon.sim$cmat), "integer")
    
    Decon.sim.SingleDelta = simulateContaminatedMatrix(K=10, delta=1)
    expect_equivalent( object=colSums(Decon.sim.SingleDelta$rmat) + colSums(Decon.sim.SingleDelta$cmat), expected=Decon.sim.SingleDelta$N.by.C)

    Decon.sim.KTooLarge = simulateContaminatedMatrix(K=101, C=10)
    expect_equal( unique(Decon.sim.KTooLarge$z), 1:ncol(Decon.sim.KTooLarge$eta) )

    } )

  # DecontX 
  test_that( desc = "Testing DecontX", {
    expect_equal( ncol( Decon.sim$rmat ) + ncol( Decon.sim2$rmat ) ,  ncol(  batch_DecontX$res.list$est.rmat )  ) 
    expect_equal( length( batch_DecontX$res.list$est.conp) , ncol(  batch_DecontX$res.list$est.rmat )  )  
    expect_equal( batch_DecontX.bg$method, "background" ) 
  } )


  # DecontXoneBatch
  test_that( desc = "Testing DecontXoneBatch", {
    expect_equal( model_DecontXoneBatch$res.list$est.conp  , 1 - colSums(model_DecontXoneBatch$res.list$est.rmat) /  colSums( Decon.sim$rmat + Decon.sim$cmat) )
    expect_equal( model_DecontXoneBatch$res.list$theta, (model_DecontXoneBatch$run.params$delta + colSums(model_DecontXoneBatch$res.list$est.rmat) ) / ( 2*model_DecontXoneBatch$run.params$delta + colSums(Decon.sim$cmat + Decon.sim$rmat) ) )
    expect_error( DecontXoneBatch(counts=Decon.sim$rmat+Decon.sim$cmat, z=Decon.sim$z, beta=-1), "'beta' should be a single positive value.")
    expect_error( DecontXoneBatch(counts=Decon.sim$rmat+Decon.sim$cmat, z=Decon.sim$z, beta=c(1,1) ), "'beta' should be a single positive value.")
    expect_error( DecontXoneBatch(counts=Decon.sim$rmat+Decon.sim$cmat, z=Decon.sim$z, delta=-1), "'delta' should be a single positive value.")
    expect_error( DecontXoneBatch(counts=Decon.sim$rmat+Decon.sim$cmat, z=Decon.sim$z, delta=c(1,1) ), "'delta' should be a single positive value.")
    expect_error( DecontXoneBatch(counts=Decon.sim$rmat+Decon.sim$cmat, z=c(Decon.sim$z, 1) ), "'z' must be of the same length as the number of cells in the 'counts' matrix.")   
    expect_error( DecontXoneBatch(counts=Decon.sim$rmat+Decon.sim$cmat, z=rep(1, ncol(Decon.sim$rmat)) ), "'z' must have at least 2 different values.") 

    counts.NA = Decon.sim$rmat + Decon.sim$cmat
    counts.NA[1,1] = NA
    expect_error( DecontXoneBatch(counts=counts.NA, z=Decon.sim$z), "Missing value in 'counts' matrix.") 
    } )

  test_that( desc = " Testing DecontXoneBatch using background distribution", {
    expect_equal( model_DecontXoneBatchbg$res.list$est.conp, 1- colSums(model_DecontXoneBatchbg$res.list$est.rmat) / Decon.sim$N.by.C  ) 
    } )     


  # logLikelihood
  test_that( desc = "Testing logLikelihood.DecontXoneBatch", {
    z.process = processCellLabels(Decon.sim$z, num.cells=ncol(Decon.sim$rmat) )
    expect_equal( decon.calcLL(counts=Decon.sim$cmat+Decon.sim$rmat, z=z.process  ,  theta=model_DecontXoneBatch$res.list$theta, eta=model_DecontXoneBatch$res.list$est.ConDist, phi=model_DecontXoneBatch$res.list$est.GeneDist ), model_DecontXoneBatch$res.list$logLikelihood[ model_DecontXoneBatch$run.params$iteration  ] )

    cellDist.model.bg = normalizeCounts( model_DecontXoneBatchbg$res.list$est.rmat, normalize="proportion", pseudocount.normalize= model_DecontXoneBatchbg$run.params$beta) 
    bgDist.model.bg = rowSums( Decon.sim$rmat+ Decon.sim$cmat) / sum( Decon.sim$N.by.C) 
    bgDist.model.bg = matrix( rep(bgDist.model.bg, length(Decon.sim$N.by.C)   ), ncol= length(Decon.sim$N.by.C)  )
    expect_equal( bg.calcLL( counts=Decon.sim$cmat+Decon.sim$rmat, theta=model_DecontXoneBatchbg$res.list$theta, cellDist= cellDist.model.bg, bgDist= bgDist.model.bg), model_DecontXoneBatchbg$res.list$logLikelihood[ model_DecontXoneBatchbg$run.params$iteration  ] )
    } )

  # decontamination EM updates
  test_that( desc = "Testing decontamination EM updates", {
    z.process = processCellLabels(Decon.sim$z, num.cells=ncol(Decon.sim$rmat) )
    expect_equal( cD.calcEMDecontamination( counts=Decon.sim$cmat+Decon.sim$rmat, z=z.process, K=length(unique(Decon.sim$z)), theta=model_DecontXoneBatch.iter1$res.list$theta, phi=model_DecontXoneBatch.iter1$res.list$est.GeneDist, eta=model_DecontXoneBatch.iter1$res.list$est.ConDist, beta=model_DecontXoneBatch.iter1$run.params$beta, delta=model_DecontXoneBatch.iter1$run.params$delta)$theta,   model_DecontXoneBatch$res.list$theta )
    } )


