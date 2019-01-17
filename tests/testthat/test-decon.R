  # DecontX 
  library(celda) 
  context("Testing Deconx")


  Decon.sim = simulateObservedMatrix(K=10, delta=c(1,5))
  model_DecontX = DecontX( Decon.sim$rmat + Decon.sim$cmat, z= Decon.sim$z, max.iter=2, seed=1234567) 
  model_DecontX.iter1 = DecontX( Decon.sim$rmat + Decon.sim$cmat, z= Decon.sim$z, max.iter=1, seed=1234567)

  model_DecontXbg = DecontX( Decon.sim$rmat+Decon.sim$cmat, max.iter=2, seed=1234567) 
  model_DecontXbg.iter1 = DecontX( Decon.sim$rmat+Decon.sim$cmat, max.iter=1, seed=1234567) 

  # simulateObservedMatrix
  test_that( desc = "Testing simulateObservedMatrix", { 
    expect_equivalent( object= colSums(Decon.sim$rmat) + colSums(Decon.sim$cmat),  expected=Decon.sim$N.by.C)
    expect_equal( object=dim(Decon.sim$phi), expected=dim(Decon.sim$eta) ) 
    expect_equal( typeof(Decon.sim$rmat), "integer")
    expect_equal( typeof(Decon.sim$cmat), "integer")
    
    Decon.sim.SingleDelta = simulateObservedMatrix(K=10, delta=1)
    expect_equivalent( object=colSums(Decon.sim.SingleDelta$rmat) + colSums(Decon.sim.SingleDelta$cmat), expected=Decon.sim.SingleDelta$N.by.C)

    Decon.sim.KTooLarge = simulateObservedMatrix(K=101, C=10)
    expect_equal( unique(Decon.sim.KTooLarge$z), 1:ncol(Decon.sim.KTooLarge$eta) )

    } )

  # DecontX
  test_that( desc = "Testing DecontX", {
    expect_equal( model_DecontX$res.list$est.conp  , 1 - colSums(model_DecontX$res.list$est.rmat) /  colSums( Decon.sim$rmat + Decon.sim$cmat) )
    expect_equal( dim(model_DecontX$res.list$est.GeneDist),  dim( model_DecontX$res.list$est.ConDist) ) 
    expect_equal( model_DecontX$res.list$theta, (model_DecontX$run.params$delta + colSums(model_DecontX$res.list$est.rmat) ) / ( 2*model_DecontX$run.params$delta + colSums(Decon.sim$cmat + Decon.sim$rmat) ) )
    expect_error( DecontX(omat=Decon.sim$rmat+Decon.sim$cmat, z=Decon.sim$z, beta=-1), "'beta' should be a single positive value.")
    expect_error( DecontX(omat=Decon.sim$rmat+Decon.sim$cmat, z=Decon.sim$z, beta=c(1,1) ), "'beta' should be a single positive value.")
    expect_error( DecontX(omat=Decon.sim$rmat+Decon.sim$cmat, z=Decon.sim$z, delta=-1), "'delta' should be a single positive value.")
    expect_error( DecontX(omat=Decon.sim$rmat+Decon.sim$cmat, z=Decon.sim$z, delta=c(1,1) ), "'delta' should be a single positive value.")
    expect_error( DecontX(omat=Decon.sim$rmat+Decon.sim$cmat, z=c(Decon.sim$z, 1) ), "'z' must be of the same length as the number of cells in the 'counts' matrix.")   
    expect_error( DecontX(omat=Decon.sim$rmat+Decon.sim$cmat, z=rep(1, ncol(Decon.sim$rmat)) ), "'z' must have at least 2 different values.") 

    omat.NA = Decon.sim$rmat + Decon.sim$cmat
    omat.NA[1,1] = NA
    expect_error( DecontX(omat=omat.NA, z=Decon.sim$z), "Missing value in 'omat' matrix.") 
    } )

  test_that( desc = " Testing DecontX using background distribution", {
    expect_equal( model_DecontXbg$res.list$est.conp, 1- colSums(model_DecontXbg$res.list$est.rmat) / Decon.sim$N.by.C  ) 
    } )     


  # logLikelihood
  test_that( desc = "Testing logLikelihood.DecontX", {
    z.process = processCellLabels(Decon.sim$z, num.cells=ncol(Decon.sim$rmat) )
    expect_equal( decon.calcLL(omat=Decon.sim$cmat+Decon.sim$rmat, z=z.process  ,  theta=model_DecontX$res.list$theta, eta=model_DecontX$res.list$est.ConDist, phi=model_DecontX$res.list$est.GeneDist ), model_DecontX$res.list$logLikelihood[ model_DecontX$run.params$iteration  ] )

    cellDist.model.bg = normalizeCounts( model_DecontXbg$res.list$est.rmat, normalize="proportion", pseudocount.normalize= model_DecontXbg$run.params$beta) 
    bgDist.model.bg = rowSums( Decon.sim$rmat+ Decon.sim$cmat) / sum( Decon.sim$N.by.C) 
    bgDist.model.bg = matrix( rep(bgDist.model.bg, length(Decon.sim$N.by.C)   ), ncol= length(Decon.sim$N.by.C)  )
    expect_equal( bg.calcLL( omat=Decon.sim$cmat+Decon.sim$rmat, theta=model_DecontXbg$res.list$theta, cellDist= cellDist.model.bg, bgDist= bgDist.model.bg), model_DecontXbg$res.list$logLikelihood[ model_DecontXbg$run.params$iteration  ] )
    } )

  # decontamination EM updates
  test_that( desc = "Testing decontamination EM updates", {
    z.process = processCellLabels(Decon.sim$z, num.cells=ncol(Decon.sim$rmat) )
    expect_equal( cD.calcEMDecontamination( omat=Decon.sim$cmat+Decon.sim$rmat, z=z.process, K=length(unique(Decon.sim$z)), theta=model_DecontX.iter1$res.list$theta, phi=model_DecontX.iter1$res.list$est.GeneDist, eta=model_DecontX.iter1$res.list$est.ConDist, beta=model_DecontX.iter1$run.params$beta, delta=model_DecontX.iter1$run.params$delta)$theta,   model_DecontX$res.list$theta )
    } )


