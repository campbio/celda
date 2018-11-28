  # DeconX 
  library(celda) 
  context("Testing Deconx")


  Decon.sim = simulateObservedMatrix(K=10)
  model_DeconX = DeconX( Decon.sim$rmat + Decon.sim$cmat, z= Decon.sim$z, max.iter=2, seed=1234567) 
  model_DeconX.iter1 = DeconX( Decon.sim$rmat + Decon.sim$cmat, z= Decon.sim$z, max.iter=1, seed=1234567)

  # simulateObservedMatrix
  test_that( desc = "Testing simulateObservedMatrix", { 
    expect_equivalent( object= colSums(Decon.sim$rmat) + colSums(Decon.sim$cmat),  expected=Decon.sim$N.by.C)
    expect_equal( object=dim(Decon.sim$phi), expected=dim(Decon.sim$eta) ) 
    expect_equal( typeof(Decon.sim$rmat), "integer")
    expect_equal( typeof(Decon.sim$cmat), "integer")
    } )

  # DeconX
  test_that( desc = "Testing DeconX", {
    expect_equal( model_DeconX$res.list$est.conp  , 1 - colSums(model_DeconX$res.list$est.rmat) /  colSums( Decon.sim$rmat + Decon.sim$cmat) )
    expect_equal( dim(model_DeconX$res.list$est.GeneDist),  dim( model_DeconX$res.list$est.ConDist) ) 
    expect_equal( model_DeconX$res.list$theta, (model_DeconX$run.params$delta + colSums(model_DeconX$res.list$est.rmat) ) / ( 2*model_DeconX$run.params$delta + colSums(Decon.sim$cmat + Decon.sim$rmat) ) )
    } )


  # logLikelihood
  test_that( desc = "Testing logLikelihood.DeconX", {
    z.process = processCellLabels(Decon.sim$z, num.cells=ncol(Decon.sim$rmat) )
    expect_equal( decon.calcLL(omat=Decon.sim$cmat+Decon.sim$rmat, z=z.process  ,  theta=model_DeconX$res.list$theta, eta=t(model_DeconX$res.list$est.ConDist), phi=t(model_DeconX$res.list$est.GeneDist) ), model_DeconX$res.list$logLikelihood[ model_DeconX$run.params$iteration  ] )
    } )

  # decontamination EM updates
  test_that( desc = "Testing decontamination EM updates", {
    z.process = processCellLabels(Decon.sim$z, num.cells=ncol(Decon.sim$rmat) )
    expect_equal( cD.calcEMDecontamination( omat=Decon.sim$cmat+Decon.sim$rmat, z=z.process, K=length(unique(Decon.sim$z)), theta=model_DeconX.iter1$res.list$theta, phi=model_DeconX.iter1$res.list$est.GeneDist, eta=model_DeconX.iter1$res.list$est.ConDist, beta=model_DeconX.iter1$run.params$beta, delta=model_DeconX.iter1$run.params$delta)$theta,   model_DeconX$res.list$theta )
    } )
