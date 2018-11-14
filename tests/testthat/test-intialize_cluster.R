test_that(desc = "Test 'split' initialization for all models", {
  sim.res = simulateCells(model="celda_CG")
  model_CG <- celda_CG(sim.res$counts, sim.res$sample.label, K=sim.res$K, L=sim.res$L, 
                       initialize="split", max.iter=5, split.on.last=FALSE, 
                       split.on.iter=-1)
  expect_true(is(model_CG, "celda_CG"))
  model_G <- celda_G(sim.res$counts, L=sim.res$L, initialize="split", max.iter=5, split.on.last=FALSE, split.on.iter=-1)
  expect_true(is(model_G, "celda_G"))
  model_C <- celda_C(sim.res$counts, sim.res$sample.label, K=sim.res$K, initialize="split", max.iter=5, split.on.last=FALSE, split.on.iter=-1)
  expect_true(is(model_C, "celda_C"))
})


test_that(desc = "Testing initialize.cluster for random initialization", {
  
  ## Completely random
  z = initialize.cluster(10, 100)
  expect_true(length(z) == 100 & length(unique(z) == 10))
  expect_error(z <- initialize.cluster(100, 10))    
  
  ## With all values initialized
  init.z = rep(1:10, each=10)
  z = initialize.cluster(10, 100, initial=init.z)
  expect_true(all(init.z == z))
  expect_error(z <- initialize.cluster(10, 100, initial=init.z[1:99]))
  expect_error(z <- initialize.cluster(11, 100, initial=init.z))
  expect_error(z <- initialize.cluster(10, 99, initial=init.z))    

  
  ## With only a few values initialized
  fixed.z = rep(NA, 100)
  fixed.z[1:10] = 1
  z = initialize.cluster(10, 100, fixed=fixed.z)
  expect_true(all(z[1:10] == 1) & length(z) == 100 & length(unique(z)) == 10)
  expect_error(z <- initialize.cluster(10, 100, fixed=fixed.z[1:99]))    
  fixed.z[1:10] = 11
  expect_error(z <- initialize.cluster(10, 100, fixed=fixed.z))    
})



test_that(desc = "Testing additional recursive.splitY error checking", {
  L = 10
  sim.res = simulateCells(model="celda_G", L=L)
  z = initialize.cluster(10, ncol(sim.res$counts))
  y = recursive.splitY(sim.res$counts, sim.res$L, beta=1, delta=1, gamma=1, z=z, K=15, K.subclusters=15, min.feature=3, max.cells=50, seed=12345)
  expect_true(length(y) == nrow(sim.res$counts) & length(unique(y)) == L)

  z[z == 1] = 2
  z[1] = 1
  y = recursive.splitY(sim.res$counts, sim.res$L, beta=1, delta=1, gamma=1, z=z, K=15, K.subclusters=15, min.feature=3, max.cells=50, seed=12345)
  expect_true(length(y) == nrow(sim.res$counts) & length(unique(y)) == L)  
  y = recursive.splitY(sim.res$counts, sim.res$L, beta=1, delta=1, gamma=1, z=NULL, K=NULL, K.subclusters=NULL, min.feature=3, max.cells=1000, seed=12345)
  expect_true(length(y) == nrow(sim.res$counts) & length(unique(y)) == L)  
})
