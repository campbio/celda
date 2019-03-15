test_that(desc = "Test 'random' initialization for all models", {
  sim.res = simulateCells(model="celda_CG")
  model_CG <- celda_CG(sim.res$counts, sim.res$sample.label, K=sim.res$K, L=sim.res$L, 
                       z.initialize="random", y.initialize="random", max.iter=5, split.on.last=FALSE, 
                       split.on.iter=-1)
  expect_true(is(model_CG, "celda_CG"))
  model_G <- celda_G(sim.res$counts, L=sim.res$L, y.initialize="random", max.iter=5, split.on.last=FALSE, split.on.iter=-1)
  expect_true(is(model_G, "celda_G"))
  model_C <- celda_C(sim.res$counts, sim.res$sample.label, K=sim.res$K, z.initialize="random", max.iter=5, split.on.last=FALSE, split.on.iter=-1)
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


