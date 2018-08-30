# celda_CG.R
library(celda)
library(testthat)
context("Testing model performance functions")

#Loading pre-made simulatedcells/celda objects
load("../celdaCGsim.rda")
load("../celdaCG.rda")
load("../celdaCsim.rda")
load("../celdaC.rda")
load("../celdaGsim.rda")
load("../celdaG.rda")

test_that(desc = "Testing calculatePerplexityWithResampling for celda_CG", {
  expect_equal(is.null(celdaCG.res$perplexity), TRUE)
  expect_error(plotGridSearchPerplexity(celdaCG.res))

  celdaCG.res = calculatePerplexityWithResampling(celdaCG.sim$counts, celdaCG.res, resample=2)
  expect_equal(is.null(celdaCG.res$perplexity), FALSE)
  expect_is(celdaCG.res, "celda_list")
  expect_is(celdaCG.res, "celda_CG")
  expect_error(calculatePerplexityWithResampling(celdaCG.sim$counts, celdaCG.res, resample="2"))
  expect_error(calculatePerplexityWithResampling(celdaCG.sim$counts, "celdaCG.res", resample=2))
  
  plot.obj = plotGridSearchPerplexity(celdaCG.res)
  expect_is(plot.obj, "ggplot")

  model_CG = filterCeldaList(celdaCG.res, K = 5, L = 3)[[1]]
  expect_error(celdaCG.res <- calculatePerplexityWithResampling(celdaCG.sim$counts, model_CG, resample=2))
  expect_error(celdaCG.res <- calculatePerplexityWithResampling(celdaCG.sim$counts, celdaCG.res, resample='a'))

  celdaC.res = calculatePerplexityWithResampling(celdaC.sim$counts, celdaC.res, resample=2)
  expect_error(plotGridSearchPerplexity.celda_CG(celdaC.res))
  
  res <- calculatePerplexity.celda_CG(celdaCG.sim$counts, model_CG)
  res2 <- calculatePerplexity.celda_CG(celdaCG.sim$counts, model_CG, new.counts = celdaCG.sim$counts + 1)
  
  expect_error(res <- calculatePerplexity.celda_CG(celdaCG.sim$counts, model_CG, new.counts = celdaCG.sim$counts[-1,]))  
})


test_that(desc = "Testing calculatePerplexityWithResampling for celda_C", {
  expect_equal(is.null(celdaC.res$perplexity), TRUE)
  expect_error(plotGridSearchPerplexity(celdaC.res))

  celdaC.res = calculatePerplexityWithResampling(celdaC.sim$counts, celdaC.res, resample=2)
  expect_equal(is.null(celdaC.res$perplexity), FALSE)
  expect_is(celdaC.res, "celda_list")
  expect_is(celdaC.res, "celda_C")
  expect_error(calculatePerplexityWithResampling(celdaC.sim$counts, celdaC.res, resample="2"))
  expect_error(calculatePerplexityWithResampling(celdaC.sim$counts, "celdaC.res", resample=2))
  
  plot.obj = plotGridSearchPerplexity(celdaC.res)
  expect_is(plot.obj, "ggplot")

  celdaG.res = calculatePerplexityWithResampling(celdaG.sim$counts, celdaG.res, resample=2)
  expect_error(plotGridSearchPerplexity.celda_C(celdaG.res))
  
  model_C <- filterCeldaList(celdaC.res, K = 5)[[1]]
  res <- calculatePerplexity.celda_C(celdaC.sim$counts, model_C)
  res2 <- calculatePerplexity.celda_C(celdaC.sim$counts, model_C, new.counts = celdaC.sim$counts + 1)
  
  expect_error(res <- calculatePerplexity.celda_C(celdaC.sim$counts, model_C, new.counts = celdaC.sim$counts[-1,]))
})


test_that(desc = "Testing calculatePerplexityWithResampling for celda_G", {
  expect_equal(is.null(celdaG.res$perplexity), TRUE)
  expect_error(plotGridSearchPerplexity(celdaG.res))

  celdaG.res = calculatePerplexityWithResampling(celdaG.sim$counts, celdaG.res, resample=2)
  expect_equal(is.null(celdaG.res$perplexity), FALSE)
  expect_is(celdaG.res, "celda_list")
  expect_is(celdaG.res, "celda_G")
  expect_error(calculatePerplexityWithResampling(celdaG.sim$counts, celdaG.res, resample="2"))
  expect_error(calculatePerplexityWithResampling(celdaG.sim$counts, "celdaG.res", resample=2))
  
  plot.obj = plotGridSearchPerplexity(celdaG.res)
  expect_is(plot.obj, "ggplot")

  celdaC.res = calculatePerplexityWithResampling(celdaC.sim$counts, celdaC.res, resample=2)
  expect_error(plotGridSearchPerplexity.celda_G(celdaC.res))

  model_G <- filterCeldaList(celdaG.res, L = 5)[[1]]
  res <- calculatePerplexity.celda_G(celdaG.sim$counts, model_G)
  res2 <- calculatePerplexity.celda_G(celdaG.sim$counts, model_G, new.counts = celdaG.sim$counts + 1)
  
  expect_error(res <- calculatePerplexity.celda_G(celdaG.sim$counts, model_G, new.counts = celdaG.sim$counts[-1,]))
})



