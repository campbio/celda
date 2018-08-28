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
  celdaCG.res = calculatePerplexityWithResampling(celdaCG.sim$counts, celdaCG.res, resample=2)
  expect_equal(is.null(celdaCG.res$perplexity), FALSE)
  expect_is(celdaCG.res, "celda_list")
  expect_is(celdaCG.res, "celda_CG")

  plot.obj = plotGridSearchPerplexity(celdaCG.res)
  expect_is(plot.obj, "ggplot")
})


test_that(desc = "Testing calculatePerplexityWithResampling for celda_C", {
  expect_equal(is.null(celdaC.res$perplexity), TRUE)
  celdaC.res = calculatePerplexityWithResampling(celdaC.sim$counts, celdaC.res, resample=2)
  expect_equal(is.null(celdaC.res$perplexity), FALSE)
  expect_is(celdaC.res, "celda_list")
  expect_is(celdaC.res, "celda_C")

  plot.obj = plotGridSearchPerplexity(celdaC.res)
  expect_is(plot.obj, "ggplot")
})


test_that(desc = "Testing calculatePerplexityWithResampling for celda_G", {
  expect_equal(is.null(celdaG.res$perplexity), TRUE)
  celdaG.res = calculatePerplexityWithResampling(celdaG.sim$counts, celdaG.res, resample=2)
  expect_equal(is.null(celdaG.res$perplexity), FALSE)
  expect_is(celdaG.res, "celda_list")
  expect_is(celdaG.res, "celda_G")

  plot.obj = plotGridSearchPerplexity(celdaG.res)
  expect_is(plot.obj, "ggplot")
})


