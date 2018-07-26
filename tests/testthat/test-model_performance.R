# celda_CG.R
library(celda)
library(testthat)
context("Testing model performance functions")

#Loading pre-made simulatedcells/celda objects
load("../celdaCGsim.rda")
load("../celdaCG.rda")

test_that(desc = "Testing calculatePerplexityWithResampling", {
  perplexity.results = calculatePerplexityWithResampling(celdaCG.res, celdaCG.sim$counts, resample=2, validate.counts=FALSE)
  expect_named(perplexity.results, c("perplexity.info", "plot"))
  expect_is(perplexity.results$plot, "ggplot")
})


test_that(desc = "Testing visualizePerplexityByKL", {
  perplexity.results = calculatePerplexityWithResampling(celdaCG.res, celdaCG.sim$counts, resample=2, validate.counts=FALSE)
  plot.obj = visualizePerplexity(perplexity.results$perplexity.info)
  expect_is(plot.obj, "ggplot")
})
