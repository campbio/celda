# celda_CG.R
library(celda)
library(testthat)
library(Rtsne)
context("Testing model performance functions")

#Loading pre-made simulatedcells/celda objects
load("../celdaCGsim.rda")
load("../celdaCG.rda")
model_CG = getModel(celdaCG.res, K = 5, L = 3)[[1]]
factorized <- factorizeMatrix(model_CG, celdaCG.sim$counts)
counts.matrix <- celdaCG.sim$counts


test_that(desc = "Testing calculatePerplexityWithResampling", {
  perplexity.results = calculatePerplexityWithResampling(celdaCG.res, celdaCG.sim$counts, resample=2)
  expect_named(perplexity.results, c("perplexity.info", "plot"))
  expect_is(perplexity.results$plot, "ggplot")
})


test_that(desc = "Testing visualizePerplexityByKL", {
  perplexity.results = calculatePerplexityWithResampling(celdaCG.res, celdaCG.sim$counts, resample=2)
  plot.obj = visualizePerplexityByKL(perplexity.results$perplexity.info)
  expect_is(plot.obj, "ggplot")
})
