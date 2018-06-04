# celda_C.R
library(celda)
library(Rtsne)
context("Testing celda Visualization Experiments")

load("../celdaCGsim.rda")
load("../celdaCG.rda")
model_CG = getModel(celdaCG.res, K = 5, L=3)[[1]]
factorized <- factorizeMatrix(celda.mod = model_CG, counts = celdaCG.sim$counts)
counts.matrix <- celdaCG.sim$counts

test_that(desc = "Testing that absoluteProbabilityHeatmap runs",{
  plot.obj = absoluteProbabilityHeatmap(counts=counts.matrix, celda.mod=model_CG)
  expect_true(!is.null(plot.obj))
})
