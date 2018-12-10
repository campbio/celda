library(celda)
context("Testing miscellaneous celda functions")

celdaCG.sim = simulateCells("celda_CG", K=5, L=10)
model_CG = celda_CG(counts=celdaCG.sim$counts, sample.label=celdaCG.sim$sample.label, K=celdaCG.sim$K, L=celdaCG.sim$L, algorithm="EM", verbose=FALSE)
factorized <- factorizeMatrix(celda.mod = model_CG, counts = celdaCG.sim$counts)

test_that(desc = "Testing compareCountMatrix with numeric matrix input", {
  # Case from GitHub issue #137
  counts = celdaCG.sim$counts
  storage.mode(counts) = "numeric"
  expect_true(compareCountMatrix(counts, model_CG, error.on.mismatch = TRUE))
})

