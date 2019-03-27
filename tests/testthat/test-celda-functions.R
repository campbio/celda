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

test_that(desc = "Testing appendCeldaList", {
  eg.list = celda::celda.CG.grid.search.res
  expect_error(appendCeldaList(eg.list, matrix(0)),
               "Both parameters to appendCeldaList must be of class celdaList.")
  modified.eg.list = eg.list
  modified.eg.list@count.checksum = "abcd12345"
  expect_warning(appendCeldaList(eg.list, modified.eg.list),
                 "Provided lists have different count.checksums and may have been generated from different count matrices. Using checksum from first list...")
  expect_equal(length(eg.list@res.list)*2,
               length(appendCeldaList(eg.list, eg.list)@res.list))
})
