# celda_CG.R
library(celda)
context("Testing celda_CG")

celdacg <- simulateCells.celda_CG(K=3, L=5)
celdaCG.res <- celda(counts=celdacg$counts, nchains=1, K = 2:4, L = 4:6, 
                     ncores=1, model="celda_CG")

test_that("visualizeModelPerformance Returns a Plot",{
  expect_equal(TRUE, all(!is.na(visualizeModelPerformance(celdaCG.res, celdacg$counts))))
})
