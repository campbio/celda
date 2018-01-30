# celda_C.R
library(celda)
context("Testing celda_C")

celdac <- simulateCells.celda_C(K=10)
celdaC.res <- celda(counts=celdac$counts, model="celda_C", nchains=1, K=10)

test_that("visualizeModelPerformance Returns a Plot",{
	expect_equal(TRUE, all(!is.na(visualizeModelPerformance(celdaC.res, celdac$counts))))
})
