# celda_C.R
library(celda)
context("Testing celda_C")

celdac <- simulateCells.celda_C(K=10)
celdaC.res <- celda(counts=celdac$counts, model="celda_C", nchains=1, K=10)

test_that("CheckingVisualizeModelPerformace",{
	expect_equal(TRUE, all(!is.na(visualizeModelPerformance(celdaC.res))))
	})
#test_that("finalClusterAssignment.celda_C",{
#  expect_equal(celdaC.res$res.list[[1]]$z, finalClusterAssignment(celdaC.res$res.list[[1]]))
#})


