# celda_CG.R
library(celda)
context("Testing celda_CG")

celdacg <- simulateCells.celda_CG(K=10,L=10)
celdaCG.res <- celda(counts=celdacg$counts, model="celda_CG", nchains=1, K=10, L=10)

test_that("CheckingVisualizeModelPerformace",{
        expect_equal(TRUE, all(!is.na(visualize_model_performance(celdaCG.res))))
        })

