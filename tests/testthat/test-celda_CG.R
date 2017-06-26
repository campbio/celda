# celda_CG.R
library(celda)
context("Testing celda_CG")

celdacg <- simulateCells.celda_CG(K = 3, L = 5, G = 200)
celdaCG.res <- celda_CG(counts=celdacg$counts, nchains=1, K = 3, L = 5, ncores=1)

test_that("CheckingVisualizeModelPerformace",{
        expect_equal(TRUE, all(!is.na(visualize_model_performance(celdaCG.res))))
        })

