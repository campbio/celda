# with_seed
library(celda)
context("Testing seed setting behavior")

set.seed(123)
celdaCGSim1 <- simulateCells("celda_CG", K = 5, L = 10, seed = 1234)
celdaCGSim2 <- simulateCells("celda_CG", K = 5, L = 10, seed = NULL)
celdaCGSim3 <- simulateCells("celda_CG", K = 5, L = 10, seed = 123)
celdaCGSim4 <- simulateCells("celda_CG", K = 5, L = 10, seed = 1234)
celdaCGSim5 <- simulateCells("celda_CG", K = 5, L = 10, seed = 1234)
celdaCGSim6 <- simulateCells("celda_CG", K = 5, L = 10, seed = 12345)
celdaCGSim7 <- simulateCells("celda_CG", K = 5, L = 10, seed = NULL)

set.seed(123)
celdaCGSim8 <- simulateCells("celda_CG", K = 5, L = 10, seed = NULL)

test_that(desc = "Testing seed setting behavior in count matrix simulation", {
    expect_equal(celdaCGSim1, celdaCGSim4)
    expect_equal(celdaCGSim1, celdaCGSim5)
    #expect_equal(celdaCGSim2, celdaCGSim3)
    expect_equal(celdaCGSim2, celdaCGSim8)

    # expect_false(isTRUE(all.equal(celdaCGSim1, celdaCGSim2)))
    # expect_false(isTRUE(all.equal(celdaCGSim1, celdaCGSim3)))
    # expect_false(isTRUE(all.equal(celdaCGSim1, celdaCGSim6)))
    # expect_false(isTRUE(all.equal(celdaCGSim2, celdaCGSim7)))
})
