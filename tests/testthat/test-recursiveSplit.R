library(celda)
context("Testing recursiveSplits")

# celdaCGSim <- simulateCells("celda_CG", K = 5, L = 10)
# 
# test_that(desc = "Testing recursiveSplitModule", {
#     moduleSplit <- recursiveSplitModule(
#         celdaCGSim$counts,
#         initialL = 8,
#         maxL = 15,
#         tempK = 5,
#         zInit = NULL)
#     expect_true(is(moduleSplit, "celdaList"))
#     moduleSplit <- recursiveSplitModule(
#         celdaCGSim$counts,
#         initialL = 8,
#         maxL = 15,
#         tempK = NULL,
#         zInit = celdaCGSim$z)
#     expect_true(is(moduleSplit, "celdaList"))
#     moduleSplit <- recursiveSplitModule(
#         celdaCGSim$counts,
#         initialL = 8,
#         maxL = 15,
#         tempK = NULL,
#         zInit = NULL)
#     expect_true(is(moduleSplit, "celdaList"))
#     plotGridSearchPerplexity(moduleSplit)
# })
# 
# test_that(desc = "Testing recursiveSplitCell", {
#     cellSplit <- recursiveSplitCell(
#         celdaCGSim$counts,
#         initialK = 3,
#         maxK = 8,
#         tempL = 20,
#         yInit = NULL)
#     expect_true(is(cellSplit, "celdaList"))
#     cellSplit <- recursiveSplitCell(
#         celdaCGSim$counts,
#         initialK = 3,
#         maxK = 8,
#         tempL = NULL,
#         yInit = celdaCGSim$y)
#     expect_true(is(cellSplit, "celdaList"))
#     cellSplit <- recursiveSplitCell(
#         celdaCGSim$counts,
#         initialK = 3,
#         maxK = 8,
#         tempL = NULL,
#         yInit = NULL)
#     expect_true(is(cellSplit, "celdaList"))
#     plotGridSearchPerplexity(cellSplit)
# })
