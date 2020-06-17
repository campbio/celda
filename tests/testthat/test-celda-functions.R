# library(celda)
# context("Testing miscellaneous celda functions")
#
# data(celdaCGGridSearchRes)
# celdaCGSim <- simulateCells("celda_CG", K = 5, L = 10)
# modelCG <- celda_CG(
#     counts = celdaCGSim$counts,
#     sampleLabel = celdaCGSim$sampleLabel,
#     K = celdaCGSim$K,
#     L = celdaCGSim$L,
#     algorithm = "EM",
#     verbose = FALSE,
#     nchains = 1)
#
# factorized <- factorizeMatrix(celdaMod = modelCG, counts = celdaCGSim$counts)
#
# test_that(desc = "Testing compareCountMatrix with numeric matrix input", {
#     # Case from GitHub issue #137
#     counts <- celdaCGSim$counts
#     storage.mode(counts) <- "numeric"
#     expect_true(compareCountMatrix(counts, modelCG, errorOnMismatch = TRUE))
# })
#
# test_that(desc = "Testing appendCeldaList", {
#     expect_error(
#         appendCeldaList(celdaCGGridSearchRes, matrix(0)),
#         "Both parameters to appendCeldaList must be of class celdaList."
#     )
#     modifiedEgList <- celdaCGGridSearchRes
#     modifiedEgList@countChecksum <- "abcd12345"
#     expect_warning(
#         appendCeldaList(celdaCGGridSearchRes, modifiedEgList),
#         paste0("Provided lists have different countChecksums and may have",
#         " been generated from different count matrices. Using checksum",
#         " from first list..."))
#     expect_equal(length(celdaCGGridSearchRes@resList) * 2,
#         length(appendCeldaList(celdaCGGridSearchRes,
#             celdaCGGridSearchRes)@resList))
# })
