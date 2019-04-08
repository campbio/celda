library(celda)
context("Testing miscellaneous celda functions")

celdaCGSim <- simulateCells("celda_CG", K = 5, L = 10)
modelCG <- celda_CG(
    counts = celdaCGSim$counts,
    sample.label = celdaCGSim$sample.label,
    K = celdaCGSim$K,
    L = celdaCGSim$L,
    algorithm = "EM",
    verbose = FALSE
)
factorized <-
    factorizeMatrix(celdaMod = modelCG, counts = celdaCGSim$counts)

test_that(desc = "Testing compareCountMatrix with numeric matrix input", {
    # Case from GitHub issue #137
    counts <- celdaCGSim$counts
    storage.mode(counts) <- "numeric"
    expect_true(compareCountMatrix(counts, modelCG, errorOnMismatch = TRUE))
})

test_that(desc = "Testing appendCeldaList", {
    egList <- celda::celdaCGGridSearchRes
    expect_error(
        appendCeldaList(egList, matrix(0)),
        "Both parameters to appendCeldaList must be of class celdaList."
    )
    modifiedEgList <- egList
    modifiedEgList@count.checksum <- "abcd12345"
    expect_warning(
        appendCeldaList(egList, modifiedEgList),
        paste0("Provided lists have different count.checksums and may have",
        "been generated from different count matrices. Using checksum",
        " from first list..."))
    expect_equal(length(egList@resList) * 2,
        length(appendCeldaList(egList, egList)@resList))
})
