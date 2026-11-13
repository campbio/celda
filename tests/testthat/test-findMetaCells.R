# findMetaCells
library(celda)
context("Testing findMetaCells")

sceSim <- simulateCells("celda_CG", S = 2, K = 4, L = 6, G = 200,
    CRange = c(60, 80), seed = 3)
sampleLabel <- as.character(sceSim$celda_sample_label)
logfile <- tempfile(fileext = ".log")
sceMeta <- findMetaCells(sceSim,
    initK = 2,
    maxK = 8,
    tempL = 10,
    minCellMeta = 10,
    sampleLabel = sampleLabel,
    logfile = logfile)

test_that(desc = "Testing findMetaCells output", {
    expect_true(is(sceMeta, "SingleCellExperiment"))
    expect_setequal(colnames(sceMeta), colnames(sceSim))
    expect_true("metacell" %in% colnames(colData(sceMeta)))
    expect_true(all(grepl("^Sample_[12]_K[0-9]+$", sceMeta$metacell)))
    expect_true(all(table(sceMeta$metacell) >= 10))
})

test_that(desc = "Testing findMetaCells keeps meta-cells within a sample", {
    origSample <- sampleLabel[match(colnames(sceMeta), colnames(sceSim))]
    expect_equal(as.character(sceMeta$celda_sample_label), origSample)
    expect_equal(sub("_K[0-9]+$", "", sceMeta$metacell), origSample)
    samplesPerMeta <- tapply(origSample, sceMeta$metacell,
        function(x) length(unique(x)))
    expect_true(all(samplesPerMeta == 1))
})

test_that(desc = "Testing findMetaCells is reproducible with a seed", {
    sceMeta2 <- findMetaCells(sceSim,
        initK = 2,
        maxK = 8,
        tempL = 10,
        minCellMeta = 10,
        sampleLabel = sampleLabel,
        verbose = FALSE)
    expect_identical(sceMeta2$metacell, sceMeta$metacell)
})

test_that(desc = "Testing findMetaCells logfile", {
    logLines <- readLines(logfile)
    expect_true(any(grepl("Meta-cell aggregation starts", logLines)))
    expect_equal(sum(grepl("^Sample:", logLines)), 2)
    expect_equal(sum(grepl("Starting recursive cell population",
        logLines)), 2)
    expect_true(any(grepl("Completed meta-cell aggregation", logLines)))
})

test_that(desc = "Testing findMetaCells with sampleLabel = NULL", {
    sceOne <- sceSim[, sampleLabel == "Sample_1"]
    res <- findMetaCells(sceOne,
        initK = 2,
        maxK = 8,
        tempL = 10,
        verbose = FALSE)
    expect_true(all(grepl("^Sample_1_K[0-9]+$", res$metacell)))
    sceOne$sample <- "tumor1"
    expect_warning(findMetaCells(sceOne,
        initK = 2,
        maxK = 8,
        tempL = 10,
        verbose = FALSE),
        "all cells are treated as one sample")
})

test_that(desc = "Testing findMetaCells sampleLabel error checking", {
    expect_error(findMetaCells(sceSim,
        sampleLabel = sampleLabel[-1],
        verbose = FALSE),
        "must be the same length")
    naLabel <- sampleLabel
    naLabel[seq(5)] <- NA
    expect_error(findMetaCells(sceSim,
        sampleLabel = naLabel,
        verbose = FALSE),
        "contains NA values")
})

test_that(desc = "Testing findMetaCells skips samples below minCellMeta", {
    keep <- c(which(sampleLabel == "Sample_1"),
        which(sampleLabel == "Sample_2")[seq(5)])
    expect_message(res <- findMetaCells(sceSim[, keep],
        initK = 2,
        maxK = 8,
        tempL = 10,
        sampleLabel = sampleLabel[keep],
        verbose = FALSE),
        "Sample 'Sample_2' has 5 cells, fewer than minCellMeta \\(10\\)")
    expect_equal(unique(as.character(res$celda_sample_label)), "Sample_1")
    expect_equal(ncol(res), sum(sampleLabel == "Sample_1"))
    expect_error(suppressMessages(findMetaCells(sceSim[, seq(12)],
        sampleLabel = rep(c("a", "b"), 6),
        verbose = FALSE)),
        "No sample has at least minCellMeta")
})

test_that(desc = "Testing findMetaCells argument checking", {
    expect_error(findMetaCells(sceSim, kneePoint = 3, verbose = FALSE),
        "'kneePoint' must be 1 or 2")
    expect_error(findMetaCells(sceSim, initK = 0, verbose = FALSE),
        "Value of K cannot be zero")
    expect_error(findMetaCells(sceSim, maxK = 0, verbose = FALSE),
        "Value of K cannot be zero")
    expect_error(findMetaCells(sceSim, tempL = NULL, verbose = FALSE),
        "'tempL' must be a positive integer")
    expect_error(findMetaCells(sceSim, tempL = 0, verbose = FALSE),
        "'tempL' must be a positive integer")
})

test_that(desc = "Testing .findKnee", {
    K <- seq(2, 20)
    P <- 20 + 80 * exp(-0.3 * (K - 2))
    expect_equal(.findKnee(K, P, nPoints = 1), 8)
    expect_equal(.findKnee(K, P, nPoints = 2), 12)
    expect_equal(.findKnee(K, rep(5, length(K))), 2)
})
