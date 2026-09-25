library(celda)
context("Testing DecontX wrapper functions")

# As of celda 1.23.0 the DecontX implementation lives in the standalone
# 'decontX' package and celda's functions are thin wrappers. These tests verify
# that the wrappers delegate correctly and preserve the expected output
# contract. They are skipped when 'decontX' is not installed.

test_that(desc = "simulateContamination wrapper returns expected structure", {
    skip_if_not_installed("decontX")
    deconSim <- simulateContamination(K = 10, delta = c(1, 5))
    expect_equivalent(object = colSums(deconSim$observedCounts),
        expected = deconSim$NByC)
    expect_equal(object = dim(deconSim$phi),
        expected = dim(deconSim$eta))
    expect_equal(typeof(deconSim$observedCounts), "integer")
    expect_warning(simulateContamination(K = 101, C = 10))
    expect_error(simulateContamination(K = 3, G = 2, numMarkers = 10))
})

test_that(desc = "decontX wrapper works on a counts matrix", {
    skip_if_not_installed("decontX")
    s <- simulateContamination()
    res <- decontX(s$observedCounts, verbose = FALSE)
    expect_true(!is.null(res$decontXcounts))
    expect_true(!is.null(res$contamination))
    expect_equal(length(res$contamination), ncol(s$observedCounts))

    # Plot wrappers return ggplot objects
    expect_s3_class(plotDecontXMarkerPercentage(s$observedCounts,
        z = res$z, markers = s$markers), "ggplot")
    expect_s3_class(plotDecontXMarkerPercentage(res$decontXcounts,
        z = res$z, markers = s$markers), "ggplot")
    expect_s3_class(plotDecontXMarkerExpression(s$observedCounts,
        s$markers[[1]], z = s$z), "ggplot")
    expect_s3_class(plotDecontXContamination(res), "ggplot")

    # Background input is passed through to decontX
    b <- s$observedCounts[, seq_len(5)]
    colnames(b) <- paste(colnames(b), "_", sep = "")
    res <- decontX(s$observedCounts, background = b, verbose = FALSE)
    expect_true(!is.null(res$decontXcounts))
})

test_that(desc = "decontX wrapper works on a SingleCellExperiment", {
    skip_if_not_installed("decontX")
    s <- simulateContamination()
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = s$observedCounts))
    sce <- decontX(sce, verbose = FALSE)

    # The wrapper preserves the decontX output contract on the SCE
    expect_true(!is.null(sce$decontX_contamination))
    expect_true(!is.null(sce$decontX_clusters))
    expect_true(!is.null(decontXcounts(sce)))
    expect_true("decontX_UMAP" %in%
        SingleCellExperiment::reducedDimNames(sce))
    expect_true(!is.null(S4Vectors::metadata(sce)$decontX))

    expect_s3_class(plotDecontXContamination(sce), "ggplot")
    expect_s3_class(plotDecontXMarkerPercentage(sce, z = s$z,
        markers = s$markers, assayName = "decontXcounts"), "ggplot")
    expect_s3_class(plotDecontXMarkerExpression(sce, s$markers[[1]]), "ggplot")

    newz <- paste0("X", s$z)
    sce$newz2 <- newz
    expect_s3_class(plotDecontXMarkerPercentage(sce, z = "newz2",
        markers = s$markers, assayName = "decontXcounts"), "ggplot")

    sce <- decontX(sce, estimateDelta = FALSE, verbose = FALSE)
    expect_true(!is.null(sce$decontX_contamination))

    # Background input is passed through to decontX
    bg <- sce[, seq_len(5)]
    colnames(bg) <- paste(colnames(bg), "_", sep = "")
    sce <- decontX(sce, background = bg, verbose = FALSE)
    expect_true(!is.null(sce$decontX_contamination))
})

test_that(desc = "decontX wrapper propagates input-validation errors", {
    skip_if_not_installed("decontX")
    deconSim <- simulateContamination(K = 10, delta = c(1, 5))
    expect_error(decontX(x = deconSim$observedCounts,
        z = deconSim$z, delta = c(1, -1), verbose = FALSE))
    expect_error(decontX(x = deconSim$observedCounts,
        z = deconSim$z, delta = c(1, 1, 1), verbose = FALSE))
    expect_error(decontX(x = deconSim$observedCounts,
        z = c(deconSim$z, 1), verbose = FALSE))
    countsNA <- deconSim$observedCounts
    countsNA[1, 1] <- NA
    expect_error(decontX(countsNA, z = deconSim$z, verbose = FALSE))
})
