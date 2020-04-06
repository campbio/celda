library(celda)
context("Testing DecontX functions")

deconSim <- simulateContamination(K = 10, delta = c(1, 5))
modelDecontXoneBatch <- decontX(deconSim$observedCounts,
        z = deconSim$z,
        maxIter = 2)

deconSim2 <- simulateContamination(K = 10, delta = c(1, 5))
batchDecontX <- decontX(cbind(deconSim$observedCounts,
    deconSim2$observedCounts),
        z = c(deconSim$z, deconSim2$z),
        batch = rep(seq(2), each = ncol(deconSim$observedCounts)),
        maxIter = 2)

test_that(desc = "Testing simulateContamination", {
    expect_equivalent(object = colSums(deconSim$observedCounts),
        expected = deconSim$NByC)
    expect_equal(object = dim(deconSim$phi),
        expected = dim(deconSim$eta))
    expect_equal(typeof(deconSim$observedCounts), "integer")
    expect_warning(simulateContamination(K = 101, C = 10))
    expect_error(simulateContamination(K = 3, G = 2, numMarkers = 10))
})

## DecontX
test_that(desc = "Testing DecontX on counts matrix", {
  s <- simulateContamination()
  res <- decontX(s$observedCounts)
  p <- plotDecontXMarkerPercentage(s$observedCounts,
                                   z = res$z,
                                   markers = s$markers)
  p <- plotDecontXMarkerPercentage(res$decontXcounts,
                                   z = res$z,
                                   markers = s$markers)
  p <- plotDecontXMarkerExpression(s$observedCounts,
                                   s$markers[[1]],
                                   z = s$z)
  p <- plotDecontXContamination(res)
})

test_that(desc = "Testing DecontX on SCE", {
  s <- simulateContamination()
  sce <- SingleCellExperiment::SingleCellExperiment(
                               list(counts = s$observedCounts))
  sce <- decontX(sce)
  p <- plotDecontXContamination(sce)
  p <- plotDecontXMarkerPercentage(sce,
                                   z = s$z,
                                   markers = s$markers,
                                   assayName = "decontXcounts")
  p <- plotDecontXMarkerExpression(sce, s$markers[[1]])
  newz <- paste0("X", s$z)
  sce$newz2 <- newz
  p <- plotDecontXMarkerPercentage(sce,
                                   z = "newz2",
                                   markers = s$markers,
                                   assayName = "decontXcounts")
  sce <- decontX(sce, estimateDelta = FALSE)
})

## .decontXoneBatch
test_that(desc = "Testing .decontXoneBatch", {
    expect_error(decontX(x = deconSim$observedCounts,
        z = deconSim$z, delta = c(1, -1)))
    expect_error(decontX(x = deconSim$observedCounts,
        z = deconSim$z, delta = c(1, 1, 1)))
    expect_error(decontX(x = deconSim$observedCounts,
        z = c(deconSim$z, 1)),
        paste0("'z' must be of the same length as the number of cells in the",
            " 'counts' matrix."))
    expect_error(.decontXoneBatch(counts = deconSim$observedCounts,
        z = rep(1, ncol(
            deconSim$observedCounts))),
        "No need to decontaminate when only one cluster is in the dataset.")
    countsNA <- deconSim$observedCounts
    countsNA[1, 1] <- NA
    expect_error(.decontXoneBatch(counts = countsNA, z = deconSim$z),
        "Missing value in 'counts' matrix.")
})
