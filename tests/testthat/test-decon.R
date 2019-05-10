## decontXoneBatch
library(celda)
context("Testing Deconx")

deconSim <- simulateContaminatedMatrix(K = 10, delta = c(1, 5))
modelDecontXoneBatch <- .decontXoneBatch(deconSim$observedCounts,
        z = deconSim$z,
        maxIter = 2)
modelDecontXoneBatchIter1 <- .decontXoneBatch(deconSim$observedCounts,
        z = deconSim$z,
        maxIter = 1)
modelDecontXoneBatchbg <- decontX(deconSim$observedCounts,
        maxIter = 2)

deconSim2 <- simulateContaminatedMatrix(K = 10, delta = 5)
batchDecontX <- decontX(cbind(deconSim$observedCounts,
    deconSim2$observedCounts),
        z = c(deconSim$z, deconSim2$z),
        batch = rep(seq(2), each = ncol(deconSim$observedCounts)),
        maxIter = 2)
batchDecontXBg <- decontX(cbind(deconSim$observedCounts,
    deconSim2$observedCounts),
        batch = rep(seq(2), each = ncol(deconSim$observedCounts)),
        maxIter = 2)

## simulateContaminatedMatrix
test_that(desc = "Testing simulateContaminatedMatrix", {
    expect_equivalent(object = colSums(deconSim$observedCounts),
        expected = deconSim$NByC)
    expect_equal(object = dim(deconSim$phi),
        expected = dim(deconSim$eta))
    expect_equal(typeof(deconSim$observedCounts), "integer")
    expect_warning(simulateContaminatedMatrix(K = 101, C = 10))
    #expect_equal(unique(deconSimKTooLarge$z), seq(ncol(deconSimKTooLarge$eta)))
})

## DecontX
test_that(desc = "Testing DecontX", {
    expect_equal(ncol(deconSim$observedCounts) + ncol(deconSim2$observedCounts),
        ncol(batchDecontX$resList$estNativeCounts))
    # expect_equal(length(batchDecontX$resList$estConp) ,
    # ncol(batchDecontX$resList$estNativeCounts))
    expect_equal(batchDecontXBg$method, "background")
})

## .decontXoneBatch
test_that(desc = "Testing .decontXoneBatch", {
    expect_equal(modelDecontXoneBatch$resList$estConp,
        1 - colSums(modelDecontXoneBatch$resList$estNativeCounts) /
            colSums(deconSim$observedCounts))
    expect_error(.decontXoneBatch(counts = deconSim$observedCounts,
        z = deconSim$z,
        delta = -1),
        "'delta' should be a single positive value.")
    expect_error(.decontXoneBatch(counts = deconSim$observedCounts,
        z = deconSim$z,
        delta = c(1, 1)),
        "'delta' should be a single positive value.")
    expect_error(.decontXoneBatch(counts = deconSim$observedCounts,
        z = c(deconSim$z, 1)),
        paste0("'z' must be of the same length as the number of cells in the",
            " 'counts' matrix."))
    expect_error(.decontXoneBatch(counts = deconSim$observedCounts,
        z = rep(1, ncol(
            deconSim$observedCounts))),
        "'z' must have at least 2 different values.")
    countsNA <- deconSim$observedCounts
    countsNA[1, 1] <- NA
    expect_error(.decontXoneBatch(counts = countsNA, z = deconSim$z),
        "Missing value in 'counts' matrix.")
})

test_that(desc = "Testing .decontXoneBatch using background distribution", {
    expect_equal(
        modelDecontXoneBatchbg$resList$estConp,
        1 - colSums(modelDecontXoneBatchbg$resList$estNativeCounts) /
            deconSim$NByC)
})

## logLikelihood
#test_that(desc = "Testing logLikelihood.DecontXoneBatch", {
    # z.process = processCellLabels(deconSim$z,
    # num.cells=ncol(deconSim$observedCounts) )
    # expect_equal( decon.calcLL(counts=deconSim$observedCounts, z=z.process  ,
    #    theta=modelDecontXoneBatch$resList$theta,
    # eta=modelDecontXoneBatch$resList$est.ConDist,
    # phi=modelDecontXoneBatch$resList$est.GeneDist ),
    # modelDecontXoneBatch$resList$logLikelihood[
    # modelDecontXoneBatch$runParams$iteration  ] )

    #cellDistModelBg <- normalizeCounts(
    #    modelDecontXoneBatchbg$resList$estNativeCounts,
    #    normalize = "proportion",
    #    pseudocountNormalize = 1e-20)
    #bgDistModelBg <- rowSums(deconSim$observedCounts) / sum(deconSim$NByC)
    #bgDistModelBg <- matrix(rep(bgDistModelBg,
    #    length(deconSim$NByC)), ncol = length(deconSim$NByC))
    #expect_equal(.bgCalcLL(counts = deconSim$observedCounts,
    #    theta = modelDecontXoneBatchbg$resList$theta,
    #    cellDist = cellDistModelBg,
    #    bgDist = bgDistModelBg),
    #    modelDecontXoneBatchbg$resList$logLikelihood[
    #        modelDecontXoneBatchbg$runParams$iteration])
#})

## decontamination EM updates
# test_that( desc = "Testing decontamination EM updates", {
#    z.process = processCellLabels(deconSim$z,
# num.cells=ncol(deconSim$observedCounts) )
#    expect_equal( cD.calcEMDecontamination( counts=deconSim$observedCounts,
# z=z.process, K=length(unique(deconSim$z)),
#        theta=modelDecontXoneBatchIter1$resList$theta,
# phi=modelDecontXoneBatchIter1$resList$est.GeneDist,
# eta=modelDecontXoneBatchIter1$resList$est.ConDist,
#        beta=modelDecontXoneBatchIter1$runParams$beta,
# delta=modelDecontXoneBatchIter1$runParams$delta)$theta,
# modelDecontXoneBatch$resList$theta )
# } )
