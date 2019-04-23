# celda_CG
library(celda)
context("Testing celda_CG")

celdaCGSim <- simulateCells("celda_CG", K = 5, L = 10)
modelCG <- celda_CG(counts = celdaCGSim$counts,
    sampleLabel = celdaCGSim$sampleLabel,
    K = celdaCGSim$K,
    L = celdaCGSim$L,
    nchains = 2,
    algorithm = "EM",
    verbose = FALSE)
factorized <- factorizeMatrix(celdaMod = modelCG, counts = celdaCGSim$counts)

# celda_CG
test_that(desc = "Testing simulation and celda_CG model", {
    expect_equal(typeof(celdaCGSim$counts), "integer")
    expect_true(all(sweep(factorized$counts$cell,
        2,
        colSums(celdaCGSim$counts),
        "/") == factorized$proportions$cell))
    expect_equal(celdaCGSim$K,
        ncol(factorized$proportions$cellPopulation))
    expect_equal(celdaCGSim$L,
        nrow(factorized$proportions$cellPopulation))
    # expect_true(all(is.numeric(
    #     logLikelihoodHistory(celdaMod = modelCG))))
    # expect_equal(max(logLikelihoodHistory(celdaMod = modelCG)),
    #     bestLogLikelihood(modelCG))

    # GitHub #347
    numericCounts <- celdaCGSim$counts
    storage.mode(numericCounts) <- "numeric"
    expect_true(is(celda_CG(counts = celdaCGSim$counts,
        sampleLabel = celdaCGSim$sampleLabel,
        K = celdaCGSim$K,
        L = celdaCGSim$L,
        algorithm = "EM",
        verbose = FALSE,
        maxIter = 1,
        nchains = 1),
        "celda_CG"))
})

# Cluster probabilities
test_that(desc = "Testing clusterProbability with celda_CG", {
    clustProb <- clusterProbability(counts = celdaCGSim$counts, modelCG)
    expect_true(all(round(rowSums(clustProb$zProbability), 10) == 1) &
            nrow(clustProb$zProbability) == ncol(celdaCGSim$counts))
    expect_true(all(round(rowSums(clustProb$yProbability), 10) == 1) &
            nrow(clustProb$yProbability) == nrow(celdaCGSim$counts))

    clustProb <- clusterProbability(counts = celdaCGSim$counts,
        modelCG, log = TRUE)
    expect_true(all(round(rowSums(.normalizeLogProbs(clustProb$zProbability)),
        10) == 1) & nrow(clustProb$zProbability) == ncol(celdaCGSim$counts))
    expect_true(all(round(rowSums(.normalizeLogProbs(clustProb$yProbability)),
        10) == 1) & nrow(clustProb$yProbability) == nrow(celdaCGSim$counts))
})

test_that(desc = paste0("Testing simulateCells.celda_CG error checking with",
    " low gamma"), {
        expect_error(simulateCells(model = "celda_CG", gamma = 0.000001))
    })

test_that(desc = paste0("Testing simulateCells.celda_CG, make sure all genes",
    " expressed"), {
        simCellsLow <- simulateCells(model = "celda_CG",
            G = 1000,
            C = 300,
            CRange = c(1, 100),
            NRange = c(1, 100))
        expect_true(all(rowSums(simCellsLow$counts) > 0))
    })

# test_that(desc = "Testing celdaGridSearch with celda_CG", {
#     expect_error(celdaGridSearch(counts = celdaCGSim$counts,
#         model = "celda_CG",
#         maxIter = 1,
#         nchains = 1,
#         paramsTest = list(K = 5, L = 10, M = c(3, 4)),
#         paramsFixed = list(sampleLabel = celdaCGSim$sampleLabel),
#         bestOnly = FALSE),
#         paste0("The following elements in 'paramsTest' are not arguments of",
#             " 'celda_CG': M"))
#
#     expect_error(celdaGridSearch(counts = celdaCGSim$counts,
#         model = "celda_CG",
#         nchains = 1,
#         maxIter = 1,
#         paramsTest = list(K = 5, L = 10, sampleLabel = "Sample"),
#         paramsFixed = list(sampleLabel = celdaCGSim$sampleLabel)),
#         paste0("Setting parameters such as 'zInit', 'yInit', and",
#             " 'sampleLabel' in 'paramsTest' is not currently supported."))
#
#     expect_error(celdaGridSearch(counts = celdaCGSim$counts,
#         model = "celda_CG",
#         nchains = 1,
#         maxIter = 1,
#         paramsTest = list(),
#         paramsFixed = list(sampleLabel = celdaCGSim$sampleLabel)),
#         paste0("The following arguments are not in 'paramsTest' or",
#             " 'paramsFixed' but are required for 'celda_CG': K,L"))
#
#     expect_error(celdaGridSearch(counts = celdaCGSim$counts,
#         model = "celda_CG",
#         nchains = 1,
#         maxIter = 1,
#         paramsTest = list(K = c(4, 5), L = c(9, 10)),
#         paramsFixed = list(sampleLabel = celdaCGSim$sampleLabel,
#             xxx = "xxx")),
#         paste0("The following elements in 'paramsFixed' are not arguments",
#             " of 'celda_CG': xxx"))
#
#     expect_warning(celdaGridSearch(counts = celdaCGSim$counts,
#         model = "celda_CG",
#         maxIter = 1,
#         perplexity = FALSE,
#         paramsTest = list(K = 5, L = 10, nchains = 2),
#         paramsFixed = list(zInitialize = "random", yInitialize = "random")),
#         paste0("Parameter 'nchains' should not be used within the",
#             " paramsTest list"))
#
#     celdaCGres <- celdaGridSearch(counts = celdaCGSim$counts,
#         model = "celda_CG",
#         nchains = 2,
#         paramsTest = list(K = 5, L = 10),
#         paramsFixed = list(sampleLabel = celdaCGSim$sampleLabel,
#             zInitialize = "random",
#             yInitialize = "random"),
#         maxIter = 1,
#         verbose = FALSE,
#         bestOnly = FALSE,
#         perplexity = FALSE)
#     expect_true(is(celdaCGres, "celdaList"))
#     expect_error(plotGridSearchPerplexity(celdaCGres))
#     expect_equal(names(runParams(celdaCGres)),
#         c("index", "chain", "K", "L", "log_likelihood"))
#
#     celdaCGres <- resamplePerplexity(celdaCGSim$counts,
#         celdaCGres, resample = 2)
#     expect_is(celdaCGres, "celdaList")
#     expect_error(resamplePerplexity(celdaCGSim$counts,
#         celdaCGres, resample = "2"))
#     expect_error(resamplePerplexity(celdaCGSim$counts,
#         "celdaCGres", resample = 2))
#
#     plotObj <- plotGridSearchPerplexity(celdaCGres)
#     expect_is(plotObj, "ggplot")
#
#     expect_error(subsetCeldaList(celdaList = "celdaList"),
#         "celdaList parameter was not of class celdaList.")
#     expect_error(subsetCeldaList(celdaCGres, params = list(K = 7, L = 11)))
#     expect_error(subsetCeldaList(celdaCGres, params = list(K = 5, M = 10)))
#
#     celdaCGresK5L10 <- subsetCeldaList(celdaCGres,
#          params = list(K = 5, L = 10))
#     modelCG <- selectBestModel(celdaCGresK5L10)
#
#     expect_error(selectBestModel("celdaList"),
#         "celdaList parameter was not of class celdaList.")
#     expect_error(celdaCGres <- resamplePerplexity(celdaCGSim$counts,
#         modelCG, resample = 2))
#     expect_error(celdaCGres <- resamplePerplexity(celdaCGSim$counts,
#         celdaCGres, resample = "a"))
#
#     celdaCGresIndex1 <- subsetCeldaList(celdaCGres, params = list(index = 1))
#     expect_true(all(is(celdaCGresIndex1, "celda_CG") &&
#             !is(celdaCGresIndex1, "celdaList")))
#
#     res <- perplexity(celdaCGSim$counts, modelCG)
#     res2 <- perplexity(celdaCGSim$counts,
#         modelCG, newCounts = celdaCGSim$counts + 1)
#
#     expect_error(res <- perplexity(celdaCGSim$counts,
#         modelCG, newCounts = celdaCGSim$counts[-1, ]))
# })

# Ensure logLikelihood calculates the expected values
test_that(desc = "Testing logLikelihood.celda_CG", {
    expect_lt(logLikelihood(model = "celda_CG",
        y = celdaCGSim$y,
        z = celdaCGSim$z,
        delta = 1,
        gamma = 1,
        beta = 1,
        alpha = 1,
        K = celdaCGSim$K,
        L = celdaCGSim$L,
        s = celdaCGSim$sampleLabel,
        counts = celdaCGSim$counts),
        0)

    fakeZ <- celdaCGSim$z
    fakeZ[1] <- celdaCGSim$K + 1
    expect_error(logLikelihood(model = "celda_CG",
        y = celdaCGSim$y,
        z = fakeZ,
        delta = 1,
        gamma = 1,
        beta = 1,
        alpha = 1,
        K = celdaCGSim$K,
        L = celdaCGSim$L,
        s = celdaCGSim$sampleLabel,
        counts = celdaCGSim$counts),
        "An entry in z contains a value greater than the provided K.")

    fakeY <- celdaCGSim$y
    fakeY[1] <- celdaCGSim$L + 1
    expect_error(logLikelihood(model = "celda_CG",
        y = fakeY,
        z = celdaCGSim$z,
        delta = 1,
        gamma = 1,
        beta = 1,
        alpha = 1,
        K = celdaCGSim$K,
        L = celdaCGSim$L,
        s = celdaCGSim$sampleLabel,
        counts = celdaCGSim$counts),
        "An entry in y contains a value greater than the provided L.")
})

# normalizeCounts
test_that(desc = "Testing normalizeCounts with celda_CG", {
    normCounts <- normalizeCounts(celdaCGSim$counts)
    expect_equal(dim(normCounts), dim(celdaCGSim$counts))
    expect_equal(rownames(normCounts), rownames(celdaCGSim$counts))
    expect_equal(colnames(normCounts), colnames(celdaCGSim$counts))
    expect_error(normalizeCounts(celdaCGSim$counts,
        transformationFun = "scale"),
        "'transformationFun' needs to be of class 'function'")
    expect_error(normalizeCounts(celdaCGSim$counts, scaleFun = "scale"),
        "'scaleFun' needs to be of class 'function'")
})

# recodeClusterY
test_that(desc = "Testing recodeClusterY with celda_CG", {
    expect_error(recodeClusterY(celdaMod = modelCG,
        from = NULL,
        to = ""))
    expect_error(recodeClusterY(celdaMod = modelCG,
        from = c(1, 2, 3, 4, 5),
        to = c(1, 2, 4, 3, 6)))
    expect_error(recodeClusterY(celdaMod = modelCG,
        from = c(1, 2, 3, 4, 6),
        to = c(1, 2, 4, 3, 5)))
    newRecoded <- recodeClusterY(celdaMod = modelCG,
        from = c(1, 2, 3, 4, 5),
        to = c(3, 2, 1, 4, 5))
    expect_equal(modelCG@clusters$y == 1, newRecoded@clusters$y == 3)
})

# recodeClusterZ
test_that(desc = "Testing recodeClusterZ with celda_CG", {
    expect_error(recodeClusterZ(celdaMod = modelCG,
        from = NULL,
        to = ""))
    expect_error(recodeClusterZ(celdaMod = modelCG,
        from = c(1, 2, 3, 4, 5),
        to = c(1, 2, 3, 4, 6)))
    expect_error(recodeClusterZ(celdaMod = modelCG,
        from = c(1, 2, 3, 4, 6),
        to = c(1, 2, 3, 4, 5)))
    newRecoded <- recodeClusterZ(celdaMod = modelCG,
        from = c(1, 2, 3, 4, 5),
        to = c(5, 4, 3, 2, 1))
    expect_equal(modelCG@clusters$z == 1, newRecoded@clusters$z == 5)
})

# compareCountMatrix
test_that(desc = "Testing CompareCountMatrix with celda_CG", {
    expect_true(compareCountMatrix(counts = celdaCGSim$counts,
        celdaMod = modelCG))
    lessFeatures <- celdaCGSim$counts[seq(50), ]
    expect_error(compareCountMatrix(counts = lessFeatures, celdaMod = modelCG),
        paste0("The provided celda object was generated from a counts matrix",
            " with a different number of features than the one provided."))
    lessCells <- celdaCGSim$counts[, seq(100)]
    expect_error(compareCountMatrix(counts = lessCells, celdaMod = modelCG),
        paste0("The provided celda object was generated from a counts matrix",
            " with a different number of cells than the one provided."))
    countsMatrixError <- matrix(data = 1,
        nrow = nrow(celdaCGSim$counts),
        ncol = ncol(celdaCGSim$counts))
    expect_false(compareCountMatrix(counts = countsMatrixError,
        celdaMod = modelCG,
        errorOnMismatch = FALSE))
    expect_error(compareCountMatrix(counts = countsMatrixError,
        celdaMod = modelCG,
        errorOnMismatch = TRUE))
})

# topRank
test_that(desc = "Testing topRank with celda_CG", {
    topRank <- topRank(matrix = factorized$proportions$module,
        n = 1000,
        threshold = NULL)
    expect_equal(names(topRank), c("index", "names"))
    topRank <- topRank(matrix = factorized$proportions$module, n = 1000)
    expect_equal(nrow(celdaCGSim$counts), sum(sapply(topRank$names, length)))
    expect_equal(names(topRank), c("index", "names"))
})

# plotHeatmap
# test_that(desc = "Testing plotHeatmap with celda_CG", {
#     expect_error(plotHeatmap(counts = celdaCGSim$counts,
#         z = modelCG@clusters$y,
#         y = modelCG@clusters$y),
#         "Length of z must match number of columns in counts matrix")
#     expect_error(plotHeatmap(counts = celdaCGSim$counts,
#         z = modelCG@clusters$z,
#         y = modelCG@clusters$z),
#         "Length of y must match number of rows in counts matrix")
#     expect_error(plotHeatmap(counts = celdaCGSim$counts,
#         z = modelCG@clusters$z,
#         y = modelCG@clusters$y,
#         scaleRow = modelCG),
#         "'scaleRow' needs to be of class 'function'")
#     expect_error(plotHeatmap(counts = celdaCGSim$counts,
#         z = modelCG@clusters$z,
#         y = modelCG@clusters$y,
#         trim = 3),
#         paste0("'trim' should be a 2 element vector specifying the lower",
#             " and upper boundaries"))
#     expect_equal(names(plotHeatmap(counts = celdaCGSim$counts,
#         z = modelCG@clusters$z,
#         y = modelCG@clusters$y,
#         cellIx = seq(10))),
#         c("treeRow", "treeCol", "gtable"))
#     expect_equal(names(plotHeatmap(counts = celdaCGSim$counts,
#         z = NULL,
#         y = modelCG@clusters$y,
#         cellIx = seq(10),
#         colorScheme = "sequential")),
#         c("treeRow", "treeCol", "gtable"))
#     expect_equal(names(plotHeatmap(counts = celdaCGSim$counts,
#         z = modelCG@clusters$z,
#         y = modelCG@clusters$y,
#         featureIx = seq(10))),
#         c("treeRow", "treeCol", "gtable"))
#     expect_equal(names(plotHeatmap(counts = celdaCGSim$counts,
#         z = modelCG@clusters$z,
#         y = NULL,
#         featureIx = seq(10))),
#         c("treeRow", "treeCol", "gtable"))
#     expect_equal(names(plotHeatmap(counts = celdaCGSim$counts,
#         z = modelCG@clusters$z,
#         y = modelCG@clusters$y,
#         cellIx = seq(10),
#         colorScheme = "sequential",
#         annotationColor = "default")),
#         c("treeRow", "treeCol", "gtable"))
# })

# plotHeatmap with annotations
# test_that(desc = "Testing plotHeatmap with annotations", {
#     annotCell <- as.data.frame(c(rep(x = 1,
#         times = ncol(celdaCGSim$counts) - 100),
#         rep(x = 2, 100)))
#     annotFeature <- as.data.frame(c(rep(x = 1,
#         times = nrow(celdaCGSim$counts) - 100),
#         rep(x = 2, 100)))
#
#     rownames(annotCell) <- colnames(celdaCGSim$counts)
#     colnames(annotCell) <- "cell"
#     rownames(annotFeature) <- rownames(celdaCGSim$counts)
#     colnames(annotFeature) <- "feature"
#     expect_equal(names(plotHeatmap(celdaMod = modelCG,
#         counts = celdaCGSim$counts,
#         annotationCell = annotCell,
#         annotationFeature = annotFeature,
#         z = modelCG@clusters$z,
#         y = modelCG@clusters$y)),
#         c("treeRow", "treeCol", "gtable"))
#
#     rownames(annotCell) <- NULL
#     rownames(annotFeature) <- NULL
#     expect_equal(names(plotHeatmap(celdaMod = modelCG,
#         counts = celdaCGSim$counts,
#         annotationCell = as.matrix(annotCell),
#         annotationFeature = as.matrix(annotFeature),
#         z = modelCG@clusters$z,
#         y = modelCG@clusters$y)),
#         c("treeRow", "treeCol", "gtable"))
# })

# celdaHeatmap
test_that(desc = "Testing celdaHeatmap with celda_CG", {
    expect_equal(names(celdaHeatmap(celdaMod = modelCG,
        counts = celdaCGSim$counts)),
        c("treeRow", "treeCol", "gtable"))
})

# moduleHeatmap
test_that(desc = "Checking moduleHeatmap to see if it runs", {
    expect_equal(names(moduleHeatmap(celdaCGSim$counts,
        celdaMod = modelCG,
        featureModule = c(2, 3),
        topCells = 500)),
        c("treeRow", "treeCol", "gtable"))
    expect_equal(names(moduleHeatmap(celdaCGSim$counts,
        celdaMod = modelCG,
        topFeatures = 15,
        topCells = 15,
        normalize = FALSE)),
        c("treeRow", "treeCol", "gtable"))
    expect_equal(names(moduleHeatmap(celdaCGSim$counts,
        celdaMod = modelCG,
        topFeatures = 15,
        topCells = NULL,
        normalize = FALSE)),
        c("treeRow", "treeCol", "gtable"))
    expect_error(moduleHeatmap(counts = "counts",
        celdaMod = modelCG,
        featureModule = c(2, 3)),
        "'counts' should be a numeric count matrix")
    expect_error(moduleHeatmap(counts = celdaCGSim$counts,
        celdaMod = "model",
        featureModule = c(2, 3)),
        "'celdaMod' should be an object of class celda_G or celda_CG")
})

# celdaProbabilityMap
# test_that(desc = "Testing celdaProbabiltyMap.celda_CG for sample level", {
#     plotObj <- celdaProbabilityMap(counts = celdaCGSim$counts,
#         celdaMod = modelCG,
#         level = "sample")
#     expect_true(!is.null(plotObj))
#
#     ## Without a sample label
#     modelCG2 <- celda_CG(celdaCGSim$counts,
#         sampleLabel = NULL,
#         K = celdaCGSim$K,
#         L = celdaCGSim$L,
#         maxIter = 5,
#         nchain = 1)
#     plotObj <- celdaProbabilityMap(counts = celdaCGSim$counts,
#         celdaMod = modelCG2,
#         level = "sample")
#     expect_true(!is.null(plotObj))
# })

test_that(desc = "Testing celdaProbabiltyMap.celda_CG for cellPopulation", {
    plotObj <- celdaProbabilityMap(counts = celdaCGSim$counts,
        celdaMod = modelCG,
        level = "cellPopulation")
    expect_true(!is.null(plotObj))
})

# differentialExpression
# test_that(desc = "Testing differentialExpression for celda_CG", {
#     expect_equal(class(diffExpK1 <- differentialExpression(
#         counts = celdaCGSim$counts,
#         celdaMod = modelCG,
#         c1 = 3,
#         log2fcThreshold = 0.5,
#         onlyPos = TRUE)),
#         c("data.table", "data.frame"))
#     expect_equal(class(diffExpK1 <- differentialExpression(
#         counts = celdaCGSim$counts,
#         celdaMod = modelCG,
#         c1 = c(2, 3),
#         c2 = 4,
#         log2fcThreshold = 0.5)),
#         c("data.table", "data.frame"))
#     expect_error(differentialExpression(counts = "counts",
#         celdaMod = modelCG,
#         c1 = 3,
#         log2fcThreshold = 0.5),
#         "'counts' should be a numeric count matrix")
#     expect_error(differentialExpression(counts = celdaCGSim$counts,
#         celdaMod = NULL,
#         c1 = 3),
#         "'celdaMod' should be an object of class celda_C or celda_CG")
#     expect_error(differentialExpression(counts = celdaCGSim$counts,
#         celdaMod = modelCG,
#         c1 = NULL,
#         log2fcThreshold = 0.5,
#         onlyPos = TRUE))
# })

# plotDimReduce
test_that(desc = "Testing plotDimReduce* with celda_CG", {
    celdaTsne <- celdaTsne(counts = celdaCGSim$counts,
        maxIter = 50,
        celdaMod = modelCG,
        maxCells = 500)
    expect_equal(names(plotDimReduceCluster(dim1 = celdaTsne[, 1],
        dim2 = celdaTsne[, 2],
        cluster = as.factor(modelCG@clusters$z),
        specificClusters = c(1, 2, 3))),
        c("data",
            "layers",
            "scales",
            "mapping",
            "theme",
            "coordinates",
            "facet",
            "plot_env",
            "labels",
            "guides"))
    expect_equal(names(plotDimReduceCluster(
        dim1 = celdaTsne[, 1],
        dim2 = celdaTsne[, 2],
        cluster = as.factor(modelCG@clusters$z),
        specificClusters = c(1, 2, 3),
        labelClusters = TRUE)),
        c("data",
            "layers",
            "scales",
            "mapping",
            "theme",
            "coordinates",
            "facet",
            "plot_env",
            "labels",
            "guides"))
    expect_equal(names(plotDimReduceModule(
        dim1 = celdaTsne[, 1],
        dim2 = celdaTsne[, 2],
        counts = celdaCGSim$counts,
        celdaMod = modelCG,
        modules = c(1, 2))),
        c("data",
            "layers",
            "scales",
            "mapping",
            "theme",
            "coordinates",
            "facet",
            "plot_env",
            "labels"))
    expect_equal(names(plotDimReduceModule(dim1 = celdaTsne[, 1],
        dim2 = celdaTsne[, 2],
        counts = celdaCGSim$counts,
        celdaMod = modelCG,
        modules = c(1, 2),
        rescale = FALSE)),
        c("data",
            "layers",
            "scales",
            "mapping",
            "theme",
            "coordinates",
            "facet",
            "plot_env",
            "labels"))
    expect_equal(names(plotDimReduceFeature(dim1 = celdaTsne[, 1],
        dim2 = celdaTsne[, 2],
        counts = celdaCGSim$counts,
        features = c("Gene_99"),
        exactMatch = TRUE)),
        c("data",
            "layers",
            "scales",
            "mapping",
            "theme",
            "coordinates",
            "facet",
            "plot_env",
            "labels"))
    expect_equal(names(plotDimReduceFeature(
        dim1 = celdaTsne[, 1],
        dim2 = celdaTsne[, 2],
        counts = celdaCGSim$counts,
        features = c("Gene_99"),
        exactMatch = FALSE)),
        c("data",
            "layers",
            "scales",
            "mapping",
            "theme",
            "coordinates",
            "facet",
            "plot_env",
            "labels"))
    expect_error(plotDimReduceModule(dim1 = celdaTsne[, 1],
        dim2 = celdaTsne[, 2],
        counts = celdaCGSim$counts,
        celdaMod = modelCG,
        modules = c(11, 12)))
    expect_error(plotDimReduceFeature(dim1 = celdaTsne[, 1],
        dim2 = celdaTsne[, 2],
        counts = celdaCGSim$counts,
        features = NULL,
        exactMatch = TRUE))
    expect_error(plotDimReduceFeature(dim1 = celdaTsne[, 1],
        dim2 = celdaTsne[, 2],
        counts = celdaCGSim$counts,
        features = c("Gene_99"),
        trim = 2,
        exactMatch = TRUE))
    expect_error(plotDimReduceFeature(dim1 = celdaTsne[, 1],
        dim2 = celdaTsne[, 2],
        counts = celdaCGSim$counts,
        features = c("Nonexistent_Gene"),
        exactMatch = TRUE))

    # Check cases when there are some or all features not present in the counts
    # matrix
    expect_error(plotDimReduceFeature(dim1 = celdaTsne[, 1],
        dim2 = celdaTsne[, 2],
        counts = celdaCGSim$counts,
        features = c("Nonexistent_Gene"),
        exactMatch = TRUE))
    expect_warning(plotDimReduceFeature(dim1 = celdaTsne[, 1],
        dim2 = celdaTsne[, 2],
        counts = celdaCGSim$counts,
        features = c("Gene_99", "Nonexistent_Gene"),
        exactMatch = TRUE))
    expect_warning(plotDimReduceFeature(dim1 = celdaTsne[, 1],
        dim2 = celdaTsne[, 2],
        counts = celdaCGSim$counts,
        features = c("Gene_99", "Nonexistent_Gene"),
        exactMatch = FALSE))
})

# celdaTsne
test_that(desc = paste0("Testing celdaTsne with celda_CG when model class",
    " is changed, should error"), {
        modelX <- modelCG
        class(modelX) <- "celda_X"
        expect_error(celdaTsne(counts = celdaCGSim$counts, celdaMod = modelX),
            "unable to find")
    })

# test_that(desc = "Testing celdaTsne.celda_CG with all cells", {
#     tsne <- celdaTsne(counts = celdaCGSim$counts,
#         celdaMod = modelCG,
#         maxCells = length(modelCG@clusters$z))
#     plotObj <- plotDimReduceCluster(tsne[, 1], tsne[, 2], modelCG@clusters$z)
#     expect_true(ncol(tsne) == 2 & nrow(tsne) == length(modelCG@clusters$z))
#     expect_true(!is.null(plotObj))
#
#     tsne <- celdaTsne(counts = celdaCGSim$counts,
#         celdaMod = modelCG,
#         maxCells = ncol(celdaCGSim$counts),
#         modules = c(1, 2))
#     expect_error(tsne <- celdaTsne(counts = celdaCGSim$counts,
#         celdaMod = modelCG,
#         maxCells = ncol(celdaCGSim$counts),
#         modules = seq(1000, 1005)))
# })

test_that(desc = "Testing celdaTsne.celda_CG with subset of cells", {
    expect_success(expect_error(tsne <- celdaTsne(counts = celdaCGSim$counts,
        celdaMod = modelCG,
        maxCells = 50,
        minClusterSize = 50)))
    tsne <- celdaTsne(counts = celdaCGSim$counts,
        celdaMod = modelCG,
        maxCells = 100,
        minClusterSize = 10)
    plotObj <- plotDimReduceCluster(tsne[, 1], tsne[, 2], modelCG@clusters$z)
    expect_true(ncol(tsne) == 2 & nrow(tsne) == length(modelCG@clusters$z) &&
            sum(!is.na(tsne[, 1])) == 100)
    expect_true(!is.null(plotObj))
})

# celdaUmap
test_that(desc = paste0("Testing celdaUmap with celda_CG when model class is",
    " changed, should error"), {
        modelX <- modelCG
        class(modelX) <- "celda_X"
        expect_error(celdaUmap(counts = celdaCGSim$counts, celdaMod = modelX),
            "unable to find")
    })

# test_that(desc = "Testing celdaUmap.celda_CG with all cells", {
#     umap <- celdaUmap(counts = celdaCGSim$counts,
#         celdaMod = modelCG,
#         maxCells = length(modelCG@clusters$z))
#     plotObj <- plotDimReduceCluster(umap[, 1], umap[, 2], modelCG@clusters$z)
#     expect_true(ncol(umap) == 2 & nrow(umap) == length(modelCG@clusters$z))
#     expect_true(!is.null(plotObj))
#
#     umap <- celdaUmap(counts = celdaCGSim$counts,
#         celdaMod = modelCG,
#         maxCells = ncol(celdaCGSim$counts),
#         modules = c(1, 2))
#     expect_error(umap <- celdaUmap(counts = celdaCGSim$counts,
#         celdaMod = modelCG,
#         maxCells = ncol(celdaCGSim$counts),
#         modules = seq(1000, 1005)))
# })

test_that(desc = "Testing celdaUmap.celda_CG with subset of cells", {
    # expect_success(expect_error(umap <- celdaUmap(counts = celdaCGSim$counts,
    #     celdaMod = modelCG,
    #     maxCells = 50,
    #     minClusterSize = 50)))
    umap <- celdaUmap(counts = celdaCGSim$counts,
        celdaMod = modelCG,
        maxCells = 100,
        minClusterSize = 10)
    plotObj <- plotDimReduceCluster(umap[, 1], umap[, 2], modelCG@clusters$z)
    expect_true(ncol(umap) == 2 & nrow(umap) == length(modelCG@clusters$z) &&
            sum(!is.na(umap[, 1])) == 100)
    expect_true(!is.null(plotObj))
})

# featureModuleLookup
test_that(desc = "Testing featureModuleLookup with celda_CG", {
    res <- featureModuleLookup(celdaCGSim$counts, modelCG, "Gene_1")
    expect_true(res == modelCG@clusters$y[1])

    res <- featureModuleLookup(celdaCGSim$counts,
        modelCG, "Gene_2", exactMatch = FALSE)
    expect_true(length(res) == 11)

    res <- featureModuleLookup(celdaCGSim$counts, modelCG, "XXXXXXX")
    expect_true(grepl("No feature", res))
})

# .cCGSplitZ/.cCGSplitZ
test_that(desc = "Testing .cCGSplitZ and .cCGSplitY", {
    r <- simulateCells("celda_CG",
            S = 1,
            G = 100,
            CRange = c(50, 100),
            K = 2,
            L = 2)
    modelCG <- celda_CG(r$counts,
            K = r$K,
            L = r$L,
            maxIter = 5,
            nchain = 1)
    probs <- clusterProbability(r$counts, modelCG, log = TRUE)

    dc <- .cCGDecomposeCounts(r$counts, r$sampleLabel, r$z, r$y, r$K, r$L)
    res <- .cCGSplitZ(r$counts,
            dc$mCPByS,
            dc$nTSByC,
            dc$nTSByCP,
            dc$nByG,
            dc$nByTS,
            dc$nGByTS,
            as.integer(r$sampleLabel),
            z = r$z,
            K = r$K,
            L = r$L,
            nS = dc$nS,
            nG = dc$nG,
            alpha = 1,
            beta = 1,
            delta = 1,
            gamma = 1,
            zProb = probs$zProbability,
            minCell = 1000)
    expect_true(grepl("Cluster sizes too small", res$message))
    res <- .cCGSplitY(r$counts,
            r$y,
            dc$mCPByS,
            dc$nGByCP,
            dc$nTSByC,
            dc$nTSByCP,
            dc$nByG,
            dc$nByTS,
            dc$nGByTS,
            dc$nCP,
            s = as.integer(r$sampleLabel),
            z = r$z,
            K = r$K,
            L = r$L,
            nS = dc$nS,
            nG = dc$nG,
            alpha = 1,
            beta = 1,
            delta = 1,
            gamma = 1,
            yProb = probs$yProbability,
            minCell = 1000)
    expect_true(grepl("Cluster sizes too small", res$message))

    ## Testing KSubclusters parameter
    res <- .cCGSplitY(r$counts,
            r$y,
            dc$mCPByS,
            dc$nGByCP,
            dc$nTSByC,
            dc$nTSByCP,
            dc$nByG,
            dc$nByTS,
            dc$nGByTS,
            dc$nCP,
            s = as.integer(r$sampleLabel),
            z = r$z,
            K = r$K,
            L = r$L,
            nS = dc$nS,
            nG = dc$nG,
            alpha = 1,
            beta = 1,
            delta = 1,
            gamma = 1,
            yProb = probs$yProbability,
            KSubclusters = 1000)
    expect_true(length(res$y) == nrow(r$counts))
})

test_that(desc = "Testing perplexity.celda_CG", {
    expect_true(is.numeric(perplexity(celdaCGSim$counts, modelCG)))

    class(modelCG) <- c("celda_C")
    expect_error(perplexity.celda_CG(celdaCGSim$counts, modelCG),
        "could not find function \"perplexity.celda_CG\"")
})

# test_that(desc = "Testing featureModuleTable", {
#     table <- featureModuleTable(celdaCGSim$counts, modelCG, outputFile = NULL)
#     expect_equal(ncol(table), 10)
# })

test_that(desc = "Testing violinPlot", {
    violin <- violinPlot(counts = celdaCGSim$counts,
        celdaMod = modelCG,
        features = "Gene_1")
    expect_equal(names(violin),
        c("data",
            "layers",
            "scales",
            "mapping",
            "theme",
            "coordinates",
            "facet",
            "plot_env",
            "labels"))
})

# miscellaneous fxns
# functions used internally
test_that(desc = "Invoking error from distinctColors function", {
    expect_error(distinctColors(n = 3, hues = "xx"),
        "Only color names listed in the 'color' function can be used in 'hues'")
})

test_that(desc = "Invoking error from sample labels function", {
    expect_error(.processSampleLabels("Sample_1", ncol(celdaCGSim$counts)),
        paste0("'sampleLabel' must be the same length as the number of columns",
            " in the 'counts' matrix."))
})

test_that(desc = "Invoking error from .logMessages function", {
    expect_error(.logMessages(date(), logfile = 5))
})

test_that(desc = paste0("miscellaneous distance fxns that are not directly",
    " used within celda, but will be tested"), {
        x <- data.frame(x = seq(2, 4), y = seq(1, 3))
        expect_equal(class(.hellingerDist(x)), "dist")
        expect_equal(class(.spearmanDist(x)), "dist")
    })
