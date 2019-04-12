# celda_G
library(celda)
context("Testing celda_G")

celdaGSim <- simulateCells("celda_G", L = 5, G = 100)
modelG <- celda_G(counts = celdaGSim$counts,
        L = celdaGSim$L,
        verbose = FALSE)
factorized <- factorizeMatrix(counts = celdaGSim$counts, celdaMod = modelG)

test_that(desc = "Testing celda_G model with numeric input matrix", {
    # Github Issue #347
    numericCounts <- celdaGSim$counts
    storage.mode(numericCounts) <- "numeric"
    expect_true(is(celda_G(counts = numericCounts,
        L = celdaGSim$L,
        maxIter = 1,
        verbose = FALSE),
        "celda_G"))
})

test_that(desc = "Testing clusterProbability with celda_G", {
    expect_true(ncol(clusterProbability(celdaGSim$counts,
        modelG)$yProbability) == celdaGSim$L)
})

test_that(desc = paste0("Testing simulateCells.celda_G error checking with",
    " low gamma"), {
    expect_error(simulateCells(model = "celda_G", gamma = 0.000001))
})

test_that(desc = paste0("Testing simulateCells.celda_G, make sure all genes",
    " expressed"), {
        simCellsLow <- simulateCells(model = "celda_G",
            G = 1000,
            C = 300,
            NRange = c(1, 10))
        expect_true(all(rowSums(simCellsLow$counts) > 0))
    })

test_that(desc = "Testing LogLikelihood functions", {
    expect_true(all(is.numeric(logLikelihoodHistory(celdaMod = modelG))))
    expect_equal(max(logLikelihoodHistory(celdaMod = modelG)),
        bestLogLikelihood(modelG))
})

# test_that(desc = "Testing celdaGridSearch with celda_G", {
#     celdaGRes <- celdaGridSearch(counts = celdaGSim$counts,
#         model = "celda_G",
#         nchains = 2,
#         paramsTest = list(L = c(5, 10)),
#         maxIter = 1,
#         verbose = FALSE,
#         bestOnly = FALSE,
#         perplexity = FALSE)
#
#     expect_error(celdaGridSearch(counts = celdaGSim$counts,
#         model = "celda_G",
#         paramsTest = list(L = 10, M = 4),
#         bestOnly = FALSE),
#         paste0("The following elements in 'paramsTest' are not arguments of",
#             " 'celda_G': M"))
#
#     expect_error(celdaGridSearch(counts = celdaGSim$counts,
#         model = "celda_G",
#         nchains = 1,
#         paramsTest = list(L = c(4, 5), yInit = 10)),
#         paste0("Setting parameters such as 'z.init', 'yInit', and",
#             " 'sample.label' in 'paramsTest' is not currently supported."))
#
#     expect_error(celdaGridSearch(counts = celdaGSim$counts,
#         model = "celda_G",
#         nchains = 1,
#         paramsTest = list()),
#         paste0("The following arguments are not in 'paramsTest' or",
#             " 'paramsFixed' but are required for 'celda_G': L"))
#
#     expect_error(celdaGridSearch(counts = celdaGSim$counts,
#         model = "celda_G",
#         nchains = 1,
#         paramsTest = list(L = 10),
#         paramsFixed = list(xxx = "xxx")),
#         paste0("The following elements in 'paramsFixed' are not arguments",
#             " of 'celda_G': xxx"))
#
#     expect_warning(celdaGridSearch(counts = celdaGSim$counts,
#         model = "celda_G",
#         paramsTest = list(L = c(5, 6), nchains = 2)),
#         paste0("Parameter 'nchains' should not be used within the",
#             " paramsTest list"))
#
#     expect_true(is(celdaGRes, "celdaList"))
#     expect_error(plotGridSearchPerplexity(celdaGRes))
#     expect_equal(names(runParams(celdaGRes)),
#         c("index", "chain", "L", "log_likelihood"))
#
#     celdaGRes <- resamplePerplexity(celdaGSim$counts, celdaGRes, resample = 2)
#     expect_equal(is.null(celdaGRes@perplexity), FALSE)
#     expect_is(celdaGRes, "celdaList")
#     expect_error(resamplePerplexity(celdaGSim$counts,
#         celdaGRes, resample = "2"))
#     expect_error(resamplePerplexity(celdaGSim$counts,
#         "celdaGRes", resample = 2))
#
#     plotObj <- plotGridSearchPerplexity(celdaGRes)
#     expect_is(plotObj, "ggplot")
#
#     celdaCRes <- celdaGridSearch(counts = celdaGSim$counts,
#         model = "celda_C",
#         nchains = 2,
#         paramsTest = list(K = c(5, 10)),
#         maxIter = 1,
#         verbose = FALSE,
#         bestOnly = TRUE)
#     expect_error(plotGridSearchPerplexity.celda_G(celdaCRes))
#
#     celdaGResIndex1 <- subsetCeldaList(celdaGRes, params = list(index = 1))
#     expect_true(all(is(celdaGResIndex1, "celda_G") &&
#             !is(celdaGResIndex1, "celdaList")))
#
#     expect_error(subsetCeldaList(celdaGRes, params = list(L = 11)))
#     expect_error(subsetCeldaList(celdaGRes, params = list(L = 5, M = 10)))
#
#     celdaGResL5 <- subsetCeldaList(celdaGRes, params = list(L = 5))
#     modelG <- selectBestModel(celdaGResL5)
#     res <- perplexity(celdaGSim$counts, modelG)
#     res2 <- perplexity(celdaGSim$counts,
#         modelG, newCounts = celdaGSim$counts + 1)
#
#     expect_error(res <- perplexity(celdaGSim$counts, modelG,
#         newCounts = celdaGSim$counts[-1, ]))
# })

# logLikelihood
test_that(desc = "Testing logLikelihood.celda_G", {
    expect_lt(logLikelihood(model = "celda_G",
        counts = celdaGSim$counts,
        y = celdaGSim$y,
        L = celdaGSim$L,
        delta = 1,
        gamma = 1,
        beta = 1),
        0)

    fakeY <- celdaGSim$y
    fakeY[1] <- celdaGSim$L + 1
    expect_error(logLikelihood(model = "celda_G",
        counts = celdaGSim$counts,
        y = fakeY,
        L = celdaGSim$L,
        delta = 1,
        gamma = 1,
        beta = 1),
        "An entry in y contains a value greater than the provided L.")
})

# normalizeCounts
test_that(desc = paste0("Making sure normalizeCounts doesn't change",
    " dimensions of counts matrix"), {
        normCounts <- normalizeCounts(celdaGSim$counts)
        expect_equal(dim(normCounts), dim(celdaGSim$counts))
        expect_equal(rownames(normCounts), rownames(celdaGSim$counts))
        expect_equal(colnames(normCounts), colnames(celdaGSim$counts))
        expect_error(normalizeCounts(celdaGSim$counts,
            transformationFun = "scale"),
            "'transformationFun' needs to be of class 'function'")
        expect_error(normalizeCounts(celdaGSim$counts, scaleFun = "scale"),
            "'scaleFun' needs to be of class 'function'")
    })

# recodeClusterY
test_that(desc = "Testing recodeClusterY with celda_G", {
    expect_error(recodeClusterZ(celdaMod = modelG,
        from = c(1, 2, 3, 4, 5),
        to = c(5, 4, 3, 2, 1)))
    expect_error(recodeClusterY(celdaMod = modelG,
        from = NULL,
        to =))
    expect_error(recodeClusterY(celdaMod = modelG,
        from = c(1, 2, 3, 4, 5),
        to = c(1, 2, 3, 4, 6)))
    expect_error(recodeClusterY(celdaMod = modelG,
        from = c(1, 2, 3, 4, 6),
        to = c(1, 2, 3, 4, 5)))
    newRecoded <- recodeClusterY(celdaMod = modelG,
        from = c(1, 2, 3, 4, 5),
        to = c(5, 4, 3, 2, 1))
    expect_equal(modelG@clusters$y == 1, newRecoded@clusters$y == 5)
})

# compareCountMatrix
test_that(desc = "Testing CompareCountMatrix with celda_G", {
    expect_true(compareCountMatrix(counts = celdaGSim$counts,
        celdaMod = modelG))
    lessFeatures <- celdaGSim$counts[1:50, ]
    expect_error(compareCountMatrix(counts = lessFeatures, celdaMod = modelG),
        paste0("The provided celda object was generated from a counts matrix",
            " with a different number of features than the one provided."))

    countsMatrixError <- matrix(data = 1,
        nrow = nrow(celdaGSim$counts),
        ncol = ncol(celdaGSim$counts))
    expect_false(compareCountMatrix(counts = countsMatrixError,
        celdaMod = modelG,
        errorOnMismatch = FALSE))
    expect_error(compareCountMatrix(counts = countsMatrixError,
        celdaMod = modelG,
        errorOnMismatch = TRUE))
})

# topRank
test_that(desc = "Testing topRank function with celda_G", {
    topRank <- topRank(matrix = factorized$proportions$module,
        n = 1000,
        threshold = NULL)
    expect_equal(names(topRank),
        c("index", "names"))
    expect_equal(names(topRank(matrix = factorized$proportions$module)),
        c("index", "names"))
})

# plotHeatmap
test_that(desc = "Testing plotHeatmap with celda_G", {
    expect_error(plotHeatmap(counts = celdaGSim$counts, y = modelG@params$L),
        "Length of y must match number of rows in counts matrix")
    expect_error(plotHeatmap(counts = celdaGSim$counts,
        y = modelG@clusters$y,
        scaleRow = "scale"),
        "'scaleRow' needs to be of class 'function'")
    expect_error(plotHeatmap(counts = celdaGSim$counts,
        y = modelG@clusters$y,
        trim = 3),
        paste0("'trim' should be a 2 element vector specifying the lower",
            " and upper boundaries"))
})

# test_that(desc = "Testing plotHeatmap with celda_G, including annotations", {
#     annot <- as.data.frame(c(rep(x = 1,
#         times = nrow(celdaGSim$counts) - 100), rep(x = 2, 100)))
#     rownames(annot) <- rownames(celdaGSim$counts)
#     colnames(annot) <- "label"
#
#     expect_equal(names(plotHeatmap(celdaMod = modelG,
#         counts = celdaGSim$counts,
#         annotationFeature = annot,
#         y = modelG@clusters$y)),
#         c("treeRow", "treeCol", "gtable"))
#
#     rownames(annot) <- NULL
#     expect_equal(names(plotHeatmap(celdaMod = modelG,
#         counts = celdaGSim$counts,
#         annotationFeature = as.matrix(annot),
#         y = modelG@clusters$y)),
#         c("treeRow", "treeCol", "gtable"))
#
#     rownames(annot) <- rev(rownames(celdaGSim$counts))
#     expect_error(plotHeatmap(celdaMod = modelG,
#         counts = celdaGSim$counts,
#         annotationFeature = annot,
#         y = modelG@clusters$y),
#         paste0("Row names of 'annotationFeature' are different than the row",
#             " names of 'counts'"))
# })

# celdaHeatmap
test_that(desc = "Testing celdaHeatmap with celda_G", {
    expect_equal(names(celdaHeatmap(celdaMod = modelG,
        counts = celdaGSim$counts)),
        c("treeRow", "treeCol", "gtable"))
})

# moduleHeatmap
test_that(desc = "Testing moduleHeatmap with celda_G", {
    expect_equal(names(moduleHeatmap(celdaGSim$counts,
        celdaMod = modelG,
        topCells = 300,
        featureModule = c(1, 2))),
        c("treeRow", "treeCol", "gtable"))
    expect_equal(names(moduleHeatmap(celdaGSim$counts,
        celdaMod = modelG,
        topFeatures = 15,
        topCells = 15,
        normalize = FALSE)),
        c("treeRow", "treeCol", "gtable"))

    expect_error(moduleHeatmap("counts", celdaMod = modelG),
        "'counts' should be a numeric count matrix")
    expect_error(moduleHeatmap(celdaGSim$counts, celdaMod = "modelG"),
        "'celdaMod' should be an object of class celda_G or celda_CG")
})

# plotDimReduceModule
test_that(desc = "Testing plotDimReduceModule with celda_G", {
    celdaTsne <- celdaTsne(counts = celdaGSim$counts,
        maxIter = 50,
        celdaMod = modelG)
    expect_equal(names(plotDimReduceModule(
        dim1 = celdaTsne[, 1],
        dim2 = celdaTsne[, 2],
        counts = celdaGSim$counts,
        celdaMod = modelG)),
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
        counts = celdaGSim$counts,
        celdaMod = modelG,
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
    expect_error(plotDimReduceModule(dim1 = celdaTsne[, 1],
        dim2 = celdaTsne[, 2],
        counts = celdaGSim$counts,
        celdaMod = modelG,
        modules = c(11, 12)))
})

# celdaTsne
test_that(desc = paste0("Testing celdaTsne with celda_G when model class is",
    " changed, should error"), {
        modelX <- modelG
        class(modelX) <- "celda_X"
        expect_error(celdaTsne(counts = celdaGSim$counts, celdaMod = modelX),
            "unable to find")
    })

# test_that(desc = "Testing celdaTsne with celda_C including all cells", {
#     tsne <- celdaTsne(counts = celdaGSim$counts,
#         celdaMod = modelG,
#         maxCells = ncol(celdaGSim$counts))
#     plotObj <- plotDimReduceCluster(tsne[, 1],
#         tsne[, 2], rep(1, ncol(celdaGSim$counts)))
#     expect_true(ncol(tsne) == 2 & nrow(tsne) == ncol(celdaGSim$counts))
#     expect_true(!is.null(plotObj))
#
#     tsne <- celdaTsne(counts = celdaGSim$counts,
#         celdaMod = modelG,
#         maxCells = ncol(celdaGSim$counts),
#         modules = c(1, 2))
#     expect_error(tsne <- celdaTsne(counts = celdaGSim$counts,
#         celdaMod = modelG,
#         maxCells = ncol(celdaGSim$counts),
#         modules = seq(1000, 1005)))
# })

test_that(desc = "Testing celdaTsne with celda_G including a subset of cells", {
    tsne <- celdaTsne(counts = celdaGSim$counts,
        celdaMod = modelG,
        maxCells = 100)
    plotObj <- plotDimReduceCluster(tsne[, 1],
        tsne[, 2], rep(1, ncol(celdaGSim$counts)))
    expect_true(ncol(tsne) == 2 & nrow(tsne) == ncol(celdaGSim$counts) &&
            sum(!is.na(tsne[, 1])) == 100)
    expect_true(!is.null(plotObj))
})

# celdaUmap
test_that(desc = paste0("Testing celdaUmap with celda_G when model class is",
    " changed, should error"), {
        modelX <- modelG
        class(modelX) <- "celda_X"
        expect_error(celdaUmap(counts = celdaGSim$counts, celdaMod = modelX),
            "unable to find")
    })

test_that(desc = "Testing celdaUmap with celda_C including all cells", {
    umap <- celdaUmap(counts = celdaGSim$counts,
        celdaMod = modelG,
        maxCells = ncol(celdaGSim$counts))
    plotObj <- plotDimReduceCluster(umap[, 1],
        umap[, 2], rep(1, ncol(celdaGSim$counts)))
    expect_true(ncol(umap) == 2 & nrow(umap) == ncol(celdaGSim$counts))
    expect_true(!is.null(plotObj))

    umap <- celdaUmap(counts = celdaGSim$counts,
        celdaMod = modelG,
        maxCells = ncol(celdaGSim$counts),
        modules = c(1, 2))
    expect_error(umap <- celdaUmap(counts = celdaGSim$counts,
        celdaMod = modelG,
        maxCells = ncol(celdaGSim$counts),
        modules = seq(1000, 1005)))
})

test_that(desc = "Testing celdaUmap with celda_G including a subset of cells", {
    umap <- celdaUmap(counts = celdaGSim$counts,
        celdaMod = modelG,
        maxCells = 100)
    plotObj <- plotDimReduceCluster(umap[, 1],
        umap[, 2], rep(1, ncol(celdaGSim$counts)))
    expect_true(ncol(umap) == 2 & nrow(umap) == ncol(celdaGSim$counts) &&
            sum(!is.na(umap[, 1])) == 100)
    expect_true(!is.null(plotObj))
})

# featureModuleLookup
test_that(desc = "Testing featureModuleLookup with celda_G", {
    res <- featureModuleLookup(celdaGSim$counts, modelG, "Gene_1")
    expect_true(res == modelG@clusters$y[1])
    res <- featureModuleLookup(celdaGSim$counts,
        modelG, "Gene_2", exactMatch = FALSE)
    expect_true(length(res) == 11)
    res <- featureModuleLookup(celdaGSim$counts, modelG, "XXXXXXX")
    expect_true(grepl("No feature", res))
})

# .cGSplitY
test_that(desc = "Testing error checking for .cGSplitY", {
    r <- simulateCells("celda_G",
        C = 100,
        G = 100,
        L = 2)
    dc <- .cGDecomposeCounts(r$counts, r$y, r$L)
    res <- .cGSplitY(r$counts,
        r$y,
        dc$nTSByC,
        dc$nByTS,
        dc$nByG,
        dc$nGByTS,
        dc$nM,
        dc$nG,
        r$L,
        beta = 1,
        delta = 1,
        gamma = 1,
        yProb = NULL,
        minFeature = 1000)
    expect_true(grepl("Cluster sizes too small", res$message))
})

test_that(desc = "Testing perplexity.celda_G", {
    expect_true(is.numeric(perplexity(celdaGSim$counts, modelG)))
    class(modelG) <- c("celda_C")
    expect_error(perplexity.celda_G(celdaGSim$counts, modelG),
        "could not find function \"perplexity.celda_G\"")
})
