# celda_C
library(celda)
context("Testing celda_C")

celdaCCim <- simulateCells("celda_C", K = 10)
modelC <- celda_C(counts = celdaCCim$counts,
    sampleLabel = celdaCCim$sampleLabel,
    K = celdaCCim$K,
    algorithm = "EM",
    verbose = FALSE)
factorized <- factorizeMatrix(counts = celdaCCim$counts, celdaMod = modelC)

# celda_C
test_that(desc = "Testing simulation and celda_C model", {
    expect_equal(typeof(celdaCCim$counts), "integer")
    expect_true(all(sweep(factorized$counts$sample,
        2,
        colSums(factorized$counts$sample),
        "/") == factorized$proportions$sample))
    expect_true(ncol(factorized$proportions$module) == modelC@params$K)
    expect_true(all(is.numeric(logLikelihoodHistory(celdaMod = modelC))))
    expect_equal(max(logLikelihoodHistory(celdaMod = modelC)),
        bestLogLikelihood(modelC))

    # GitHub #347
    numericCounts <- celdaCCim$counts
    storage.mode(numericCounts) <- "numeric"
    expect_true(is(celda_C(counts = celdaCCim$counts,
        sampleLabel = celdaCCim$sampleLabel,
        K = celdaCCim$K,
        algorithm = "EM",
        verbose = FALSE,
        maxIter = 2),
        "celda_C"))
})

# clusterProbability
test_that(desc = "Testing clusterProbability with celda_C", {
    expect_true(all(round(rowSums(
        clusterProbability(modelC, counts = celdaCCim$counts)[[1]]
    )) == 1))
})

# celdaGridSearch and perplexity calculations
test_that(desc = "Testing celdaGridSearch with celda_C", {
    celdaCRes <- celdaGridSearch(counts = celdaCCim$counts,
        model = "celda_C",
        nchains = 2,
        paramsTest = list(K = c(5, 6)),
        paramsFixed = list(sampleLabel = celdaCCim$sampleLabel),
        maxIter = 2,
        verbose = FALSE,
        bestOnly = FALSE,
        perplexity = FALSE)
    expect_error(celdaGridSearch(counts = celdaCCim$counts,
        model = "celda_C",
        paramsTest = list(K = c(4, 5), M = c(3, 4)),
        paramsFixed = list(sampleLabel = celdaCCim$sampleLabel),
        bestOnly = FALSE),
        paste0("The following elements in 'paramsTest' are not arguments of",
            " 'celda_C': M"))

    expect_error(celdaGridSearch(counts = celdaCCim$counts,
        model = "celda_C",
        nchains = 1,
        maxIter = 1,
        paramsTest = list(K = c(4, 5), sampleLabel = "Sample"),
        paramsFixed = list(sampleLabel = celdaCCim$sampleLabel)),
        paste0("Setting parameters such as 'z.init', 'y.init', and",
            "'sampleLabel' in 'paramsTest' is not currently supported."))

    expect_error(celdaGridSearch(counts = celdaCCim$counts,
        model = "celda_C",
        nchains = 1,
        maxIter = 1,
        paramsTest = list(),
        paramsFixed = list(sampleLabel = celdaCCim$sampleLabel)),
        paste0("The following arguments are not in 'paramsTest' or",
            " 'paramsFixed' but are required for 'celda_C': K"))

    expect_error(celdaGridSearch(counts = celdaCCim$counts,
        model = "celda_C",
        nchains = 1,
        maxIter = 1,
        paramsTest = list(K = c(9, 10)),
        paramsFixed = list(sampleLabel = celdaCCim$sampleLabel,
            xxx = "xxx")),
        paste0("The following elements in 'paramsFixed' are not arguments",
            "of 'celda_C': xxx"))

    expect_true(class(celdaCRes)[1] == "celdaList")
    expect_equal(names(runParams(celdaCRes)),
        c("index", "chain", "K", "log_likelihood"))
    expect_error(plotGridSearchPerplexity(celdaCRes))

    celdaCRes <- resamplePerplexity(celdaCCim$counts, celdaCRes, resample = 2)
    expect_equal(is.null(celdaCRes@perplexity), FALSE)
    expect_true(is(celdaCRes, "celdaList"))
    expect_error(resamplePerplexity(celdaCCim$counts,
        celdaCRes, resample = "2"))
    expect_error(resamplePerplexity(celdaCCim$counts,
        "celdaCRes", resample = 2))

    plotObj <- plotGridSearchPerplexity(celdaCRes)
    expect_is(plotObj, "ggplot")

    celdaCResIndex1 <- subsetCeldaList(celdaCRes, params = list(index = 1))
    expect_true(all(is(celdaCResIndex1, "celda_C") &&
            !is(celdaCResIndex1, "celda_list")))

    expect_error(subsetCeldaList(celdaCRes, params = list(K = 11)))
    expect_error(subsetCeldaList(celdaCRes, params = list(K = 5, M = 10)))

    celdaCResK5 <- subsetCeldaList(celdaCRes, params = list(K = 5))
    modelC2 <- selectBestModel(celdaCResK5)
    res <- perplexity(celdaCCim$counts, modelC)
    res2 <- perplexity(celdaCCim$counts, modelC,
        newCounts = celdaCCim$counts + 1)

    expect_error(res <- perplexity.celda_C(celdaCCim$counts,
        modelC,
        newCounts = celdaCCim$counts[-1, ]))
})

# logLikelihood
test_that(desc = "Testing logLikelihood.celda_C", {
    expect_lt(logLikelihood(model = "celda_C",
        counts = celdaCCim$counts,
        z = celdaCCim$z,
        K = celdaCCim$K,
        alpha = 1,
        beta = 1,
        sampleLabel = celdaCCim$sampleLabel),
        0)

    fakeZ <- celdaCCim$z
    fakeZ[1] <- celdaCCim$K + 1
    expect_error(logLikelihood(model = "celda_C",
        z = fakeZ,
        counts = celdaCCim$counts,
        K = celdaCCim$K,
        alpha = 1,
        beta = 1,
        sampleLabel = celdaCCim$sampleLabel),
        "An entry in z contains a value greater than the provided K.")
})

# Gibbs sampling
test_that(desc = "Testing celda_C with Gibbs sampling", {
    res <- celda_C(counts = celdaCCim$counts,
        sampleLabel = celdaCCim$sampleLabel,
        K = celdaCCim$K,
        algorithm = "Gibbs",
        maxIter = 5,
        nchain = 1)
    expect_is(res, "celda_C")
})


# normalizeCounts
test_that(desc = paste0("Making sure normalizeCounts doesn't change dimensions",
    " of counts matrix"), {
        normCounts <- normalizeCounts(celdaCCim$counts)
        expect_equal(dim(normCounts), dim(celdaCCim$counts))
        expect_equal(rownames(normCounts), rownames(celdaCCim$counts))
        expect_equal(colnames(normCounts), colnames(celdaCCim$counts))
        expect_error(normalizeCounts(celdaCCim$counts,
            transformation.fun = "scale"),
            "'transformation.fun' needs to be of class 'function'")
        expect_error(normalizeCounts(celdaCCim$counts, scale.fun = "scale"),
            "'scale.fun' needs to be of class 'function'")
    })

# recodeClusterZ
test_that(desc = "Testing recodeClusterZ with celda_C", {
    expect_error(recodeClusterY(celdaMod = modelC,
        from = c(1, 2, 3, 4, 5),
        to = c(5, 4, 3, 2, 1)))
    expect_error(recodeClusterZ(celdaMod = modelC,
        from = NULL,
        to = ))
    expect_error(recodeClusterZ(celdaMod = modelC,
        from = c(1, 2, 3, 4, 5),
        to = c(1, 2, 3, 4, 6)))
    expect_error(recodeClusterZ(celdaMod = modelC,
        from = c(1, 2, 3, 4, 6),
        to = c(1, 2, 3, 4, 5)))
    newRecoded <- recodeClusterZ(celdaMod = modelC,
        from = c(1, 2, 3, 4, 5),
        to = c(5, 4, 3, 2, 1))
    expect_equal(modelC@clusters$z == 1, newRecoded@clusters$z == 5)
})

# compareCountMatrix
test_that(desc = "Testing CompareCountMatrix with celda_C", {
    expect_true(compareCountMatrix(counts = celdaCCim$counts,
        celdaMod = modelC))

    lessCells <- celdaCCim$counts[, seq(100)]
    expect_error(compareCountMatrix(counts = lessCells, celdaMod = modelC),
        paste0("The provided celda object was generated from a counts matrix",
            " with a different number of cells than the one provided."))

    countsMatrixError <- matrix(data = 1,
        nrow = nrow(celdaCCim$counts),
        ncol = ncol(celdaCCim$counts))
    expect_false(compareCountMatrix(counts = countsMatrixError,
        celdaMod = modelC,
        errorOnMismatch = FALSE))
    expect_error(compareCountMatrix(counts = countsMatrixError,
        celdaMod = modelC,
        errorOnMismatch = TRUE))
})

# topRank
test_that(desc = "Checking topRank to see if it runs without errors", {
    topRank <- topRank(matrix = factorized$proportions$module, threshold = NULL)
    expect_equal(names(topRank), c("index", "names"))
    topRank <- topRank(matrix = factorized$proportions$module, n = 1000)
})

# plotHeatmap
test_that(desc = "Testing plotHeatmap with celda_C", {
    expect_error(plotHeatmap(counts = celdaCCim$counts, z = modelC@params$K),
        "Length of z must match number of columns in counts matrix")
    expect_error(plotHeatmap(
        counts = celdaCCim$counts,
        z = modelC@clusters$z,
        scaleRow = modelC),
        "'scaleRow' needs to be of class 'function'")
    expect_error(plotHeatmap(counts = celdaCCim$counts,
        z = modelC@clusters$z,
        trim = 3),
        paste0("'trim' should be a 2 element vector specifying the lower and",
            " upper boundaries"))
})


# plotHeatmap with annotationCell
test_that(desc = "Testing plotHeatmap with celda_C, including annotations", {
    annot <- as.data.frame(c(rep(x = 1,
        times = ncol(celdaCCim$counts) - 100),
        rep(x = 2, 100)))

    rownames(annot) <- colnames(celdaCCim$counts)
    expect_equal(names(plotHeatmap(
        celdaMod = modelC,
        counts = celdaCCim$counts,
        annotationCell = annot,
        z = modelC@clusters$z)),
        c("tree_row", "tree_col", "gtable"))

    rownames(annot) <- NULL
    expect_equal(names(plotHeatmap(celdaMod = modelC,
        counts = celdaCCim$counts,
        annotation.feature = as.matrix(annot),
        z = modelC@clusters$z)),
        c("tree_row", "tree_col", "gtable"))

    rownames(annot) <- rev(colnames(celdaCCim$counts))
    expect_error(plotHeatmap(celdaMod = modelC,
        counts = celdaCCim$counts,
        annotationCell = annot,
        z = modelC@clusters$z),
        paste0("Row names of 'annotationCell' are different than the column",
            " names of 'counts'")
    )
})


# celdaHeatmap
test_that(desc = "Testing celdaHeatmap with celda_C", {
    expect_equal(names(celdaHeatmap(celdaMod = modelC,
        counts = celdaCCim$counts)),
        c("tree_row", "tree_col", "gtable"))
})


# celdaProbabilityMap
test_that(desc = "Testing celdaProbabiltyMap.celda_C for sample level", {
    plotObj <- celdaProbabilityMap(counts = celdaCCim$counts,
        celdaMod = modelC,
        level = "sample")
    expect_true(!is.null(plotObj))

    ## Without a sample label
    modelC <- celda_C(celdaCCim$counts,
        sampleLabel = NULL,
        K = celdaCCim$K,
        maxIter = 5,
        nchain = 1)
    plotObj <- celdaProbabilityMap(counts = celdaCCim$counts,
        celdaMod = modelC,
        level = "sample")
    expect_true(!is.null(plotObj))
})


# differentialExpression
test_that(desc = "Testing differentialExpression with celda_C", {
    diffexp_K1 <- differentialExpression(counts = celdaCCim$counts,
        celdaMod = modelC,
        c1 = 1)
    expect_equal(class(diffexp_K1), c("data.table", "data.frame"))
    expect_equal(class(diffExp_K1 <- differentialExpression(
        counts = celdaCCim$counts,
        celdaMod = modelC,
        c1 = c(2, 3),
        c2 = 4,
        log2fcThreshold = 0.5)),
        c("data.table", "data.frame"))
    expect_error(differentialExpression(counts = "counts",
        celdaMod = modelC,
        c1 = 3,
        log2fcThreshold = 0.5),
        "'counts' should be a numeric count matrix")
    expect_error(differentialExpression(
        counts = celdaCCim$counts,
        celdaMod = NULL,
        c1 = 3),
        "'celdaMod' should be an object of class celda_C or celda_CG")
    expect_error(differentialExpression(counts = celdaCCim$counts,
        celdaMod = modelC,
        c1 = NULL,
        log2fcThreshold = 0.5,
        only.pos = TRUE))
})


test_that(desc = paste0("Testing celdaTsne with celda_C when model class is",
    "changed, should error"), {
        modelX <- modelC
        class(modelX) <- "celda_X"
        expect_error(celdaTsne(counts = celdaCCim$counts,
            celdaMod = modelX,
            maxCells = length(modelC@clusters$z),
            minClusterSize = 10),
            "unable to find")
    })

test_that(desc = "Testing celdaTsne with celda_C including all cells", {
    tsne <- celdaTsne(counts = celdaCCim$counts,
        celdaMod = modelC,
        maxCells = length(modelC@clusters$z),
        minClusterSize = 10)
    plotObj <- plotDimReduceCluster(tsne[, 1],
        tsne[, 2],
        modelC@clusters$z,
        labelClusters = TRUE)
    expect_true(ncol(tsne) == 2 & nrow(tsne) == length(modelC@clusters$z))
    expect_true(!is.null(plotObj))
})


test_that(desc = "Testing celdaTsne with celda_C including a subset of cells", {
    expect_success(expect_error(tsne <- celdaTsne(counts = celdaCCim$counts,
        celdaMod = modelC,
        maxCells = 50,
        minClusterSize = 50)))
    tsne <- celdaTsne(counts = celdaCCim$counts,
        celdaMod = modelC,
        maxCells = 100,
        minClusterSize = 10)
    plotObj <- plotDimReduceCluster(tsne[, 1], tsne[, 2], modelC@clusters$z)
    expect_true(ncol(tsne) == 2 &
            nrow(tsne) == length(modelC@clusters$z) &&
            sum(!is.na(tsne[, 1])) == 100)
    expect_true(!is.null(plotObj))
})

test_that(desc = paste0("Testing celdaUmap with celda_C when model class is",
    " changed, should error"), {
        modelX <- modelC
        class(modelX) <- "celda_X"
        expect_error(celdaUmap(counts = celdaCCim$counts,
            celdaMod = modelX,
            maxCells = length(modelC@clusters$z),
            minClusterSize = 10),
            "unable to find")
    })

test_that(desc = "Testing celdaUmap with celda_C including all cells", {
    umap <- celdaUmap(counts = celdaCCim$counts,
        celdaMod = modelC,
        maxCells = length(modelC@clusters$z),
        minClusterSize = 10)
    plotObj <- plotDimReduceCluster(umap[, 1], umap[, 2], modelC@clusters$z)
    expect_true(ncol(umap) == 2 &
            nrow(umap) == length(modelC@clusters$z))
    expect_true(!is.null(plotObj))
})

test_that(desc = "Testing celdaUmap with celda_C including a subset of cells", {
    expect_success(expect_error(umap <- celdaUmap(
        counts = celdaCCim$counts,
        celdaMod = modelC,
        maxCells = 50,
        minClusterSize = 50)))
    umap <- celdaUmap(counts = celdaCCim$counts,
        celdaMod = modelC,
        maxCells = 100,
        minClusterSize = 10)
    plotObj <- plotDimReduceCluster(umap[, 1], umap[, 2], modelC@clusters$z)
    expect_true(ncol(umap) == 2 &
            nrow(umap) == length(modelC@clusters$z) &&
            sum(!is.na(umap[, 1])) == 100)
    expect_true(!is.null(plotObj))
})


# featureModuleLookup
test_that(desc = "Testing featureModuleLookup with celda_C", {
    expect_error(featureModuleLookup(celdaCCim$counts, modelC, "test_feat"))
})


# .cCSplitZ
test_that(desc = "Testing error checking for .cCSplitZ", {
    r <- simulateCells("celda_C",
        S = 1,
        CRange = c(50, 100),
        K = 2)
    dc <- .cCDecomposeCounts(r$counts, r$sampleLabel, r$z, r$K)
    res <- .cCSplitZ(r$counts,
        dc$mCPByS,
        dc$nGByCP,
        dc$nCP,
        s = as.integer(r$sampleLabel),
        z = r$z,
        K = r$K,
        nS = dc$nS,
        nG = dc$nG,
        alpha = 1,
        beta = 1,
        zProb = NULL,
        min.cell = 1000)
    expect_true(grepl("Cluster sizes too small", res$message))
})


test_that(desc = "Testing perplexity.celda_C", {
    expect_true(is.numeric(perplexity(celdaCCim$counts, modelC)))

    class(modelC) <- c("celda_CG")
    expect_error(perplexity.celda_C(celdaCCim$counts, modelC),
        "could not find function \"perplexity.celda_C\"")
})
