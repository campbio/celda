# celda_CG
library(celda)
context("Testing celda_CG")

K <- 5
L <- 10
celdaCGSim <- simulateCells("celda_CG", K = K, L = L)
counts(celdaCGSim) <- as(counts(celdaCGSim), "dgCMatrix")
celdaCGSim <- selectFeatures(celdaCGSim, 0, 0)
modelCG <- celda_CG(celdaCGSim,
    sampleLabel = celdaCGSim$celda_sample_label,
    K = K,
    L = L,
    zInitialize = "random",
    yInitialize = "random",
    nchains = 1,
    algorithm = "EM",
    verbose = FALSE)
factorized <- factorizeMatrix(modelCG)

# celda_CG
test_that(desc = "Testing simulation and celda_CG model", {
    expect_true(all(sweep(factorized$counts$cell,
        2,
        colSums(counts(altExp(celdaCGSim))),
        "/") == factorized$proportions$cell))
    expect_equal(K,
        ncol(factorized$proportions$cellPopulation))
    expect_equal(L,
        nrow(factorized$proportions$cellPopulation))
})

# Cluster probabilities
test_that(desc = "Testing clusterProbability with celda_CG", {
    clustProb <- clusterProbability(modelCG)
    expect_true(all(round(rowSums(clustProb$zProbability), 10) == 1) &
            nrow(clustProb$zProbability) == ncol(modelCG))
    expect_true(all(round(rowSums(clustProb$yProbability), 10) == 1) &
            nrow(clustProb$yProbability) == nrow(modelCG))
    clustProb <- clusterProbability(modelCG, log = TRUE)
})

test_that(desc = paste0("Testing simulateCells celda_CG error checking with",
    " low gamma"), {
        expect_warning(simulateCells(model = "celda_CG", gamma = 0.1))
    })

test_that(desc = paste0("Testing simulateCells celda_CG, make sure all genes",
    " expressed"), {
        simCellsLow <- simulateCells(model = "celda_CG",
            G = 1000,
            C = 300,
            CRange = c(1, 100),
            NRange = c(1, 100))
        expect_true(all(rowSums(counts(simCellsLow)) > 0))
    })


# Ensure logLikelihood calculates the expected values
test_that(desc = "Testing logLikelihood functions for celda_CG", {
    expect_lt(logLikelihood(modelCG), 0)

    fakeZ <- as.integer(celdaClusters(modelCG))
    fakeZ[1] <- K + 1
    expect_error(.logLikelihoodcelda_CG(
      K = K,
      L = L,
      y = as.integer(celdaModules(modelCG)),
      z = fakeZ,
      delta = 1,
      gamma = 1,
      beta = 1,
      alpha = 1,
      s = modelCG$celda_sample_label,
      counts = counts(altExp(modelCG))),
        paste0("Assigned value of cell cluster greater than the total ",
        "number of cell clusters!"))

    fakeY <- as.integer(celdaModules(modelCG))
    fakeY[1] <- L + 1
    expect_error(.logLikelihoodcelda_CG(
        y = fakeY,
        z = as.integer(celdaClusters(modelCG)),
        delta = 1,
        gamma = 1,
        beta = 1,
        alpha = 1,
        K = K,
        L = L,
        s = modelCG$celda_sample_label,
        counts = counts(altExp(modelCG))),
        paste0("Assigned value of feature module greater than the total ",
        "number of feature modules!"))
})

# normalizeCounts
test_that(desc = "Testing normalizeCounts with celda_CG", {
    normCounts <- normalizeCounts(counts(celdaCGSim))
    expect_equal(dim(normCounts), dim(counts(celdaCGSim)))
    expect_equal(rownames(normCounts), rownames(counts(celdaCGSim)))
    expect_equal(colnames(normCounts), colnames(counts(celdaCGSim)))
    expect_error(normalizeCounts(counts(celdaCGSim),
        transformationFun = "scale"),
        "'transformationFun' needs to be of class 'function'")
    expect_error(normalizeCounts(counts(celdaCGSim), scaleFun = "scale"),
        "'scaleFun' needs to be of class 'function'")
})

# recodeClusterY
test_that(desc = "Testing recodeClusterY with celda_CG", {
    expect_error(recodeClusterY(modelCG,
        from = NULL,
        to = ""))
    expect_error(recodeClusterY(modelCG,
        from = c(1, 2, 3, 4, 5),
        to = c(1, 2, 4, 3, 6)))
    expect_error(recodeClusterY(modelCG,
        from = c(1, 2, 3, 4, 6),
        to = c(1, 2, 4, 3, 5)))
    newRecoded <- recodeClusterY(modelCG,
        from = c(1, 2, 3, 4, 5),
        to = c(3, 2, 1, 4, 5))
    expect_equal(celdaModules(modelCG) == 1, celdaModules(newRecoded) == 3)
})

# recodeClusterZ
test_that(desc = "Testing recodeClusterZ with celda_CG", {
    expect_error(recodeClusterZ(modelCG,
        from = NULL,
        to = ""))
    expect_error(recodeClusterZ(modelCG,
        from = c(1, 2, 3, 4, 5),
        to = c(1, 2, 3, 4, 6)))
    expect_error(recodeClusterZ(modelCG,
        from = c(1, 2, 3, 4, 6),
        to = c(1, 2, 3, 4, 5)))
    newRecoded <- recodeClusterZ(modelCG,
        from = c(1, 2, 3, 4, 5),
        to = c(5, 4, 3, 2, 1))
    expect_equal(celdaClusters(modelCG) == 1, celdaClusters(newRecoded) == 5)
})


# topRank
test_that(desc = "Testing topRank with celda_CG", {
    topRank <- topRank(matrix = factorized$proportions$module,
        n = 1000,
        threshold = NULL)
    expect_equal(names(topRank), c("index", "names"))
})

# celdaHeatmap
test_that(desc = "Testing celdaHeatmap with celda_CG", {
  plt <- celdaHeatmap(modelCG)
  expect_equal(class(plt), c("gtable", "gTree", "grob", "gDesc"))
})

# moduleHeatmap
test_that(desc = "Testing moduleHeatmap with celda_CG", {
    plt <- moduleHeatmap(modelCG, featureModule = c(2, 3),
                         topCells = 10, topFeatures = 10)
    expect_is(plt, "list")
})

# celdaProbabiltyMap
test_that(desc = "Testing celdaProbabiltyMap", {
  plotObj <- celdaProbabilityMap(modelCG)
  plotObj <- celdaProbabilityMap(modelCG, level = "cellPopulation")
  expect_true(!is.null(plotObj))
})


test_that(desc = "Testing celdaUmap and celdaTsne with celda_CG", {
  modelCG <- celdaUmap(modelCG, maxCells = 100, minClusterSize = 10)
  modelCG <- celdaTsne(modelCG, maxCells = 100, minClusterSize = 10)
  plotObj <- plotDimReduceCluster(modelCG, "celda_UMAP")
  expect_true(!is.null(plotObj))
})

# featureModuleLookup
test_that(desc = "Testing featureModuleLookup with celda_CG", {
    res <- featureModuleLookup(modelCG, "Gene_1")
    expect_true(res == celdaModules(modelCG)[1])

    expect_error(featureModuleLookup(modelCG, "XXXXXXX"))
})


test_that(desc = "Testing perplexity of celda_CG", {
    expect_true(is.numeric(perplexity(modelCG)))
})

test_that(desc = "Testing featureModuleTable", {
    table <- featureModuleTable(modelCG,  outputFile = NULL)
    expect_equal(ncol(table), 10)
})

test_that(desc = "Testing plotCeldaViolin with celda_CG", {
    violin <- plotCeldaViolin(modelCG, features = "Gene_1")
    expect_is(violin, "ggplot")
})
