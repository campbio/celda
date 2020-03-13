## ---- eval= FALSE-------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly=TRUE)){
#      install.packages("BiocManager")}
#  BiocManager::install("celda")

## ---- eval = TRUE, message = FALSE--------------------------------------------
library(celda)

## ---- eval = FALSE------------------------------------------------------------
#  help(package = celda)

## -----------------------------------------------------------------------------
simCounts <- simulateCells("celda_CG",
    K = 5, L = 10, S = 5, G = 200, CRange = c(30, 50))

## -----------------------------------------------------------------------------
dim(simCounts$counts)

## -----------------------------------------------------------------------------
table(simCounts$z)

## -----------------------------------------------------------------------------
table(simCounts$y)

## -----------------------------------------------------------------------------
table(simCounts$sampleLabel)

## ---- warning = FALSE, message = FALSE----------------------------------------
celdaModel <- celda_CG(counts = simCounts$counts,
    K = 5, L = 10, verbose = FALSE)

## -----------------------------------------------------------------------------
table(clusters(celdaModel)$z, simCounts$z)
table(clusters(celdaModel)$y, simCounts$y)

## -----------------------------------------------------------------------------
factorized <- factorizeMatrix(counts = simCounts$counts, celdaMod = celdaModel)
names(factorized)

## -----------------------------------------------------------------------------
dim(factorized$proportions$cell)
head(factorized$proportions$cell[, seq(3)], 5)

## -----------------------------------------------------------------------------
cellPop <- factorized$proportions$cellPopulation
dim(cellPop)
head(cellPop, 5)

## -----------------------------------------------------------------------------
dim(factorized$proportions$module)
head(factorized$proportions$module, 5)

## -----------------------------------------------------------------------------
topGenes <- topRank(matrix = factorized$proportions$module,
    n = 10, threshold = NULL)

## -----------------------------------------------------------------------------
topGenes$names$L1

## ---- eval = TRUE, fig.width = 7, fig.height = 7------------------------------
celdaHeatmap(counts = simCounts$counts, celdaMod = celdaModel, nfeatures = 10)

## -----------------------------------------------------------------------------
tsne <- celdaTsne(counts = simCounts$counts, celdaMod = celdaModel)

## ---- eval = TRUE, fig.width = 7, fig.height = 7------------------------------
plotDimReduceCluster(dim1 = tsne[, 1],
    dim2 = tsne[, 2],
    cluster = clusters(celdaModel)$z)

plotDimReduceModule(dim1 = tsne[, 1],
    dim2 = tsne[, 2],
    celdaMod = celdaModel,
    counts = simCounts$counts,
    rescale = TRUE)

plotDimReduceFeature(dim1 = tsne[, 1],
    dim2 = tsne[, 2],
    counts = simCounts$counts,
    features = "Gene_1")

## ---- eval = TRUE, fig.width = 7, fig.height = 7------------------------------
celdaProbabilityMap(counts = simCounts$counts, celdaMod = celdaModel)

## ---- eval = TRUE, fig.width = 7, fig.height = 7------------------------------
moduleHeatmap(counts = simCounts$counts, celdaMod = celdaModel,
    featureModule = 1, topCells = 100)

## -----------------------------------------------------------------------------
genes <- c(topGenes$names$L1, topGenes$names$L10)
geneIx <- which(rownames(simCounts$counts) %in% genes)
normCounts <- normalizeCounts(counts = simCounts$counts, scaleFun = scale)

## ---- fig.width = 8, fig.height = 8-------------------------------------------
plotHeatmap(counts = normCounts,
    z = clusters(celdaModel)$z,
    y = clusters(celdaModel)$y,
    featureIx = geneIx,
    showNamesFeature = TRUE)

## ---- message=FALSE-----------------------------------------------------------
diffexpClust1 <- differentialExpression(counts = simCounts$counts,
    celdaMod = celdaModel,
    c1 = 1,
    c2 = NULL)

head(diffexpClust1, 5)

## ---- message=FALSE-----------------------------------------------------------
diffexpClust1vs2 <- differentialExpression(
    counts = simCounts$counts,
    celdaMod = celdaModel,
    c1 = 1,
    c2 = 2)

diffexpClust1vs2 <- diffexpClust1vs2[diffexpClust1vs2$FDR < 0.05 &
    abs(diffexpClust1vs2$Log2_FC) > 2, ]
head(diffexpClust1vs2, 5)

## -----------------------------------------------------------------------------
diffexpGeneIx <- which(rownames(simCounts$counts) %in% diffexpClust1vs2$Gene)

normCounts <- normalizeCounts(counts = simCounts$counts, scaleFun = scale)

## -----------------------------------------------------------------------------
plotHeatmap(counts = normCounts[, clusters(celdaModel)$z %in% c(1, 2)],
    clusterCell = TRUE,
    z = clusters(celdaModel)$z[clusters(celdaModel)$z %in% c(1, 2)],
    y = clusters(celdaModel)$y,
    featureIx = diffexpGeneIx,
    showNamesFeature = TRUE)

## ---- message = FALSE---------------------------------------------------------
moduleSplit <- recursiveSplitModule(counts = simCounts$counts,
    initialL = 2, maxL = 15)

## -----------------------------------------------------------------------------
plotGridSearchPerplexity(celdaList = moduleSplit)

## ---- message = FALSE---------------------------------------------------------
moduleSplitSelect <- subsetCeldaList(moduleSplit, params = list(L = 10))

cellSplit <- recursiveSplitCell(counts = simCounts$counts,
    initialK = 3,
    maxK = 12,
    yInit = clusters(moduleSplitSelect)$y)

## -----------------------------------------------------------------------------
plotGridSearchPerplexity(celdaList = cellSplit)

## ---- eval = TRUE-------------------------------------------------------------
celdaModel <- subsetCeldaList(celdaList = cellSplit,
    params = list(K = 5, L = 10))

## ---- eval = TRUE, message = FALSE--------------------------------------------
cgs <- celdaGridSearch(simCounts$counts,
    paramsTest = list(K = seq(4, 6), L = seq(9, 11)),
    cores = 1,
    model = "celda_CG",
    nchains = 2,
    maxIter = 100,
    verbose = FALSE,
    bestOnly = TRUE)

## ---- eval = TRUE-------------------------------------------------------------
cgs <- resamplePerplexity(counts = simCounts$counts,
    celdaList = cgs, resample = 5)

## ---- eval = TRUE, fig.width = 8, fig.height = 8, warning = FALSE, message = FALSE----
plotGridSearchPerplexity(celdaList = cgs)

## ---- eval = TRUE-------------------------------------------------------------
celdaModel <- subsetCeldaList(celdaList = cgs, params = list(K = 5, L = 10))

## ---- message=FALSE-----------------------------------------------------------
cgs <- celdaGridSearch(simCounts$counts,
    paramsTest = list(K = seq(4, 6), L = seq(9, 11)),
    cores = 1,
    model = "celda_CG",
    nchains = 2,
    maxIter = 100,
    verbose = FALSE,
    bestOnly = FALSE)

cgs <- resamplePerplexity(counts = simCounts$counts,
    celdaList = cgs,
    resample = 2)

cgsK5L10 <- subsetCeldaList(celdaList = cgs, params = list(K = 5, L = 10))

celdaModel1 <- selectBestModel(celdaList = cgsK5L10)

## -----------------------------------------------------------------------------
featureModuleLookup(counts = simCounts$counts, celdaMod = celdaModel,
    feature = c("Gene_99"))

## -----------------------------------------------------------------------------
celdaModelZRecoded <- recodeClusterZ(celdaMod = celdaModel,
    from = c(1, 2, 3, 4, 5), to = c(2, 1, 3, 4, 5))

## -----------------------------------------------------------------------------
table(clusters(celdaModel)$z, clusters(celdaModelZRecoded)$z)

## -----------------------------------------------------------------------------
sessionInfo()

