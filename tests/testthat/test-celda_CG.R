#celda_CG
library(celda)
context("Testing celda_CG")

celdaCG.sim = simulateCells("celda_CG", K=5, L=10)
model_CG = celda_CG(counts=celdaCG.sim$counts, sample.label=celdaCG.sim$sample.label, K=celdaCG.sim$K, L=celdaCG.sim$L, algorithm="EM", verbose=FALSE)
factorized <- factorizeMatrix(celda.mod = model_CG, counts = celdaCG.sim$counts)

# celda_CG
test_that(desc = "Testing simulation and celda_CG model", {
  expect_equal(typeof(celdaCG.sim$counts), "integer")
  expect_true(all(sweep(factorized$counts$cell.states, 2, colSums(celdaCG.sim$counts), "/") == factorized$proportions$cell.states))
  expect_equal(celdaCG.sim$K, ncol(factorized$proportions$population.states))
  expect_equal(celdaCG.sim$L, nrow(factorized$proportions$population.states))                        
})

# Cluster probabilities
test_that(desc = "Testing clusterProbability with celda_CG", {
  clust.prob = clusterProbability(counts = celdaCG.sim$counts, model_CG)
  expect_true(all(round(rowSums(clust.prob$z.probability), 10) == 1) & nrow(clust.prob$z.probability) == ncol(celdaCG.sim$counts))
  expect_true(all(round(rowSums(clust.prob$y.probability), 10) == 1) & nrow(clust.prob$y.probability) == nrow(celdaCG.sim$counts))

  clust.prob = clusterProbability(counts = celdaCG.sim$counts, model_CG, log=TRUE)
  expect_true(all(round(rowSums(normalizeLogProbs(clust.prob$z.probability)), 10) == 1) & nrow(clust.prob$z.probability) == ncol(celdaCG.sim$counts))
  expect_true(all(round(rowSums(normalizeLogProbs(clust.prob$y.probability)), 10) == 1) & nrow(clust.prob$y.probability) == nrow(celdaCG.sim$counts))
})

test_that(desc = "Testing simulateCells.celda_G error checking with low gamma", {
  expect_error(simulateCells(model = "celda_CG", gamma=0.000001))
})

test_that(desc = "Testing celdaGridSearch with celda_CG", {
  celdaCG.res <- celdaGridSearch(counts=celdaCG.sim$counts, model="celda_CG", nchains = 2, params.test=list(K=4:5, L=9:10), params.fixed=list(sample.label=celdaCG.sim$sample.label), max.iter = 10, verbose = FALSE, best.only=FALSE)
  expect_true(all(class(celdaCG.res) == c("celda_list", "celda_CG")))

  expect_equal(is.null(celdaCG.res$perplexity), TRUE)
  expect_error(plotGridSearchPerplexity(celdaCG.res))

  celdaCG.res = calculatePerplexityWithResampling(celdaCG.sim$counts, celdaCG.res, resample=2)
  expect_equal(is.null(celdaCG.res$perplexity), FALSE)
  expect_is(celdaCG.res, "celda_list")
  expect_is(celdaCG.res, "celda_CG")
  expect_error(calculatePerplexityWithResampling(celdaCG.sim$counts, celdaCG.res, resample="2"))
  expect_error(calculatePerplexityWithResampling(celdaCG.sim$counts, "celdaCG.res", resample=2))
  
  plot.obj = plotGridSearchPerplexity(celdaCG.res)
  expect_is(plot.obj, "ggplot")

  celdaCG.res.K5.L10 = subsetCeldaList(celdaCG.res, params=list(K = 5, L = 10))
  model_CG = selectBestModel(celdaCG.res.K5.L10)
  expect_error(celdaCG.res <- calculatePerplexityWithResampling(celdaCG.sim$counts, model_CG, resample=2))
  expect_error(celdaCG.res <- calculatePerplexityWithResampling(celdaCG.sim$counts, celdaCG.res, resample='a'))

  celdaC.res = celdaGridSearch(counts=celdaCG.sim$counts, model="celda_C", nchains = 1, params.test=list(K=4:5), params.fixed=list(sample.label=celdaCG.sim$sample.label), max.iter = 10, verbose = FALSE, best.only=TRUE)
  expect_error(plotGridSearchPerplexity.celda_CG(celdaC.res))
  
  res <- calculatePerplexity.celda_CG(celdaCG.sim$counts, model_CG)
  res2 <- calculatePerplexity.celda_CG(celdaCG.sim$counts, model_CG, new.counts = celdaCG.sim$counts + 1)
  
  expect_error(res <- calculatePerplexity.celda_CG(celdaCG.sim$counts, model_CG, new.counts = celdaCG.sim$counts[-1,]))    
})

# Ensure logLikelihood calculates the expected values
test_that(desc = "Testing logLikelihood.celda_CG", {
  expect_lt(logLikelihood(model="celda_CG",
                                         y = celdaCG.sim$y, z = celdaCG.sim$z,
                                         delta = 1, gamma = 1,  beta = 1, 
                                         alpha = 1, K = celdaCG.sim$K, L = celdaCG.sim$L, 
                                         s = celdaCG.sim$sample.label, 
                                         counts=celdaCG.sim$counts),0)
})

# normalizeCounts
test_that(desc = "Testing normalizeCounts with celda_CG", {
  norm.counts <- normalizeCounts(celdaCG.sim$counts)
  expect_equal(dim(norm.counts), dim(celdaCG.sim$counts))
  expect_equal(rownames(norm.counts), rownames(celdaCG.sim$counts))
  expect_equal(colnames(norm.counts), colnames(celdaCG.sim$counts))
})


# recodeClusterY
test_that(desc = "Testing recodeClusterY with celda_CG", {
  expect_error(recodeClusterY(celda.mod = model_CG, from = NULL, to = ))
  expect_error(recodeClusterY(celda.mod = model_CG, from = c(1,2,3,4,5), to = c(1,2,4,3,6)))
  expect_error(recodeClusterY(celda.mod = model_CG, from = c(1,2,3,4,6), to = c(1,2,4,3,5)))
  new.recoded <- recodeClusterY(celda.mod = model_CG, from = c(1,2,3,4,5), to = c(3,2,1,4,5))
  expect_equal(model_CG$y == 1, new.recoded$y == 3)
})

# recodeClusterZ
test_that(desc = "Testing recodeClusterZ with celda_CG", {
  expect_error(recodeClusterZ(celda.mod = model_CG, from = NULL, to = ))
  expect_error(recodeClusterZ(celda.mod = model_CG, from = c(1,2,3,4,5), to = c(1,2,3,4,6)))
  expect_error(recodeClusterZ(celda.mod = model_CG, from = c(1,2,3,4,6), to = c(1,2,3,4,5)))  
  new.recoded <- recodeClusterZ(celda.mod = model_CG, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
  expect_equal(model_CG$z == 1, new.recoded$z == 5)
})

# compareCountMatrix
test_that(desc = "Testing CompareCountMatrix with celda_CG", {
  expect_true(compareCountMatrix(counts = celdaCG.sim$counts, celda.mod = model_CG))
})

# topRank
test_that(desc = "Testing topRank with celda_CG", {
  top.rank <- topRank(matrix = factorized$proportions$gene.states, n = 1000, threshold = NULL)
  expect_equal(names(top.rank), c("index","names"))
  top.rank <- topRank(matrix = factorized$proportions$gene.states, n = 1000)
  expect_equal(nrow(celdaCG.sim$counts), sum(sapply(top.rank$names,length)))
  expect_equal(names(top.rank), c("index","names"))
})

# plotHeatmap
test_that(desc = "Testing plotHeatmap with celda_CG", {
  expect_error(plotHeatmap(counts = celdaCG.sim$counts, z = model_CG$y, y = model_CG$y), "Length of z must match number of columns in counts matrix")
  expect_error(plotHeatmap(counts = celdaCG.sim$counts, z = model_CG$z, y = model_CG$z), "Length of y must match number of rows in counts matrix")
  expect_error(plotHeatmap(counts = celdaCG.sim$counts, z = model_CG$z, y = model_CG$y, scale.row = model_CG), "'scale.row' needs to be of class 'function'")
  expect_error(plotHeatmap(counts = celdaCG.sim$counts, z = model_CG$z, y = model_CG$y, trim = 3), "'trim' should be a 2 element vector specifying the lower and upper boundaries")
  expect_equal(names(plotHeatmap(counts = celdaCG.sim$counts, z = model_CG$z, y = model_CG$y, cell.ix = 1:10)), c("tree_row", "tree_col", "gtable"))
  expect_equal(names(plotHeatmap(counts = celdaCG.sim$counts, z = NULL, y = model_CG$y, cell.ix = 1:10)), c("tree_row", "tree_col", "gtable"))
  expect_equal(names(plotHeatmap(counts = celdaCG.sim$counts, z = model_CG$z, y = model_CG$y, feature.ix = 1:10)), c("tree_row", "tree_col", "gtable"))
  expect_equal(names(plotHeatmap(counts = celdaCG.sim$counts, z = model_CG$z, y = NULL, feature.ix = 1:10)), c("tree_row", "tree_col", "gtable"))
})

# celdaHeatmap
test_that(desc = "Testing celdaHeatmap with celda_CG", {
  expect_equal(names(celdaHeatmap(celda.mod = model_CG, counts = celdaCG.sim$counts)),
               c("tree_row", "tree_col", "gtable"))
})

# moduleHeatmap
test_that(desc = "Checking moduleHeatmap to see if it runs", {
  expect_equal(names(moduleHeatmap(celdaCG.sim$counts, celda.mod = model_CG, feature.module = c(2,3))), c("tree_row","tree_col","gtable"))
  expect_error(moduleHeatmap(counts = "counts", celda.mod = model_CG, feature.module = c(2,3)),"'counts' should be a numeric count matrix")
  expect_error(moduleHeatmap(counts = celdaCG.sim$counts, celda.mod = "model", feature.module = c(2,3)), "'celda.mod' should be an object of class celda_G or celda_CG")               
})

# celdaProbabilityMap
test_that(desc = "Testing celdaProbabiltyMap.celda_CG for sample level",{
  plot.obj = celdaProbabilityMap(counts=celdaCG.sim$counts, celda.mod=model_CG, level="sample")
  expect_true(!is.null(plot.obj))
  
  ## Without a sample label
  model_CG_2 = celda_CG(celdaCG.sim$counts, sample.label=NULL, K=celdaCG.sim$K, L=celdaCG.sim$L, max.iter=5, nchain=1)
  plot.obj = celdaProbabilityMap(counts=celdaCG.sim$counts, celda.mod=model_CG_2, level="sample")
  expect_true(!is.null(plot.obj))  
})

test_that(desc = "Testing celdaProbabiltyMap.celda_CG for cell.population",{
  plot.obj = celdaProbabilityMap(counts=celdaCG.sim$counts, celda.mod=model_CG, level="cell.population")
  expect_true(!is.null(plot.obj))
})

# differentialExpression
test_that(desc = "Testing differentialExpression for celda_CG", {
  expect_equal(class(diffExp_K1 <- differentialExpression(counts = celdaCG.sim$counts, celda.mod = model_CG, c1 = 3, log2fc.threshold = 0.5)),
		c("data.table", "data.frame"))
  expect_equal(class(diffExp_K1 <- differentialExpression(counts = celdaCG.sim$counts, celda.mod = model_CG, c1 = 3, c2 = 4, log2fc.threshold = 0.5)),
               c("data.table", "data.frame"))		
  expect_error(differentialExpression(counts = "counts", celda.mod = model_CG, c1 = 3, log2fc.threshold = 0.5),"'counts' should be a numeric count matrix")
  expect_error(differentialExpression(counts = celdaCG.sim$counts, celda.mod = NULL, c1 = 3), "'celda.mod' should be an object of class celda_C or celda_CG")               
})

# plotDimReduce
test_that(desc = "Testing plotDimReduce* with celda_CG", {
  celda.tsne <- celdaTsne(counts = celdaCG.sim$counts, max.iter = 50, celda.mod = model_CG, max.cells = 500)
  expect_equal(names(plotDimReduceCluster(dim1 = celda.tsne[,1], dim2 = celda.tsne[,2],cluster = as.factor(model_CG$z),specific_clusters = c(1,2,3))),
               c("data", "layers", "scales", "mapping", "theme", "coordinates", "facet", "plot_env", "labels", "guides"))
  expect_equal(names(plotDimReduceState(dim1 = celda.tsne[,1], dim2 = celda.tsne[,2], counts = celdaCG.sim$counts, celda.mod = model_CG, modules = c("L1","L2"))),
               c("data", "layers", "scales", "mapping", "theme", "coordinates", "facet", "plot_env", "labels"))  
  expect_equal(names(plotDimReduceGene(dim1 = celda.tsne[,1],dim2 = celda.tsne[,2],counts = celdaCG.sim$counts,features = c("Gene_99"), exact.match = TRUE)),
               c("data", "layers", "scales", "mapping", "theme", "coordinates", "facet", "plot_env", "labels"))  
  expect_equal(names(plotDimReduceGene(dim1 = celda.tsne[,1],dim2 = celda.tsne[,2],counts = celdaCG.sim$counts,features = c("Gene_99"), exact.match = FALSE)),
               c("data", "layers", "scales", "mapping", "theme", "coordinates", "facet", "plot_env", "labels"))  
  expect_error(plotDimReduceState(dim1 = celda.tsne[,1], dim2 = celda.tsne[,2], matrix = factorized$proportions$cell.states, distance = "char"))
})

# celdaTsne
test_that(desc = "Testing celdaTsne.celda_CG with all cells",{
  tsne = celdaTsne(counts=celdaCG.sim$counts, celda.mod=model_CG, max.cells=length(model_CG$z))
  plot.obj = plotDimReduceCluster(tsne[,1], tsne[,2], model_CG$z)
  expect_true(ncol(tsne) == 2 & nrow(tsne) == length(model_CG$z))
  expect_true(!is.null(plot.obj))
  
  tsne = celdaTsne(counts=celdaCG.sim$counts, celda.mod=model_CG, max.cells=ncol(celdaCG.sim$counts), modules=1:2)
  expect_error(tsne <- celdaTsne(counts=celdaCG.sim$counts, celda.mod=model_CG, max.cells=ncol(celdaCG.sim$counts), modules=1000:1005))
})

test_that(desc = "Testing celdaTsne.celda_CG with subset of cells",{
  expect_success(expect_error(tsne <- celdaTsne(counts=celdaCG.sim$counts, celda.mod=model_CG, max.cells=50, min.cluster.size=50)))
  tsne = celdaTsne(counts=celdaCG.sim$counts, celda.mod=model_CG, max.cells=100, min.cluster.size=10)
  plot.obj = plotDimReduceCluster(tsne[,1], tsne[,2], model_CG$z)
  expect_true(ncol(tsne) == 2 & nrow(tsne) == length(model_CG$z) && sum(!is.na(tsne[,1])) == 100)
  expect_true(!is.null(plot.obj))
})

# featureModuleLookup
test_that(desc = "Testing featureModuleLookup with celda_CG", {
  res = featureModuleLookup(celdaCG.sim$counts, model_CG, "Gene_1")
  expect_true(res == model_CG$y[1])
  
  res = featureModuleLookup(celdaCG.sim$counts, model_CG, "XXXXXXX")
  expect_true(grepl("No feature", res))
})


# cCG.splitZ/cCG.splitZ
test_that(desc = "Testing cCG.splitZ and cCG.splitY", {
  r = simulateCells("celda_CG", S=1, G=100, C.Range=c(50,100), K=2, L=2)
  model_CG = celda_CG(r$counts, K=r$K, L=r$L, max.iter=5, nchain=1)
  probs = clusterProbability(r$counts, model_CG, log=TRUE)
    
  dc = cCG.decomposeCounts(r$counts, r$sample.label, r$z, r$y, r$K, r$L)
  res = cCG.splitZ(r$counts, dc$m.CP.by.S, dc$n.TS.by.C, dc$n.TS.by.CP, dc$n.by.G, dc$n.by.TS, dc$nG.by.TS, as.integer(r$sample.label), z=r$z, K=r$K, L=r$L, nS=dc$nS, nG=dc$nG, alpha=1, beta=1, delta=1, gamma=1,z.prob=probs$z.probability, min.cell=1000)
  expect_true(grepl("Cluster sizes too small", res$message))
  res = cCG.splitY(r$counts, r$y, dc$m.CP.by.S, dc$n.G.by.CP, dc$n.TS.by.C, dc$n.TS.by.CP, dc$n.by.G, dc$n.by.TS, dc$nG.by.TS, dc$n.CP, s=as.integer(r$sample.label), z=r$z, K=r$K, L=r$L, nS=dc$nS, nG=dc$nG, alpha=1, beta=1, delta=1, gamma=1, y.prob=probs$y.probability, min.cell=1000)
  expect_true(grepl("Cluster sizes too small", res$message))
  
  ## Testing K.subclusters parameter
  res = cCG.splitY(r$counts, r$y, dc$m.CP.by.S, dc$n.G.by.CP, dc$n.TS.by.C, dc$n.TS.by.CP, dc$n.by.G, dc$n.by.TS, dc$nG.by.TS, dc$n.CP, s=as.integer(r$sample.label), z=r$z, K=r$K, L=r$L, nS=dc$nS, nG=dc$nG, alpha=1, beta=1, delta=1, gamma=1, y.prob=probs$y.probability, K.subclusters=1000)
  expect_true(length(res$y) == nrow(r$counts))  
})

