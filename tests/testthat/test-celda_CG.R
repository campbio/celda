#celda_CG
library(celda)
context("Testing celda_CG")

celdaCG.sim = simulateCells("celda_CG", K=5, L=10)
model_CG = celda_CG(counts=celdaCG.sim$counts, 
                    sample.label=celdaCG.sim$sample.label, 
                    K=celdaCG.sim$K, L=celdaCG.sim$L, 
                    max.iter = 1, nchains = 2,
                    algorithm="EM", verbose=FALSE)
factorized <- factorizeMatrix(celda.mod = model_CG, counts = celdaCG.sim$counts)

# celda_CG
test_that(desc = "Testing simulation and celda_CG model", {
  expect_equal(typeof(celdaCG.sim$counts), "integer")
  expect_true(all(sweep(factorized$counts$cell, 2, colSums(celdaCG.sim$counts), "/") == factorized$proportions$cell))
  expect_equal(celdaCG.sim$K, ncol(factorized$proportions$cell.population))
  expect_equal(celdaCG.sim$L, nrow(factorized$proportions$cell.population))   
  expect_true(all(is.numeric(logLikelihoodHistory(celda.mod = model_CG))))
  expect_equal(max(logLikelihoodHistory(celda.mod = model_CG)), bestLogLikelihood(model_CG))
  
  # GitHub #347
  numeric.counts = celdaCG.sim$counts
  storage.mode(numeric.counts) = "numeric"
  expect_true(is(celda_CG(counts=celdaCG.sim$counts, 
                          sample.label=celdaCG.sim$sample.label, 
                          K=celdaCG.sim$K,  L=celdaCG.sim$L, 
                          algorithm="EM", verbose=FALSE,
                          max.iter=1, nchains=1),
               "celda_CG"))
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

test_that(desc = "Testing simulateCells.celda_CG error checking with low gamma", {
  expect_error(simulateCells(model = "celda_CG", gamma=0.000001))
})

test_that(desc = "Testing simulateCells.celda_CG, make sure all genes expressed", {
  sim.cells.low <- simulateCells(model = "celda_CG", G = 1000, C = 300,C.Range= c(0,100), N.Range = c(0,100))
  expect_true(all(rowSums(sim.cells.low$counts) > 0))
})

test_that(desc = "Testing celdaGridSearch with celda_CG", {  
  expect_error(celdaGridSearch(counts=celdaCG.sim$counts, 
                               model="celda_CG", 
                               max.iter=1, nchains=1,
                               params.test=list(K=5, L=10, M = 3:4), 
                               params.fixed=list(sample.label=celdaCG.sim$sample.label),
                               best.only=FALSE),
               "The following elements in 'params.test' are not arguments of 'celda_CG': M")
  
  expect_error(celdaGridSearch(counts=celdaCG.sim$counts, 
                               model="celda_CG", nchains = 1, 
                               max.iter=1,
                               params.test=list(K=5, L=10, 
                                                sample.label = "Sample"), 
                               params.fixed=list(sample.label=celdaCG.sim$sample.label)),
               "Setting parameters such as 'z.init', 'y.init', and 'sample.label' in 'params.test' is not currently supported.")
  
  expect_error(celdaGridSearch(counts=celdaCG.sim$counts, 
                               model="celda_CG", nchains = 1, max.iter=1,
                               params.test=list(), 
                               params.fixed=list(sample.label=celdaCG.sim$sample.label)),
               "The following arguments are not in 'params.test' or 'params.fixed' but are required for 'celda_CG': K,L")
  
  expect_error(celdaGridSearch(counts=celdaCG.sim$counts, model="celda_CG", 
                               nchains = 1, max.iter=1, params.test=list(K=4:5, L=9:10), 
                               params.fixed=list(sample.label=celdaCG.sim$sample.label, xxx = "xxx")),
               "The following elements in 'params.fixed' are not arguments of 'celda_CG': xxx")
  
  expect_warning(celdaGridSearch(counts=celdaCG.sim$counts, 
                                 model="celda_CG", max.iter=1,
                                 params.test=list(K=5, L=10, nchains = 2)),
                                 "Parameter 'nchains' should not be used within the params.test list")
  
  
  celdaCG.res <- celdaGridSearch(counts=celdaCG.sim$counts, 
                                 model="celda_CG", nchains = 2, 
                                 params.test=list(K=5, L=10), 
                                 params.fixed=list(sample.label=celdaCG.sim$sample.label), 
                                 max.iter = 1, verbose = FALSE, 
                                 best.only=FALSE)
  expect_true(is(celdaCG.res, "celdaList"))
  expect_error(plotGridSearchPerplexity(celdaCG.res))
  expect_equal(names(runParams(celdaCG.res)), 
               c("index","chain","K","L","log_likelihood"))
  
  celdaCG.res = resamplePerplexity(celdaCG.sim$counts, celdaCG.res, resample=2)
  expect_is(celdaCG.res, "celdaList")
  expect_error(resamplePerplexity(celdaCG.sim$counts, celdaCG.res, resample="2"))
  expect_error(resamplePerplexity(celdaCG.sim$counts, "celdaCG.res", resample=2))
  
  plot.obj = plotGridSearchPerplexity(celdaCG.res)
  expect_is(plot.obj, "ggplot")
  
  expect_error(subsetCeldaList(celda.list = "celda_list"),
               "celda.list parameter was not of class celdaList.")
  expect_error(subsetCeldaList(celdaCG.res, params = list(K = 7, L = 11)))
  expect_error(subsetCeldaList(celdaCG.res, params = list(K = 5, M = 10)))
                 
  celdaCG.res.K5.L10 = subsetCeldaList(celdaCG.res, params=list(K = 5, L = 10))
  model_CG = selectBestModel(celdaCG.res.K5.L10)
  
  expect_error(selectBestModel("celda_list"),
               "celda.list parameter was not of class celdaList.")
  expect_error(celdaCG.res <- resamplePerplexity(celdaCG.sim$counts, model_CG, resample=2))
  expect_error(celdaCG.res <- resamplePerplexity(celdaCG.sim$counts, celdaCG.res, resample='a'))
  
  celdaCG.res.index1 = subsetCeldaList(celdaCG.res, params=list(index = 1))
  expect_true(all(is(celdaCG.res.index1, "celda_CG") && !is(celdaCG.res.index1, "celda_list")))
  
  res <- perplexity(celdaCG.sim$counts, model_CG)
  res2 <- perplexity(celdaCG.sim$counts, model_CG, new.counts = celdaCG.sim$counts + 1)
  
  expect_error(res <- perplexity(celdaCG.sim$counts, model_CG, new.counts = celdaCG.sim$counts[-1,]))    
})

# Ensure logLikelihood calculates the expected values
test_that(desc = "Testing logLikelihood.celda_CG", {
  expect_lt(logLikelihood(model="celda_CG",
                          y = celdaCG.sim$y, z = celdaCG.sim$z,
                          delta = 1, gamma = 1,  beta = 1, 
                          alpha = 1, K = celdaCG.sim$K, L = celdaCG.sim$L, 
                          s = celdaCG.sim$sample.label, 
                          counts=celdaCG.sim$counts),
            0)
  
  fake.z = celdaCG.sim$z
  fake.z[1] = celdaCG.sim$K + 1
  expect_error(logLikelihood(model="celda_CG",
                             y = celdaCG.sim$y, z = fake.z,
                             delta = 1, gamma = 1,  beta = 1, 
                             alpha = 1, K = celdaCG.sim$K, L = celdaCG.sim$L, 
                             s = celdaCG.sim$sample.label, 
                             counts=celdaCG.sim$counts),
                             "An entry in z contains a value greater than the provided K.")
  
  fake.y = celdaCG.sim$y
  fake.y[1] = celdaCG.sim$L + 1
  expect_error(logLikelihood(model="celda_CG",
                             y = fake.y, z = celdaCG.sim$z,
                             delta = 1, gamma = 1,  beta = 1, 
                             alpha = 1, K = celdaCG.sim$K, L = celdaCG.sim$L, 
                             s = celdaCG.sim$sample.label, 
                             counts=celdaCG.sim$counts),
                             "An entry in y contains a value greater than the provided L.")
})

# normalizeCounts
test_that(desc = "Testing normalizeCounts with celda_CG", {
  norm.counts <- normalizeCounts(celdaCG.sim$counts)
  expect_equal(dim(norm.counts), dim(celdaCG.sim$counts))
  expect_equal(rownames(norm.counts), rownames(celdaCG.sim$counts))
  expect_equal(colnames(norm.counts), colnames(celdaCG.sim$counts))
  expect_error(normalizeCounts(celdaCG.sim$counts, transformation.fun = "scale"), 
               "'transformation.fun' needs to be of class 'function'")
  expect_error(normalizeCounts(celdaCG.sim$counts, scale.fun = "scale"), 
               "'scale.fun' needs to be of class 'function'")
})


# recodeClusterY
test_that(desc = "Testing recodeClusterY with celda_CG", {
  expect_error(recodeClusterY(celda.mod = model_CG, from = NULL, to = ))
  expect_error(recodeClusterY(celda.mod = model_CG, from = c(1,2,3,4,5), to = c(1,2,4,3,6)))
  expect_error(recodeClusterY(celda.mod = model_CG, from = c(1,2,3,4,6), to = c(1,2,4,3,5)))
  new.recoded <- recodeClusterY(celda.mod = model_CG, from = c(1,2,3,4,5), to = c(3,2,1,4,5))
  expect_equal(model_CG@clusters$y == 1, new.recoded@clusters$y == 3)
})

# recodeClusterZ
test_that(desc = "Testing recodeClusterZ with celda_CG", {
  expect_error(recodeClusterZ(celda.mod = model_CG, from = NULL, to = ))
  expect_error(recodeClusterZ(celda.mod = model_CG, from = c(1,2,3,4,5), to = c(1,2,3,4,6)))
  expect_error(recodeClusterZ(celda.mod = model_CG, from = c(1,2,3,4,6), to = c(1,2,3,4,5)))  
  new.recoded <- recodeClusterZ(celda.mod = model_CG, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
  expect_equal(model_CG@clusters$z == 1, new.recoded@clusters$z == 5)
})

# compareCountMatrix
test_that(desc = "Testing CompareCountMatrix with celda_CG", {
  expect_true(compareCountMatrix(counts = celdaCG.sim$counts, celda.mod = model_CG))
  
  less.features <- celdaCG.sim$counts[1:50,]
  expect_error(compareCountMatrix(counts = less.features, celda.mod = model_CG),
               "There was a mismatch between the provided count matrix and the count matrix used to generate the provided celda result.")
  less.cells <- celdaCG.sim$counts[,1:100]
  expect_error(compareCountMatrix(counts = less.cells, celda.mod = model_CG),
               "There was a mismatch between the provided count matrix and the count matrix used to generate the provided celda result.")
  
  counts.matrix.error <- matrix(data = 1, nrow = nrow(celdaCG.sim$counts), ncol = ncol(celdaCG.sim$counts))
  expect_false(compareCountMatrix(counts = counts.matrix.error, celda.mod = model_CG, error.on.mismatch = FALSE))
  expect_error(compareCountMatrix(counts = counts.matrix.error, celda.mod = model_CG, error.on.mismatch = TRUE))
})

# topRank
test_that(desc = "Testing topRank with celda_CG", {
  top.rank <- topRank(matrix = factorized$proportions$module, n = 1000, threshold = NULL)
  expect_equal(names(top.rank), c("index","names"))
  top.rank <- topRank(matrix = factorized$proportions$module, n = 1000)
  expect_equal(nrow(celdaCG.sim$counts), sum(sapply(top.rank$names,length)))
  expect_equal(names(top.rank), c("index","names"))
})

# plotHeatmap
test_that(desc = "Testing plotHeatmap with celda_CG", {
  expect_error(plotHeatmap(counts = celdaCG.sim$counts, z = model_CG@clusters$y, y = model_CG@clusters$y), "Length of z must match number of columns in counts matrix")
  expect_error(plotHeatmap(counts = celdaCG.sim$counts, z = model_CG@clusters$z, y = model_CG@clusters$z), "Length of y must match number of rows in counts matrix")
  expect_error(plotHeatmap(counts = celdaCG.sim$counts, z = model_CG@clusters$z, y = model_CG@clusters$y, scale.row = model_CG), "'scale.row' needs to be of class 'function'")
  expect_error(plotHeatmap(counts = celdaCG.sim$counts, z = model_CG@clusters$z, y = model_CG@clusters$y, trim = 3), "'trim' should be a 2 element vector specifying the lower and upper boundaries")
  expect_equal(names(plotHeatmap(counts = celdaCG.sim$counts, z = model_CG@clusters$z, y = model_CG@clusters$y, cell.ix = 1:10)), c("tree_row", "tree_col", "gtable"))
  expect_equal(names(plotHeatmap(counts = celdaCG.sim$counts, z = NULL, y = model_CG@clusters$y, cell.ix = 1:10, color.scheme = "sequential")), c("tree_row", "tree_col", "gtable"))
  expect_equal(names(plotHeatmap(counts = celdaCG.sim$counts, z = model_CG@clusters$z, y = model_CG@clusters$y, feature.ix = 1:10)), c("tree_row", "tree_col", "gtable"))
  expect_equal(names(plotHeatmap(counts = celdaCG.sim$counts, z = model_CG@clusters$z, y = NULL, feature.ix = 1:10)), c("tree_row", "tree_col", "gtable"))
  expect_equal(names(plotHeatmap(counts = celdaCG.sim$counts, z = model_CG@clusters$z, y = model_CG@clusters$y, cell.ix = 1:10, color.scheme = "sequential", annotation.color = "default")), c("tree_row", "tree_col", "gtable"))
})

#plotHeatmap with annotations
test_that(desc = "Testing plotHeatmap with annotations", {
  annot.cell <- as.data.frame(c(rep(x = 1, times = ncol(celdaCG.sim$counts) - 100),rep(x = 2, 100)))
  annot.feature <- as.data.frame(c(rep(x = 1, times = nrow(celdaCG.sim$counts) - 100),rep(x = 2, 100)))
  
  rownames(annot.cell) <- colnames(celdaCG.sim$counts)
  colnames(annot.cell) <- "cell"
  rownames(annot.feature) <- rownames(celdaCG.sim$counts)
  colnames(annot.feature) <- "feature"
  expect_equal(names(plotHeatmap(celda.mod = model_CG, counts = celdaCG.sim$counts, annotation.cell = annot.cell, annotation.feature = annot.feature, z = model_CG@clusters$z, y = model_CG@clusters$y)),
               c("tree_row", "tree_col", "gtable"))
  
  rownames(annot.cell) <- NULL
  rownames(annot.feature) <- NULL
  expect_equal(names(plotHeatmap(celda.mod = model_CG, counts = celdaCG.sim$counts, annotation.cell = as.matrix(annot.cell), annotation.feature = as.matrix(annot.feature), z = model_CG@clusters$z, y = model_CG@clusters$y)),
               c("tree_row", "tree_col", "gtable"))
  
})

# celdaHeatmap
test_that(desc = "Testing celdaHeatmap with celda_CG", {
  expect_equal(names(celdaHeatmap(celda.mod = model_CG, counts = celdaCG.sim$counts)),
               c("tree_row", "tree_col", "gtable"))
})

# moduleHeatmap
test_that(desc = "Checking moduleHeatmap to see if it runs", {
  expect_equal(names(moduleHeatmap(celdaCG.sim$counts, celda.mod = model_CG, feature.module = c(2,3), top.cells = 500)), c("tree_row","tree_col","gtable"))
  expect_equal(names(moduleHeatmap(celdaCG.sim$counts, celda.mod = model_CG, top.features = 15, top.cells = 15, normalize = FALSE)), c("tree_row","tree_col","gtable"))
  expect_equal(names(moduleHeatmap(celdaCG.sim$counts, celda.mod = model_CG, top.features = 15, top.cells = NULL, normalize = FALSE)), c("tree_row","tree_col","gtable"))
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
  expect_equal(class(diffExp_K1 <- differentialExpression(counts = celdaCG.sim$counts, celda.mod = model_CG, c1 = 3, log2fc.threshold = 0.5, only.pos = TRUE)),
		c("data.table", "data.frame"))
  expect_equal(class(diffExp_K1 <- differentialExpression(counts = celdaCG.sim$counts, celda.mod = model_CG, c1 = 2:3, c2 = 4, log2fc.threshold = 0.5)),
               c("data.table", "data.frame"))		
  expect_error(differentialExpression(counts = "counts", celda.mod = model_CG, c1 = 3, log2fc.threshold = 0.5),"'counts' should be a numeric count matrix")
  expect_error(differentialExpression(counts = celdaCG.sim$counts, celda.mod = NULL, c1 = 3), "'celda.mod' should be an object of class celda_C or celda_CG")               
  expect_error(differentialExpression(counts = celdaCG.sim$counts, celda.mod = model_CG, c1 = NULL, log2fc.threshold = 0.5, only.pos = TRUE))
})

# plotDimReduce
test_that(desc = "Testing plotDimReduce* with celda_CG", {
  celda.tsne <- celdaTsne(counts = celdaCG.sim$counts, max.iter = 50, celda.mod = model_CG, max.cells = 500)
  expect_equal(names(plotDimReduceCluster(dim1 = celda.tsne[,1], dim2 = celda.tsne[,2],cluster = as.factor(model_CG@clusters$z),specific_clusters = c(1,2,3))),
               c("data", "layers", "scales", "mapping", "theme", "coordinates", "facet", "plot_env", "labels", "guides"))
  expect_equal(names(plotDimReduceModule(dim1 = celda.tsne[,1], dim2 = celda.tsne[,2], counts = celdaCG.sim$counts, celda.mod = model_CG, modules = c("L1", "L2"))),
               c("data", "layers", "scales", "mapping", "theme", "coordinates", "facet", "plot_env", "labels"))  
  expect_equal(names(plotDimReduceModule(dim1 = celda.tsne[,1], dim2 = celda.tsne[,2], counts = celdaCG.sim$counts, celda.mod = model_CG, modules = c("L1", "L2"), rescale = FALSE)),
               c("data", "layers", "scales", "mapping", "theme", "coordinates", "facet", "plot_env", "labels"))    
  expect_equal(names(plotDimReduceFeature(dim1 = celda.tsne[,1],dim2 = celda.tsne[,2],counts = celdaCG.sim$counts,features = c("Gene_99"), exact.match = TRUE)),
               c("data", "layers", "scales", "mapping", "theme", "coordinates", "facet", "plot_env", "labels"))  
  expect_equal(names(plotDimReduceFeature(dim1 = celda.tsne[,1],dim2 = celda.tsne[,2],counts = celdaCG.sim$counts,features = c("Gene_99"), exact.match = FALSE)),
               c("data", "layers", "scales", "mapping", "theme", "coordinates", "facet", "plot_env", "labels"))  
  expect_error(plotDimReduceModule(dim1 = celda.tsne[,1], dim2 = celda.tsne[,2], counts = celdaCG.sim$counts, celda.mod = model_CG, modules = c(11,12)))
  expect_error(plotDimReduceFeature(dim1 = celda.tsne[,1],dim2 = celda.tsne[,2],counts = celdaCG.sim$counts,features = NULL, exact.match = TRUE))
  expect_error(plotDimReduceFeature(dim1 = celda.tsne[,1],dim2 = celda.tsne[,2],counts = celdaCG.sim$counts,features = c("Gene_99"), trim = 2, exact.match = TRUE))
  expect_error(plotDimReduceFeature(dim1 = celda.tsne[,1],dim2 = celda.tsne[,2],counts = celdaCG.sim$counts,features = c("Nonexistent_Gene"), exact.match = TRUE))
  
  # Check cases when there are some or all features not present in the counts matrix
  expect_error(plotDimReduceFeature(dim1 = celda.tsne[,1],dim2 = celda.tsne[,2],counts = celdaCG.sim$counts,features = c("Nonexistent_Gene"), exact.match = TRUE))
  expect_warning(plotDimReduceFeature(dim1 = celda.tsne[,1],
                                      dim2 = celda.tsne[,2],
                                      counts = celdaCG.sim$counts,
                                      features = c("Gene_99","Nonexistent_Gene"), 
                                      exact.match = TRUE))
})

# celdaTsne
test_that(desc = "Testing celdaTsne with celda_CG when model class is changed, should error",{
  model_X <- model_CG
  class(model_X) <- "celda_X"
  expect_error(celdaTsne(counts=celdaCG.sim$counts, celda.mod=model_X),
               "unable to find an inherited method for function 'celdaTsne' for signature '\"celda_X\"'")
})

test_that(desc = "Testing celdaTsne.celda_CG with all cells",{
  tsne = celdaTsne(counts=celdaCG.sim$counts, celda.mod=model_CG, max.cells=length(model_CG@clusters$z))
  plot.obj = plotDimReduceCluster(tsne[,1], tsne[,2], model_CG@clusters$z)
  expect_true(ncol(tsne) == 2 & nrow(tsne) == length(model_CG@clusters$z))
  expect_true(!is.null(plot.obj))
  
  tsne = celdaTsne(counts=celdaCG.sim$counts, celda.mod=model_CG, max.cells=ncol(celdaCG.sim$counts), modules=1:2)
  expect_error(tsne <- celdaTsne(counts=celdaCG.sim$counts, celda.mod=model_CG, max.cells=ncol(celdaCG.sim$counts), modules=1000:1005))
})

test_that(desc = "Testing celdaTsne.celda_CG with subset of cells",{
  expect_success(expect_error(tsne <- celdaTsne(counts=celdaCG.sim$counts, celda.mod=model_CG, max.cells=50, min.cluster.size=50)))
  tsne = celdaTsne(counts=celdaCG.sim$counts, celda.mod=model_CG, max.cells=100, min.cluster.size=10)
  plot.obj = plotDimReduceCluster(tsne[,1], tsne[,2], model_CG@clusters$z)
  expect_true(ncol(tsne) == 2 & nrow(tsne) == length(model_CG@clusters$z) && sum(!is.na(tsne[,1])) == 100)
  expect_true(!is.null(plot.obj))
})

# featureModuleLookup
test_that(desc = "Testing featureModuleLookup with celda_CG", {
  res = featureModuleLookup(celdaCG.sim$counts, model_CG, "Gene_1")
  expect_true(res == model_CG@clusters$y[1])

  res = featureModuleLookup(celdaCG.sim$counts, model_CG, "Gene_2", exact.match = FALSE)
  expect_true(length(res) == 11)
  
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

test_that(desc = "Testing perplexity.celda_CG", {
  expect_true(is.numeric(perplexity(celdaCG.sim$counts, model_CG)))
  
  class(model_CG) = c("celda_C")
  expect_error(perplexity.celda_CG(celdaCG.sim$counts, model_CG),
               "could not find function \"perplexity.celda_CG\"")
})

#miscellaneous fxns

#functions used internally
test_that(desc = "Invoking error from distinct_colors function",{
  expect_error(distinct_colors(n = 3, hues = "xx"), "Only color names listed in the 'color' function can be used in 'hues'")
})

test_that(desc = "Invoking error from distinct_colors function",{
  expect_error(processSampleLabels("Sample_1", ncol(celdaCG.sim$counts)), "'sample.label' must be the same length as the number of columns in the 'counts' matrix.")
})

test_that(desc = "Invoking error from logMessages function", {
  expect_error(logMessages(date(), logfile = 5))
})

test_that(desc = "miscellaneous distance fxns that are not directly used within celda, but will be tested", {
  x = data.frame(x = 2:4, y=  1:3)
  expect_equal(class(hellingerDist(x)), "dist")
  expect_equal(class(spearmanDist(x)), "dist")
})

