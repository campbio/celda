#celda_G
library(celda)
context("Testing celda_G")

celdaG.sim = simulateCells("celda_G", L=5, G=100)
model_G = celda_G(counts=celdaG.sim$counts, L=celdaG.sim$L, max.iter=10, verbose=FALSE)
factorized = factorizeMatrix(counts=celdaG.sim$counts, celda.mod = model_G)  

test_that(desc = "Testing clusterProbability with celda_G", {
  expect_true(ncol(clusterProbability(celdaG.sim$counts, model_G)$y.probability) == celdaG.sim$L)
})

test_that(desc = "Testing simulateCells.celda_G error checking with low gamma", {
  expect_error(simulateCells(model = "celda_G", gamma=0.000001))
})

test_that(desc = "Testing simulateCells.celda_G, make sure all genes expressed", {
  sim.cells.low <- simulateCells(model = "celda_G", G = 1000, C = 300, N.Range = c(0,10))
  expect_true(all(rowSums(sim.cells.low$counts) > 0))
})

test_that(desc = "Testing LogLikelihood functions", {
  expect_true(all(is.numeric(completeLogLikelihood(celda.mod = model_G))))
  expect_equal(max(completeLogLikelihood(celda.mod = model_G)), finalLogLikelihood(model_G))
})

test_that(desc = "Testing celdaGridSearch with celda_G", {
  celdaG.res <- celdaGridSearch(counts = celdaG.sim$counts, model = "celda_G", nchains = 2, params.test=list(L=c(5,10)), max.iter = 10, verbose = FALSE, best.only=FALSE)
  expect_true(all(class(celdaG.res) == c("celda_list", "celda_G")))
  expect_equal(is.null(celdaG.res$perplexity), TRUE)
  expect_error(plotGridSearchPerplexity(celdaG.res))
  expect_equal(names(runParams(celda.list = celdaG.res)), c("index","chain","L","log_likelihood"))
  
  celdaG.res = calculatePerplexityWithResampling(celdaG.sim$counts, celdaG.res, resample=2)
  expect_equal(is.null(celdaG.res$perplexity), FALSE)
  expect_is(celdaG.res, "celda_list")
  expect_is(celdaG.res, "celda_G")
  expect_error(calculatePerplexityWithResampling(celdaG.sim$counts, celdaG.res, resample="2"))
  expect_error(calculatePerplexityWithResampling(celdaG.sim$counts, "celdaG.res", resample=2))
  
  plot.obj = plotGridSearchPerplexity(celdaG.res)
  expect_is(plot.obj, "ggplot")

  celdaC.res = celdaGridSearch(counts = celdaG.sim$counts, model = "celda_C", nchains = 1, params.test=list(K=c(5,10)), max.iter = 10, verbose = FALSE, best.only=TRUE)
  expect_error(plotGridSearchPerplexity.celda_G(celdaC.res))

  celdaG.res.L5 <- subsetCeldaList(celdaG.res, params=list(L = 5))
  model_G = selectBestModel(celdaG.res.L5)
  res <- calculatePerplexity.celda_G(celdaG.sim$counts, model_G)
  res2 <- calculatePerplexity.celda_G(celdaG.sim$counts, model_G, new.counts = celdaG.sim$counts + 1)
  
  expect_error(res <- calculatePerplexity.celda_G(celdaG.sim$counts, model_G, new.counts = celdaG.sim$counts[-1,]))
})

# calculateLoglikFromVariables
test_that(desc = "Testing calculateLoglikFromVariables.celda_G", {
  expect_lt(calculateLoglikFromVariables(model = "celda_G", 
                                         counts = celdaG.sim$counts, 
                                         y = celdaG.sim$y, L = celdaG.sim$L, delta = 1, 
                                         gamma = 1, beta = 1),0)
})

# normalizeCounts
test_that(desc = "Making sure normalizeCounts doesn't change dimensions of counts matrix", {
  norm.counts <- normalizeCounts(celdaG.sim$counts)
  expect_equal(dim(norm.counts), dim(celdaG.sim$counts))
  expect_equal(rownames(norm.counts), rownames(celdaG.sim$counts))
  expect_equal(colnames(norm.counts), colnames(celdaG.sim$counts))
  expect_error(normalizeCounts(celdaG.sim$counts, transformation.fun = "scale"), 
               "'transformation.fun' needs to be of class 'function'")
  expect_error(normalizeCounts(celdaG.sim$counts, scale.fun = "scale"), 
               "'scale.fun' needs to be of class 'function'")
})

# recodeClusterY
test_that(desc = "Testing recodeClusterY with celda_G", {
  expect_error(recodeClusterZ(celda.mod = model_G, from = c(1,2,3,4,5), to = c(5,4,3,2,1)))
  expect_error(recodeClusterY(celda.mod = model_G, from = NULL, to = ))
  expect_error(recodeClusterY(celda.mod = model_G, from = c(1,2,3,4,5), to = c(1,2,3,4,6)))
  expect_error(recodeClusterY(celda.mod = model_G, from = c(1,2,3,4,6), to = c(1,2,3,4,5)))  
  new.recoded = recodeClusterY(celda.mod = model_G, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
  expect_equal(model_G$y == 1, new.recoded$y == 5)
})

# compareCountMatrix
test_that(desc = "Testing CompareCountMatrix with celda_G", {
  expect_true(compareCountMatrix(counts = celdaG.sim$counts, celda.mod = model_G))
  less.features <- celdaG.sim$counts[1:50,]
  expect_error(compareCountMatrix(counts = less.features, celda.mod = model_G),
               "The provided celda object was generated from a counts matrix with a different number of genes than the one provided.")
  less.cells <- celdaG.sim$counts[,1:100]
  expect_error(compareCountMatrix(counts = less.cells, celda.mod = model_G),
               "The provided celda object was generated from a counts matrix with a different number of cells than the one provided.")
  
  counts.matrix.error <- matrix(x = 1, nrow = nrow(celdaG.sim$counts), ncol = ncol(celdaG.sim$counts))
  expect_false(compareCountMatrix(counts = counts.matrix.error, celda.mod = model_G, error.on.mismatch = FALSE))
  expect_error(compareCountMatrix(counts = counts.matrix.error, celda.mod = model_G, error.on.mismatch = TRUE))
  
})

# topRank
test_that(desc = "Testing topRank function with celda_G", {
  top.rank <- topRank(matrix = factorized$proportions$gene.states, n = 1000, threshold = NULL)
  expect_equal(names(top.rank),
               c("index","names"))
  expect_equal(names(topRank(matrix = factorized$proportions$gene.states)),
               c("index","names"))
})

# plotHeatmap
test_that(desc = "Testing plotHeatmap with celda_G", {
  expect_error(plotHeatmap(counts = celdaG.sim$counts, y = model_G$L), "Length of y must match number of rows in counts matrix")
  expect_error(plotHeatmap(counts = celdaG.sim$counts, y = model_G$y, scale.row = "scale"), "'scale.row' needs to be of class 'function'")
  expect_error(plotHeatmap(counts = celdaG.sim$counts, y = model_G$y, trim = 3), "'trim' should be a 2 element vector specifying the lower and upper boundaries")
})

test_that(desc = "Testing plotHeatmap with celda_G, including annotations",{
  annot <- as.data.frame(c(rep(x = 1, times = nrow(celdaG.sim$counts) - 100),rep(x = 2, 100)))
  rownames(annot) <- rownames(celdaG.sim$counts)
  colnames(annot) <- "label"
  
  expect_equal(names(plotHeatmap(celda.mod = model_G, counts = celdaG.sim$counts, annotation.feature = annot, y = model_G$y)),
               c("tree_row", "tree_col", "gtable"))
  
  rownames(annot) <- NULL
  expect_equal(names(plotHeatmap(celda.mod = model_G, counts = celdaG.sim$counts, annotation.feature = as.matrix(annot), y = model_G$y)),
               c("tree_row", "tree_col", "gtable"))

  rownames(annot) <- rev(rownames(celdaG.sim$counts))
  expect_error(plotHeatmap(celda.mod = model_G, counts = celdaG.sim$counts, annotation.feature = annot, y = model_G$y),
               "Row names of 'annotation.feature' are different than the row names of 'counts'")
})

# celdaHeatmap
test_that(desc = "Testing celdaHeatmap with celda_G",{
  expect_equal(names(celdaHeatmap(celda.mod = model_G, counts = celdaG.sim$counts)),
               c("tree_row","tree_col","gtable"))
})

# moduleHeatmap
test_that(desc = "Testing moduleHeatmap with celda_G",{
  expect_equal(names(moduleHeatmap(celdaG.sim$counts, celda.mod = model_G, top.cells = 300, feature.module = c(1,2))),
               c("tree_row","tree_col","gtable"))
  expect_equal(names(moduleHeatmap(celdaG.sim$counts, celda.mod = model_G, top.features = 15, top.cells = 15, normalize = FALSE)),
               c("tree_row","tree_col","gtable"))

  expect_error(moduleHeatmap("counts", celda.mod = model_G), "'counts' should be a numeric count matrix")
  expect_error(moduleHeatmap(celdaG.sim$counts, celda.mod = "model_G"), "'celda.mod' should be an object of class celda_G or celda_CG")
})

# plotDimReduceState
test_that(desc = "Testing plotDimReduceState with celda_G", {
  celda.tsne <- celdaTsne(counts = celdaG.sim$counts,max.iter = 50,celda.mod = model_G)
  expect_equal(names(plotDimReduceState(dim1 = celda.tsne[,1], dim2 = celda.tsne[,2], counts = celdaG.sim$counts, celda.mod = model_G)),
               c("data", "layers", "scales", "mapping", "theme", "coordinates", "facet", "plot_env", "labels"))
  expect_equal(names(plotDimReduceState(dim1 = celda.tsne[,1], dim2 = celda.tsne[,2], counts = celdaG.sim$counts, celda.mod = model_G, modules = c("L1","L2"), rescale = F)),
               c("data", "layers", "scales", "mapping", "theme", "coordinates", "facet", "plot_env", "labels"))    
  expect_error(plotDimReduceState(dim1 = celda.tsne[,1], dim2 = celda.tsne[,2], counts = celdaG.sim$counts, celda.mod = model_G, modules = c("L11","L12")))
})

# celdaTsne
test_that(desc = "Testing celdaTsne with celda_G when model class is changed, should error",{
  model_X <- model_G
  class(model_X) <- "celda_X"
  expect_error(celdaTsne(counts=celdaG.sim$counts, celda.mod=model_X),
               "celda.mod argument is not of class celda_C, celda_G or celda_CG")
})

test_that(desc = "Testing celdaTsne with celda_C including all cells",{
  tsne = celdaTsne(counts=celdaG.sim$counts, celda.mod=model_G, max.cells=ncol(celdaG.sim$counts))
  plot.obj = plotDimReduceCluster(tsne[,1], tsne[,2], rep(1,ncol(celdaG.sim$counts)))
  expect_true(ncol(tsne) == 2 & nrow(tsne) == ncol(celdaG.sim$counts))
  expect_true(!is.null(plot.obj))
  
  tsne = celdaTsne(counts=celdaG.sim$counts, celda.mod=model_G, max.cells=ncol(celdaG.sim$counts), modules=1:2)
  expect_error(tsne <- celdaTsne(counts=celdaG.sim$counts, celda.mod=model_G, max.cells=ncol(celdaG.sim$counts), modules=1000:1005))
})

test_that(desc = "Testing celdaTsne with celda_G including a subset of cells",{
  tsne = celdaTsne(counts=celdaG.sim$counts, celda.mod=model_G, max.cells=100)
  plot.obj = plotDimReduceCluster(tsne[,1], tsne[,2], rep(1, ncol(celdaG.sim$counts)))
    expect_true(ncol(tsne) == 2 & nrow(tsne) == ncol(celdaG.sim$counts) && sum(!is.na(tsne[,1])) == 100)
  expect_true(!is.null(plot.obj))
})

# featureModuleLookup
test_that(desc = "Testing featureModuleLookup with celda_G", {
  res = featureModuleLookup(celdaG.sim$counts, model_G, "Gene_1")
  expect_true(res == model_G$y[1])
  res = featureModuleLookup(celdaG.sim$counts, model_G, "Gene_2", exact.match = FALSE)
  expect_true(length(res) == 11)
  res = featureModuleLookup(celdaG.sim$counts, model_G, "XXXXXXX")
  expect_true(grepl("No feature", res))
})

# cG.splitY
test_that(desc = "Testing error checking for cG.splitY", {
  r = simulateCells("celda_G", C=100, G=100, L=2)
  dc = cG.decomposeCounts(r$counts, r$y, r$L)
  res = cG.splitY(r$counts, r$y, dc$n.TS.by.C, dc$n.by.TS, dc$n.by.G, dc$nG.by.TS, dc$nM, dc$nG, r$L, beta=1, delta=1, gamma=1, y.prob=NULL, min.feature=1000)
  expect_true(grepl("Cluster sizes too small", res$message))
})


