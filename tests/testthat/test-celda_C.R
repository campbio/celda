#celda_C
library(celda)
context("Testing celda_C")

load("../celdaCsim.rda")
load("../celdaC.rda")
counts.matrix <- celdaC.sim$counts
model_C = filterCeldaList(celdaC.res, K = 5, chain = 1)[[1]]
factorized = factorizeMatrix(counts=counts.matrix, celda.mod = model_C)

#Checking pre-loaded data
#Making sure filterCeldaList if functioning correctly
test_that(desc = "Sanity checking filterCeldaList", {
  expect_equal(celdaC.res$content.type, class(model_C))
  expect_equal(celdaC.res$content.type, 
               class(filterCeldaList(celdaC.res, K = 5, index = 2)[[1]]))
})

#selectBestModel
test_that(desc = "Checking selectBestModel", {
  expect_error(selectBestModel(model_C, K = 5), "Provided object is not of class celda_list")
})

#distinct_colors#
test_that(desc = "Checking distinct_colors", {
  expect_equal(distinct_colors(2), c("#FF4D4D", "#4DFFFF"))
})

test_that(desc = "Checking clusterProbability, celdaC", {
  expect_true(all(rowSums(clusterProbability(model_C, counts = counts.matrix)[[1]]) == 1))
})


#simulateCells
test_that(desc = "simulateCells.celda_C returns correctly typed output", {
  sim.res = simulateCells(model = "celda_C")
  expect_equal(typeof(sim.res$counts), "integer")
})

#celda_C.R#
#celdaGridSearch
test_that(desc = "Checking celda_C to see if it runs without errors", {
  celdaC.res <- celdaGridSearch(counts = celdaC.sim$counts, model = "celda_C",  nchains = 2, K.to.test = c(5,10), max.iter = 15, verbose = F)
  expect_true(class(celdaC.res)[1] == "celda_list")  # Only best chain is returned
})

#Ensure calculateLoglikFromVariables calculates the expected values
test_that(desc = "calculateLoglikFromVariables.celda_C returns correct output for various params", {
  expect_lt(calculateLoglikFromVariables(celda.mod="celda_C", 
                                         counts = celdaC.sim$counts,z = celdaC.sim$z,
                                         K = celdaC.sim$K, alpha = 1, beta = 1,
                                         sample.label = celdaC.sim$sample.label), 0)
})


# Gibbs sampling
test_that(desc = "Test Gibbs sampling for celda_C", {
  res <- celda_C(celdaC.sim$counts, sample.label=celdaC.sim$sample.label, K=celdaC.sim$K, algorithm="Gibbs", max.iter=5, nchain=1)
  expect_is(res, "celda_C")
})


#normalizeCounts
test_that(desc = "Making sure normalizeCounts doesn't change dimensions of counts matrix", {
  norm.counts <- normalizeCounts(counts.matrix)
  expect_equal(dim(norm.counts),dim(counts.matrix))
  expect_equal(rownames(norm.counts),rownames(counts.matrix))
  expect_equal(colnames(norm.counts),colnames(counts.matrix))
})

#feature_selection.R#
#topRank
test_that(desc = "Checking topRank to see if it runs without errors", {
  factorized = factorizeMatrix(counts=counts.matrix, celda.mod = model_C)
  top.rank <- topRank(matrix = factorized$proportions$gene.states, threshold = NULL)
  expect_equal(names(top.rank),
               c("index","names"))
  top.rank <- topRank(matrix = factorized$proportions$gene.states, n = 1000)
  expect_equal(names(top.rank),
               c("index","names"))
})

#plotHeatmap, testing for errors
test_that(desc = "plotHeatmap, testing for errors", {
  expect_error(plotHeatmap(counts = celdaC.sim$counts, z = model_C$K), "Length of z must match number of columns in counts matrix")
  expect_error(plotHeatmap(counts = celdaC.sim$counts, z = model_C$z, scale.row = model_C), "'scale.row' needs to be of class 'function'")
  expect_error(plotHeatmap(counts = celdaC.sim$counts, z = model_C$z, trim = 3), "'trim' should be a 2 element vector specifying the lower and upper boundaries")
})


#plotHeatmap, annotation.cell
test_that(desc = "plotHeatmap, testing for annotation.cell",{
  annot <- as.data.frame(c(rep(x = 1, times = ncol(celdaC.sim$counts) - 100),rep(x = 2, 100)))
  rownames(annot) <- colnames(celdaC.sim$counts)
  
  expect_equal(names(plotHeatmap(celda.mod = model_C, counts = celdaC.sim$counts, annotation.cell = annot, z = model_C$z)),
               c("tree_row", "tree_col", "gtable"))
  
  rownames(annot) <- rev(colnames(celdaC.sim$counts))
  expect_error(plotHeatmap(celda.mod = model_C, counts = celdaC.sim$counts, annotation.cell = annot, z = model_C$z),
               "Row names of 'annotation.cell' are different than the column names of 'counts'")
})

#celdaHeatmap#
test_that(desc = "Checking celdaHeatmap to see if it runs without errors", {
  expect_equal(names(celdaHeatmap(celda.mod = model_C, counts = celdaC.sim$counts)),
               c("tree_row", "tree_col", "gtable"))
})

# celdaProbabilityMap
test_that(desc = "Testing celdaProbabiltyMap.celda_C for sample",{
  plot.obj = celdaProbabilityMap(counts=counts.matrix, celda.mod=model_C, level="sample")
  expect_true(!is.null(plot.obj))
  
  model_C = celda_C(celdaC.sim$counts, sample.label=celdaC.sim$sample.label, K=celdaC.sim$K, max.iter=5, nchain=1)
  plot.obj = celdaProbabilityMap(counts=counts.matrix, celda.mod=model_C, level="sample")
  expect_true(!is.null(plot.obj))  
})

#differentialExpression#
test_that(desc = "Checking differentialExpression", {
  diffexp_K1 <- differentialExpression(counts = counts.matrix, 
                                         celda.mod = model_C, c1 = 1)
  expect_equal(class(diffexp_K1), c("data.table", "data.frame"))
})

# celdaTsne
test_that(desc = "Testing celdaTsne.celda_C with all cells",{
  tsne = celdaTsne(counts=counts.matrix, celda.mod=model_C, max.cells=length(model_C$z), min.cluster.size=50)
  plot.obj = plotDimReduceCluster(tsne[,1], tsne[,2], model_C$z)
  expect_true(ncol(tsne) == 2 & nrow(tsne) == length(model_C$z))
  expect_true(!is.null(plot.obj))
})

# celdaTsne
test_that(desc = "Testing celdaTsne.celda_C with subset of cells",{
  expect_success(expect_error(tsne <- celdaTsne(counts=counts.matrix, celda.mod=model_C, max.cells=50, min.cluster.size=50)))
  tsne <- celdaTsne(counts=counts.matrix, celda.mod=model_C, max.cells=100, min.cluster.size=10)
  plot.obj = plotDimReduceCluster(tsne[,1], tsne[,2], model_C$z)
  expect_true(ncol(tsne) == 2 & nrow(tsne) == length(model_C$z) && sum(!is.na(tsne[,1])) == 100)
  expect_true(!is.null(plot.obj))
})

# featureModuleLookup
test_that(desc = "Testing featureModuleLookup() fails for celda_C models", {
  expect_error(featureModuleLookup(counts.matrix, model_C, "test_feat"))
})

# cC.splitZ
test_that(desc = "Testing cC.splitZ", {
  r = simulateCells("celda_C", S=1, C.Range=c(50,100), K=2)
  dc = cC.decomposeCounts(r$counts, r$sample.label, r$z, r$K)
  res = cC.splitZ(r$counts, dc$m.CP.by.S, dc$n.G.by.CP, dc$n.CP, s=as.integer(r$sample.label), z=r$z, K=r$K, nS=dc$nS, nG=dc$nG, alpha=1, beta=1, z.prob=NULL, min.cell=1000)
  expect_true(grepl("Cluster sizes too small", res$message))
})
