#celda_CG
library(celda)
library(Rtsne)
library(SummarizedExperiment)
context("Testing celda_CG")

#Loading pre-made simulatedcells/celda objects
load("../celdaCGsim.rda")
load("../celdaCG.rda")
model_CG = filterCeldaList(celdaCG.res, K = 5, L = 3)[[1]]
factorized <- factorizeMatrix(celda.mod = model_CG, counts = celdaCG.sim$counts)
counts.matrix <- celdaCG.sim$counts

#Making sure filterCeldaList if functioning correctly
test_that(desc = "Sanity checking filterCeldaList", {
  expect_equal(celdaCG.res$content.type, class(model_CG))
})

test_that(desc = "Checking clusterProbability, celdaCG", {
  clust.prob = clusterProbability(model_CG, counts = counts.matrix)
  expect_true(length(clust.prob) == 2 && ncol(clust.prob[[1]]) == 5)
})

#Making sure relationship of counts vs proportions is correct in factorize matrix
test_that(desc = "Checking factorize matrix, counts vs proportions", {
  expect_equal(TRUE, all(factorized$counts$sample.states/sum(factorized$counts$sample.states) 
                        == factorized$proportions$sample.states))
})

#Checking dimension of factorize matrix
test_that(desc = "Checking factorize matrix dimension size", {
  expect_equal(5, ncol(factorized$proportions$population.states))
  expect_equal(3, nrow(factorized$proportions$population.states))
})

test_that(desc = "simulateCells.celda_CG returns correctly typed output", {
  sim.res = simulateCells(model="celda_CG")
  expect_equal(typeof(sim.res$counts), "integer")
})

#celda_CG.R#
test_that(desc = "Making sure celda_CG runs without crashing", {
  celdacg <- simulateCells(K.to.test = 5, L = 4, model = "celda_CG", S = 3)
  celdaCG.res <- celdaGridSearch(counts = celdacg$counts, model = "celda_CG", nchains = 2, K.to.test = 5, L = c(3,5), max.iter = 15, verbose = FALSE)
  expect_equal(length(celdaCG.res$res.list[[1]]$z), ncol(celdacg$counts))
  expect_equal(length(celdaCG.res$res.list[[1]]$y), nrow(celdacg$counts)) 
})

# Ensure calculateLoglikFromVariables calculates the expected values
test_that(desc = "calculateLoglikFromVariables.celda_CG returns correct output for various params", {
  expect_lt(calculateLoglikFromVariables(model="celda_CG",
                                         y = celdaCG.sim$y, z = celdaCG.sim$z, 
                                         delta = 1, gamma = 1,  beta = 1, 
                                         alpha = 1, K = 5, L = 3, 
                                         s = celdaCG.sim$sample.label, 
                                         counts=celdaCG.sim$counts, celda.mod = model_CG),0)
})

#normalizeCounts
test_that(desc = "Making sure normalizeCounts doesn't change dimensions", {
  norm.counts <- normalizeCounts(counts.matrix)
  expect_equal(dim(norm.counts), dim(counts.matrix))
  expect_equal(rownames(norm.counts), rownames(counts.matrix))
  expect_equal(colnames(norm.counts), colnames(counts.matrix))
})


#recodeClusterY
test_that(desc = "Checking to see if recodeClusterY gives/doesn't give error", {
  expect_error(recodeClusterY(celda.mod = model_CG, from = NULL, to = ))
  expect_error(recodeClusterY(celda.mod = model_CG, from = c(1,2,3), to = c(1,2,4)))
  new.recoded <- recodeClusterY(celda.mod = model_CG, from = c(1,2,3), to = c(3,2,1))
  expect_equal(model_CG$y == 1, new.recoded$y == 3)
})

#recodeClusterZ
test_that(desc = "Checking to see if recodeClusterZ gives/doesn't give error", {
  expect_error(recodeClusterZ(celda.mod = model_CG, from = NULL, to = ))
  expect_error(recodeClusterZ(celda.mod = model_CG, from = c(1,2,3,4,5), to = c(1,2,3,4,6)))
  new.recoded <- recodeClusterZ(celda.mod = model_CG, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
  expect_equal(model_CG$z == 1, new.recoded$z == 5)
})

#compareCountMatrix
test_that(desc = "Checking CompareCountMatrix", {
  expect_true(compareCountMatrix(counts = celdaCG.sim$counts, celda.mod = model_CG))
})

##feature_selection.R##
#topRank
test_that(desc = "Checking topRank", {
  top.rank <- topRank(matrix = factorized$proportions$gene.states, n = 1000)
  # TODO: find a better way to validate lengths of topRank names
  expect_equal(nrow(counts.matrix),
               sum(sapply(top.rank$names,length)))
  expect_equal(names(top.rank),
               c("index","names"))
})

#celdaHeatmap#
test_that(desc = "Checking celdaHeatmap", {
  expect_equal(names(celdaHeatmap(celda.mod = model_CG, counts = celdaCG.sim$counts)),
               c("tree_row", "tree_col", "gtable"))
})

#moduleHeatmap
test_that(desc = "Checking moduleHeatmap to see if it runs", {
  expect_equal(names(moduleHeatmap(celdaCG.sim$counts, celda.mod = model_CG, feature.module = c(2,3))),
               c("tree_row","tree_col","gtable"))
})

#celdaProbabilityMap
test_that(desc = "Testing celdaProbabiltyMap.celda_CG for sample",{
  plot.obj = celdaProbabilityMap(counts=counts.matrix, celda.mod=model_CG, level="sample")
  expect_true(!is.null(plot.obj))
})

test_that(desc = "Testing celdaProbabiltyMap.celda_CG for cell.population",{
  plot.obj = celdaProbabilityMap(counts=counts.matrix, celda.mod=model_CG, level="cell.population")
  expect_true(!is.null(plot.obj))
})

#differentialExpression
test_that(desc = "Checking differentialExpression", {
 expect_equal(class(diffExp_K1 <- differentialExpression(counts = counts.matrix, celda.mod = model_CG, c1 = 3, log2fc.threshold = 0.5)),
		c("data.table", "data.frame"))
})

#differentialExpression, compare cluster1 and cluster2
test_that(desc = "Checking differentialExpression", {
  expect_equal(class(diffExp_K1 <- differentialExpression(counts = counts.matrix, celda.mod = model_CG, c1 = 3, c2 = 4, log2fc.threshold = 0.5)),
               c("data.table", "data.frame"))
})

#plotDimReduceCluster,State,Gene
test_that(desc = "Checking plotDimReduceCluster to see if it runs", {
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
  tsne = celdaTsne(counts=counts.matrix, celda.mod=model_CG, max.cells=length(model_CG$z))
  plot.obj = plotDimReduceCluster(tsne[,1], tsne[,2], model_CG$z)
  expect_true(ncol(tsne) == 2 & nrow(tsne) == length(model_CG$z))
  expect_true(!is.null(plot.obj))
})

# celdaTsne
test_that(desc = "Testing celdaTsne.celda_CG with subset of cells",{
  expect_success(expect_error(tsne <- celdaTsne(counts=counts.matrix, celda.mod=model_CG, max.cells=50, min.cluster.size=50)))
  tsne = celdaTsne(counts=counts.matrix, celda.mod=model_CG, max.cells=100, min.cluster.size=10)
  plot.obj = plotDimReduceCluster(tsne[,1], tsne[,2], model_CG$z)
  expect_true(ncol(tsne) == 2 & nrow(tsne) == length(model_CG$z) && sum(!is.na(tsne[,1])) == 100)
  expect_true(!is.null(plot.obj))
})


# featureModuleLookup
test_that(desc = "Testing featureModuleLookup() roundtrip", {
  res = featureModuleLookup(counts.matrix, model_CG, "Gene_1")
  expect_true(res == 1)
})
