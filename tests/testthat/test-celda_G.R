#celda_G
library(celda)
context("Testing celda_G")

load("../celdaGsim.rda")
load("../celdaG.rda")
model_G = filterCeldaList(celdaG.res, L = 5)[[1]]
factorized <- factorizeMatrix(celda.mod = model_G, counts = celdaG.sim$counts)
counts.matrix <- celdaG.sim$counts

#Checking pre-loaded data
#Making sure filterCeldaList if functioning correctly
test_that(desc = "Sanity checking filterCeldaList", {
  expect_equal(celdaG.res$content.type, class(model_G))
})

#Making sure relationship of counts vs proportions is correct in factorize matrix
test_that(desc = "Checking factorize matrix, counts vs proportions", {
  expect_equal(TRUE,all(factorized$counts$sample.states/sum(factorized$counts$sample.states) 
                        == factorized$proportions$sample.states))
})

#Checking dimension of factorize matrix
test_that(desc = "Checking factorize matrix dimension size", {
  expect_equal(5, ncol(factorized$proportions$gene.states))  
})

test_that(desc = "Checking clusterProbability, celdaG", {
  expect_true(ncol(clusterProbability(model_G, counts = counts.matrix)[[1]]) == 5)
})


test_that(desc = "Checking getL", {
  expect_equal(5, getL(celda.mod = model_G))
})

test_that(desc = "simulateCells.celda_G returns correctly typed output", {
  sim.res = simulateCells(model = "celda_G")
  expect_equal(typeof(sim.res$counts), "integer")
})

#celda_G.R#
test_that(desc = "Making sure celda_G runs without errors", {
  celdaG.res <- celdaGridSearch(counts = celdaG.sim$counts, model = "celda_G", nchains = 2, L = c(5,10), max.iter = 15, verbose = F)
  expect_true(all(class(celdaG.res) == c("celda_list", "celda_G")))  # Only best chain returned by default
})

# Ensure calculateLoglikFromVariables calculates the expected values
test_that(desc = "calculateLoglikFromVariables.celda_G returns correct output for various params", {
  expect_lt(calculateLoglikFromVariables(celda.mod = "celda_G", 
                                         counts = celdaG.sim$counts, 
                                         y = celdaG.sim$y, L = celdaG.sim$L, delta = 1, 
                                         gamma = 1, beta = 1),0)
})

#normalizeCounts
test_that(desc = "Making sure normalizeCounts doesn't change dimensions of counts matrix", {
  norm.counts <- normalizeCounts(counts.matrix)
  expect_equal(dim(norm.counts), dim(counts.matrix))
  expect_equal(rownames(norm.counts), rownames(counts.matrix))
  expect_equal(colnames(norm.counts), colnames(counts.matrix))
})

#recodeClusterY
test_that(desc = "Checking recodeClusterY gives/doesn't give error", {
  expect_error(recodeClusterY(celda.mod = model_G, from = NULL, to = ))
  expect_error(recodeClusterY(celda.mod = model_G, from = c(1,2,3,4,5), to = c(1,2,3,4,6)))
  new.recoded = recodeClusterY(celda.mod = model_G, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
  expect_equal(model_G$y == 1, new.recoded$y == 5)
})

#compareCountMatrix
test_that(desc = "Checking CompareCountMatrix", {
  expect_true(compareCountMatrix(counts = celdaG.sim$counts, celda.mod = model_G))
})

##feature_selection.R##
#topRank
test_that(desc = "Checking topRank function", {
  top.rank <- topRank(matrix = factorized$proportions$gene.states, n = 1000, threshold = NULL)
  expect_equal(names(top.rank),
               c("index","names"))
  expect_equal(names(topRank(matrix = factorized$proportions$gene.states)),
               c("index","names"))
})

#feature_selection.R#

###celdaHeatmap###
test_that(desc = "Checking celdaHeatmap output",{
  expect_equal(names(celdaHeatmap(celda.mod = model_G, counts = celdaG.sim$counts)),
               c("tree_row","tree_col","gtable"))
})

#moduleHeatmap
test_that(desc = "Checking moduleHeatmap to see if it runs",{
  expect_equal(names(moduleHeatmap(celdaG.sim$counts, celda.mod = model_G)),
               c("tree_row","tree_col","gtable"))
})

#plotDimReduceState
test_that(desc = "Checking plotDimReduceState", {
  celda.tsne <- celdaTsne(counts = celdaG.sim$counts,max.iter = 50,celda.mod = model_G)
  expect_equal(names(plotDimReduceState(dim1 = celda.tsne[,1], dim2 = celda.tsne[,2], counts = celdaG.sim$counts, celda.mod = model_G)),
               c("data", "layers", "scales", "mapping", "theme", "coordinates", "facet", "plot_env", "labels"))
})

# celdaTsne
test_that(desc = "Testing celdaTsne.celda_G with all cells",{
  tsne = celdaTsne(counts=counts.matrix, celda.mod=model_G, max.cells=ncol(counts.matrix))
  plot.obj = plotDimReduceCluster(tsne[,1], tsne[,2], rep(1,ncol(counts.matrix)))
  expect_true(ncol(tsne) == 2 & nrow(tsne) == ncol(counts.matrix))
  expect_true(!is.null(plot.obj))
})

# celdaTsne
test_that(desc = "Testing celdaTsne.celda_G with subset of cells",{
  tsne = celdaTsne(counts=counts.matrix, celda.mod=model_G, max.cells=100)
  plot.obj = plotDimReduceCluster(tsne[,1], tsne[,2], rep(1, ncol(counts.matrix)))
    expect_true(ncol(tsne) == 2 & nrow(tsne) == ncol(counts.matrix) && sum(!is.na(tsne[,1])) == 100)
  expect_true(!is.null(plot.obj))
})

# featureModuleLookup
test_that(desc = "Testing featureModuleLookup() roundtrip", {
  res = featureModuleLookup(counts.matrix, model_G, "Gene_1")
  expect_true(res == 5)
})

# cG.splitY
test_that(desc = "Testing cG.splitY", {
  r = simulateCells("celda_G", C=100, G=100, L=2)
  dc = cG.decomposeCounts(r$counts, r$y, r$L)
  res = cG.splitY(r$counts, r$y, dc$n.TS.by.C, dc$n.by.TS, dc$n.by.G, dc$nG.by.TS, dc$nM, dc$nG, r$L, beta=1, delta=1, gamma=1, y.prob=NULL, min.feature=1000)
  expect_true(grepl("Cluster sizes too small", res$message))
})


