#celda_C
library(celda)
context("Testing celda_C")

load("../celdaCsim.rda")
load("../celdaC.rda")
counts.matrix <- celdaC.sim$counts
model_C = filterCeldaList(celdaC.res, K = 5)[[1]]
factorized = factorizeMatrix(counts=counts.matrix, celda.mod = model_C)

#Checking pre-loaded data
#Making sure filterCeldaList if functioning correctly
test_that(desc = "Sanity checking filterCeldaList", {
  expect_equal(celdaC.res$content.type, class(model_C))
})

#distinct_colors#
test_that(desc = "Checking distinct_colors", {
  expect_equal(distinct_colors(2), c("#FF4D4D", "#4DFFFF"))
})

#Convenience functions#
test_that(desc = "Checking finalClusterAssignment, celdaC", {
  expect_true(all(finalClusterAssignment(celda.mod = model_C) <= 5))
})

test_that(desc = "Checking clusterProbability, celdaC", {
  expect_true(all(rowSums(clusterProbability(model_C, counts = counts.matrix)[[1]]) == 1))
})


test_that(desc = "Checking getK", {
  expect_equal(5, getK(celda.mod = model_C))
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
  expect_lt(calculateLoglikFromVariables.celda_C(counts = celdaC.sim$counts,z = celdaC.sim$z,
                                                 K = celdaC.sim$K, alpha = 1, beta = 1,sample.label = celdaC.sim$sample.label), 0)
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
  top.rank <- topRank(matrix = factorized$proportions$gene.states, n = 1000)
  expect_equal(names(top.rank),
               c("index","names"))
})


#celdaHeatmap#
test_that(desc = "Checking renderCeldaHeatmap to see if it runs without errors", {
  expect_equal(names(celdaHeatmap(celda.mod = model_C, counts = celdaC.sim$counts)),
               c("tree_row", "tree_col", "gtable"))
})

# celdaProbabilityMap
test_that(desc = "Testing celdaProbabiltyMap.celda_C for sample",{
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
  plot.obj = plotDrCluster(tsne[,1], tsne[,2], model_C$z)
  expect_true(ncol(tsne) == 2 & nrow(tsne) == length(model_C$z))
  expect_true(!is.null(plot.obj))
})

# celdaTsne
test_that(desc = "Testing celdaTsne.celda_C with subset of cells",{
  expect_success(expect_error(tsne <- celdaTsne(counts=counts.matrix, celda.mod=model_C, max.cells=50, min.cluster.size=50)))
  tsne <- celdaTsne(counts=counts.matrix, celda.mod=model_C, max.cells=100, min.cluster.size=10)
  plot.obj = plotDrCluster(tsne[,1], tsne[,2], model_C$z)
  expect_true(ncol(tsne) == 2 & nrow(tsne) == length(model_C$z) && sum(!is.na(tsne[,1])) == 100)
  expect_true(!is.null(plot.obj))
})

