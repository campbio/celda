#celda_C
library(celda)
context("Testing celda_C")

load("../celdaCsim.rda")
load("../celdaC.rda")
counts.matrix <- celdaC.sim$counts
model_C = filterCeldaList(celdaC.res, K = 5)[[1]]
factorized = factorizeMatrix(counts=counts.matrix, celda.mod=model_C)

##distinct_colors##
test_that(desc = "Checking distinct_colors",{
  expect_equal(distinct_colors(2), c("#FF4D4D", "#4DFFFF"))
})

###Convenience functions###

test_that(desc = "Checking finalClusterAssignment, celdaC",{
  expect_true(all(finalClusterAssignment(celda.mod = model_C) <= 5))
})

test_that(desc = "Checking clusterProbability, celdaC",{
  expect_true(all(rowSums(clusterProbability(model_C, counts = counts.matrix)[[1]]) == 1))
})


test_that(desc = "Checking getK",{
  expect_equal(5,getK(celda.mod = model_C))
})

###celdaHeatmap###
test_that(desc = "Checking renderCeldaHeatmap to see if it runs without errors",{
  expect_equal(names(celdaHeatmap(celda.mod = model_C, counts = celdaC.sim$counts)),
               c("tree_row","tree_col","kmeans","gtable"))
})

##feature_selection.R##
#topRank
test_that(desc = "Checking topRank to see if it runs without errors",{
  top.rank <- topRank(matrix = factorized$proportions$gene.states, n = 1000)
  expect_equal(names(top.rank),
               c("index","names"))
})

#differentialExpression
test_that(desc = "Checking differentialExpression",{
  diffexp_K1 <- differentialExpression(counts = counts.matrix, 
                                         celda.mod = model_C, c1 = 1)
  expect_equal(class(diffexp_K1), c("data.table","data.frame"))
})


test_that(desc = "simulateCells.celda_C returns correctly typed output", {
  sim.res = simulateCells(model="celda_C")
  expect_equal(typeof(sim.res$counts), "integer")
})

# Ensure calculateLoglikFromVariables calculates the expected values
test_that(desc = "calculateLoglikFromVariables.celda_C returns correct output for various params", {
  expect_lt(calculateLoglikFromVariables.celda_C(counts = celdaC.sim$counts,z = celdaC.sim$z,
                                                 K = celdaC.sim$K, alpha = 1, beta = 1,sample.label = celdaC.sim$sample.label 
  ),
  0)
})

##celda_C.R##
test_that(desc = "Checking celda_C to see if it runs without errors",{
  celdaC.res <- celdaGridSearch(counts = celdaC.sim$counts, celda.mod = "celda_C",  nchains = 2, K.to.test = c(5,10), max.iter = 15)
  expect_true(class(celdaC.res)[1] == "celda_list")  # Only best chain is returned
})


##celdaTsne##
test_that(desc = "testing celdaTsne", {
  celdaC.tsne <-celdaTsne(counts = celdaC.sim$counts,celda.mod = model_C, max.iter = 50, max.cells = 500)
  plotDrCluster(dim1 = celdaC.tsne[,1], dim2 = celdaC.tsne[,2], cluster = model_C$z)
  expect_equal(class(plotDrCluster(dim1 = celdaC.tsne[,1], dim2 = celdaC.tsne[,2], cluster = model_C$z)), c("gg","ggplot"))
})

