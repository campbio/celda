#celda_G
library(celda)
context("Testing celda_G")

load("../celdaGsim.rda")
load("../celdaG.rda")
model_G = getModel(celdaG.res, L = 5)[[1]]
factorized <- factorizeMatrix(celda.mod = model_G, counts = celdaG.sim$counts)
counts.matrix <- celdaG.sim$counts

#Making sure getModel if functioning correctly
test_that(desc = "Sanity checking getModel",{
  expect_equal(celdaG.res$content.type, class(model_G))
})

#Making sure relationship of counts vs proportions is correct in factorize matrix
test_that(desc = "Checking factorize matrix, counts vs proportions",{
  expect_equal(TRUE,all(factorized$counts$sample.states/sum(factorized$counts$sample.states) 
                        == factorized$proportions$sample.states))
})

#Checking dimension of factorize matrix
test_that(desc = "Checking factorize matrix dimension size",{
  expect_equal(5, ncol(factorized$proportions$gene.states))  
})


# Ensure calculateLoglikFromVariables calculates the expected values
test_that(desc = "calculateLoglikFromVariables.celda_G returns correct output for various params", {
  expect_lt(calculateLoglikFromVariables(counts = celdaG.sim$counts, y = celdaG.sim$y,
                                         L = celdaG.sim$L, delta = 1, gamma = 1, beta = 1, 
                                         model="celda_G"),
               0)
})


test_that(desc = "simulateCells.celda_G returns correctly typed output", {
  sim.res = simulateCells(model="celda_G")
  expect_equal(typeof(sim.res$counts), "integer")
})


#normalizeCounts
test_that(desc = "Making sure normalizeCounts doesn't change dimensions of counts matrix",{
  norm.counts <- normalizeCounts(counts.matrix)
  expect_equal(dim(norm.counts),dim(counts.matrix))
  expect_equal(rownames(norm.counts),rownames(counts.matrix))
  expect_equal(colnames(norm.counts),colnames(counts.matrix))
})

#recodeClusterY
test_that(desc = "Checking recodeClusterY gives/doesn't give error",{
  expect_error(recodeClusterY(celda.mod = model_G, from = NULL, to = ))
  expect_error(recodeClusterY(celda.mod = model_G, from=c(1,2,3,4,5), to = c(1,2,3,4,6)))
  new.recoded = recodeClusterY(celda.mod = model_G, from=c(1,2,3,4,5), to = c(5,4,3,2,1))
  expect_equal(model_G$y == 1, new.recoded$y == 5)
})

#compareCountMatrix
test_that(desc = "Checking CompareCountMatrix",{
  expect_true(compareCountMatrix(count.matrix = celdaG.sim$counts, celda.obj = model_G))
})

#distinct_colors
test_that(desc = "Checking distinct_colors",{
  expect_equal(distinct_colors(2), c("#FF4D4D", "#4DFFFF"))
})


###renderCeldaHeatmap###
test_that(desc = "Checking renderCeldaHeatmap output",{
  expect_equal(names(renderCeldaHeatmap(counts = celdaG.sim$counts, z = model_G$z, y = model_G$y)),
               c("tree_row","tree_col","kmeans","gtable"))
})

##feature_selection.R##
#topRank
test_that(desc = "Checking topRank function",{
  expect_equal(names(topRank(fm = factorized$proportions$gene.states)),
               c("index","names"))
})

#stateHeatmap
test_that(desc = "Checking stateHeatmap to see if it runs",{
  expect_equal(names(stateHeatmap(celdaG.sim$counts, celda.mod = model_G)),
               c("tree_row","tree_col","kmeans","gtable"))
})

##celda_G.R##
test_that(desc = "Making sure celda_G runs without errors",{
  celdaG.res <- celda(counts = celdaG.sim$counts, model = "celda_G", nchains = 2, L = c(5,10), max.iter = 15)
  expect_true(class(celdaG.res) == "celda_list")  # Only best chain returned by default
})

#plotDrState
test_that(desc = "Checking plotDrState",{
  celda.tsne <- celdaTsne(counts = celdaG.sim$counts,max.iter = 50,celda.mod=model_G)
  expect_equal(names(plotDrState(dim1 = celda.tsne[,1], dim2 = celda.tsne[,2],matrix = factorized$proportions$cell.states)),
               c("data","layers","scales","mapping","theme","coordinates","facet","plot_env","labels"))
})
