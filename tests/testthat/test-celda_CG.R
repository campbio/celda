# celda_CG.R
library(celda)
library(testthat)
library(Rtsne)
library(SummarizedExperiment)
context("Testing celda_CG")

##celda_CG.R##
test_that(desc = "Making sure celda_CG runs without crashing",{
  celdacg <- simulateCells(K = 5, L = 3, model = "celda_CG")
  celdaCG.res <- celdaGridSearch(counts = celdacg$counts, model = "celda_CG", nchains = 2,K = c(5,6), L = c(3,5), max.iter = 15)
  expect_equal(length(celdaCG.res$res.list[[1]]$z), ncol(celdacg$counts))
  expect_equal(length(celdaCG.res$res.list[[1]]$y), nrow(celdacg$counts)) 
})

#Loading pre-made simulatedcells/celda objects
load("../celdaCGsim.rda")
load("../celdaCG.rda")
model_CG = getModel(celdaCG.res, K = 5, L = 3)[[1]]
factorized <- factorizeMatrix(celda.mod = model_CG, counts = celdaCG.sim$counts)
counts.matrix <- celdaCG.sim$counts


#Making sure getModel if functioning correctly
test_that(desc = "Sanity checking getModel",{
  expect_equal(celdaCG.res$content.type, class(model_CG))
})

#Making sure relationship of counts vs proportions is correct in factorize matrix
test_that(desc = "Checking factorize matrix, counts vs proportions",{
  expect_equal(TRUE,all(factorized$counts$sample.states/sum(factorized$counts$sample.states) 
                        == factorized$proportions$sample.states))
})

#Checking dimension of factorize matrix
test_that(desc = "Checking factorize matrix dimension size",{
  expect_equal(5, ncol(factorized$proportions$population.states))
  expect_equal(3, nrow(factorized$proportions$population.states))
})


# Ensure calculateLoglikFromVariables calculates the expected values
test_that(desc = "calculateLoglikFromVariables.celda_CG returns correct output for various params", {
  expect_lt(calculateLoglikFromVariables(y = celdaCG.sim$y, z = celdaCG.sim$z, 
                                            delta = 1, gamma = 1,  beta = 1, 
                                            alpha = 1, K = 5, L = 3, model="celda_CG", 
                                            s = celdaCG.sim$sample.label, 
                                            counts=celdaCG.sim$counts),
               0)
})

test_that(desc = "simulateCells.celda_CG returns correctly typed output", {
  sim.res = simulateCells(model="celda_CG")
  expect_equal(typeof(sim.res$counts), "integer")
})

#normalizeCounts
test_that(desc = "Making sure normalizeCounts doesn't change dimensions",{
  norm.counts <- normalizeCounts(counts.matrix)
  expect_equal(dim(norm.counts),dim(counts.matrix))
  expect_equal(rownames(norm.counts),rownames(counts.matrix))
  expect_equal(colnames(norm.counts),colnames(counts.matrix))
})

#recodeClusterY
test_that(desc = "Checking to see if recodeClusterY gives/doesn't give error",{
  expect_error(recodeClusterY(celda.mod = model_CG, from = NULL, to = ))
  expect_error(recodeClusterY(celda.mod = model_CG, from = c(1,2,3), to = c(1,2,4)))
  new.recoded <- recodeClusterY(celda.mod = model_CG, from = c(1,2,3), to = c(3,2,1))
  expect_equal(model_CG$y == 1,new.recoded$y == 3)
})

#recodeClusterZ
test_that(desc = "Checking to see if recodeClusterZ gives/doesn't give error",{
  expect_error(recodeClusterZ(celda.mod = model_CG, from = NULL, to = ))
  expect_error(recodeClusterZ(celda.mod = model_CG, from = c(1,2,3,4,5), to = c(1,2,3,4,6)))
  new.recoded <- recodeClusterZ(celda.mod = model_CG, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
  expect_equal(model_CG$z == 1,new.recoded$z == 5)
})

#compareCountMatrix
test_that(desc = "Checking CompareCountMatrix",{
  expect_true(compareCountMatrix(count.matrix = celdaCG.sim$counts, celda.obj = model_CG))
})

#distinct_colors
test_that(desc = "Making sure distinct_colors gives expected output",{
  expect_equal(distinct_colors(2), c("#FF4D4D", "#4DFFFF"))
  expect_equal(distinct_colors(4), c("#FF4D4D", "#4DFFFF", "#FFC04D", "#4D4DFF"))
})


###renderCeldaHeatmap###
test_that(desc = "Checking renderCeldaHeatmap",{
  expect_equal(names(renderCeldaHeatmap(counts = celdaCG.sim$counts, z = model_CG$z, y = model_CG$y)),
               c("tree_row","tree_col","kmeans","gtable"))
})

##feature_selection.R##
#topRank
test_that(desc = "Checking topRank",{
  top.rank <- topRank(fm = factorized$proportions$gene.states, n = 1000)
  expect_equal(nrow(counts.matrix),
               sum(sapply(top.rank$names,FUN = length)))
  expect_equal(names(top.rank),
               c("index","names"))
})

#GiniPlot
test_that(desc = "Checking GiniPlot to see if it runs",{
  gini.plot <- GiniPlot(counts = celdaCG.sim$counts, celda.mod = model_CG)
  expect_equal(class(gini.plot),
               c("gg","ggplot"))
})


#stateHeatmap
test_that(desc = "Checking stateHeatmap to see if it runs",{
  expect_equal(names(stateHeatmap(celdaCG.sim$counts, celda.mod = model_CG)),
               c("tree_row","tree_col","kmeans","gtable"))
})

#diffExpBetweenCellStates
test_that(desc = "Checking diffExpBetweenCellStates",{
 expect_equal(class(diffExp_K1 <- diffExpBetweenCellStates(counts = counts.matrix, celda.mod = model_CG, c1 = 1)),
		c("data.table","data.frame"))
})

#plotDrCluster,State
test_that(desc = "Checking plotDrCluster to see if it runs",{
  celda.tsne <- celdaTsne(counts = celdaCG.sim$counts, max.iter = 50,celda.mod = model_CG)
  expect_equal(names(plotDrCluster(dim1 = celda.tsne[,1], dim2 = celda.tsne[,2],cluster = as.factor(model_CG$z))),
               c("data","layers","scales","mapping","theme","coordinates","facet","plot_env","labels","guides"))
  expect_equal(names(plotDrState(dim1 = celda.tsne[,1], dim2 = celda.tsne[,2],matrix = factorized$proportions$cell.states)),
               c("data","layers","scales","mapping","theme","coordinates","facet","plot_env","labels"))  
  expect_error(plotDrState(dim1 = celda.tsne[,1], dim2 = celda.tsne[,2], matrix = factorized$proportions$cell.states, distance = "char"))
})
