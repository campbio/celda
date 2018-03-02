# celda_CG.R
library(celda)
library(testthat)
library(Rtsne)
context("Testing celda_CG")

celdacg <- simulateCells.celda_CG(K=3, L=5)
celdaCG.res <- celda(counts=celdacg$counts, nchains=1, K = 2:4, L = 4:6, 
                     ncores=1, model="celda_CG")

##celda_CG.R##
test_that(desc = "Checking celda_CG",{
  celdacg <- simulateCells(K = 5, L = 3, model = "celda_CG")
  celdaCG.res <- celda(counts = celdacg$counts, model = "celda_CG", nchains = 2, K = 5, L = 3)
  expect_equal(length(celdaCG.res$res.list[[1]]$z),ncol(celdacg$counts))
  expect_equal(length(celdaCG.res$res.list[[1]]$y),nrow(celdacg$counts)) 
})


# celdaCG.sim <- simulateCells(K = 5, L = 3, model = "celda_CG")
# save(celdaCG.sim,file = "celdaCGsim.rda")
# celdaCG.res <- celda(counts = celdaCG.sim$counts, model = "celda_CG", nchains = 2, K = 5, L = 3)
# save(celdaCG.res, file = "celdaCG.rda")


load("../celdaCGsim.rda")
load("../celdaCG.rda")
model_CG = getModel(celdaCG.res, K = 5, L = 3)
factorized <- factorizeMatrix(model_CG, celdaCG.sim$counts)


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
test_that(desc = "Checking factorize matrix",{
  expect_equal(5, ncol(factorized$proportions$population.states))  
})
counts.matrix <- celdaCG.sim$counts


#normalizeCounts
test_that(desc = "Checking normalizeCounts",{
  norm.counts <- normalizeCounts(counts.matrix)
  expect_equal(dim(norm.counts),dim(counts.matrix))
  expect_equal(rownames(norm.counts),rownames(counts.matrix))
  expect_equal(colnames(norm.counts),colnames(counts.matrix))
})

#recodeClusterY
test_that(desc = "Checking recodeClusterY",{
  expect_error(recodeClusterY(celda.mod = model_CG, from = NULL, to = ))
})

#recodeClusterZ
test_that(desc = "Checking recodeClusterZ",{
  expect_error(recodeClusterY(celda.mod = model_CG, from = NULL, to = ))
})

#compareCountMatrix
test_that(desc = "Checking CompareCountMatrix",{
  expect_true(compareCountMatrix(count.matrix = celdaCG.sim$counts, celda.checksum = model_CG$count.checksum))
})

#distinct_colors
test_that(desc = "Checking distinct_colors",{
  expect_equal(distinct_colors(2), c("#FF9999","#99FFFF"))
})


###renderCeldaHeatmap###
test_that(desc = "Checking renderCeldaHeatmap",{
  expect_equal(names(renderCeldaHeatmap(counts = celdaCG.sim$counts, z = model_CG$z, y = model_CG$y)),
               c("tree_row","tree_col","kmeans","gtable"))
})

##feature_selection.R##
#topRank
test_that(desc = "Checking topRank",{
  expect_equal(names(topRank(fm = factorized$proportions$gene.states)),
               c("index","names"))
})

#GiniPlot
test_that(desc = "Checking GiniPlot",{
  expect_equal(class(GiniPlot(counts = celdaCG.sim$counts, celda.mod = model_CG)),
               c("gg","ggplot"))
})


#stateHeatmap
test_that(desc = "Checking stateHeatmap",{
  expect_equal(names(stateHeatmap(celdaCG.sim$counts, celda.mod = model_CG)),
               c("tree_row","tree_col","kmeans","gtable"))
})

#plotDrCluster
test_that(desc = "Checking plotDrCluster",{
  rtsne <- Rtsne::Rtsne(X = t(celdaCG.sim$counts),max_iter = 100,pca = FALSE)
  expect_equal(names(plotDrCluster(dim1 = rtsne$Y[,1], dim2 = rtsne$Y[,2],cluster = as.factor(model_CG$z))),
               c("data","layers","scales","mapping","theme","coordinates","facet","plot_env","labels","guides"))
})

#plotDrState
test_that(desc = "Checking plotDrState",{
  rtsne <- Rtsne::Rtsne(X = t(celdaCG.sim$counts),max_iter = 100,pca = FALSE)
  expect_equal(names(plotDrState(dim1 = rtsne$Y[,1], dim2 = rtsne$Y[,2],matrix = factorized$proportions$cell.states)),
               c("data","layers","scales","mapping","theme","coordinates","facet","plot_env","labels"))
})
