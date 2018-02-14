#celda_G
library(celda)
context("Testing celda_G")

load("../celdaGsim.rda")
load("../celdaG.rda")
model_G = getModel(celdaG.res, L = 5, best = "loglik")
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
  expect_true(compareCountMatrix(count.matrix = celdaG.sim$counts, celda.checksum = model_G$count.checksum))
})

#distinct_colors
test_that(desc = "Checking distinct_colors",{
  expect_equal(distinct_colors(2), c("#FF9999","#99FFFF"))
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



##celda_G.R##
test_that(desc = "Making sure celda_G runs without errors",{
  celdaG.res <- celda(counts = celdaG.sim$counts, model = "celda_G", nchains = 2, L = 5)
  expect_equal(celdaG.res$run.params$chain,c(1,2))
})

#plotDrState
test_that(desc = "Checking plotDrState",{
  rtsne <- Rtsne::Rtsne(X = t(celdaG.sim$counts),max_iter = 100,pca = FALSE)
  expect_equal(names(plotDrState(dim1 = rtsne$Y[,1], dim2 = rtsne$Y[,2],matrix = factorized$proportions$cell.states)),
               c("data","layers","scales","mapping","theme","coordinates","facet","plot_env","labels"))
})