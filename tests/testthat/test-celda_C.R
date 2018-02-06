# celda_C.R
library(celda)
library(Rtsne)
context("Testing celda_C")

#celdac <- simulateCells.celda_C(K=10)
#celdaC.res <- celda(counts=celdac$counts, model="celda_C", nchains=1, K=10)

#test_that("CheckingVisualizeModelPerformace",{
#	expect_equal(TRUE, all(!is.na(visualizeModelPerformance(celdaC.res))))
#	})
#test_that("finalClusterAssignment.celda_C",{
#  expect_equal(celdaC.res$res.list[[1]]$z, finalClusterAssignment(celdaC.res$res.list[[1]]))
#})


celdaC.sim <- simulateCells(K = 5, model = "celda_C")
celdaC.res <- celda(counts = celdaC.sim$counts, model = "celda_C", nchains = 2, K = 5)
model = getModel(celdaC.res, K = 5, best = "loglik")
factorized <- factorizeMatrix(celda.mod = model, counts = celdaC.sim$counts)

#Making sure getModel if functioning correctly
test_that(desc = "Sanity checking getModel",{
  expect_equal(celdaC.res$content.type, class(model))
})

#Making sure relationship of counts vs proportions is correct in factorize matrix
test_that(desc = "Checking factorize matrix, counts vs proportions",{
  expect_equal(TRUE,all(factorized$counts$sample.states/sum(factorized$counts$sample.states) 
                        == factorized$proportions$sample.states))
})

#Checking dimension of factorize matrix
test_that(desc = "Checking factorize matrix",{
  expect_equal(5, nrow(factorized$proportions$sample.states))  
})


#01-31-18 
#Creating unit tests for celda_functions
counts.matrix <- celdaC.sim$counts

# #cosineDist
# test_that(desc = "Checking cosineDist",{
#   expect_equal("dist",class(cosineDist(factorize.matrix$posterior$population.states)))
# })
# 
# #cosine
# test_that(desc = "Checking cosineDist",{
#   expect_equal("matrix",class(cosine(t(factorize.matrix$posterior$population.states))))
# })
# 
# #spearmanDist
# test_that(desc = "Checking spearmanDist",{
#   expect_equal("dist",class(spearmanDist(factorize.matrix$posterior$population.states)))
# })

#normalizeLogProbs


#normalizeCounts
test_that(desc = "Checking normalizeCounts",{
  norm.counts <- normalizeCounts(counts.matrix)
  expect_equal(dim(norm.counts),dim(counts.matrix))
  expect_equal(rownames(norm.counts),rownames(counts.matrix))
  expect_equal(colnames(norm.counts),colnames(counts.matrix))
})

# #reorder.label.by.size
# test_that(desc = "Checking reorder.label.by.size",{
#   expect_true(all(table((reorder.label.by.size(model$z,model$K))$new.labels)[1] >= 
#                     table((reorder.label.by.size(model$z,model$K))$new.labels)))
# })
# 
# #reorder.labels.by.size.then.counts
# test_that(desc = "Checking reorder.labels.by.size.then.counts",{
#   expect_true(all(table((reorder.labels.by.size.then.counts(counts = celdaC$counts,z = model$z,y = model$y, K = model$K, L = model$L))$new.labels)[1] >= 
#                     table((reorder.labels.by.size.then.counts(counts = celdaC$counts,z = model$z,y = model$y, K = model$K, L = model$L))$new.labels)))
# })

#recodeClusterY
test_that(desc = "Checking recodeClusterY",{
  expect_error(recodeClusterY(celda.mod = model, from = NULL, to = ))
})

#recodeClusterZ
test_that(desc = "Checking recodeClusterZ",{
  expect_error(recodeClusterY(celda.mod = model, from = NULL, to = ))
})

#compareCountMatrix
test_that(desc = "Checking CompareCountMatrix",{
  expect_true(compareCountMatrix(count.matrix = celdaC.sim$counts, celda.checksum = model$count.checksum))
})

#distinct_colors
test_that(desc = "Checking distinct_colors",{
  expect_equal(distinct_colors(2), c("#FF9999","#99FFFF"))
})

#initialize.cluster
#test_that(desc = "Checking distinct_colors",{
#  expect_equal(initialize.cluster(N = 3, len = 3), c(2,3,1))
#})


###renderCeldaHeatmap###
test_that(desc = "Checking renderCeldaHeatmap",{
  expect_equal(names(renderCeldaHeatmap(counts = celdaC.sim$counts, z = model$z, y = model$y)),
               c("tree_row","tree_col","kmeans","gtable"))
})

##feature_selection.R##
#topRank
test_that(desc = "Checking topRank",{
  expect_equal(names(topRank(fm = factorized$proportions$gene.states)),
               c("index","names"))
})

##celda_C.R##
test_that(desc = "Checking celda_C",{
  celdaC.res <- celda(counts = celdaC.sim$counts, model = "celda_C", nchains = 2, K = 5)
  expect_equal(celdaC.res$run.params$chain,c(1,2))
})




#plotDrCluster
test_that(desc = "Checking plotDrCluster",{
  rtsne <- Rtsne::Rtsne(X = t(celdaC.sim$counts),max_iter = 500,pca = FALSE)
  expect_equal(names(plotDrCluster(dim1 = rtsne$Y[,1], dim2 = rtsne$Y[,2],cluster = as.factor(model$z))),
               c("data","layers","scales","mapping","theme","coordinates","facet","plot_env","labels","guides"))
})