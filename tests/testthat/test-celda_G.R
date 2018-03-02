#celda_G
library(celda)
context("Testing celda_G")

#celdaG 
# celdaG.sim <- simulateCells(L = 5, model = "celda_G")
# save(celdaG.sim,file = "celdaGsim.rda")
# celdaG.res <- celda(counts = celdaG.sim$counts, model = "celda_G", nchains = 2, L = 5)
# save(celdaG.res, file = "celdaG.rda")


load("../celdaGsim.rda")
load("../celdaG.rda")
model_G = getModel(celdaG.res, L = 5, best = "loglik")
factorized <- factorizeMatrix(celda.mod = model_G, counts = celdaG.sim$counts)

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
test_that(desc = "Checking factorize matrix",{
  expect_equal(5, ncol(factorized$proportions$gene.states))  
})


#01-31-18 
#Creating unit tests for celda_functions
counts.matrix <- celdaG.sim$counts

# #cosineDist
# test_that(desc = "Checking cosineDist",{
#   expect_equal("dist",class(cosineDist(factorized$posterior$population.states)))
# })
# 
# #cosine
# test_that(desc = "Checking cosineDist",{
#   expect_equal("matrix",class(cosine(t(factorized$posterior$population.states))))
# })
# 
# #spearmanDist
# test_that(desc = "Checking spearmanDist",{
#   expect_equal("dist",class(spearmanDist(factorized$posterior$population.states)))
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
#   expect_true(all(table((reorder.labels.by.size.then.counts(counts = celdaG$counts,z = model$z,y = model$y, L = model$K, L = model$L))$new.labels)[1] >= 
#                     table((reorder.labels.by.size.then.counts(counts = celdaG$counts,z = model$z,y = model$y, L = model$K, L = model$L))$new.labels)))
# })

#recodeClusterY
test_that(desc = "Checking recodeClusterY",{
  expect_error(recodeClusterY(celda.mod = model_G, from = NULL, to = ))
})

#recodeClusterZ
test_that(desc = "Checking recodeClusterZ",{
  expect_error(recodeClusterY(celda.mod = model_G, from = NULL, to = ))
})

#compareCountMatrix
test_that(desc = "Checking CompareCountMatrix",{
  expect_true(compareCountMatrix(count.matrix = celdaG.sim$counts, celda.obj = model_G))
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
  expect_equal(names(renderCeldaHeatmap(counts = celdaG.sim$counts, z = model_G$z, y = model_G$y)),
               c("tree_row","tree_col","kmeans","gtable"))
})

##feature_selection.R##
#topRank
test_that(desc = "Checking topRank",{
  expect_equal(names(topRank(fm = factorized$proportions$gene.states)),
               c("index","names"))
})



##celda_G.R##
test_that(desc = "Checking celda_G",{
  celdaG.res <- celda(counts = celdaG.sim$counts, model = "celda_G", nchains = 2, L = 5)
  expect_equal(celdaG.res$run.params$chain,c(1,2))
})

#plotDrState
test_that(desc = "Checking plotDrState",{
  rtsne <- Rtsne::Rtsne(X = t(celdaG.sim$counts),max_iter = 100,pca = FALSE)
  expect_equal(names(plotDrState(dim1 = rtsne$Y[,1], dim2 = rtsne$Y[,2],matrix = factorized$proportions$cell.states)),
               c("data","layers","scales","mapping","theme","coordinates","facet","plot_env","labels"))
})
