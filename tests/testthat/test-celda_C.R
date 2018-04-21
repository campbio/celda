# celda_C.R
library(celda)
library(Rtsne)
context("Testing celda_C")

load("../celdaCsim.rda")
load("../celdaC.rda")
model_C = getModel(celdaC.res, K = 5)[[1]]
factorized <- factorizeMatrix(celda.mod = model_C, counts = celdaC.sim$counts)
counts.matrix <- celdaC.sim$counts

#Making sure getModel if functioning correctly
test_that(desc = "Sanity checking getModel",{
  expect_equal(celdaC.res$content.type, class(model_C))
})

#Making sure relationship of counts vs proportions is correct in factorize matrix
test_that(desc = "Checking factorize matrix, counts vs proportions",{
  expect_equal(TRUE,all(factorized$counts$sample.states/sum(factorized$counts$sample.states) 
                        == factorized$proportions$sample.states))
})

#Checking dimension of factorize matrix
test_that(desc = "Checking factorize matrix dimension size",{
  expect_equal(5, nrow(factorized$proportions$sample.states))  
})


# Ensure calculateLoglikFromVariables calculates the expected values
test_that(desc = "calculateLoglikFromVariables.celda_C returns correct output for various params", {
  expect_lt(calculateLoglikFromVariables(z = celdaC.sim$z, 
                                            beta = 1, alpha = 1, K = 5, model="celda_C", 
                                            s = celdaC.sim$sample.label, 
                                            counts=celdaC.sim$counts),
               0)
})


#normalizeCounts
test_that(desc = "Checking normalizeCounts doesn't change dimensions",{
  norm.counts <- normalizeCounts(counts.matrix)
  expect_equal(dim(norm.counts),dim(counts.matrix))
  expect_equal(rownames(norm.counts),rownames(counts.matrix))
  expect_equal(colnames(norm.counts),colnames(counts.matrix))
})

#recodeClusterZ
test_that(desc = "Checking recodeClusterZ gets correct order",{
  expect_error(recodeClusterZ(celda.mod = model_C, from = NULL, to = ))
  expect_error(recodeClusterZ(celda.mod = model_C, from=c(1,2,3,4,5), to = c(1,2,3,4,6)))
  new.recoded = recodeClusterZ(celda.mod = model_C, from=c(1,2,3,4,5), to = c(5,4,3,2,1))
  expect_equal(model_C$z == 1, new.recoded$z == 5)
})

#compareCountMatrix
test_that(desc = "Checking CompareCountMatrix",{
  expect_true(compareCountMatrix(count.matrix = celdaC.sim$counts, celda.obj = model_C))
})

#distinct_colors
test_that(desc = "Checking distinct_colors",{
  expect_equal(distinct_colors(2), c("#FF9999","#99FFFF"))
})

###renderCeldaHeatmap###
test_that(desc = "Checking renderCeldaHeatmap to see if it runs without errors",{
  expect_equal(names(renderCeldaHeatmap(counts = celdaC.sim$counts, z = model_C$z, y = model_C$y)),
               c("tree_row","tree_col","kmeans","gtable"))
})

##feature_selection.R##
#topRank
test_that(desc = "Checking topRank to see if it runs without errors",{
  top.rank <- topRank(fm = factorized$proportions$gene.states, n = 1000)
  expect_equal(names(top.rank),
               c("index","names"))
})

#diffExp
test_that(desc = "Checking diffExp",{
 expect_equal(class(diffexp_K1 <- diffExp(counts = counts.matrix, celda.mod = model_C, c1 = 1)),
                c("data.table","data.frame"))
})

##celda_C.R##
test_that(desc = "Checking celda_C to see if it runs without errors",{
  celdaC.res <- celda(counts = celdaC.sim$counts, model = "celda_C",  nchains = 2, K = c(5,10), max.iter = 15)
  expect_true(class(celdaC.res) == "celda_list")  # Only best chain is returned
})
