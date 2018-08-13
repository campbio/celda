#celda_C
library(celda)
context("Testing celda_C")

load("../celdaCsim.rda")
load("../celdaC.rda")
counts.matrix <- celdaC.sim$counts
model_C = getModel(celdaC.res, K = 5)[[1]]
factorized = factorizeMatrix(counts=counts.matrix, celda.mod=model_C)

distinct_colors
test_that(desc = "Checking distinct_colors",{
  expect_equal(distinct_colors(2), c("#FF4D4D", "#4DFFFF"))
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

#diffExpBetweenCellStates
test_that(desc = "Checking diffExpBetweenCellStates",{
  diffexp_K1 <- diffExpBetweenCellStates(counts = counts.matrix, 
                                         celda.mod = model_C, c1 = 1)
  expect_equal(class(diffexp_K1), c("data.table","data.frame"))
})

##celda_C.R##
test_that(desc = "Checking celda_C to see if it runs without errors",{
  celdaC.res <- celdaGridSearch(counts = celdaC.sim$counts, model = "celda_C",  nchains = 2, K = c(5,10), max.iter = 15)
  expect_true(class(celdaC.res)[1] == "celda_list")  # Only best chain is returned
})
