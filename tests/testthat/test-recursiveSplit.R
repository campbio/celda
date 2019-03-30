library(celda)
context("Testing recursiveSplits")

celdaCG.sim = simulateCells("celda_CG", K=5, L=10)

test_that(desc = "Testing recursiveSplitModule", { 
  module.split = recursiveSplitModule(celdaCG.sim$counts, initial.L = 8, max.L=15, temp.K=5, z.init=NULL)
  expect_true(is(module.split, "celdaList"))
  module.split = recursiveSplitModule(celdaCG.sim$counts, initial.L = 8, max.L=15, temp.K=NULL, z.init=celdaCG.sim$z)
  expect_true(is(module.split, "celdaList"))
  module.split = recursiveSplitModule(celdaCG.sim$counts, initial.L = 8, max.L=15, temp.K=NULL, z.init=NULL)
  expect_true(is(module.split, "celdaList"))
  plotGridSearchPerplexity(module.split)
})  

test_that(desc = "Testing recursiveSplitCell", { 
  cell.split = recursiveSplitCell(celdaCG.sim$counts, initial.K = 3, max.K=8, temp.L=20, y.init=NULL)
  expect_true(is(cell.split, "celdaList"))
  cell.split = recursiveSplitCell(celdaCG.sim$counts, initial.K = 3, max.K=8, temp.L=NULL, y.init=celdaCG.sim$y)
  expect_true(is(cell.split, "celdaList"))
  cell.split = recursiveSplitCell(celdaCG.sim$counts, initial.K = 3, max.K=8, temp.L=NULL, y.init=NULL)
  expect_true(is(cell.split, "celdaList"))
  plotGridSearchPerplexity(cell.split)
})  
