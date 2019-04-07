library(celda)
context("Testing recursiveSplits")

celdaCGSim = simulateCells("celda_CG", K=5, L=10)

test_that(desc = "Testing recursiveSplitModule", { 
    moduleSplit = recursiveSplitModule(celdaCGSim$counts, initial.L = 8, 
        max.L=15, temp.K=5, z.init=NULL)
    expect_true(is(moduleSplit, "celdaList"))
    moduleSplit = recursiveSplitModule(celdaCGSim$counts, initial.L = 8, 
        max.L=15, temp.K=NULL, z.init=celdaCGSim$z)
    expect_true(is(moduleSplit, "celdaList"))
    moduleSplit = recursiveSplitModule(celdaCGSim$counts, initial.L = 8, 
        max.L=15, temp.K=NULL, z.init=NULL)
    expect_true(is(moduleSplit, "celdaList"))
    plotGridSearchPerplexity(moduleSplit)
})  

test_that(desc = "Testing recursiveSplitCell", { 
    cellSplit = recursiveSplitCell(celdaCGSim$counts, initial.K = 3, max.K=8, 
        temp.L=20, y.init=NULL)
    expect_true(is(cellSplit, "celdaList"))
    cellSplit = recursiveSplitCell(celdaCGSim$counts, initial.K = 3, max.K=8, 
        temp.L=NULL, y.init=celdaCGSim$y)
    expect_true(is(cellSplit, "celdaList"))
    cellSplit = recursiveSplitCell(celdaCGSim$counts, initial.K = 3, max.K=8, 
        temp.L=NULL, y.init=NULL)
    expect_true(is(cellSplit, "celdaList"))
    plotGridSearchPerplexity(cellSplit)
})  
