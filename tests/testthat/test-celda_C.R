#celda_C
library(celda)
context("Testing celda_C")

celdaC.sim = simulateCells("celda_C", K=10)
model_C = celda_C(counts=celdaC.sim$counts, sample.label=celdaC.sim$sample.label, 
                  K=celdaC.sim$K, algorithm="EM", verbose=FALSE,
                  max.iter=2)
factorized = factorizeMatrix(counts=celdaC.sim$counts, celda.mod = model_C)  

# celda_C
test_that(desc = "Testing simulation and celda_C model", {
  expect_equal(typeof(celdaC.sim$counts), "integer")
  expect_true(all(sweep(factorized$counts$sample, 2, colSums(factorized$counts$sample), "/") == factorized$proportions$sample))  
  expect_true(ncol(factorized$proportions$module) == model_C@params$K)
  expect_true(all(is.numeric(logLikelihoodHistory(celda.mod = model_C))))
  expect_equal(max(logLikelihoodHistory(celda.mod = model_C)), bestLogLikelihood(model_C))
  
  # GitHub #347
  numeric.counts = celdaC.sim$counts
  storage.mode(numeric.counts) = "numeric"
  expect_true(is(celda_C(counts=celdaC.sim$counts, 
                             sample.label=celdaC.sim$sample.label, 
                             K=celdaC.sim$K, algorithm="EM", verbose=FALSE,
                             max.iter=2),
               "celda_C"))
})

# clusterProbability
test_that(desc = "Testing clusterProbability with celda_C", {
  expect_true(all(round(rowSums(clusterProbability(model_C, counts = celdaC.sim$counts)[[1]])) == 1))
})

# celdaGridSearch and perplexity calculations
test_that(desc = "Testing celdaGridSearch with celda_C", {
  celdaC.res <- celdaGridSearch(counts = celdaC.sim$counts, model = "celda_C",  
                                nchains = 2, params.test=list(K=c(5,6)), 
                                params.fixed=list(sample.label=celdaC.sim$sample.label),
                                max.iter = 2, verbose = FALSE, best.only=FALSE)
  expect_error(celdaGridSearch(counts=celdaC.sim$counts, model="celda_C", 
                               params.test=list(K=4:5, M = 3:4), 
                               params.fixed=list(sample.label=celdaC.sim$sample.label),
                               best.only=FALSE),
               "The following elements in 'params.test' are not arguments of 'celda_C': M")
  
  expect_error(celdaGridSearch(counts=celdaC.sim$counts, 
                               model="celda_C", nchains = 1, max.iter=1,
                               params.test=list(K=4:5, sample.label = "Sample"), 
                               params.fixed=list(sample.label=celdaC.sim$sample.label)),
               "Setting parameters such as 'z.init', 'y.init', and 'sample.label' in 'params.test' is not currently supported.")
  
  expect_error(celdaGridSearch(counts=celdaC.sim$counts, 
                               model="celda_C", nchains = 1, max.iter=1,
                               params.test=list(), 
                               params.fixed=list(sample.label=celdaC.sim$sample.label)),
               "The following arguments are not in 'params.test' or 'params.fixed' but are required for 'celda_C': K")
  
  expect_error(celdaGridSearch(counts=celdaC.sim$counts, model="celda_C", 
                               nchains = 1, max.iter = 1,
                               params.test=list(K=9:10), 
                               params.fixed=list(sample.label=celdaC.sim$sample.label,
                                                 xxx = "xxx")),
               "The following elements in 'params.fixed' are not arguments of 'celda_C': xxx")
  
  expect_true(class(celdaC.res)[1] == "celdaList")
  expect_equal(names(runParams(celdaC.res)), c("index","chain","K","log_likelihood"))
  expect_error(plotGridSearchPerplexity(celdaC.res))

  celdaC.res = resamplePerplexity(celdaC.sim$counts, celdaC.res, resample=2)
  expect_equal(is.null(celdaC.res@perplexity), FALSE)
  expect_true(is(celdaC.res, "celdaList"))
  expect_error(resamplePerplexity(celdaC.sim$counts, celdaC.res, resample="2"))
  expect_error(resamplePerplexity(celdaC.sim$counts, "celdaC.res", resample=2))
  
  plot.obj = plotGridSearchPerplexity(celdaC.res)
  expect_is(plot.obj, "ggplot")

  celdaC.res.index1 = subsetCeldaList(celdaC.res, params=list(index = 1))
  expect_true(all(is(celdaC.res.index1, "celda_C") && !is(celdaC.res.index1, "celda_list")))
  
  expect_error(subsetCeldaList(celdaC.res, params = list(K = 11)))
  expect_error(subsetCeldaList(celdaC.res, params = list(K = 5, M = 10)))
  
  celdaC.res.K5 <- subsetCeldaList(celdaC.res, params=list(K = 5))
  model_C_2 = selectBestModel(celdaC.res.K5)
  res <- perplexity(celdaC.sim$counts, model_C)
  res2 <- perplexity(celdaC.sim$counts, model_C, new.counts = celdaC.sim$counts + 1)
  
  expect_error(res <- perplexity.celda_C(celdaC.sim$counts, model_C, new.counts = celdaC.sim$counts[-1,]))
})

# logLikelihood
test_that(desc = "Testing logLikelihood.celda_C", {
  expect_lt(logLikelihood(model="celda_C", 
                          counts = celdaC.sim$counts,z = celdaC.sim$z,
                          K = celdaC.sim$K, alpha = 1, beta = 1,
                          sample.label = celdaC.sim$sample.label), 0)
  
  fake.z = celdaC.sim$z
  fake.z[1] = celdaC.sim$K + 1
  expect_error(logLikelihood(model="celda_C",
                             z = fake.z, counts = celdaC.sim$counts,
                             K = celdaC.sim$K, alpha = 1, beta = 1,
                             sample.label = celdaC.sim$sample.label),
                             "An entry in z contains a value greater than the provided K.")
})

# Gibbs sampling
test_that(desc = "Testing celda_C with Gibbs sampling", {
  res <- celda_C(counts=celdaC.sim$counts, sample.label=celdaC.sim$sample.label, K=celdaC.sim$K, algorithm="Gibbs", max.iter=5, nchain=1)
  expect_is(res, "celda_C")
})


# normalizeCounts
test_that(desc = "Making sure normalizeCounts doesn't change dimensions of counts matrix", {
  norm.counts <- normalizeCounts(celdaC.sim$counts)
  expect_equal(dim(norm.counts),dim(celdaC.sim$counts))
  expect_equal(rownames(norm.counts),rownames(celdaC.sim$counts))
  expect_equal(colnames(norm.counts),colnames(celdaC.sim$counts))
  expect_error(normalizeCounts(celdaC.sim$counts, transformation.fun = "scale"), 
               "'transformation.fun' needs to be of class 'function'")
  expect_error(normalizeCounts(celdaC.sim$counts, scale.fun = "scale"), 
               "'scale.fun' needs to be of class 'function'")
})

# recodeClusterZ
test_that(desc = "Testing recodeClusterZ with celda_C", {
  expect_error(recodeClusterY(celda.mod = model_C, from = c(1,2,3,4,5), to = c(5,4,3,2,1)))
  expect_error(recodeClusterZ(celda.mod = model_C, from = NULL, to = ))
  expect_error(recodeClusterZ(celda.mod = model_C, from = c(1,2,3,4,5), to = c(1,2,3,4,6)))
  expect_error(recodeClusterZ(celda.mod = model_C, from = c(1,2,3,4,6), to = c(1,2,3,4,5)))    
  new.recoded <- recodeClusterZ(celda.mod = model_C, from = c(1,2,3,4,5), to = c(5,4,3,2,1))
  expect_equal(model_C@clusters$z == 1, new.recoded@clusters$z == 5)
})

# compareCountMatrix
test_that(desc = "Testing CompareCountMatrix with celda_C", {
  expect_true(compareCountMatrix(counts = celdaC.sim$counts, celda.mod = model_C))
  
  less.cells <- celdaC.sim$counts[,1:100]
  expect_error(compareCountMatrix(counts = less.cells, celda.mod = model_C),
               "There was a mismatch between the provided count matrix and the count matrix used to generate the provided celda result.")
  
  counts.matrix.error <- matrix(data = 1, nrow = nrow(celdaC.sim$counts), ncol = ncol(celdaC.sim$counts))
  expect_false(compareCountMatrix(counts = counts.matrix.error, celda.mod = model_C, error.on.mismatch = FALSE))
  expect_error(compareCountMatrix(counts = counts.matrix.error, celda.mod = model_C, error.on.mismatch = TRUE))
  
})

# topRank
test_that(desc = "Checking topRank to see if it runs without errors", {

  top.rank <- topRank(matrix = factorized$proportions$module, threshold = NULL)
  expect_equal(names(top.rank),
               c("index","names"))
  top.rank <- topRank(matrix = factorized$proportions$module, n = 1000)
})

# plotHeatmap
test_that(desc = "Testing plotHeatmap with celda_C", {
  expect_error(plotHeatmap(counts = celdaC.sim$counts, z = model_C@params$K), "Length of z must match number of columns in counts matrix")
  expect_error(plotHeatmap(counts = celdaC.sim$counts, z = model_C@clusters$z, scale.row = model_C), "'scale.row' needs to be of class 'function'")
  expect_error(plotHeatmap(counts = celdaC.sim$counts, z = model_C@clusters$z, trim = 3), "'trim' should be a 2 element vector specifying the lower and upper boundaries")
})


# plotHeatmap with annotation.cell
test_that(desc = "Testing plotHeatmap with celda_C, including annotations",{
  annot <- as.data.frame(c(rep(x = 1, times = ncol(celdaC.sim$counts) - 100),rep(x = 2, 100)))
  
  rownames(annot) <- colnames(celdaC.sim$counts)
  expect_equal(names(plotHeatmap(celda.mod = model_C, counts = celdaC.sim$counts, annotation.cell = annot, z = model_C@clusters$z)),
               c("tree_row", "tree_col", "gtable"))
  
  rownames(annot) <- NULL
  expect_equal(names(plotHeatmap(celda.mod = model_C, counts = celdaC.sim$counts, annotation.feature = as.matrix(annot), z = model_C@clusters$z)),
               c("tree_row", "tree_col", "gtable"))
  
  rownames(annot) <- rev(colnames(celdaC.sim$counts))
  expect_error(plotHeatmap(celda.mod = model_C, counts = celdaC.sim$counts, annotation.cell = annot, z = model_C@clusters$z),
               "Row names of 'annotation.cell' are different than the column names of 'counts'")
})

# celdaHeatmap
test_that(desc = "Testing celdaHeatmap with celda_C", {
  expect_equal(names(celdaHeatmap(celda.mod = model_C, counts = celdaC.sim$counts)),
               c("tree_row", "tree_col", "gtable"))
})

# celdaProbabilityMap
test_that(desc = "Testing celdaProbabiltyMap.celda_C for sample level",{
  plot.obj = celdaProbabilityMap(counts=celdaC.sim$counts, celda.mod=model_C, level="sample")
  expect_true(!is.null(plot.obj))
  
  ## Without a sample label
  model_C = celda_C(celdaC.sim$counts, sample.label=NULL, K=celdaC.sim$K, max.iter=5, nchain=1)
  plot.obj = celdaProbabilityMap(counts=celdaC.sim$counts, celda.mod=model_C, level="sample")
  expect_true(!is.null(plot.obj))  
})

# differentialExpression
test_that(desc = "Testing differentialExpression with celda_C", {
  diffexp_K1 <- differentialExpression(counts = celdaC.sim$counts, celda.mod = model_C, c1 = 1)
  expect_equal(class(diffexp_K1), c("data.table", "data.frame"))
  expect_equal(class(diffExp_K1 <- differentialExpression(counts = celdaC.sim$counts, celda.mod = model_C, c1 = 2:3, c2 = 4, log2fc.threshold = 0.5)),
               c("data.table", "data.frame"))  	
  expect_error(differentialExpression(counts = "counts", celda.mod = model_C, c1 = 3, log2fc.threshold = 0.5),"'counts' should be a numeric count matrix")
  expect_error(differentialExpression(counts = celdaC.sim$counts, celda.mod = NULL, c1 = 3), "'celda.mod' should be an object of class celda_C or celda_CG")               
  expect_error(differentialExpression(counts = celdaC.sim$counts, celda.mod = model_C, c1 = NULL, log2fc.threshold = 0.5, only.pos = TRUE))
  
})

test_that(desc = "Testing celdaTsne with celda_C when model class is changed, should error",{
  model_X <- model_C
  class(model_X) <- "celda_X"
  expect_error(celdaTsne(counts=celdaC.sim$counts, celda.mod=model_X, max.cells=length(model_C@clusters$z), min.cluster.size=10),
               "unable to find an inherited method for function 'celdaTsne' for signature '\"celda_X\"'")
})

test_that(desc = "Testing celdaTsne with celda_C including all cells",{
  tsne = celdaTsne(counts=celdaC.sim$counts, celda.mod=model_C, max.cells=length(model_C@clusters$z), min.cluster.size=10)
  plot.obj = plotDimReduceCluster(tsne[,1], tsne[,2], model_C@clusters$z)
  expect_true(ncol(tsne) == 2 & nrow(tsne) == length(model_C@clusters$z))
  expect_true(!is.null(plot.obj))
})

test_that(desc = "Testing celdaTsne with celda_C including a subset of cells",{
  expect_success(expect_error(tsne <- celdaTsne(counts=celdaC.sim$counts, celda.mod=model_C, max.cells=50, min.cluster.size=50)))
  tsne <- celdaTsne(counts=celdaC.sim$counts, celda.mod=model_C, max.cells=100, min.cluster.size=10)
  plot.obj = plotDimReduceCluster(tsne[,1], tsne[,2], model_C@clusters$z)
  expect_true(ncol(tsne) == 2 & nrow(tsne) == length(model_C@clusters$z) && sum(!is.na(tsne[,1])) == 100)
  expect_true(!is.null(plot.obj))
})

# featureModuleLookup
test_that(desc = "Testing featureModuleLookup with celda_C", {
  expect_error(featureModuleLookup(celdaC.sim$counts, model_C, "test_feat"))
})

# cC.splitZ
test_that(desc = "Testing error checking for cC.splitZ", {
  r = simulateCells("celda_C", S=1, C.Range=c(50,100), K=2)
  dc = cC.decomposeCounts(r$counts, r$sample.label, r$z, r$K)
  res = cC.splitZ(r$counts, dc$m.CP.by.S, dc$n.G.by.CP, dc$n.CP, s=as.integer(r$sample.label), z=r$z, K=r$K, nS=dc$nS, nG=dc$nG, alpha=1, beta=1, z.prob=NULL, min.cell=1000)
  expect_true(grepl("Cluster sizes too small", res$message))
})

test_that(desc = "Testing perplexity.celda_C", {
  expect_true(is.numeric(perplexity(celdaC.sim$counts, model_C)))
  
  class(model_C) = c("celda_CG")
  expect_error(perplexity.celda_C(celdaC.sim$counts, model_C),
               "could not find function \"perplexity.celda_C\"")
})
