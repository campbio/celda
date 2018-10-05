## ----eval=FALSE------------------------------------------------------------
#  library(devtools)
#  source("https://bioconductor.org/biocLite.R")
#  install.packages("Rcpp")
#  biocLite("SummarizedExperiment")
#  biocLite("MAST")
#  install_github("compbiomed/celda")

## ----eval=FALSE------------------------------------------------------------
#  library(celda)

## ----eval=FALSE------------------------------------------------------------
#  help(package=celda)

## ----eval=T, warning = FALSE, echo = FALSE, message = FALSE----------------
library(devtools)
library(celda)
library(Matrix)
library(gtools)
library(ggplot2)
library(Rtsne)
library(reshape2)

## --------------------------------------------------------------------------
sim_counts <- simulateCells("celda_CG", K = 5, L = 10, S = 5, G = 500)

## --------------------------------------------------------------------------
dim(sim_counts$counts)

## --------------------------------------------------------------------------
table(sim_counts$z)

## --------------------------------------------------------------------------
table(sim_counts$y)

## --------------------------------------------------------------------------
table(sim_counts$sample.label)

## ---- warning = FALSE, message = FALSE-------------------------------------
celda.model <- celda_CG(counts = sim_counts$counts, K = 5, L = 10, verbose = FALSE)

## --------------------------------------------------------------------------
table(celda.model$z, sim_counts$z)
table(celda.model$y, sim_counts$y)


## ---- fig.width = 8, fig.height = 8, warning = FALSE, message = FALSE------
celdaHeatmap(counts = sim_counts$counts, celda.mod = celda.model)

## --------------------------------------------------------------------------
factorized <- factorizeMatrix(counts = sim_counts$counts, celda.mod = celda.model)
names(factorized)

## --------------------------------------------------------------------------
dim(factorized$proportions$module)
head(factorized$proportions$module)

## --------------------------------------------------------------------------
dim(factorized$proportions$cell)
factorized$proportions$cell[,1:3]

## --------------------------------------------------------------------------
cell.pop <- factorized$proportions$cell.population
dim(cell.pop)
cell.pop

## --------------------------------------------------------------------------
top.genes <- topRank(matrix = factorized$proportions$module, n = 50, threshold = NULL)

## --------------------------------------------------------------------------
top.genes$names$L1

## ---- fig.width = 8, fig.height = 8----------------------------------------
top.genes.ix <- unique(unlist(top.genes$index))
celdaHeatmap(counts = sim_counts$counts, celda.mod = celda.model, nfeatures = 10)

## ---- fig.width = 8, fig.height = 8----------------------------------------
tsne <- celdaTsne(counts = sim_counts$counts, celda.mod = celda.model)

plotDimReduceCluster(dim1 = tsne[,1], 
                     dim2 = tsne[,2], 
                     cluster = celda.model$z)

plotDimReduceModule(dim1 = tsne[,1], dim2 = tsne[,2], 
                   celda.mod = celda.model, counts = sim_counts$counts, rescale = TRUE)

plotDimReduceFeature(dim1 = tsne[,1], dim2 = tsne[,2],
                     counts = sim_counts$counts, features = "Gene_1")

## ---- fig.width = 8, fig.height = 8----------------------------------------
celdaProbabilityMap(counts = sim_counts$counts, celda.mod = celda.model)

## ---- fig.width = 8, fig.height = 8----------------------------------------
moduleHeatmap(counts = sim_counts$counts, celda.mod = celda.model, feature.module = 1:2)

## ---- fig.width = 8, fig.height = 8----------------------------------------
genes = c("Gene_1","Gene_2","Gene_3")
gene.ix = which(rownames(sim_counts$counts) %in% genes)

norm.counts <- normalizeCounts(counts = sim_counts$counts, scale.fun = scale)
plotHeatmap(counts = norm.counts, z = celda.model$z, y = celda.model$y, feature.ix = gene.ix, show.names.feature = TRUE)

## ----message=FALSE---------------------------------------------------------
diff.exp.clust1 <- differentialExpression(counts = sim_counts$counts, 
                                            celda.mod = celda.model, 
                                            c1 = 1, c2 = NULL)

head(diff.exp.clust1,10)

## ---- message=FALSE--------------------------------------------------------
diff.exp.clust1vs2 <- differentialExpression(counts = sim_counts$counts, celda.mod = celda.model, c1 = 1, c2 = 2)

diff.exp.clust1vs2 <- diff.exp.clust1vs2[diff.exp.clust1vs2$fdr < 0.25,]

## --------------------------------------------------------------------------
upreg.genes <- head(diff.exp.clust1vs2[order(diff.exp.clust1vs2$log2fc, decreasing = TRUE),],10)

## ---- fig.height = 8, fig.width = 8----------------------------------------
upreg.gene.ix = which(rownames(sim_counts$counts) %in% upreg.genes$Gene)

norm.counts <- normalizeCounts(counts = sim_counts$counts, scale.fun = scale)
plotHeatmap(counts = norm.counts, z = celda.model$z, y = celda.model$y, feature.ix = upreg.gene.ix, show.names.feature = TRUE)

## ---- fig.height = 8, fig.width = 8----------------------------------------
plotDimReduceFeature(dim1 = tsne[,1], dim2 = tsne[,2], counts = sim_counts$counts, features = upreg.genes$Gene[1:9])

## ---- message = FALSE------------------------------------------------------
cgs <- celdaGridSearch(sim_counts$counts,
                            params.test = list(K = 3:7, L = 8:12),
                            cores = 1,
                            model = "celda_CG",
                            nchains = 2,
                            max.iter = 100,
                            best.only = TRUE)

## --------------------------------------------------------------------------
head(cgs$run.params)

## --------------------------------------------------------------------------
cgs <- resamplePerplexity(counts = sim_counts$counts,
                                     celda.list = cgs,
                                     resample = 5)

## ---- fig.width = 8, fig.height = 8, warning = FALSE, message = FALSE------
plotGridSearchPerplexity(celda.list = cgs)

## --------------------------------------------------------------------------
celda.model = subsetCeldaList(celda.list = cgs, params = list(K = 5, L = 10))

## ---- message=FALSE, eval=FALSE--------------------------------------------
#  cgs <- celdaGridSearch(sim_counts$counts,
#                              params.test = list(K = 3:7, L = 8:12),
#                              cores = 1,
#                              model = "celda_CG",
#                              nchains = 2,
#                              max.iter = 100,
#                              best.only = FALSE)
#  
#  cgs <- resamplePerplexity(counts = sim_counts$counts,
#                                       celda.list = cgs,
#                                       resample = 2)
#  
#  celda.model.list = subsetCeldaList(celda.list = cgs, params = list(K = 5, L = 10))
#  
#  celda.model1 <- selectBestModel(celda.list = celda.model.list)
#  
#  celda.model2 <- subsetCeldaList(celda.list = cgs, params = list(K = 5, L = 10, chain = 2))

## --------------------------------------------------------------------------
featureModuleLookup(counts = sim_counts$counts, celda.mod = celda.model, 
                          feature = c("Gene_99"))

## --------------------------------------------------------------------------
celda.model.y.recoded <- recodeClusterY(celda.mod = celda.model, from = c(1,2,3,4,5), to = c(2,1,3,4,5))

celda.model.z.recoded <- recodeClusterZ(celda.mod = celda.model, from = c(1,2,3,4,5), to = c(2,1,3,4,5))

## --------------------------------------------------------------------------
table(celda.model$z)
table(celda.model$y)

## --------------------------------------------------------------------------
table(celda.model.z.recoded$z)
table(celda.model.y.recoded$y)

