## ----eval=FALSE------------------------------------------------------------
#  library(celda)

## ----eval=FALSE------------------------------------------------------------
#  help(package=celda)

## ----eval=T, warning = FALSE, echo = FALSE, message = FALSE----------------
library(celda)

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


## --------------------------------------------------------------------------
factorized <- factorizeMatrix(counts = sim_counts$counts, celda.mod = celda.model)
names(factorized)

## --------------------------------------------------------------------------
dim(factorized$proportions$cell)
head(factorized$proportions$cell[,1:3], 5)

## --------------------------------------------------------------------------
cell.pop <- factorized$proportions$cell.population
dim(cell.pop)
head(cell.pop, 5)

## --------------------------------------------------------------------------
dim(factorized$proportions$module)
head(factorized$proportions$module, 5)

## --------------------------------------------------------------------------
top.genes <- topRank(matrix = factorized$proportions$module, n = 10, threshold = NULL)

## --------------------------------------------------------------------------
top.genes$names$L1

## ---- eval = FALSE, fig.width = 8, fig.height = 8--------------------------
#  celdaHeatmap(counts = sim_counts$counts, celda.mod = celda.model, nfeatures = 10)

## --------------------------------------------------------------------------
tsne <- celdaTsne(counts = sim_counts$counts, celda.mod = celda.model)

## ---- eval = FALSE---------------------------------------------------------
#  plotDimReduceCluster(dim1 = tsne[,1],
#                       dim2 = tsne[,2],
#                       cluster = celda.model$z)
#  
#  plotDimReduceModule(dim1 = tsne[,1], dim2 = tsne[,2],
#                     celda.mod = celda.model, counts = sim_counts$counts, rescale = TRUE)
#  
#  plotDimReduceFeature(dim1 = tsne[,1], dim2 = tsne[,2],
#                       counts = sim_counts$counts, features = "Gene_1")

## ---- eval = FALSE, fig.width = 8, fig.height = 8--------------------------
#  celdaProbabilityMap(counts = sim_counts$counts, celda.mod = celda.model)

## ---- eval = FALSE, fig.width = 8, fig.height = 8--------------------------
#  moduleHeatmap(counts = sim_counts$counts, celda.mod = celda.model, feature.module = 1, top.cells = 100)

## --------------------------------------------------------------------------
genes = c(top.genes$names$L1, top.genes$names$L10)
gene.ix = which(rownames(sim_counts$counts) %in% genes)
norm.counts <- normalizeCounts(counts = sim_counts$counts, scale.fun = scale)

## ---- eval = FALSE, fig.width = 8, fig.height = 8--------------------------
#  plotHeatmap(counts = norm.counts, z = celda.model$z, y = celda.model$y, feature.ix = gene.ix, show.names.feature = TRUE)

## ----message=FALSE---------------------------------------------------------
diffexp.clust1 <- differentialExpression(counts = sim_counts$counts, 
                                            celda.mod = celda.model, 
                                            c1 = 1, c2 = NULL)

head(diffexp.clust1,5)

## ---- message=FALSE--------------------------------------------------------
diffexp.clust1vs2 <- differentialExpression(counts = sim_counts$counts, celda.mod = celda.model, c1 = 1, c2 = 2)

diffexp.clust1vs2 <- diffexp.clust1vs2[diffexp.clust1vs2$FDR < 0.05 & abs(diffexp.clust1vs2$Log2_FC) > 2,]

head(diffexp.clust1vs2, 5)

## --------------------------------------------------------------------------
diffexp.gene.ix = which(rownames(sim_counts$counts) %in% diffexp.clust1vs2$Gene)

norm.counts <- normalizeCounts(counts = sim_counts$counts, scale.fun = scale)

## ---- eval = FALSE---------------------------------------------------------
#  plotHeatmap(counts = norm.counts[,celda.model$z %in% c(1,2)], z = celda.model$z[celda.model$z %in% c(1,2)], y = celda.model$y, feature.ix = diffexp.gene.ix, show.names.feature = TRUE)

## ---- eval = FALSE, message = FALSE----------------------------------------
#  cgs <- celdaGridSearch(sim_counts$counts,
#                              params.test = list(K = 3:7, L = 8:12),
#                              cores = 1,
#                              model = "celda_CG",
#                              nchains = 2,
#                              max.iter = 100,
#                              best.only = TRUE)

## ---- eval = FALSE---------------------------------------------------------
#  cgs <- resamplePerplexity(counts = sim_counts$counts,
#                                       celda.list = cgs,
#                                       resample = 5)

## ---- eval = FALSE, fig.width = 8, fig.height = 8, warning = FALSE, message = FALSE----
#  plotGridSearchPerplexity(celda.list = cgs)

## ---- eval = FALSE---------------------------------------------------------
#  celda.model = subsetCeldaList(celda.list = cgs, params = list(K = 5, L = 10))

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
#  cgs.K5.L10 = subsetCeldaList(celda.list = cgs, params = list(K = 5, L = 10))
#  
#  celda.model1 <- selectBestModel(celda.list = cgs.K5.L10)

## --------------------------------------------------------------------------
featureModuleLookup(counts = sim_counts$counts, celda.mod = celda.model, 
                          feature = c("Gene_99"))

## --------------------------------------------------------------------------
celda.model.z.recoded <- recodeClusterZ(celda.mod = celda.model, from = c(1,2,3,4,5), to = c(2,1,3,4,5))

## --------------------------------------------------------------------------
table(celda.model$z, celda.model.z.recoded$z)

