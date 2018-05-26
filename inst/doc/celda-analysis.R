## ----eval=FALSE----------------------------------------------------------
#  library(devtools)
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("SummarizedExperiment")
#  biocLite("MAST")
#  install_github("compbiomed/celda")

## ----eval=FALSE----------------------------------------------------------
#  library(celda)

## ----eval=FALSE----------------------------------------------------------
#  help(package=celda)

## ----eval=T, warning = FALSE, echo = FALSE-------------------------------
library(devtools)
library(celda)
library(Matrix)
library(gtools)
library(ggplot2)
library(Rtsne)
library(reshape2)

## ------------------------------------------------------------------------
sim_counts <- simulateCells("celda_CG", K = 6, L = 6, S = 10)

## ------------------------------------------------------------------------
dim(sim_counts$counts)

## ------------------------------------------------------------------------
table(sim_counts$z)

## ------------------------------------------------------------------------
table(sim_counts$y)

## ------------------------------------------------------------------------
table(sim_counts$sample.label)

## ------------------------------------------------------------------------
table(sim_counts$z, sim_counts$sample.label)

## ---- warning = FALSE, message = FALSE-----------------------------------

celda.res <- celda(counts = sim_counts$counts, model = "celda_CG", K = 6, L = 6, max.iter = 10, cores = 1, nchains = 1)


## ------------------------------------------------------------------------
names(celda.res)

model <- celda.res$res.list[[1]]
z <- model$z
y <- model$y

table(z, sim_counts$z)
table(y, sim_counts$y)


## ---- fig.width = 7, fig.height = 7, warning = FALSE, message = FALSE----
norm.counts <- normalizeCounts(sim_counts$counts, scale.factor = 1e6)
renderCeldaHeatmap(counts = norm.counts, z=z, y=y, normalize = NULL, color_scheme = "divergent",cluster_gene = TRUE, cluster_cell = TRUE)

## ------------------------------------------------------------------------
factorized <- factorizeMatrix(model, sim_counts$count)

## ------------------------------------------------------------------------
dim(factorized$proportions$gene.states)
head(factorized$proportions$gene.states)

## ------------------------------------------------------------------------
dim(factorized$proportions$cell.states)
factorized$proportions$cell.states[1:6,1:8]

## ------------------------------------------------------------------------
pop.states <- factorized$proportions$population.states
dim(pop.states)

## ------------------------------------------------------------------------
pop.states

## ---- fig.width = 7, fig.height = 7--------------------------------------
data.pca <- prcomp(t(scale(t(factorized$proportions$cell.states))),scale = F, center = F)

plotDrCluster(dim1 = data.pca$rotation[,1], dim2 = data.pca$rotation[,2], cluster = celda.res$res.list[[1]]$z)

plotDrState(dim1 = data.pca$rotation[,1], dim2 = data.pca$rotation[,2], matrix = factorized$proportions$cell.states, rescale = TRUE)

## ---- fig.width = 7, fig.height = 7--------------------------------------
renderCeldaHeatmap(pop.states, color_scheme = "sequential", show_cellnames=T, show_genenames=T, breaks = NA, z= 1:ncol(pop.states), y = 1:nrow(pop.states))

## ---- fig.width = 7, fig.height = 7--------------------------------------
rel.states <- sweep(pop.states, 1, rowSums(pop.states), "/")
renderCeldaHeatmap(rel.states, color_scheme = "sequential", show_cellnames=T, show_genenames=T, breaks = NA, z= 1:ncol(rel.states), y = 1:nrow(rel.states))

## ---- eval = TRUE--------------------------------------------------------
celda.res.list <- celda(counts = sim_counts$counts,K = 5:7, L = 6:7, nchains = 2, max.iter = 10, model = "celda_CG", cores = 1)

## ---- fig.width = 7, fig.height = 7--------------------------------------
calc.perplexity <- calculatePerplexityWithResampling(celda.res.list, counts = sim_counts$counts, resample = TRUE)

visualizePerplexityByKL(calc.perplexity$perplexity.info)

## ---- fig.width = 7, fig.height = 7--------------------------------------
cluster.plot.sim <- gettingClusters(celda.list = celda.res.list, matrix = sim_counts$counts, iterations = 3)
cluster.plot.sim

## ---- eval = FALSE-------------------------------------------------------
#  model <- getBestModel(celda.res.list, K = 6, L = 8)

## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,fig.align = "center")
library(celda)
library(Matrix)
library(data.table)
library(pheatmap)
library(Rtsne)

## ------------------------------------------------------------------------
dim(pbmc_data)
head(rownames(pbmc_data))
head(colnames(pbmc_data))

## ------------------------------------------------------------------------
pbmc_select <- pbmc_data[rowSums(pbmc_data>4) > 4,]

## ---- eval = FALSE-------------------------------------------------------
#  pbmc_res <- celda(pbmc_select, K = seq(10:30,by = 2), L = c(20,30), cores = 1, model = "celda_CG", nchains = 4, max.iter = 100)

## ---- eval = FALSE-------------------------------------------------------
#  cluster.plot <- gettingClusters(celda.list = pbmc_res, matrix = pbmc_select, iterations = 3)

## ---- fig.width = 7, fig.height = 7--------------------------------------
cluster.plot

## ------------------------------------------------------------------------
model.pbmc <- getBestModel(pbmc_res, K = 16, L = 30)

## ---- fig.width = 7, fig.height = 7, message=FALSE-----------------------
stateHeatmap(counts = pbmc_select, celda.mod = model.pbmc, state.use = 3)

## ---- message = FALSE----------------------------------------------------
factorize.matrix <- factorizeMatrix(model.pbmc, counts=pbmc_select)
norm_pbmc <- normalizeCounts(factorize.matrix$counts$cell.states)
set.seed(123)
pbmc_tsne <- celdaTsne(counts = pbmc_select, celda.mod = model.pbmc, distance = "cosine")

## ---- fig.width = 7, fig.height = 7--------------------------------------
plotDrCluster(dim1 = pbmc_tsne[,1], dim2 = pbmc_tsne[,2], cluster = as.factor(model.pbmc$z))

plotDrState(dim1 = pbmc_tsne[,1], dim2 = pbmc_tsne[,2], matrix = factorize.matrix$proportions$cell.states, rescale = TRUE)

marker.genes <- c("ENSG00000168685_IL7R","ENSG00000198851_CD3E","ENSG00000105374_NKG7",
                 "ENSG00000203747_FCGR3A","ENSG00000090382_LYZ","ENSG00000179639_FCER1A",
                 "ENSG00000156738_MS4A1", "ENSG00000163736_PPBP")
gene.counts <- pbmc_select[marker.genes,]
plotDrGene(dim1 = pbmc_tsne[,1],dim2 = pbmc_tsne[,2], matrix = gene.counts, rescale = TRUE)

## ----message=FALSE-------------------------------------------------------
diff.exp.clust1 <- diffExpBetweenCellStates(counts = pbmc_select, celda.mod = model.pbmc, c1 = 1, c2 = NULL)

diff.exp.clust1

## ---- message=FALSE------------------------------------------------------
diff.exp.clust1vs2 <- diffExpBetweenCellStates(counts = pbmc_select, celda.mod = model.pbmc, c1 = 1, c2 = 2)

diff.exp.clust1vs2

## ------------------------------------------------------------------------
factorize.matrix <- factorizeMatrix(model.pbmc, counts=pbmc_select)
top.genes <- topRank(factorize.matrix$proportions$gene.states, n = 25)

## ------------------------------------------------------------------------
top.genes$names$L3

## ------------------------------------------------------------------------
top.genes$names$L29

## ---- fig.width = 7, fig.height = 7--------------------------------------
top.genes.ix <- unique(unlist(top.genes$index))
norm.pbmc.counts <- normalizeCounts(pbmc_select)
renderCeldaHeatmap(norm.pbmc.counts[top.genes.ix,], z = model.pbmc$z, y = model.pbmc$y[top.genes.ix], normalize = NULL, color_scheme = "divergent")

## ---- fig.width = 7, fig.height = 7--------------------------------------
gini <- GiniPlot(counts = pbmc_select, celda.mod = model.pbmc)

gini

## ---- fig.width = 7, fig.height = 7--------------------------------------
filtered.tr.states <- gini$data$Transcriptional_States[gini$data$Gini_Coefficient > 0.19]

top.genes.filtered.ix <- unique(unlist(top.genes$index[as.numeric(levels(filtered.tr.states))][as.numeric(filtered.tr.states)]))

renderCeldaHeatmap(norm.pbmc.counts[top.genes.filtered.ix,], z = model.pbmc$z, y = model.pbmc$y[top.genes.filtered.ix], normalize = NULL, color_scheme = "divergent")

