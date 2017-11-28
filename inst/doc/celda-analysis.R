## ----eval=FALSE----------------------------------------------------------
#  library(devtools)
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
sim_counts = simulateCells("celda_CG", K = 3, L = 10)

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

## ---- fig.show='hold', warning = FALSE, message = FALSE------------------

celdaCG.res = celda_CG(counts = sim_counts$counts, K = 3, L = 10, max.iter = 10, sample.label = sim_counts$sample.label)

## ------------------------------------------------------------------------
z = celdaCG.res$z
y = celdaCG.res$y

table(z, sim_counts$z)
table(y, sim_counts$y)


## ---- fig.width = 7, fig.height = 7, warning = FALSE, message = FALSE----
norm.counts <- normalizeCounts(sim_counts$counts)
renderCeldaHeatmap(counts = norm.counts, z=z, y=y, normalize = NULL, color_scheme = "divergent",cluster_gene = TRUE, cluster_cell = TRUE)

## ------------------------------------------------------------------------
factorized <- factorizeMatrix(celdaCG.res, sim_counts$count)

pop.states = factorized$proportions$population.states
head(pop.states)

## ---- fig.width = 7, fig.height = 7--------------------------------------
data.pca <- prcomp(t(scale(t(factorized$proportions$cell.states))),scale = F, center = F)

plotDrCluster(dim1 = data.pca$rotation[,1], dim2 = data.pca$rotation[,2], cluster = as.factor(celdaCG.res$z))

plotDrState(dim1 = data.pca$rotation[,1], dim2 = data.pca$rotation[,2], matrix = factorized$proportions$cell.states, rescale = TRUE)

## ---- fig.width = 7, fig.height = 7--------------------------------------
renderCeldaHeatmap(pop.states, color_scheme = "sequential", show_cellnames=T, show_genenames=T, breaks = NA, z= 1:ncol(pop.states), y = 1:nrow(pop.states))

## ---- fig.width = 7, fig.height = 7--------------------------------------
rel.states = sweep(pop.states, 1, rowSums(pop.states), "/")
renderCeldaHeatmap(rel.states, color_scheme = "sequential", show_cellnames=T, show_genenames=T, breaks = NA, z= 1:ncol(rel.states), y = 1:nrow(rel.states))

## ------------------------------------------------------------------------
celda.res.list <- celda(counts = sim_counts$counts,K = 2:4, L = 9:11, nchains = 3, max.iter = 10, model = "celda_CG")

## ------------------------------------------------------------------------
model = getModel(celda.res.list, K = 3, L = 10, best = "loglik")

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

## ---- message = FALSE, verbose = FALSE-----------------------------------
pbmc.res = celda(pbmc_select, K = 15, L = 20, cores = 1, model = "celda_CG", nchains = 10, max.iter = 100)
model.res.list<-getModel(pbmc.res, K = 15, L = 20, best = "loglik")

## ------------------------------------------------------------------------
factorize.matrix = factorizeMatrix(model.res.list, counts=pbmc_select)
top.genes <- topRank(factorize.matrix$proportions$gene.states, n = 25)

## ------------------------------------------------------------------------
top.genes$names$L6

## ------------------------------------------------------------------------
top.genes$names$L15

## ---- fig.width = 7, fig.height = 7--------------------------------------
top.genes.ix <- unique(unlist(top.genes$index))
norm.pbmc.counts <- normalizeCounts(pbmc_select)
renderCeldaHeatmap(norm.pbmc.counts[top.genes.ix,], z = model.res.list$z, y = model.res.list$y[top.genes.ix], normalize = NULL, color_scheme = "divergent")

## ---- message = FALSE----------------------------------------------------
norm_pbmc <- normalizeCounts(factorize.matrix$counts$cell.states)
set.seed(123)
tsne <- Rtsne(t(norm_pbmc), pca = FALSE, max_iter = 5000)
pbmc_tsne <- tsne$Y

## ---- fig.width = 7, fig.height = 7--------------------------------------
plotDrCluster(dim1 = pbmc_tsne[,1], dim2 = pbmc_tsne[,2], cluster = as.factor(model.res.list$z))

plotDrState(dim1 = pbmc_tsne[,1], dim2 = pbmc_tsne[,2], matrix = factorize.matrix$proportions$cell.states, rescale = TRUE)

marker.genes <- c("ENSG00000168685_IL7R","ENSG00000198851_CD3E","ENSG00000105374_NKG7","ENSG00000203747_FCGR3A","ENSG00000090382_LYZ","ENSG00000179639_FCER1A","ENSG00000156738_MS4A1", "ENSG00000163736_PPBP")
gene.counts <- pbmc_select[marker.genes,]
plotDrGene(dim1 = pbmc_tsne[,1],dim2 = pbmc_tsne[,2], matrix = gene.counts, rescale = TRUE)

