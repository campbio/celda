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

## ----eval=T, warning = FALSE, echo = FALSE---------------------------------
library(devtools)
library(celda)
library(Matrix)
library(gtools)
library(ggplot2)
library(Rtsne)
library(reshape2)

## --------------------------------------------------------------------------
sim_counts <- simulateCells("celda_CG", K = 5, L = 10, S = 10)

## --------------------------------------------------------------------------
dim(sim_counts$counts)

## --------------------------------------------------------------------------
table(sim_counts$z)

## --------------------------------------------------------------------------
table(sim_counts$y)

## --------------------------------------------------------------------------
table(sim_counts$sample.label)

## ---- warning = FALSE, message = FALSE-------------------------------------

celda.res <- celda(counts = sim_counts$counts, model = "celda_CG", K = 5, L = 10, max.iter = 10, cores = 1, nchains = 1)


## --------------------------------------------------------------------------
names(celda.res)

model <- celda.res$res.list[[1]]
z <- model$z
y <- model$y

table(z, sim_counts$z)
table(y, sim_counts$y)


## ---- fig.width = 8, fig.height = 8, warning = FALSE, message = FALSE------
norm.counts <- normalizeCounts(sim_counts$counts, scale.factor = 1e6)
renderCeldaHeatmap(counts = norm.counts, z=z, y=y, normalize = NULL, 
                   color_scheme = "divergent",cluster_gene = TRUE, 
                   cluster_cell = TRUE)

## --------------------------------------------------------------------------
factorized <- factorizeMatrix(counts = sim_counts$count, celda.mod = model)
names(factorized)

## --------------------------------------------------------------------------
dim(factorized$proportions$gene.states)
head(factorized$proportions$gene.states)

## --------------------------------------------------------------------------
dim(factorized$proportions$cell.states)
factorized$proportions$cell.states[,1:8]

## --------------------------------------------------------------------------
pop.states <- factorized$proportions$population.states
dim(pop.states)
pop.states

## ---- fig.width = 8, fig.height = 8----------------------------------------
data.pca <- prcomp(t(scale(t(factorized$proportions$cell.states))),
                   scale = F, center = F)

plotDrCluster(dim1 = data.pca$rotation[,1], dim2 = data.pca$rotation[,2], 
              cluster = celda.res$res.list[[1]]$z)

plotDrState(dim1 = data.pca$rotation[,1], dim2 = data.pca$rotation[,2],
            matrix = factorized$proportions$cell.states, rescale = TRUE)

## ---- fig.width = 8, fig.height = 8----------------------------------------
absoluteProbabilityHeatmap(counts = sim_counts$counts, celda.mod = model)

## ---- fig.width = 8, fig.height = 8----------------------------------------
relativeProbabilityHeatmap(counts = sim_counts$counts, celda.mod = model)

## ----setup, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,fig.align = "center")
library(celda)
library(Matrix)
library(data.table)
library(pheatmap)
library(Rtsne)

## --------------------------------------------------------------------------
##Original pbmc_data: 7090 genes, 2700 cells
#pbmc_select <- pbmc_data[rowSums(pbmc_data>3) > 3,]

dim(pbmc_select)
head(rownames(pbmc_select))
head(colnames(pbmc_select))

## ---- eval = FALSE---------------------------------------------------------
#  pbmc_res1 <- celda(pbmc_select,
#                     K = seq(5, 50, by = 5),
#                     L = seq(10, 50, by = 5),
#                     cores = 1,
#                     model = "celda_CG",
#                     nchains = 4,
#                     max.iter = 100)

## ---- eval=FALSE-----------------------------------------------------------
#  calc.perplexity <- calculatePerplexityWithResampling(pbmc_res1,
#                                                       counts = pbmc_select,
#                                                       resample = 1)

## ---- fig.width = 8, fig.height = 8----------------------------------------
calc.perplexity$plot

## ---- eval = FALSE---------------------------------------------------------
#  pbmc_res <- celda(pbmc_select, K = 10:20, L = seq(10, 50, by = 5),
#                    cores = 1, model = "celda_CG", nchains = 4, max.iter = 100)

## --------------------------------------------------------------------------
model.pbmc <- getBestModel(pbmc_res, K = 13, L = 50)

## ---- fig.width = 8, fig.height = 8----------------------------------------
relativeProbabilityHeatmap(counts = pbmc_select, celda.mod = model.pbmc)

## ---- message = FALSE------------------------------------------------------
factorize.matrix <- factorizeMatrix(counts=pbmc_select, celda.mod = model.pbmc)
norm_pbmc <- normalizeCounts(factorize.matrix$counts$cell.states)
set.seed(123)
pbmc_tsne <- celdaTsne(counts = pbmc_select, celda.mod = model.pbmc, distance = "hellinger")

## ---- fig.width = 8, fig.height = 8----------------------------------------
plotDrCluster(dim1 = pbmc_tsne[,1], dim2 = pbmc_tsne[,2], cluster = as.factor(model.pbmc$z))

plotDrState(dim1 = pbmc_tsne[,1], dim2 = pbmc_tsne[,2], matrix = factorize.matrix$proportions$cell.states, rescale = TRUE)

marker.genes <- c("ENSG00000168685_IL7R","ENSG00000132646_PCNA",
                  "ENSG00000105374_NKG7","ENSG00000203747_FCGR3A",
                  "ENSG00000090382_LYZ","ENSG00000153563_CD8A",
                  "ENSG00000156738_MS4A1", "ENSG00000163736_PPBP",
                  "ENSG00000101439_CST3")
gene.counts <- pbmc_select[marker.genes,]
plotDrGene(dim1 = pbmc_tsne[,1],dim2 = pbmc_tsne[,2], matrix = gene.counts, 
           rescale = TRUE)

## ----message=FALSE---------------------------------------------------------
diff.exp.clust1 <- diffExpBetweenCellStates(counts = pbmc_select, 
                                            celda.mod = model.pbmc, 
                                            c1 = 1, c2 = NULL)

head(diff.exp.clust1,10)

## ---- message=FALSE--------------------------------------------------------
diff.exp.clust1vs2 <- diffExpBetweenCellStates(counts = pbmc_select, celda.mod = model.pbmc, c1 = 2, c2 = 10)

diff.exp.clust1vs2 <- diff.exp.clust1vs2[diff.exp.clust1vs2$fdr < 0.25,]

## --------------------------------------------------------------------------
head(diff.exp.clust1vs2[order(diff.exp.clust1vs2$log2fc, decreasing = TRUE),],10)

## --------------------------------------------------------------------------
head(diff.exp.clust1vs2[order(diff.exp.clust1vs2$log2fc),],10)

## --------------------------------------------------------------------------
factorize.matrix <- factorizeMatrix(model.pbmc, counts=pbmc_select)
top.genes <- topRank(factorize.matrix$proportions$gene.states, n = 25)

## --------------------------------------------------------------------------
top.genes$names$L42

## --------------------------------------------------------------------------
top.genes$names$L4

## ---- fig.width = 8, fig.height = 8----------------------------------------
top.genes.ix <- unique(unlist(top.genes$index))
norm.pbmc.counts <- normalizeCounts(pbmc_select)
renderCeldaHeatmap(norm.pbmc.counts[top.genes.ix,], z = model.pbmc$z, 
                   y = model.pbmc$y[top.genes.ix], normalize = NULL, 
                   color_scheme = "divergent")

## ---- fig.width = 8, fig.height = 8----------------------------------------
gini <- GiniPlot(counts = pbmc_select, celda.mod = model.pbmc)

gini

## ---- fig.width = 8, fig.height = 8----------------------------------------
filtered.tr.states <- gini$data$Transcriptional_States[gini$data$Gini_Coefficient > 0.3]

top.genes.filtered.ix <- unique(unlist(top.genes$index[as.numeric(levels(filtered.tr.states))][as.numeric(filtered.tr.states)]))

renderCeldaHeatmap(norm.pbmc.counts[top.genes.filtered.ix,], z = model.pbmc$z, 
                   y = model.pbmc$y[top.genes.filtered.ix], 
                   normalize = NULL, color_scheme = "divergent")

## --------------------------------------------------------------------------
lookupGeneModule(counts = pbmc_select, model = model.pbmc, 
                          gene = c("ENSG00000203747_FCGR3A",
                                   "ENSG00000163736_PPBP"))

## ---- fig.width = 8, fig.height = 8, message=FALSE-------------------------
stateHeatmap(counts = pbmc_select, celda.mod = model.pbmc, state.use = 4)

## ----table2, echo = FALSE, message = FALSE, results='asis'-----------------
tabl <- "
| Celda Cluster | Type              | Marker     |
|---------------|:-----------------:|-----------:|
|      1        | Megakaryocytes    | PPBP       |
|      2        | CD8 T cells       | CD8A       |
|      3        | CD4 T cells       | IL7R       |
|      4        | CD4 T cells       | IL7R       |
|      5        | B cells           | MS4A1      |
|      6        | Mitochondrial     |            |
|      7        | CD4 T cells       | IL7R       |
|      8        |                   | PCNA/TUBB  |
|      9        | Dendritic cells   | CST3       |
|      10       | NK cells          | NKG7/GNLY  |
|      11       | FCGR3A+ monocytes | FCGR3A     |
|      12       | CD14+ monocytes   | LYZ/CD14   |
|      13       | CD14+ monocytes   | LYZ/CD14   |
"
cat(tabl)

