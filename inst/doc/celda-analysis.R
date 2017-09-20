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

## ------------------------------------------------------------------------
sim_counts = simulateCells("celda_CG", K=3,L=10)

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

celdaCG.res = celda_CG(counts = sim_counts$counts, K =3, L = 10, max.iter = 10, sample.label = sim_counts$sample.label)

## ------------------------------------------------------------------------
z = celdaCG.res$z
y = celdaCG.res$y

table(z, sim_counts$z)
table(y, sim_counts$y)


## ---- fig.width = 7, fig.height = 7, warning = FALSE, message = FALSE----
norm.counts <- normalizeCounts(sim_counts$counts)
render_celda_heatmap(counts=norm.counts,z=z,y=y, cluster.column = FALSE, cluster.row = FALSE)

## ------------------------------------------------------------------------
factorized <- factorizeMatrix.celda_CG(counts = sim_counts$counts, celda.obj = celdaCG.res)

pop.states = t(factorized$proportions$population.states)
head(pop.states)

## ---- fig.width = 7, fig.height = 7--------------------------------------
render_celda_heatmap(pop.states, col=colorRampPalette(c("white", "blue"))(70), scale.row=F, show_cellnames=T, show_genenames=T, cluster.column=F, cluster.row=F, breaks = NA, z= 1:3, y = 1:10)

## ---- fig.width = 7, fig.height = 7--------------------------------------
rel.states = sweep(pop.states, 1, rowSums(pop.states), "/")
render_celda_heatmap(rel.states, col=colorRampPalette(c("white","blue"))(70), scale.row=F, show_cellnames=T, show_genenames=T, cluster.column=F, cluster.row=F, breaks = NA, z= 1:3, y = 1:10)

## ------------------------------------------------------------------------
celda.res.list <- celda(counts = sim_counts$counts,K = 2:4, L = 9:11, nchains = 3, max.iter = 10, model = "celda_CG")

## ---- fig.width = 7, fig.height = 7, message = FALSE---------------------
visualize.model = visualize_model_performance(celda.res.list,method="perplexity")

visualize.model$K

## ------------------------------------------------------------------------
model = getModel(celda.res.list, K = 3, L = 10, best = "loglik")

## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,fig.align = "center")
library(celda)
library(Matrix)
library(data.table)
library(pheatmap)
library(biomaRt)

## ------------------------------------------------------------------------
load(file="../data/cmatp_data.rda")
dim(cmatp_data)
head(rownames(cmatp_data))
head(colnames(cmatp_data))

## ------------------------------------------------------------------------
cmatp_select <- cmatp_data[rowSums(cmatp_data>4) > 4,]

## ---- message = FALSE, verbose = FALSE-----------------------------------
cmatp.res = celda_CG(cmatp_select, sample.label = rep(1,ncol(cmatp_select)), K = 15, L = 20)

## ------------------------------------------------------------------------
factorize.matrix = factorizeMatrix.celda_CG(counts = cmatp_select, celda.obj = cmatp.res)
top.genes <- topRank(factorize.matrix$proportions$gene.states)

## ------------------------------------------------------------------------
top.genes$names$L16

## ------------------------------------------------------------------------
top.genes$names$L19

## ---- fig.width = 7, fig.height = 7--------------------------------------
top.genes.ix <- unique(unlist(top.genes$index))
norm.cmatp.counts <- normalizeCounts(cmatp_select)
render_celda_heatmap(norm.cmatp.counts[top.genes.ix,], z = cmatp.res$z, y = cmatp.res$y[top.genes.ix], cluster.row = FALSE, cluster.column = FALSE)

## ---- fig.width = 7, fig.height = 7--------------------------------------
cmatp.states = t(factorize.matrix$proportions$population.states)
cmatp.rel.states = sweep(cmatp.states, 1, rowSums(cmatp.states), "/")
render_celda_heatmap(cmatp.rel.states, col=colorRampPalette(c("white","blue","darkgreen","green"))(100), scale.row=F, show_cellnames=T, show_genenames=TRUE, cluster.column=FALSE, cluster.row=FALSE, breaks = NA, z= 1:15, y = 1:20)

