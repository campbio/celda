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

