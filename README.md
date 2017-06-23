[![Build Status](https://travis-ci.org/compbiomed/celda.svg?branch=master)](https://travis-ci.org/compbiomed/celda)

# celda: CEllular Latent Dirichlet Allocation

"celda" stands for "**CE**llular **L**atent **D**irichlet **A**llocation", which is a suite of Bayesian hierarchical models and supporting functions to perform gene and cell clustering for count data generated by single cell RNA-seq platforms. This algorithm is an extension of the Latent Dirichlet Allocation (LDA) topic modeling framework that has been popular in text mining applications. Celda has advantages over other clustering frameworks:

1. Celda can simultaneously cluster genes into transcriptional states and cells into subpopulations
2. Celda uses count-based Dirichlet-multinomial distributions so no additional normalization is required for 3' DGE single cell RNA-seq
3. These types of models have shown good performance with sparse data.


## Installation Instructions

To install the beta release of celda via devtools:
```
library(devtools)
install_github("compbiomed/celda@v0.1")
```
The most up-to-date (but potentially less stable) version of celda can similarly be installed with:
```
install_github("compbiomed/celda")
```

## Examples and vignettes

Vignettes are available in the package. 

An analysis example using celda with RNASeq via vignette('celda-analysis')



## New Features and announcements
The v0.1 release of celda represents a useable implementation of the various celda clustering models.
Please submit any usability issues or bugs to the issue tracker at https://github.com/compbiomed/celda

You can discuss celda, or ask the developers usage questions, in our [Google Group.](https://groups.google.com/forum/#!forum/celda-list)
