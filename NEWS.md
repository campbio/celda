# celda v1.14.2 (2023-01-19)
* Update to match Bioconductor release version

# celda v1.13.0 (2022-10-20)
* Bug fixes related to cluster labels stored as factors and plotting
* Updated sparse matrix conversion to work with Matrix v1.4-2

# celda v1.12.0 (2022-04-30)
* Update to match Bioconductor 3.15 release version

# celda v1.11.1 (2022-03-31)
* Fixes to reports
* Use smoothe splines for perplexity and RPC plots

# celda v1.11.0 (2022-03-31)
* Improvments to decontX vignette
* Added ability to subsample to speed up perplexity calculations
* Added ability to use batch parameter with the raw matrix in decontX

# celda v1.10.0 (2021-12-28)
* Update to match Bioconductor release version

# celda v1.9.3 (2021-10-04)
* Fixed bug in checking background matrix with decontX
* Switched to using Github Actions for Continuous Integration
* Fixed plotting bugs in celda results reports
* Speed up final step in decontX when creating final decontaminated matrix

# celda v1.9.2 (2021-07-19)
* Added a `NEWS.md` file to track changes to the package.
* Added new tutorials and documentation generated with pkgdown.
* Removed warnings in plotRPC functions.
* Added use of "displayName" to several functions that show feature names. 
* Minor bug fix when the input matrix was sparse and contained non-integer values.
* Several improvements to plotting functions. 

# celda v1.7.7 (2021-04-12):
* Added handling for sparse matrices

# celda v1.7.6 (2021-04-04):
* Added functions for creating HTML reports
* Fixed bug in decontX plotting

# celda v1.7.4 (2021-03-09):
* Enable input of raw/droplet matrix into decontX to estimate ambient RNA

# celda v1.1.6 (2019-07-16):
* Add multiclass decision tree

# celda v1.1.4 (2019-05-28):
* Add Alternate headings support for plotDimReduceFeature

# celda v1.1.3 (2019-05-14):
* Add multiclass decision tree (MCDT) cell cluster annotation

# celda v1.1.2 (2019-05-14):
* Fix a bug in celdaHeatmap

# celda v1.0.1 (2019-05-09):
* Default seed setting to maintain reproducibility

# celda v0.99.34 (2019-04-23):
* Minor changes to the vignettes

# celda v0.99.23 (2019-04-10):
* Remove pheatmap import

# celda v0.99.22 (2019-04-09):
* Package celda, for bi-clustering of single-cell 'omics data.

# celda v0.99.8 (2019-03-11):
* Second submission to Bioconductor

# celda v0.99.0 (2018-05-15):
* First submission to Bioconductor
