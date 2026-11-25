# celda (development version)
* New `findMetaCells()` groups the cells of each sample into meta-cells:
  small sets of similar cells found with `recursiveSplitCell()` and a
  kneedle choice of K, with meta-cells smaller than `minCellMeta` merged into
  their nearest neighbour. Results are returned in the input's cell order;
  a sample with fewer than `minCellMeta` cells becomes one meta-cell.
  Based on work by Nageen Zahra (#429).
* `celda_CG()` gains `zInitialize = "fixed"`, which keeps the cell clusters
  at `zInit` and fits only the feature modules.
* `recursiveSplitModule()` gains `useMetaCells`, which keeps `zInit` fixed
  while modules are split, so meta-cells from `findMetaCells()` can be used
  directly. `zInit` can now hold character or factor labels.
* `plotRPC()` gains `plotKnees` and `knees` to mark knee points of the mean
  perplexity curve. For `celda_CG` results they are placed on whichever of K
  or L varies.
* `retrieveFeatureIndex()` no longer matches a feature that is `NA` to a row
  named `NA`, or, with `exactMatch = FALSE`, warns that it matched multiple
  items; such features are reported as not present. A `by` that is not a
  single string now stops with a clear error (#428).
* `plotDimReduceFeature()` no longer puts labels from `headers` or
  `displayName` on the wrong panels when one of the `features` is not
  found (#433).
* `plotDimReduceFeature()` now labels the panels with `headers` when both
  `headers` and `displayName` are given, instead of silently ignoring
  `headers`. `displayName` is still used to look up `features`, and
  `headers` of the wrong length now stop with an error in this case too.
  As elsewhere, `headers` are dropped with a warning when
  `exactMatch = FALSE`. The error and warning about `headers` are reworded
  to be clearer (#434).
* `findMarkersTree()` no longer fails with "subscript out of bounds" when
  more than one cluster or metacluster has no positive marker, as can happen
  at a low `threshold` once the top-level tree has used up the features. An
  alternative split is now tried for each such class, and only added if it
  picks out that class. With metaclusters, the threshold is no longer
  lowered once no unused feature could become a marker that way:
  metaclusters still without a positive marker are merged and separated in
  their own subtree, with a warning, instead of the threshold being lowered
  indefinitely. Merged metacluster names (joined with "+") are now matched
  correctly when the dendrogram is built. The documented default of
  `reuseFeatures` is corrected to `FALSE` (#410).
* `findMarkersTree()` keeps only one-off splits that mark exactly one
  cluster. Splits that marked no cluster, or whose cut also took in half or
  more of another cluster's cells, were kept before and threw off how
  features were grouped by cluster: markers could be attributed to the wrong
  cluster, and an empty node could make `findMarkersTree()` or
  `plotDendro()` fail with "argument is of length zero" or "replacement has
  length zero". Markers can change for trees that had such splits (#440).
* With metaclusters, `findMarkersTree()` now merges the metaclusters of a
  balanced split at the top level and separates them in their own subtree,
  as intended. Before, the merged tree was checked against the tree before
  the merge, so the merge was always rejected and the threshold was lowered
  instead, which could leave a metacluster without a positive marker (#442).
  When a merge was accepted, single-cluster metaclusters were not found in
  the metacluster list, so the merged metacluster got no subtree of its own
  and its clusters were never separated (#441). Merged metaclusters are now
  named without the "(cluster)" of single-cluster metaclusters, for example
  "M1+M2", as for metaclusters merged for lack of markers. Trees, markers,
  and performance can change for metacluster runs that had a balanced split
  at the top level.
* With metaclusters, `findMarkersTree()` no longer lowers the top-level
  threshold indefinitely when `reuseFeatures = TRUE`, or when merging the
  metaclusters of a top-level balanced split can't give them a positive
  marker. In both settings the threshold is lowered only while an unused
  feature could still become a marker, then metaclusters without one are
  merged, with a warning, as with `reuseFeatures = FALSE` (#438, #455).
  Metaclusters that share the marker of a top-level balanced split now count
  as having no marker of their own. A balanced split of a single metacluster
  against several, or of every metacluster, is no longer merged (#455).
* `findMarkersTree()` keeps metacluster names that contain parentheses, such
  as "T (CD4)", in `metaclusterLabels`, `classLabels`, the `metacluster` and
  `class` columns of `rules` and `branchPoints`, and `plotDendro()`. Before,
  everything from the first "(" was treated as the added "(cluster)" (#455).
  Each metacluster's subtree in the dendrogram now gets only its own
  clusters: with 10 or more metaclusters, `M1` also took those of `M10`,
  `M11` and so on (#444). A metacluster whose subtree has only one feature
  left no longer fails with "$ operator is invalid for atomic vectors". When
  no feature is left, because all are the metacluster's top-level markers,
  the error now says so and suggests `reuseFeatures = TRUE`.
* `findMarkersTree()` with `oneoffMetric = "pairwise AUC"` no longer warns
  that pROC's `coords(transpose = TRUE)` is deprecated (#447), and
  `plotDendro()` no longer warns that ggplot2's `label.size` argument of
  `geom_label()` is deprecated (#445). Results and plots are unchanged.
  celda now requires ggplot2 4.0.0 or later and pROC 1.16.0 or later.

# celda v1.29.3 (2026-10-09)
* The `decontX` vignette is now a short pointer to the decontX package's
  vignette, which has the full, maintained DecontX tutorial; celda no longer
  keeps its own copy. The duplicate pkgdown article (`decontX_pbmc4k`) was
  removed, and the site's DecontX menu entry now links to the decontX vignette
  on Bioconductor.

# celda v1.29.2 (2026-10-09)
* Fixed the `decontX` vignette, which failed to build on Bioconductor devel
  with "'groupClusters' not found in 'z': 6". It hard-coded cluster numbers
  from decontX's heuristic clustering, which differ between decontX versions
  and platforms. The vignette now maps each cell type to the cluster where its
  marker genes are detected most often, as the decontX package's own vignette
  does.

# celda v1.29.1 (2026-09-25)
* Replaced deprecated `ggplot2::aes_string()` (and `size` for line geoms) with
  tidy-eval `aes()` / `linewidth` in the dimension-reduction, violin, and
  perplexity plotting functions, silencing ggplot2 deprecation warnings
  (#418).
* Migrated the DecontX implementation to the standalone `decontX` Bioconductor
  package. celda's `decontX()`, `decontXcounts()`, `plotDecontXContamination()`,
  `plotDecontXMarkerExpression()`, `plotDecontXMarkerPercentage()`, and
  `simulateContamination()` are now thin, backward-compatible wrappers around
  the `decontX` package (added to `Suggests`).
* Removed celda's DecontX C++ backend (`src/DecontX.cpp`) and dropped the
  DecontX-only dependencies `MCMCprecision`, `scater`, `scran`, and
  `DelayedArray` from `Imports` (`scater` retained in `Suggests` for the
  vignettes).
* Bioconductor release readiness: replaced `T`/`F` logical literals with
  `TRUE`/`FALSE` in `findMarkersTree.R` (BiocCheck WARNING); added a `URL` field
  and bumped the R dependency to `R (>= 4.6.0)` in DESCRIPTION; set
  `LazyData: false`; added labels to unlabeled vignette knitr chunks; and added
  the required `url` field to `_pkgdown.yml`.
* Declared `Seurat` in `Suggests` and namespaced the `Seurat::Idents()` and
  `Seurat::RunUMAP()` calls in `findMarkersTree()`'s Seurat-object input path,
  clearing the remaining R CMD check "no visible global function" NOTE.

# celda v1.22.1 (2024-10-29)
* Fixed issue with enrichR not being loaded

# celda v1.22.0 (2024-10-29)
* Update to match Bioconductor release version

# celda v1.18.2 (2024-04-02)
* Updated Makevar files to new CRAN standards 
* Fixed unit test causing error

# celda v1.18.1 (2023-11-05)
* Update to match Bioconductor release version
* Removed multipanelfigure as a dependency

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
