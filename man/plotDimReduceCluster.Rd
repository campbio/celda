% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_dr.R
\name{plotDimReduceCluster}
\alias{plotDimReduceCluster}
\alias{plotDimReduceCluster,SingleCellExperiment-method}
\alias{plotDimReduceCluster,vector-method}
\title{Plotting the cell labels on a dimension reduction plot}
\usage{
plotDimReduceCluster(
  x,
  reducedDimName,
  altExpName = "featureSubset",
  dim1 = NULL,
  dim2 = NULL,
  size = 0.5,
  xlab = NULL,
  ylab = NULL,
  specificClusters = NULL,
  labelClusters = FALSE,
  groupBy = NULL,
  labelSize = 3.5
)

\S4method{plotDimReduceCluster}{SingleCellExperiment}(
  x,
  reducedDimName,
  altExpName = "featureSubset",
  dim1 = 1,
  dim2 = 2,
  size = 0.5,
  xlab = NULL,
  ylab = NULL,
  specificClusters = NULL,
  labelClusters = FALSE,
  groupBy = NULL,
  labelSize = 3.5
)

\S4method{plotDimReduceCluster}{vector}(
  x,
  dim1,
  dim2,
  size = 0.5,
  xlab = "Dimension_1",
  ylab = "Dimension_2",
  specificClusters = NULL,
  labelClusters = FALSE,
  groupBy = NULL,
  labelSize = 3.5
)
}
\arguments{
\item{x}{Integer vector of cell cluster labels or a
\linkS4class{SingleCellExperiment} object
containing cluster labels for each cell in \code{"celda_cell_cluster"}
column in \code{colData(x)}.}

\item{reducedDimName}{The name of the dimension reduction slot in
\code{reducedDimNames(x)} if \code{x} is a
\linkS4class{SingleCellExperiment} object. Ignored if both \code{dim1} and
\code{dim2} are set.}

\item{altExpName}{The name for the \link{altExp} slot
to use. Default "featureSubset".}

\item{dim1}{Integer or numeric vector. If \code{reducedDimName} is supplied,
then, this will be used as an index to determine which dimension will be
plotted on the x-axis. If \code{reducedDimName} is not supplied, then this
should be a vector which will be plotted on the x-axis. Default \code{1}.}

\item{dim2}{Integer or numeric vector. If \code{reducedDimName} is supplied,
then, this will be used as an index to determine which dimension will be
plotted on the y-axis. If \code{reducedDimName} is not supplied, then this
should be a vector which will be plotted on the y-axis. Default \code{2}.}

\item{size}{Numeric. Sets size of point on plot. Default \code{0.5}.}

\item{xlab}{Character vector. Label for the x-axis. Default \code{NULL}.}

\item{ylab}{Character vector. Label for the y-axis. Default \code{NULL}.}

\item{specificClusters}{Numeric vector.
Only color cells in the specified clusters.
All other cells will be grey.
If NULL, all clusters will be colored. Default \code{NULL}.}

\item{labelClusters}{Logical. Whether the cluster labels are plotted.
Default FALSE.}

\item{groupBy}{Character vector. Contains sample labels for each cell.
If NULL, all samples will be plotted together. Default NULL.}

\item{labelSize}{Numeric. Sets size of label if labelClusters is TRUE.
Default 3.5.}
}
\value{
The plot as a ggplot object
}
\description{
Create a scatterplot for each row of a normalized
 gene expression matrix where x and y axis are from a
 data dimension reduction tool.
 The cells are colored by "celda_cell_cluster" column in
  \code{colData(altExp(x, altExpName))} if \code{x} is a
  \linkS4class{SingleCellExperiment} object, or \code{x} if \code{x} is
  a integer vector of cell cluster labels.
}
\examples{
data(sceCeldaCG)
sce <- celdaTsne(sceCeldaCG)
plotDimReduceCluster(x = sce,
  reducedDimName = "celda_tSNE",
  specificClusters = c(1, 2, 3))
library(SingleCellExperiment)
data(sceCeldaCG, celdaCGMod)
sce <- celdaTsne(sceCeldaCG)
plotDimReduceCluster(x = celdaClusters(celdaCGMod)$z,
  dim1 = reducedDim(altExp(sce), "celda_tSNE")[, 1],
  dim2 = reducedDim(altExp(sce), "celda_tSNE")[, 2],
  specificClusters = c(1, 2, 3))
}
