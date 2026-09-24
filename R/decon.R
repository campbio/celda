#' @title Contamination estimation with decontX
#'
#' @description Identifies contamination from factors such as ambient RNA
#' in single cell genomic datasets.
#'
#' As of celda 1.23.0 the DecontX algorithm lives in the standalone
#' \pkg{decontX} package. The functions here are thin wrappers that call
#' \code{\link[decontX]{decontX}} so the implementation is maintained in a
#' single place. The \pkg{decontX} package must be installed
#' (\code{BiocManager::install("decontX")}); the wrappers accept and return
#' exactly what \code{\link[decontX]{decontX}} does.
#'
#' @name decontX
#'
#' @param x A numeric matrix of counts or a \linkS4class{SingleCellExperiment}
#' with the matrix located in the assay slot under \code{assayName}.
#' This object should only contain filtered cells after cell calling.
#' @param ... Additional arguments passed to \code{\link[decontX]{decontX}}
#' (for example \code{assayName}, \code{z}, \code{batch}, \code{background},
#' \code{maxIter}, \code{delta}, \code{estimateDelta}, \code{varGenes},
#' \code{dbscanEps}, \code{seed}, \code{verbose}). See
#' \code{\link[decontX]{decontX}} for the full argument documentation.
#'
#' @return If \code{x} is a matrix-like object, a list is returned with the
#' decontaminated matrix (\code{decontXcounts}), per-cell contamination
#' (\code{contamination}), estimated parameters (\code{estimates}), cluster
#' labels (\code{z}), and the run parameters (\code{runParams}).
#'
#' If \code{x} is a \linkS4class{SingleCellExperiment}, the decontaminated
#' counts are stored as the \code{decontXcounts} assay (accessed with
#' \code{decontXcounts(x)}); contamination and cluster labels are stored in
#' \code{colData(x)}; \code{estimates} and \code{runParams} are stored in
#' \code{metadata(x)$decontX}; and the UMAPs used to generate cluster labels
#' are stored in the \code{reducedDims} slot.
#'
#' @author Shiyi Yang, Yuan Yin, Joshua Campbell
#'
#' @seealso \code{\link[decontX]{decontX}} in the \pkg{decontX} package for the
#' full implementation and argument documentation.
#'
#' @examplesIf requireNamespace("decontX", quietly = TRUE)
#' # Generate matrix with contamination
#' s <- simulateContamination(seed = 12345)
#'
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(list(counts = s$observedCounts))
#' sce <- decontX(sce)
#'
#' # Plot contamination on UMAP
#' plotDecontXContamination(sce)
#'
#' # Plot decontX cluster labels
#' umap <- reducedDim(sce)
#' plotDimReduceCluster(x = sce$decontX_clusters,
#'     dim1 = umap[, 1], dim2 = umap[, 2], )
#'
#' # Plot percentage of marker genes detected
#' # in each cell cluster before decontamination
#' s$markers
#' plotDecontXMarkerPercentage(sce, markers = s$markers, assayName = "counts")
#'
#' # Plot percentage of marker genes detected
#' # in each cell cluster after contamination
#' plotDecontXMarkerPercentage(sce, markers = s$markers,
#'                             assayName = "decontXcounts")
#'
#' # Plot percentage of marker genes detected in each cell
#' # comparing original and decontaminated counts side-by-side
#' plotDecontXMarkerPercentage(sce, markers = s$markers,
#'                             assayName = c("counts", "decontXcounts"))
#'
#' # Plot raw counts of indiviual markers genes before
#' # and after decontamination
#' plotDecontXMarkerExpression(sce, unlist(s$markers))
NULL


# Ensure the decontX package (which now owns the implementation) is available.
.requireDecontX <- function() {
  if (!requireNamespace("decontX", quietly = TRUE)) {
    stop("The 'decontX' package is required for this function. Install it ",
         "with BiocManager::install('decontX').", call. = FALSE)
  }
  invisible(TRUE)
}


#' @export
#' @rdname decontX
setGeneric("decontX", function(x, ...) standardGeneric("decontX"))


#########################
# Setting up S4 methods #
#########################


#' @export
#' @rdname decontX
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom Matrix dgCMatrix
setMethod("decontX", "SingleCellExperiment", function(x, ...) {
  .requireDecontX()
  decontX::decontX(x, ...)
})

#' @export
#' @rdname decontX
setMethod("decontX", "ANY", function(x, ...) {
  .requireDecontX()
  decontX::decontX(x, ...)
})


## Copied from SingleCellExperiment Package

GET_FUN <- function(exprs_values, ...) {
  (exprs_values) # To ensure evaluation
  function(object, ...) {
    SummarizedExperiment::assay(object, i = exprs_values, ...)
  }
}

SET_FUN <- function(exprs_values, ...) {
  (exprs_values) # To ensure evaluation
  function(object, ..., value) {
    SummarizedExperiment::assay(object, i = exprs_values, ...) <- value
    object
  }
}



#' @title Get or set decontaminated counts matrix
#'
#' @description Gets or sets the decontaminated counts matrix from a
#' a \linkS4class{SingleCellExperiment} object.
#' @name decontXcounts
#' @param object A \linkS4class{SingleCellExperiment} object.
#' @param value A matrix to save as an assay called \code{decontXcounts}
#' @param ... For the generic, further arguments to pass to each method.
#' @return If getting, the assay from \code{object} with the name
#' \code{decontXcounts} will be returned. If setting, a
#' \linkS4class{SingleCellExperiment} object will be returned with
#' \code{decontXcounts} listed in the \code{assay} slot.
#' @seealso \code{\link{assay}} and \code{\link{assay<-}}
NULL

#' @export
#' @rdname decontXcounts
setGeneric("decontXcounts", function(object, ...) {
  standardGeneric("decontXcounts")
})

#' @export
#' @rdname decontXcounts
setGeneric("decontXcounts<-", function(object, ..., value) {
  standardGeneric("decontXcounts<-")
})

#' @export
#' @rdname decontXcounts
setMethod("decontXcounts", "SingleCellExperiment", GET_FUN("decontXcounts"))

#' @export
#' @rdname decontXcounts
setMethod(
  "decontXcounts<-", c("SingleCellExperiment", "ANY"),
  SET_FUN("decontXcounts")
)




#########################
# Simulating Data       #
#########################

#' @title Simulate contaminated count matrix
#' @description This function generates a list containing two count matrices --
#'  one for real expression, the other one for contamination, as well as other
#'  parameters used in the simulation which can be useful for running
#'  decontamination.
#'
#'  As of celda 1.23.0 this is a thin wrapper around
#'  \code{\link[decontX]{simulateContamination}} in the \pkg{decontX} package,
#'  which must be installed.
#' @param ... Arguments passed to
#'  \code{\link[decontX]{simulateContamination}} (for example \code{C},
#'  \code{G}, \code{K}, \code{NRange}, \code{beta}, \code{delta},
#'  \code{numMarkers}, \code{seed}).
#' @return A list containing the \code{nativeMatrix} (real expression),
#' \code{observedMatrix} (real expression + contamination), as well as other
#' parameters used in the simulation.
#' @author Shiyi Yang, Yuan Yin, Joshua Campbell
#' @seealso \code{\link[decontX]{simulateContamination}}
#' @examplesIf requireNamespace("decontX", quietly = TRUE)
#' contaminationSim <- simulateContamination(K = 3, delta = c(1, 10))
#' @export
simulateContamination <- function(...) {
  .requireDecontX()
  decontX::simulateContamination(...)
}
