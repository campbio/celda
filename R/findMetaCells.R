#' @title Find Meta-Cells per Sample
#' @description Groups the cells of each sample into meta-cells: small sets
#'  of very similar cells that can be used in place of single cells to make
#'  the data less sparse and downstream models faster. Each sample is handled
#'  on its own, so a meta-cell never mixes cells from two samples.
#' @param inSCE A
#'  \link[SingleCellExperiment:SingleCellExperiment-class]{SingleCellExperiment}
#'  object with raw counts in the \code{"counts"} assay.
#' @param initK Integer. Smallest number of cell clusters to try in each
#'  sample. Default \code{5}.
#' @param maxK Integer. Largest number of cell clusters to try. If NULL, it is
#'  set per sample to the square root of the number of cells, kept between
#'  10 and 100. Default NULL.
#' @param minCellMeta Integer. Minimum number of cells in a meta-cell. Smaller
#'  clusters are merged into their nearest neighbour, so cell states rarer
#'  than this are absorbed. Default \code{10}.
#' @param tempL Integer. Number of temporary gene modules used while
#'  splitting cells. Cannot be NULL, because these modules are also used to
#'  find the nearest neighbour when merging small meta-cells. Default
#'  \code{50}.
#' @param sampleLabel Vector or factor. Denotes the sample label for each cell
#'  (column) in \code{inSCE}, e.g. \code{inSCE$sample}. Meta-cells are built
#'  within each sample. If NULL, all cells are treated as one sample,
#'  \code{"Sample_1"}. Default NULL.
#' @param kneePoint Integer, 1 or 2. Which knee of the perplexity curve sets
#'  K: \code{1} takes the first (sharpest) bend, \code{2} the next bend after
#'  it, which gives more, smaller meta-cells. Default \code{2}.
#' @param seed Integer passed to \code{recursiveSplitCell} so results can be
#'  reproduced. Default \code{12345}.
#' @param verbose Logical. Whether to print log messages. Default TRUE.
#' @param logfile Character. File to write the log to. It is started fresh
#'  for each run and keeps the messages of every sample, including each
#'  sample's \code{recursiveSplitCell} log. If NULL, messages are printed to
#'  the console. Default NULL.
#' @returns The input
#'  \link[SingleCellExperiment:SingleCellExperiment-class]{SingleCellExperiment},
#'  with cells in the same order and a new colData column \code{metacell},
#'  labelled \code{<sample>_K<n>}. The labels can be passed to
#'  \link{recursiveSplitModule} as \code{zInit}. A sample with fewer than
#'  \code{minCellMeta} cells becomes a single meta-cell,
#'  \code{<sample>_K1}, with a message.
#' @export
#'
#' @examples
#' library(celda)
#' # Simulated counts with 2 samples
#' sce <- simulateCells("celda_CG",
#'   S = 2, K = 5, L = 10, G = 200,
#'   CRange = c(150, 200), seed = 12345
#' )
#' sce <- findMetaCells(sce,
#'   initK = 3, maxK = 8, tempL = 20,
#'   sampleLabel = sce$celda_sample_label,
#'   minCellMeta = 10, verbose = FALSE
#' )
#' table(sce$metacell)
findMetaCells <- function(
  inSCE,
  initK = 5,
  maxK = NULL,
  minCellMeta = 10,
  tempL = 50,
  sampleLabel = NULL,
  kneePoint = 2,
  seed = 12345,
  logfile = NULL,
  verbose = TRUE
) {
  # Check the arguments before any slow work starts
  if (!kneePoint %in% c(1L, 2L)) {
    stop("'kneePoint' must be 1 or 2.", call. = FALSE)
  }
  if (isTRUE(initK == 0) || (!is.null(maxK) && isTRUE(maxK == 0))) {
    stop("Value of K cannot be zero.", call. = FALSE)
  }
  validTempL <- !is.null(tempL) && is.numeric(tempL) &&
    length(tempL) == 1 && tempL >= 1
  if (!validTempL) {
    stop("'tempL' must be a positive integer.", call. = FALSE)
  }
  # Same rules as other celda functions: one label per cell, NULL = one sample
  hasSampleColumn <- "sample" %in%
    colnames(SummarizedExperiment::colData(inSCE))
  if (is.null(sampleLabel) && hasSampleColumn) {
    warning("'sampleLabel' is NULL, so all cells are treated as one sample. ",
      "colData(inSCE)$sample is not used; pass ",
      "sampleLabel = inSCE$sample to build meta-cells per sample.",
      call. = FALSE
    )
  }
  if (anyNA(sampleLabel)) {
    stop("'sampleLabel' contains NA values.", call. = FALSE)
  }
  sampleLabel <- .processSampleLabels(sampleLabel,
    numCells = ncol(inSCE)
  )
  samples <- levels(sampleLabel)
  startTime <- Sys.time()
  # append = FALSE: every run starts a new log file
  .logMessages("Meta-cell aggregation starts!\nSystem time: ",
    startTime,
    logfile = logfile,
    append = FALSE,
    verbose = verbose
  )
  # One label per input cell, filled in sample by sample, so the result
  # keeps the input's cell order
  metacell <- rep(NA_character_, ncol(inSCE))
  for (s in samples) {
    idx <- which(sampleLabel == s)
    sce <- inSCE[, idx]
    nCells <- ncol(sce)
    # A sample smaller than one meta-cell becomes a single meta-cell.
    # verbose = TRUE so the note is always printed (or written to the log)
    if (nCells < minCellMeta) {
      .logMessages("Sample '", s, "' has ", nCells,
        " cells, fewer than minCellMeta (", minCellMeta,
        "), so its cells form one meta-cell.",
        sep = "",
        logfile = logfile,
        append = TRUE,
        verbose = TRUE
      )
      metacell[idx] <- paste0(s, "_K1")
      next
    }

    # Upper limit for K: sqrt(cells), at least 10 and at most 100,
    # unless the user gave one
    sampleMaxK <- if (is.null(maxK)) {
      pmin(pmax(round(sqrt(nCells)), 10), 100)
    } else {
      maxK
    }
    # Keep genes seen in this sample (>= 3 counts in >= 3 cells); they are
    # stored in the "featureSubset" altExp
    sce <- selectFeatures(sce,
      minCount = 3, minCell = 3,
      useAssay = "counts"
    )
    .logMessages("Sample: ",
      s,
      logfile = logfile,
      append = TRUE,
      verbose = verbose
    )
    # ------------------------------------
    # 1) recursive cell split to explore K
    # ------------------------------------
    .logMessages("running recursiveSplitCell from K = ",
      initK, " to K = ", sampleMaxK,
      logfile = logfile,
      append = TRUE,
      verbose = verbose
    )
    # recursiveSplitCell starts its own log with append = FALSE, which would
    # erase this log. Give it a temporary file and append that to the log.
    rscLog <- if (is.null(logfile)) NULL else tempfile(fileext = ".log")
    sce <- recursiveSplitCell(
      sce,
      initialK = initK,
      maxK = sampleMaxK,
      tempL = tempL,
      sampleLabel = rep(s, nCells),
      useAssay = "counts",
      seed = seed,
      logfile = rscLog,
      verbose = verbose
    )
    if (!is.null(rscLog) && file.exists(rscLog)) {
      file.append(logfile, rscLog)
      unlink(rscLog)
    }

    # ------------------------
    # 2) perplexity across K's
    # ------------------------
    # One model per K was fitted above; read back the perplexity of each
    altExpObj <- SingleCellExperiment::altExp(sce, "featureSubset")
    gridSearch <- S4Vectors::metadata(altExpObj)$celda_grid_search
    perplexDf <- runParams(gridSearch)[, c("K", "mean_perplexity")]

    # ------------------------------------------
    # 3) run kneedle algorithm to find optimal k
    # ------------------------------------------
    # Perplexity falls as K grows; the knee is where adding clusters stops
    # paying off
    bestK <- .findKnee(perplexDf$K, perplexDf$mean_perplexity,
      nPoints = kneePoint
    )
    # A knee at maxK usually means maxK was too small for this sample
    if (bestK == max(perplexDf$K)) {
      warning("Sample '", s, "': selected K is at maxK (", bestK,
        "); the upper bound is binding, not the curve.",
        call. = FALSE
      )
    }
    # ---------------------------------------------
    # 4) Add updated metacell labels to the colData
    # ---------------------------------------------
    # Take the model with K = bestK, then merge clusters below minCellMeta
    sce <- subsetCeldaList(sce, list(K = bestK))
    sce <- .updateMetaCellLabel(sce, sampleName = s, minCellMeta = minCellMeta)
    metacell[idx] <- sce$metacell
  }

  elapsed <- difftime(Sys.time(), startTime, units = "auto")
  .logMessages("Completed meta-cell aggregation. Total time: ",
    format(round(elapsed, 2)),
    logfile = logfile,
    append = TRUE,
    verbose = verbose
  )
  SummarizedExperiment::colData(inSCE)$metacell <- metacell
  inSCE
}

# Kneedle algorithm: find the "knee" (point of strongest bend) of a falling
# curve P(K).
#  1. Scale K and P to 0-1 so the units do not matter.
#  2. Draw a straight line from the first point to the last point.
#  3. The knee is the K where the curve lies furthest below that line
#     (measured straight down).
# For nPoints > 1, repeat on the part of the curve from the knee to the end,
# so nPoints = 2 returns the second knee.
.findKnee <- function(values, perplexity, nPoints = 1L) {
  # A flat curve has no knee; return the first value
  if (min(values) == max(values) || min(perplexity) == max(perplexity)) {
    return(values[1])
  }
  x <- (values - min(values)) / (max(values) - min(values))
  y <- (perplexity - min(perplexity)) / (max(perplexity) - min(perplexity))

  x1 <- x[1]
  y1 <- y[1]
  xl <- x[length(x)]
  yl <- y[length(y)]

  # Straight line from the first to the last point, and the gap below it
  strLine <- y1 + (((yl - y1) / (xl - x1)) * (x - x1))
  i <- which.max(strLine - y)
  if (nPoints == 1L) {
    values[i]
  } else {
    rest <- seq(i, length(values))
    .findKnee(values[rest], perplexity[rest], nPoints = nPoints - 1L)
  }
}

# Turn celda cell clusters into meta-cell labels and make sure every
# meta-cell has at least minCellMeta cells. The smallest cluster is merged
# into the nearest cluster, one at a time, until none are too small.
# "Nearest" is measured on the temporary gene modules (tempY): each cell is
# scored per module (log CPM, scaled), and clusters are compared by the
# Euclidean distance between their average scores.
.updateMetaCellLabel <- function(sce, sampleName, minCellMeta = 10) {
  colDat <- SingleCellExperiment::colData(sce)
  colDat$metacell <- as.character(celdaClusters(sce))

  smallCluster <- names(which(table(colDat$metacell) < minCellMeta))
  if (length(smallCluster) > 0) {
    counts <- SummarizedExperiment::assay(
      SingleCellExperiment::altExp(sce, "featureSubset"), "counts"
    )

    # Temporary gene modules saved by the local recursiveSplitCell
    tempY <- S4Vectors::metadata(
      SingleCellExperiment::altExp(sce, "featureSubset")
    )[["tempY"]]

    # Module-by-cell counts -> log CPM -> scaled cell-by-module scores
    countsY <- .rowSumByGroup(
      counts,
      tempY,
      length(unique(tempY))
    )

    countsY <- normalizeCounts(
      countsY,
      scaleFactor = 1e6,
      transformationFun = log1p
    )

    countsY <- scale(t(countsY)) # cell-by-module matrix
    colnames(countsY) <- paste0("L", seq_len(ncol(countsY)))
    rownames(countsY) <- colnames(sce)

    while (length(smallCluster) > 0) {
      # Average module scores of each current cluster
      clusters <- unique(colDat$metacell)
      centroidK <- matrix(
        nrow = length(clusters),
        ncol = ncol(countsY)
      )
      rownames(centroidK) <- clusters
      colnames(centroidK) <- colnames(countsY)

      for (i in seq_along(clusters)) {
        cellsInK <- which(colDat$metacell == clusters[i])
        centroidK[i, ] <- colMeans(countsY[cellsInK, , drop = FALSE])
      }

      # Merge the smallest too-small cluster into its nearest neighbour
      clustDist <- as.matrix(stats::dist(x = centroidK, method = "euclidean"))
      s <- smallCluster[which.min(table(colDat$metacell)[smallCluster])]

      dvec <- clustDist[s, ]
      dvec <- dvec[names(dvec) != s]

      # Only one cluster left: nothing to merge into
      if (length(dvec) == 0) break
      newCluster <- names(which.min(dvec))

      colDat$metacell[colDat$metacell == s] <- newCluster
      smallCluster <- names(which(table(colDat$metacell) < minCellMeta))
    }
  }
  # Renumber the surviving clusters 1..n and name them <sample>_K<n>
  prevLabel <- unique(colDat$metacell)
  newLabel <- seq_along(prevLabel)
  labelMap <- stats::setNames(newLabel, prevLabel)

  SummarizedExperiment::colData(sce)$metacell <- paste0(
    sampleName, "_K",
    labelMap[as.character(colDat$metacell)]
  )

  sce
}
