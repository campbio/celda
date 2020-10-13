.sampleLl <- function(llProbs) {
  probsSub <- exp(llProbs - max(llProbs))
  probsNorm <- probsSub / sum(probsSub)
  probsSelect <- sample.int(
    length(probsNorm),
    size = 1L,
    replace = TRUE,
    prob = probsNorm
  )
  return(probsSelect)
}


.cosineDist <- function(x) {
  x <- t(x)
  y <- (1 - .cosine(x)) / 2
  return(stats::as.dist(y))
}


.cosine <- function(x) {
  y <- x %*% t(x) / (sqrt(rowSums(x^2) %*% t(rowSums(x^2))))
  return(y)
}


.spearmanDist <- function(x) {
  y <- (1 - stats::cor(x, method = "spearman")) / 2
  return(stats::as.dist(y))
}


.hellingerDist <- function(x) {
  y <- stats::dist(t(sqrt(x)), method = "euclidean") * 1 / sqrt(2)
  return(y)
}


.normalizeLogProbs <- function(llProbs) {
  llProbs <- exp(sweep(llProbs, 1, base::apply(llProbs, 1, max), "-"))
  probs <- sweep(llProbs, 1, rowSums(llProbs), "/")
  return(probs)
}


#' @title Normalization of count data
#' @description Performs normalization, transformation, and/or scaling of a
#'  counts matrix
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells.
#' @param normalize Character.
#'  Divides counts by the library sizes for each cell. One of 'proportion',
#'  'cpm', 'median', or 'mean'. 'proportion' uses the total counts for each
#'  cell as the library size. 'cpm' divides the library size of each cell by
#'  one million to produce counts per million. 'median' divides the library
#'  size of each cell by the median library size across all cells. 'mean'
#'  divides the library size of each cell by the mean library size across all
#'  cells.
#' @param scaleFactor Numeric. Sets the scale factor for cell-level
#'  normalization. This scale factor is multiplied to each cell after the
#'  library size of each cell had been adjusted in \code{normalize}. Default
#'  \code{NULL} which means no scale factor is applied.
#' @param transformationFun Function. Applys a transformation such as
#'  \link{sqrt}, \link{log}, \link{log2}, \link{log10}, or \link{log1p}.
#'  If NULL, no transformation will be applied. Occurs after normalization.
#'  Default NULL.
#' @param scaleFun Function. Scales the rows of the normalized and transformed
#'  count matrix. For example, 'scale' can be used to z-score normalize the
#'  rows. Default NULL.
#' @param pseudocountNormalize Numeric. Add a pseudocount to counts before
#'  normalization. Default 0.
#' @param pseudocountTransform Numeric. Add a pseudocount to normalized counts
#'  before applying the transformation function. Adding a pseudocount
#'  can be useful before applying a log transformation. Default  0.
#' @return Numeric Matrix. A normalized matrix.
#' @examples
#' data(celdaCGSim)
#' normalizedCounts <- normalizeCounts(celdaCGSim$counts, "proportion",
#'   pseudocountNormalize = 1)
#' @export
normalizeCounts <- function(counts,
                            normalize = c("proportion", "cpm",
                              "median", "mean"),
                            scaleFactor = NULL,
                            transformationFun = NULL,
                            scaleFun = NULL,
                            pseudocountNormalize = 0,
                            pseudocountTransform = 0) {

  normalize <- match.arg(normalize)
  counts <- as.matrix(counts)

  if (!is.null(transformationFun) &&
    !is.function(transformationFun)) {
    stop("'transformationFun' needs to be of class 'function'")
  }
  if (!is.null(scaleFun) && !is.function(scaleFun)) {
    stop("'scaleFun' needs to be of class 'function'")
  }
  # Perform normalization
  if (normalize == "proportion") {
    norm <- fastNormProp(counts, pseudocountNormalize)
  } else {
    counts <- counts + pseudocountNormalize
    cs <- .colSums(counts, nrow(counts), ncol(counts))
    norm <- switch(
      normalize,
      "proportion" = sweep(counts, 2, cs, "/"),
      "cpm" = sweep(counts, 2, cs / 1e6, "/"),
      "median" = sweep(counts, 2, cs / stats::median(cs), "/"),
      "mean" = sweep(counts, 2, cs / mean(cs), "/")
    )
  }

  if (!is.null(scaleFactor)) {
      norm <- norm * scaleFactor
  }

  if (!is.null(transformationFun)) {
    norm <- do.call(
      transformationFun,
      list(norm + pseudocountTransform)
    )
  }
  if (!is.null(scaleFun)) {
    norm <- t(base::apply(norm, 1, scaleFun))
  }
  colnames(norm) <- colnames(counts)
  rownames(norm) <- rownames(counts)
  return(norm)
}


#' @title Recode cell cluster labels
#' @description Recode cell subpopulaton clusters using a mapping in the
#'  \code{from} and \code{to} arguments.
#' @param sce \linkS4class{SingleCellExperiment} object returned from
#'  \link{celda_C} or \link{celda_CG}. Must contain column
#'  \code{celda_cell_cluster} in
#'  \code{\link[SummarizedExperiment]{colData}(altExp(sce, altExpName))}.
#' @param from Numeric vector. Unique values in the range of
#'  \code{seq(celdaClusters(sce, altExpName = altExpName))} that correspond to
#'  the original cluster
#'  labels in \code{sce}.
#' @param to Numeric vector. Unique values in the range of
#'  \code{seq(celdaClusters(sce, altExpName = altExpName))} that correspond to
#'  the new cluster labels.
#' @param altExpName The name for the \link[SingleCellExperiment]{altExp} slot
#'  to use. Default "featureSubset".
#' @return \linkS4class{SingleCellExperiment} object with recoded cell
#'  cluster labels.
#' @examples
#' data(sceCeldaCG)
#' sceReorderedZ <- recodeClusterZ(sceCeldaCG, c(1, 3), c(3, 1))
#' @importFrom plyr mapvalues
#' @export
recodeClusterZ <- function(sce, from, to, altExpName = "featureSubset") {
    if (length(setdiff(from, to)) != 0) {
        stop("All values in 'from' must have a mapping in 'to'")
    }
    if (is.null(celdaClusters(sce, altExpName = altExpName))) {
        stop("Provided 'sce' argument does not have a 'celda_cell_cluster'",
            " column in 'colData(altExp(sce, altExpName))'")
    }
    celdaClusters(sce, altExpName = altExpName) <- plyr::mapvalues(
        celdaClusters(sce, altExpName = altExpName), from, to)
    return(sce)
}


# for deprecated celda model objects
.recodeClusterZ <- function(celdaMod, from, to) {
    if (length(setdiff(from, to)) != 0) {
        stop("All values in 'from' must have a mapping in 'to'")
    }
    if (is.null(celdaClusters(celdaMod)$z)) {
        stop("Provided celdaMod argument does not have a z attribute")
    }
    celdaMod@clusters$z <- plyr::mapvalues(celdaClusters(celdaMod)$z,
        from, to)
    return(celdaMod)
}


#' @title Recode feature module labels
#' @description Recode feature module clusters using a mapping in the
#'  \code{from} and \code{to} arguments.
#' @param sce \linkS4class{SingleCellExperiment} object returned from
#'  \link{celda_G} or \link{celda_CG}. Must contain column
#'  \code{celda_feature_module} in
#'  \code{\link[SummarizedExperiment]{rowData}(altExp(sce, altExpName))}.
#' @param from Numeric vector. Unique values in the range of
#'  \code{seq(celdaModules(sce))} that correspond to the original module labels
#'  in \code{sce}.
#' @param to Numeric vector. Unique values in the range of
#'  \code{seq(celdaModules(sce))} that correspond to the new module labels.
#' @param altExpName The name for the \link[SingleCellExperiment]{altExp} slot
#'  to use. Default "featureSubset".
#' @return @return \linkS4class{SingleCellExperiment} object with recoded
#'  feature module labels.
#' @examples
#' data(sceCeldaCG)
#' sceReorderedY <- recodeClusterY(sceCeldaCG, c(1, 3), c(3, 1))
#' @export
recodeClusterY <- function(sce, from, to, altExpName = "featureSubset") {
    if (length(setdiff(from, to)) != 0) {
        stop("All values in 'from' must have a mapping in 'to'")
    }
    if (is.null(celdaModules(sce, altExpName = altExpName))) {
        stop("Provided 'sce' argument does not have a 'celda_feature_module'",
            " column in 'rowData(altExp(sce, altExpName))'")
    }
    celdaModules(sce, altExpName = altExpName) <- plyr::mapvalues(
        celdaModules(sce, altExpName = altExpName), from, to)
    return(sce)
}


# for deprecated celda model objects
.recodeClusterY <- function(celdaMod, from, to) {
    if (length(setdiff(from, to)) != 0) {
        stop("All values in 'from' must have a mapping in 'to'")
    }
    if (is.null(celdaClusters(celdaMod)$y)) {
        stop("Provided celdaMod argument does not have a y attribute")
    }
    celdaMod@clusters$y <- plyr::mapvalues(celdaClusters(celdaMod)$y,
        from, to)
    return(celdaMod)
}


#' @title Check count matrix consistency
#' @description Checks if the counts matrix is the same one used to generate
#'  the celda model object by comparing dimensions and MD5 checksum.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells.
#' @param celdaMod A \code{celdaModel} or \code{celdaList} object.
#' @param errorOnMismatch Logical. Whether to throw an error in the event of
#'  a mismatch. Default TRUE.
#' @param ... Ignored. Placeholder to prevent check warning.
#' @return Returns TRUE if provided count matrix matches the one used in the
#'  celda object and/or \code{errorOnMismatch = FALSE}, FALSE otherwise.
#' @export
setGeneric("compareCountMatrix", function(counts, celdaMod, ...) {
    standardGeneric("compareCountMatrix")})


#' @rdname compareCountMatrix
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' compareCountMatrix(celdaCGSim$counts, celdaCGMod, errorOnMismatch = FALSE)
#' @export
setMethod("compareCountMatrix",
    signature(celdaMod = "celdaModel"),
    function(counts, celdaMod, errorOnMismatch = TRUE) {

        if ("y" %in% names(celdaClusters(celdaMod))) {
            if (nrow(counts) != length(celdaClusters(celdaMod)$y)) {
                stop(
                    "The provided celda object was generated from a counts",
                    " matrix with a different number of features than the one",
                    " provided."
                )
            }
        }

        if ("z" %in% names(celdaClusters(celdaMod))) {
            if (ncol(counts) != length(celdaClusters(celdaMod)$z)) {
                stop(
                    "The provided celda object was generated from a counts",
                    " matrix with a different number of cells than the one",
                    " provided."
                )
            }
        }
        celdaChecksum <- params(celdaMod)$countChecksum

        counts <- .processCounts(counts)
        # Checksums are generated in celdaGridSearch and model after processing
        count.md5 <- .createCountChecksum(counts)
        res <- isTRUE(count.md5 == celdaChecksum)
        if (res) {
            return(TRUE)
        }
        if (!res && errorOnMismatch) {
            stop(
                "There was a mismatch between the provided count matrix and",
                " the count matrix used to generate the provided celda result."
            )
        } else if (!res && !errorOnMismatch) {
            warning("There was a mismatch between the provided count matrix",
                " and the count matrix used to generate the provided celda",
                " result.")
            return(FALSE)
        }
    }
)


#' @rdname compareCountMatrix
#' @examples
#' data(celdaCGSim, celdaCGGridSearchRes)
#' compareCountMatrix(celdaCGSim$counts, celdaCGGridSearchRes,
#'     errorOnMismatch = FALSE)
#' @export
setMethod("compareCountMatrix",
    signature(celdaMod = "celdaList"),
    function(counts, celdaMod, errorOnMismatch = TRUE) {

        if ("y" %in% names(celdaMod@resList[[1]]@clusters)) {
            if (nrow(counts) != length(celdaMod@resList[[1]]@clusters$y)) {
                stop(
                    "The provided celda object was generated from a counts",
                    " matrix with a different number of features than the one",
                    " provided."
                )
            }
        }

        if ("z" %in% names(celdaMod@resList[[1]]@clusters)) {
            if (ncol(counts) != length(celdaMod@resList[[1]]@clusters$z)) {
                stop(
                    "The provided celda object was generated from a counts",
                    " matrix with a different number of cells than the one",
                    " provided."
                )
            }
        }
        celdaChecksum <- celdaMod@countChecksum

        counts <- .processCounts(counts)
        # Checksums are generated in celdaGridSearch and model after processing
        count.md5 <- .createCountChecksum(counts)
        res <- isTRUE(count.md5 == celdaChecksum)
        if (res) {
            return(TRUE)
        }
        if (!res && errorOnMismatch) {
            stop(
                "There was a mismatch between the provided count matrix and",
                " the count matrix used to generate the provided celda result."
            )
        } else if (!res && !errorOnMismatch) {
            warning("There was a mismatch between the provided count matrix",
                " and the count matrix used to generate the provided celda",
                " result.")
            return(FALSE)
        }
    }
)


.logMessages <- function(...,
                         sep = " ",
                         logfile = NULL,
                         append = FALSE,
                         verbose = TRUE) {
  if (isTRUE(verbose)) {
    if (!is.null(logfile)) {
      if (!is.character(logfile) || length(logfile) > 1) {
        stop(
          "The log file parameter needs to be a single character",
          " string."
        )
      }
      cat(paste(..., "\n", sep = sep),
        file = logfile,
        append = append
      )
    } else {
      message(paste(..., sep = sep))
    }
  }
}


#' @title Create a color palette
#' @description Generate a palette of `n` distinct colors.
#' @param n Integer. Number of colors to generate.
#' @param hues Character vector. Colors available from `colors()`. These will
#'  be used as the base colors for the clustering scheme in HSV. Different
#'  saturations and values will be generated for each hue. Default c("red",
#'  "cyan", "orange", "blue", "yellow", "purple", "green", "magenta").
#' @param saturationRange Numeric vector. A vector of length 2 denoting the
#'  saturation for HSV. Values must be in [0,1]. Default: c(0.25, 1).
#' @param valueRange Numeric vector. A vector of length 2 denoting the range
#'  of values for HSV. Values must be in [0,1]. Default: `c(0.5, 1)`.
#' @return A vector of distinct colors that have been converted to HEX from HSV.
#' @examples
#' colorPal <- distinctColors(6) # can be used in plotting functions
#' @importFrom grDevices colors
#' @importFrom grDevices rgb2hsv
#' @importFrom grDevices hsv
#' @export
distinctColors <- function(n,
                           hues = c(
                             "red",
                             "cyan",
                             "orange",
                             "blue",
                             "yellow",
                             "purple",
                             "green",
                             "magenta"
                           ),
                           saturationRange = c(0.7, 1),
                           valueRange = c(0.7, 1)) {
  if (!(all(hues %in% grDevices::colors()))) {
    stop(
      "Only color names listed in the 'color' function can be used in",
      " 'hues'"
    )
  }
  # Convert R colors to RGB and then to HSV color format
  huesHsv <- grDevices::rgb2hsv(grDevices::col2rgb(hues))
  # Calculate all combination of saturation/value pairs
  # Note that low saturation with low value (i.e. high darkness) is too dark
  # for all hues. Likewise, high saturation with high value (i.e. low
  # darkness) is hard to distinguish. Therefore, saturation and value are
  # set to be anticorrelated
  numVs <- ceiling(n / length(hues))
  s <- seq(
    from = saturationRange[1],
    to = saturationRange[2],
    length = numVs
  )
  v <- seq(
    from = valueRange[2],
    to = valueRange[1],
    length = numVs
  )
  # Create all combination of hues with saturation/value pairs
  list <- lapply(seq(numVs), function(x) {
    rbind(huesHsv[1, ], s[x], v[x])
  })
  newHsv <- do.call(cbind, list)
  # Convert to hex
  col <- grDevices::hsv(newHsv[1, ], newHsv[2, ], newHsv[3, ])
  return(col[seq(n)])
}


.processCounts <- function(counts) {
  counts <- as.matrix(counts)
  if (typeof(counts) != "integer") {
    counts <- round(counts)
    storage.mode(counts) <- "integer"
  }
  return(counts)
}


# Perform some simple checks on the counts matrix, to ensure celda modeling
# expectations are met
.validateCounts <- function(counts) {
  # And each row/column of the count matrix must have at least one count
  countRowSum <- rowSums(counts)
  countColSum <- colSums(counts)
  if (sum(countRowSum == 0) > 0 | sum(countColSum == 0) > 0) {
    stop(
      "Each row and column of the count matrix must have at least",
      " one count"
    )
  }
}


# Wrapper function, creates checksum for matrix.
# Feature names, cell names are not taken into account.
#' @importFrom digest digest
.createCountChecksum <- function(counts) {
  rownames(counts) <- NULL
  colnames(counts) <- NULL
  countChecksum <- digest::digest(counts, algo = "md5")
  return(countChecksum)
}


# Generate n random deviates from the Dirichlet function with shape parameters
# alpha. Adapted from gtools v3.5
.rdirichlet <- function(n, alpha) {
  l <- length(alpha)
  x <- matrix(stats::rgamma(l * n, alpha),
    ncol = l,
    byrow = TRUE
  )
  # Check for case where all sampled entries are zero due to round off
  # One entry will be randomly chosen to be one
  isZero <- rowSums(x) == 0
  assignment <- sample(seq(l), size = sum(isZero), replace = TRUE)
  x[cbind(which(isZero), assignment)] <- 1
  # Normalize
  sm <- x %*% rep(1, l)
  y <- x / as.vector(sm)
  return(y)
}


# Make sure provided sample labels are the right type,
# or generate some if none were provided
.processSampleLabels <- function(sampleLabel, numCells) {
  if (is.null(sampleLabel)) {
    sampleLabel <- as.factor(rep("Sample_1", numCells))
  } else {
    if (length(sampleLabel) != numCells) {
      stop(
        "'sampleLabel' must be the same length as the number of",
        " columns in the 'counts' matrix."
      )
    }
  }

  if (!is.factor(sampleLabel)) {
    sampleLabel <- as.factor(sampleLabel)
  }
  return(sampleLabel)
}


#' @title Output a feature module table
#' @description Creates a table that contains the list of features in
#'  each feature module.
#' @param sce A \linkS4class{SingleCellExperiment} object returned by
#'  \link{celda_G}, or \link{celda_CG}, with the matrix
#'  located in the \code{useAssay} assay slot.
#'  Rows represent features and columns represent cells.
#' @param useAssay A string specifying which \link[SummarizedExperiment]{assay}
#'  slot to use. Default "counts".
#' @param altExpName The name for the \link[SingleCellExperiment]{altExp} slot
#'  to use. Default "featureSubset".
#' @param outputFile File name for feature module table. If NULL, file will
#'  not be created. Default NULL.
#' @return Matrix. Contains a list of features per each column (feature module)
#' @examples
#' data(sceCeldaCG)
#' featureModuleTable(sceCeldaCG)
#' @importFrom stringi stri_list2matrix
#' @export
featureModuleTable <- function(sce, useAssay = "counts",
    altExpName = "featureSubset", outputFile = NULL) {

  factorizeMatrix <- factorizeMatrix(sce,
      useAssay = useAssay,
      altExpName = altExpName,
      type = "proportion")
  allGenes <- topRank(factorizeMatrix$proportions$module, n = nrow(sce))
  res <- as.data.frame(stringi::stri_list2matrix(allGenes$names))
  res <- apply(res, c(1, 2), function(x) {
    if (is.na(x)) {
      return("")
    } else {
      return(x)
    }
  })
  colnames(res) <- gsub(
    pattern = "V",
    replacement = "L",
    x = colnames(res)
  )
  if (is.null(outputFile)) {
    return(res)
  } else {
    utils::write.table(
      res,
      file = outputFile,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
}


#' @title Retrieve row index for a set of features
#' @description This will return indices of features among the rownames
#' or rowData of a data.frame, matrix, or a \linkS4class{SummarizedExperiment}
#' object including a \linkS4class{SingleCellExperiment}.
#' Partial matching (i.e. grepping) can be used by setting
#' \code{exactMatch = FALSE}.
#' @param features Character vector of feature names to find in the rows of
#' \code{x}.
#' @param x A data.frame, matrix, or \linkS4class{SingleCellExperiment}
#' object to search.
#' @param by Character. Where to search for features in \code{x}. If set to
#' \code{"rownames"} then the features will be searched for among
#' \code{rownames(x)}. If \code{x} inherits from class
#' \linkS4class{SummarizedExperiment}, then \code{by} can be one of the
#' fields in the row annotation data.frame (i.e. one of
#' \code{colnames(rowData(x))}).
#' @param exactMatch Boolean. Whether to only identify exact matches
#' or to identify partial matches using \code{\link{grep}}.
#' @param removeNA Boolean. If set to \code{FALSE}, features not found in
#' \code{x} will be given \code{NA} and the returned vector will be the same
#' length as \code{features}. If set to \code{TRUE}, then the \code{NA}
#' values will be removed from the returned vector. Default \code{FALSE}.
#' @return A vector of row indices for the matching features in \code{x}.
#' @author Yusuke Koga, Joshua Campbell
#' @seealso '\link[scater]{retrieveFeatureInfo}' from package \code{'scater'}
#' and \code{link{regex}} for how to use regular expressions when
#' \code{exactMatch = FALSE}.
#' @examples
#' data(celdaCGSim)
#' retrieveFeatureIndex(c("Gene_1", "Gene_5"), celdaCGSim$counts)
#' retrieveFeatureIndex(c("Gene_1", "Gene_5"), celdaCGSim$counts,
#'                                             exactMatch = FALSE)
#' @export
retrieveFeatureIndex <- function(features,
                                 x,
                                 by = "rownames",
                                 exactMatch = TRUE,
                                 removeNA = FALSE) {

  # Extract vector to search through
  if (by == "rownames") {
    if (is.null(rownames(x))) {
      stop("'rownames' of 'x' are 'NULL'. Please set 'rownames' or change",
           " 'by' to search a different column in 'x'.")
    }
    search <- rownames(x)
  } else if (length(ncol(x)) > 0) {
    if (inherits(x, "SummarizedExperiment")) {
      if (!(by %in% colnames(rowData(x)))) {
        stop("'by' is not a column in 'rowData(x)'.")
      }
      search <- rowData(x)[, by]
    } else {
      if (!(by %in% colnames(x))) {
        stop("'by' is not a column in 'x'.")
      }
      search <- x[, by]
    }
  } else {
    search <- as.character(x)
  }

  # Match each element of 'pattern' in vector 'search'
  if (!isTRUE(exactMatch)) {
    featuresIndices <- rep(NA, length(features))
    for (i in seq_along(features)) {
      g <- grep(features[i], search)
      if (length(g) == 1) {
        featuresIndices[i] <- g
      } else if (length(g) > 1) {
        warning(
          "Feature '", features[i], "' matched multiple items in '",
          by, "': ", paste(search[g], collapse = ","),
          ". Only the first match will be selected."
        )
        featuresIndices[i] <- g[1]
      }
    }
  } else {
    featuresIndices <- match(features, search)
  }

  if (sum(is.na(featuresIndices)) > 0) {
    if (sum(is.na(featuresIndices)) == length(features)) {
      if (isTRUE(exactMatch)) {
        stop(
          "None of the provided features had matching",
          " items in '", by, "' within 'x'. ",
          "Check the spelling or try setting",
          " 'exactMatch = FALSE'."
        )
      } else {
        stop(
          "None of the provided features had matching",
          " items in '", by, "' within 'x'. ",
          "Check the spelling and make sure 'by' is set",
          " to the appropriate place in 'x'."
        )
      }
    }
    warning(
      "The following features were not present in 'x': ",
      paste(features[which(is.na(featuresIndices))],
        collapse = ","
      )
    )
  }

  if (isTRUE(removeNA)) {
    featuresIndices <- featuresIndices[!is.na(featuresIndices)]
  }
  return(featuresIndices)
}
