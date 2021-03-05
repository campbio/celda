#' @title Contamination estimation with decontX
#'
#' @description Identifies contamination from factors such as ambient RNA
#' in single cell genomic datasets.
#'
#' @name decontX
#'
#' @param x A numeric matrix of counts or a \linkS4class{SingleCellExperiment}
#' with the matrix located in the assay slot under \code{assayName}.
#' Cells in each batch will be subsetted and converted to a sparse matrix
#' of class \code{dgCMatrix} from package \link{Matrix} before analysis. 
#' If background_idx is not provided, this object should only contain filtered
#' cells after cell calling. Empty cell barcodes 
#' (low expression droplets before cell calling) are not needed. If background_idx
#' is provided, this objet should contain all barcodes, including both cell and
#' empty cell barcodes.
#' @param assayName Character. Name of the assay to use if \code{x} is a
#' \linkS4class{SingleCellExperiment}.
#' @param z Numeric or character vector. Cell cluster labels. If NULL,
#' PCA will be used to reduce the dimensionality of the dataset initially,
#' '\link[uwot]{umap}' from the 'uwot' package
#' will be used to further reduce the dataset to 2 dimenions and
#' the '\link[dbscan]{dbscan}' function from the 'dbscan' package
#' will be used to identify clusters of broad cell types. Default NULL.
#' @param batch Numeric or character vector. Batch labels for cells.
#' If batch labels are supplied, DecontX is run on cells from each
#' batch separately. Cells run in different channels or assays
#' should be considered different batches. Default NULL.
#' @param background_idx Numeric or caracter vector. Index of barcodes 
#' (i.e. column id of counts matrix) that are empty droplets. 
#' @param maxIter Integer. Maximum iterations of the EM algorithm. Default 500.
#' @param convergence Numeric. The EM algorithm will be stopped if the maximum
#' difference in the contamination estimates between the previous and
#' current iterations is less than this. Default 0.001.
#' @param iterLogLik Integer. Calculate log likelihood every \code{iterLogLik}
#' iteration. Default 10.
#' @param delta Numeric Vector of length 2. Concentration parameters for
#' the Dirichlet prior for the contamination in each cell. The first element
#' is the prior for the native counts while the second element is the prior for
#' the contamination counts. These essentially act as pseudocounts for the
#' native and contamination in each cell. If \code{estimateDelta = TRUE},
#' this is only used to produce a random sample of proportions for an initial
#' value of contamination in each cell. Then
#' \code{\link[MCMCprecision]{fit_dirichlet}} is used to update
#' \code{delta} in each iteration.
#' If \code{estimateDelta = FALSE}, then \code{delta} is fixed with these
#' values for the entire inference procedure. Fixing \code{delta} and
#' setting a high number in the second element will force \code{decontX}
#' to be more aggressive and estimate higher levels of contamination at
#' the expense of potentially removing native expression.
#' Default \code{c(10, 10)}.
#' @param estimateDelta Boolean. Whether to update \code{delta} at each
#' iteration.
#' @param varGenes Integer. The number of variable genes to use in
#' dimensionality reduction before clustering. Variability is calcualted using
#' \code{\link[scran]{modelGeneVar}} function from the 'scran' package.
#' Used only when z is not provided. Default 5000.
#' @param dbscanEps Numeric. The clustering resolution parameter
#' used in '\link[dbscan]{dbscan}' to estimate broad cell clusters.
#' Used only when z is not provided. Default 1.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @param logfile Character. Messages will be redirected to a file named
#'  `logfile`. If NULL, messages will be printed to stdout.  Default NULL.
#' @param verbose Logical. Whether to print log messages. Default TRUE.
#' @param ... For the generic, further arguments to pass to each method.
#'
#' @return If \code{x} is a matrix-like object, a list will be returned
#' with the following items:
#' \describe{
#' \item{\code{decontXcounts}:}{The decontaminated matrix. Values obtained
#' from the variational inference procedure may be non-integer. However,
#' integer counts can be obtained by rounding,
#' e.g. \code{round(decontXcounts)}.}
#' \item{\code{contamination}:}{Percentage of contamination in each cell.}
#' \item{\code{estimates}:}{List of estimated parameters for each batch. If z
#' was not supplied, then the UMAP coordinates used to generated cell
#' cluster labels will also be stored here.}
#' \item{\code{z}:}{Cell population/cluster labels used for analysis.}
#' \item{\code{runParams}:}{List of arguments used in the function call.}
#' }
#'
#' If \code{x} is a \linkS4class{SingleCellExperiment}, then the decontaminated
#' counts will be stored as an assay and can be accessed with
#' \code{decontXcounts(x)}. The contamination values and cluster labels
#' will be stored in \code{colData(x)}. \code{estimates} and \code{runParams}
#' will be stored in \code{metadata(x)$decontX}. The UMAPs used to generated
#' cell cluster labels will be stored in
#' \code{reducedDims} slot in \code{x}.
#'
#' @examples
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

#' @export
#' @rdname decontX
setGeneric("decontX", function(x, ...) standardGeneric("decontX"))


#########################
# Setting up S4 methods #
#########################


#' @export
#' @rdname decontX
setMethod("decontX", "SingleCellExperiment", function(x,
                                                      assayName = "counts",
                                                      z = NULL,
                                                      batch = NULL,
                                                      background_idx = NULL,
                                                      maxIter = 500,
                                                      delta = c(10, 10),
                                                      estimateDelta = TRUE,
                                                      convergence = 0.001,
                                                      iterLogLik = 10,
                                                      varGenes = 5000,
                                                      dbscanEps = 1,
                                                      seed = 12345,
                                                      logfile = NULL,
                                                      verbose = TRUE) {
  
  counts_background <- NULL
  if (!is.null(background_idx)) {
     counts_background <- x[, background_idx]
     counts_background <- SummarizedExperiment::assay(counts_background, i = assayName)
     x <- x[, -background_idx]
  }

  mat <- SummarizedExperiment::assay(x, i = assayName)
  
  result <- .decontX(
    counts = mat,
    z = z,
    batch = batch,
    counts_background = counts_background,
    maxIter = maxIter,
    convergence = convergence,
    iterLogLik = iterLogLik,
    delta = delta,
    estimateDelta = estimateDelta,
    varGenes = varGenes,
    dbscanEps = dbscanEps,
    seed = seed,
    logfile = logfile,
    verbose = verbose
  )

  ## Add results into column annotation
  colData(x)$decontX_contamination <- result$contamination
  colData(x)$decontX_clusters <- result$z

  ## Put estimated UMAPs into SCE
  batchIndex <- unique(result$runParams$batch)
  if (length(batchIndex) > 1) {
    for (i in batchIndex) {

      ## Each individual UMAP will only be for one batch so need
      ## to put NAs in for cells in other batches
      tempUMAP <- matrix(NA, ncol = 2, nrow = ncol(mat))
      tempUMAP[result$runParams$batch == i, ] <- result$estimates[[i]]$UMAP
      colnames(tempUMAP) <- c("UMAP_1", "UMAP_2")
      rownames(tempUMAP) <- colnames(mat)

      SingleCellExperiment::reducedDim(
        x,
        paste("decontX", i, "UMAP", sep = "_")
      ) <- tempUMAP
    }
  } else {
    SingleCellExperiment::reducedDim(x, "decontX_UMAP") <-
      result$estimates[[batchIndex]]$UMAP
  }

  ## Save the rest of the result object into metadata
  decontXcounts(x) <- result$decontXcounts
  result$decontXcounts <- NULL
  S4Vectors::metadata(x)$decontX <- result

  return(x)
})

#' @export
#' @rdname decontX
setMethod("decontX", "ANY", function(x,
                                     z = NULL,
                                     batch = NULL,
                                     background_idx = NULL,
                                     maxIter = 500,
                                     delta = c(10, 10),
                                     estimateDelta = TRUE,
                                     convergence = 0.001,
                                     iterLogLik = 10,
                                     varGenes = 5000,
                                     dbscanEps = 1,
                                     seed = 12345,
                                     logfile = NULL,
                                     verbose = TRUE) {

  counts_background <- NULL
  if (!is.null(background_idx)) {
     counts_background <- x[, background_idx]
     counts_background <- SummarizedExperiment::assay(counts_background, i = assayName)
     x <- x[, -background_idx]
  }

  .decontX(
    counts = x,
    z = z,
    batch = batch,
    counts_background = counts_background,
    maxIter = maxIter,
    convergence = convergence,
    iterLogLik = iterLogLik,
    delta = delta,
    estimateDelta = estimateDelta,
    varGenes = varGenes,
    dbscanEps = dbscanEps,
    seed = seed,
    logfile = logfile,
    verbose = verbose
  )
})


## Copied from SingleCellExperiment Package

GET_FUN <- function(exprs_values, ...) {
  (exprs_values) # To ensure evaluation
  function(object, ...) {
    assay(object, i = exprs_values, ...)
  }
}

SET_FUN <- function(exprs_values, ...) {
  (exprs_values) # To ensure evaluation
  function(object, ..., value) {
    assay(object, i = exprs_values, ...) <- value
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
setMethod("decontXcounts", "SingleCellExperiment", GET_FUN("decontXcounts"))

#' @export
setReplaceMethod(
  "decontXcounts", c("SingleCellExperiment", "ANY"),
  SET_FUN("decontXcounts")
)




##########################
# Core Decontx Functions #
##########################

.decontX <- function(counts,
                     z = NULL,
                     batch = NULL,
                     counts_background = NULL,
                     maxIter = 200,
                     convergence = 0.001,
                     iterLogLik = 10,
                     delta = c(10, 10),
                     estimateDelta = TRUE,
                     varGenes = NULL,
                     dbscanEps = NULL,
                     seed = 12345,
                     logfile = NULL,
                     verbose = TRUE) {
  startTime <- Sys.time()
  .logMessages(paste(rep("-", 50), collapse = ""),
    logfile = logfile,
    append = TRUE,
    verbose = verbose
  )
  .logMessages("Starting DecontX",
    logfile = logfile,
    append = TRUE,
    verbose = verbose
  )
  .logMessages(paste(rep("-", 50), collapse = ""),
    logfile = logfile,
    append = TRUE,
    verbose = verbose
  )

  runParams <- list(
    z = z,
    batch = batch,
    maxIter = maxIter,
    delta = delta,
    estimateDelta = estimateDelta,
    convergence = convergence,
    varGenes = varGenes,
    dbscanEps = dbscanEps,
    logfile = logfile,
    verbose = verbose
  )


  totalGenes <- nrow(counts)
  totalCells <- ncol(counts)
  geneNames <- rownames(counts)
  nC <- ncol(counts)
  allCellNames <- colnames(counts)

  ## Set up final deconaminated matrix
  estRmat <- Matrix::Matrix(
    data = 0,
    ncol = totalCells,
    nrow = totalGenes,
    sparse = TRUE,
    dimnames = list(geneNames, allCellNames)
  )

  ## Generate batch labels if none were supplied
  if (is.null(batch)) {
    batch <- rep("all_cells", nC)
  }
  runParams$batch <- batch
  batchIndex <- unique(batch)

  ## Set result lists upfront for all cells from different batches
  logLikelihood <- c()
  estConp <- rep(NA, nC)
  returnZ <- rep(NA, nC)
  resBatch <- list()

  ## Cycle through each sample/batch and run DecontX
  for (bat in batchIndex) {
    if (length(batchIndex) == 1) {
      .logMessages(
        date(),
        ".. Analyzing all cells",
        logfile = logfile,
        append = TRUE,
        verbose = verbose
      )
    } else {
      .logMessages(
        date(),
        " .. Analyzing cells in batch '",
        bat, "'",
        sep = "",
        logfile = logfile,
        append = TRUE,
        verbose = verbose
      )
    }

    zBat <- NULL
    countsBat <- counts[, batch == bat]

    ## Convert to sparse matrix
    if (!inherits(countsBat, "dgCMatrix")) {
      .logMessages(
        date(),
        ".... Converting to sparse matrix",
        logfile = logfile,
        append = TRUE,
        verbose = verbose
      )
      countsBat <- methods::as(countsBat, "dgCMatrix")
    }

    if (!is.null(z)) {
      zBat <- z[batch == bat]
    }
    if (is.null(seed)) {
      res <- .decontXoneBatch(
        counts = countsBat,
        z = zBat,
        batch = bat,
        counts_background = counts_background,
        maxIter = maxIter,
        delta = delta,
        estimateDelta = estimateDelta,
        convergence = convergence,
        iterLogLik = iterLogLik,
        logfile = logfile,
        verbose = verbose,
        varGenes = varGenes,
        dbscanEps = dbscanEps,
        seed = seed
      )
    } else {
      withr::with_seed(
        seed,
        res <- .decontXoneBatch(
          counts = countsBat,
          z = zBat,
          batch = bat,
          counts_background = counts_background,
          maxIter = maxIter,
          delta = delta,
          estimateDelta = estimateDelta,
          convergence = convergence,
          iterLogLik = iterLogLik,
          logfile = logfile,
          verbose = verbose,
          varGenes = varGenes,
          dbscanEps = dbscanEps,
          seed = seed
        )
      )
    }
    estRmat.temp <- calculateNativeMatrix(
      counts = countsBat,
      theta = res$theta,
      eta = res$eta,
      phi = res$phi,
      z = as.integer(res$z),
      pseudocount = 1e-20
    )
    estRmat[seq(nrow(counts)), which(batch == bat)] <- estRmat.temp
    dimnames(estRmat) <- list(geneNames, allCellNames)

    resBatch[[bat]] <- list(
      z = res$z,
      phi = res$phi,
      eta = res$eta,
      delta = res$delta,
      theta = res$theta,
      contamination = res$contamination,
      logLikelihood = res$logLikelihood,
      UMAP = res$UMAP,
      z = res$z,
      iteration = res$iteration
    )

    estConp[batch == bat] <- res$contamination
    if (length(batchIndex) > 1) {
      returnZ[batch == bat] <- paste0(bat, "-", res$z)
    } else {
      returnZ[batch == bat] <- res$z
    }
  }
  names(resBatch) <- batchIndex

  returnResult <- list(
    "runParams" = runParams,
    "estimates" = resBatch,
    "decontXcounts" = estRmat,
    "contamination" = estConp,
    "z" = returnZ
  )

  ## Try to convert class of new matrix to class of original matrix
  if (inherits(counts, "dgCMatrix")) {
    .logMessages(
      date(),
      ".. Finalizing decontaminated matrix",
      logfile = logfile,
      append = TRUE,
      verbose = verbose
    )
  }

  if (inherits(counts, c("DelayedMatrix", "DelayedArray"))) {

    ## Determine class of seed in DelayedArray
    seed.class <- unique(DelayedArray::seedApply(counts, class))[[1]]
    if (seed.class == "HDF5ArraySeed") {
      returnResult$decontXcounts <-
        methods::as(returnResult$decontXcounts, "HDF5Matrix")
    } else {
      if (isTRUE(methods::canCoerce(returnResult$decontXcounts, seed.class))) {
        returnResult$decontXcounts <-
          methods::as(returnResult$decontXcounts, seed.class)
      }
    }
    returnResult$decontXcounts <-
      DelayedArray::DelayedArray(returnResult$decontXcounts)
  } else {
    try({
      if (methods::canCoerce(returnResult$decontXcounts, class(counts))) {
          returnResult$decontXcounts <-
            methods::as(returnResult$decontXcounts, class(counts))
        }
      },
      silent = TRUE
    )
  }

  endTime <- Sys.time()
  .logMessages(paste(rep("-", 50), collapse = ""),
    logfile = logfile,
    append = TRUE,
    verbose = verbose
  )
  .logMessages("Completed DecontX. Total time:",
    format(difftime(endTime, startTime)),
    logfile = logfile,
    append = TRUE,
    verbose = verbose
  )
  .logMessages(paste(rep("-", 50), collapse = ""),
    logfile = logfile,
    append = TRUE,
    verbose = verbose
  )

  return(returnResult)
}


# This function updates decontamination for one batch
# seed passed to this function is to be furhter passed to
# function .decontxInitializeZ()
.decontXoneBatch <- function(counts,
                             z = NULL,
                             batch = NULL,
                             counts_background = NULL,
                             maxIter = 200,
                             delta = c(10, 10),
                             estimateDelta = TRUE,
                             convergence = 0.01,
                             iterLogLik = 10,
                             logfile = NULL,
                             verbose = TRUE,
                             varGenes = NULL,
                             dbscanEps = NULL,
                             seed = 12345) {
  .checkCountsDecon(counts)
  .checkDelta(delta)

  # nG <- nrow(counts)
  nC <- ncol(counts)
  deconMethod <- "clustering"

  ## Generating UMAP and cell cluster labels if none are provided
  umap <- NULL
  if (is.null(z)) {
    m <- ".... Generating UMAP and estimating cell types"
    estimateCellTypes <- TRUE
  } else {
    m <- ".... Generating UMAP"
    estimateCellTypes <- FALSE
  }
  .logMessages(
    date(),
    m,
    logfile = logfile,
    append = TRUE,
    verbose = verbose
  )

  varGenes <- .processvarGenes(varGenes)
  dbscanEps <- .processdbscanEps(dbscanEps)

  celda.init <- .decontxInitializeZ(
    object = counts,
    varGenes = varGenes,
    dbscanEps = dbscanEps,
    estimateCellTypes = estimateCellTypes,
    seed = seed
  )
  if (is.null(z)) {
    z <- celda.init$z
  }
  umap <- celda.init$umap
  colnames(umap) <- c(
    "DecontX_UMAP_1",
    "DecontX_UMAP_2"
  )
  rownames(umap) <- colnames(counts)

  z <- .processCellLabels(z, numCells = nC)
  K <- length(unique(z))

  iter <- 1L
  numIterWithoutImprovement <- 0L
  stopIter <- 3L

  .logMessages(
    date(),
    ".... Estimating contamination",
    logfile = logfile,
    append = TRUE,
    verbose = verbose
  )

  if (deconMethod == "clustering") {
    ## Initialization
    theta <- stats::rbeta(
      n = nC,
      shape1 = delta[1],
      shape2 = delta[2]
    )

    nextDecon <- decontXInitialize(
      counts = counts,
      theta = theta,
      z = z,
      pseudocount = 1e-20
    )
    phi <- nextDecon$phi
    eta <- nextDecon$eta

    # if counts_background is not null, use empirical dist. to replace eta
    if (!is.null(counts_background)) {
      # Add pseudocount to each gene in eta
       eta_tilda <- Matrix::rowSums(counts_background) + 1e-20
       eta <- eta_tilda/sum(eta_tilda)
       
       # Make eta a matrix same dimension as phi
       eta <- matrix(eta, length(eta), dim(phi)[2])
    }
    

    ll <- c()
    llRound <- decontXLogLik(
      counts = counts,
      z = z,
      phi = phi,
      eta = eta,
      theta = theta,
      pseudocount = 1e-20
    )

    ## EM updates
    theta.previous <- theta
    converged <- FALSE
    counts.colsums <- Matrix::colSums(counts)
    while (iter <= maxIter & !isTRUE(converged) &
      numIterWithoutImprovement <= stopIter) {
        if (is.null(counts_background)) {
          nextDecon <- decontXEM(
            counts = counts,
            counts_colsums = counts.colsums,
            phi = phi,
            estimate_eta = TRUE,
            eta = eta,
            theta = theta,
            z = z,
            estimate_delta = isTRUE(estimateDelta),
            delta = delta,
            pseudocount = 1e-20
          )
        } else {
           nextDecon <- decontXEM(
            counts = counts,
            counts_colsums = counts.colsums,
            phi = phi,
            estimate_eta = FALSE,
            eta = eta,
            theta = theta,
            z = z,
            estimate_delta = isTRUE(estimateDelta),
            delta = delta,
            pseudocount = 1e-20
          )
        }
      

      theta <- nextDecon$theta
      phi <- nextDecon$phi
      eta <- nextDecon$eta
      delta <- nextDecon$delta

      max.divergence <- max(abs(theta.previous - theta))
      if (max.divergence < convergence) {
        converged <- TRUE
      }
      theta.previous <- theta

      ## Calculate likelihood and check for convergence
      if (iter %% iterLogLik == 0 || converged) {
        llTemp <- decontXLogLik(
          counts = counts,
          z = z,
          phi = phi,
          eta = eta,
          theta = theta,
          pseudocount = 1e-20
        )

        ll <- c(ll, llTemp)

        .logMessages(date(),
          "...... Completed iteration:",
          iter,
          "| converge:",
          signif(max.divergence, 4),
          logfile = logfile,
          append = TRUE,
          verbose = verbose
        )
      }

      iter <- iter + 1L
    }
  }

  resConp <- nextDecon$contamination
  names(resConp) <- colnames(counts)

  return(list(
    "logLikelihood" = ll,
    "contamination" = resConp,
    "theta" = theta,
    "delta" = delta,
    "phi" = phi,
    "eta" = eta,
    "UMAP" = umap,
    "iteration" = iter - 1L,
    "z" = z
  ))
}





# This function calculates the log-likelihood
#
# counts Numeric/Integer matrix. Observed count matrix, rows represent features
# and columns represent cells
# z Integer vector. Cell population labels
# phi Numeric matrix. Rows represent features and columns represent cell
# populations
# eta Numeric matrix. Rows represent features and columns represent cell
# populations
# theta Numeric vector. Proportion of truely expressed transcripts
.deconCalcLL <- function(counts, z, phi, eta, theta) {
  # ll = sum( t(counts) * log( (1-conP )*geneDist[z,] + conP * conDist[z, ] +
  # 1e-20 ) )  # when dist_mat are K x G matrices
  ll <- sum(Matrix::t(counts) * log(theta * t(phi)[z, ] +
    (1 - theta) * t(eta)[z, ] + 1e-20))
  return(ll)
}

# DEPRECATED. This is not used, but is kept as it might be useful in the future
# This function calculates the log-likelihood of background distribution
# decontamination
# bgDist Numeric matrix. Rows represent feature and columns are the times that
# the background-distribution has been replicated.
.bgCalcLL <- function(counts, globalZ, cbZ, phi, eta, theta) {
  # ll <- sum(t(counts) * log(theta * t(cellDist) +
  #        (1 - theta) * t(bgDist) + 1e-20))
  ll <- sum(t(counts) * log(theta * t(phi)[cbZ, ] +
    (1 - theta) * t(eta)[globalZ, ] + 1e-20))
  return(ll)
}


# This function updates decontamination
#  phi Numeric matrix. Rows represent features and columns represent cell
# populations
#  eta Numeric matrix. Rows represent features and columns represent cell
# populations
#  theta Numeric vector. Proportion of truely expressed transctripts
#' @importFrom MCMCprecision fit_dirichlet
.cDCalcEMDecontamination <- function(counts,
                                     phi,
                                     eta,
                                     theta,
                                     z,
                                     K,
                                     delta) {
  ## Notes: use fix-point iteration to update prior for theta, no need
  ## to feed delta anymore

  logPr <- log(t(phi)[z, ] + 1e-20) + log(theta + 1e-20)
  logPc <- log(t(eta)[z, ] + 1e-20) + log(1 - theta + 1e-20)
  Pr.e <- exp(logPr)
  Pc.e <- exp(logPc)
  Pr <- Pr.e / (Pr.e + Pc.e)

  estRmat <- t(Pr) * counts
  rnGByK <- .colSumByGroupNumeric(estRmat, z, K)
  cnGByK <- rowSums(rnGByK) - rnGByK

  counts.cs <- colSums(counts)
  estRmat.cs <- colSums(estRmat)
  estRmat.cs.n <- estRmat.cs / counts.cs
  estCmat.cs.n <- 1 - estRmat.cs.n
  temp <- cbind(estRmat.cs.n, estCmat.cs.n)
  deltaV2 <- MCMCprecision::fit_dirichlet(temp)$alpha

  ## Update parameters
  theta <-
    (estRmat.cs + deltaV2[1]) / (counts.cs + sum(deltaV2))
  phi <- normalizeCounts(rnGByK,
    normalize = "proportion",
    pseudocountNormalize = 1e-20
  )
  eta <- normalizeCounts(cnGByK,
    normalize = "proportion",
    pseudocountNormalize = 1e-20
  )

  return(list(
    "estRmat" = estRmat,
    "theta" = theta,
    "phi" = phi,
    "eta" = eta,
    "delta" = deltaV2
  ))
}

# DEPRECATED. This is not used, but is kept as it might be useful in the
# feature.
# This function updates decontamination using background distribution
.cDCalcEMbgDecontamination <-
  function(counts, globalZ, cbZ, trZ, phi, eta, theta) {
    logPr <- log(t(phi)[cbZ, ] + 1e-20) + log(theta + 1e-20)
    logPc <-
      log(t(eta)[globalZ, ] + 1e-20) + log(1 - theta + 1e-20)

    Pr <- exp(logPr) / (exp(logPr) + exp(logPc))
    Pc <- 1 - Pr
    deltaV2 <-
      MCMCprecision::fit_dirichlet(matrix(c(Pr, Pc), ncol = 2))$alpha

    estRmat <- t(Pr) * counts
    phiUnnormalized <-
      .colSumByGroupNumeric(estRmat, cbZ, max(cbZ))
    etaUnnormalized <-
      rowSums(phiUnnormalized) - .colSumByGroupNumeric(
        phiUnnormalized,
        trZ, max(trZ)
      )

    ## Update paramters
    theta <-
      (colSums(estRmat) + deltaV2[1]) / (colSums(counts) + sum(deltaV2))
    phi <-
      normalizeCounts(phiUnnormalized,
        normalize = "proportion",
        pseudocountNormalize = 1e-20
      )
    eta <-
      normalizeCounts(etaUnnormalized,
        normalize = "proportion",
        pseudocountNormalize = 1e-20
      )

    return(list(
      "estRmat" = estRmat,
      "theta" = theta,
      "phi" = phi,
      "eta" = eta,
      "delta" = deltaV2
    ))
  }

## Make sure provided count matrix is the right type
.checkCountsDecon <- function(counts) {
  if (sum(is.na(counts)) > 0) {
    stop("Missing value in 'counts' matrix.")
  }
  if (is.null(dim(counts))) {
    stop("At least 2 genes need to have non-zero expressions.")
  }
}


## Make sure provided cell labels are the right type
#' @importFrom plyr mapvalues
.processCellLabels <- function(z, numCells) {
  if (length(z) != numCells) {
    stop(
      "'z' must be of the same length as the number of cells in the",
      " 'counts' matrix."
    )
  }
  if (length(unique(z)) < 2) {
    stop(
      "No need to decontaminate when only one cluster",
      " is in the dataset."
    ) # Even though
    # everything runs smoothly when length(unique(z)) == 1, result is not
    # trustful
  }
  if (!is.factor(z)) {
    z <- plyr::mapvalues(z, unique(z), seq(length(unique(z))))
    z <- as.factor(z)
  }
  return(z)
}


## Add two (veried-length) vectors of logLikelihood
addLogLikelihood <- function(llA, llB) {
  lengthA <- length(llA)
  lengthB <- length(llB)

  if (lengthA >= lengthB) {
    llB <- c(llB, rep(llB[lengthB], lengthA - lengthB))
    ll <- llA + llB
  } else {
    llA <- c(llA, rep(llA[lengthA], lengthB - lengthA))
    ll <- llA + llB
  }

  return(ll)
}

.decontxInitializeZ <- function(object,
                                varGenes = 2000,
                                dbscanEps = 1,
                                estimateCellTypes = TRUE,
                                seed = 12345) {
  if (!is(object, "SingleCellExperiment")) {
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = object)
    )
  }
  sce <- scater::logNormCounts(sce, log = TRUE)

  if (nrow(sce) <= varGenes) {
    topVariableGenes <- seq_len(nrow(sce))
  } else if (nrow(sce) > varGenes) {
    sce.var <- scran::modelGeneVar(sce)
    topVariableGenes <- order(sce.var$bio,
      decreasing = TRUE
    )[seq(varGenes)]
  }
  sce <- sce[topVariableGenes, ]

  if (!is.null(seed)) {
    with_seed(
      seed,
      resUmap <- scater::calculateUMAP(sce, n_threads = 1)
    )
  } else {
    resUmap <- scater::calculateUMAP(sce, n_threads = 1)
  }

  z <- NULL
  if (isTRUE(estimateCellTypes)) {
    # Find clusters with dbSCAN
    totalClusters <- 1
    iter <- 1
    while (totalClusters <= 1 & dbscanEps > 0 & iter < 10) {
      resDbscan <- dbscan::dbscan(resUmap, dbscanEps)
      dbscanEps <- dbscanEps - (0.25 * dbscanEps)
      totalClusters <- length(unique(resDbscan$cluster))
      iter <- iter + 1
    }

    # If dbscan was not able to get more than 2 clusters,
    # use kmeans to force 2 clusters as a last resort
    if (totalClusters == 1) {
      cl <- stats::kmeans(t(SingleCellExperiment::logcounts(sce)), 2)
      z <- cl$cluster
    } else {
      z <- resDbscan$cluster
    }
  }

  return(list(
    "z" = z,
    "umap" = resUmap
  ))
}


## Initialization of cell labels for DecontX when they are not given
.decontxInitializeZ_prevous <-
  function(object, # object is either a sce object or a count matrix
           varGenes = 5000,
           dbscanEps = 1.0,
           verbose = TRUE,
           seed = 12345,
           logfile = NULL) {
    if (!is(object, "SingleCellExperiment")) {
      sce <- SingleCellExperiment::SingleCellExperiment(
        assays =
          list(counts = object)
      )
    }

    sce <- sce[Matrix::rowSums(SingleCellExperiment::counts(sce)) > 0, ]
    sce <- scater::logNormCounts(sce, log = TRUE)
    # sce <- scater::normalize(sce)


    if (nrow(sce) <= varGenes) {
      topVariableGenes <- seq_len(nrow(sce))
    } else if (nrow(sce) > varGenes) {
      sce.var <- scran::modelGeneVar(sce)
      topVariableGenes <- order(sce.var$bio,
        decreasing = TRUE
      )[seq(varGenes)]
    }
    countsFiltered <- as.matrix(SingleCellExperiment::counts(
      sce[topVariableGenes, ]
    ))
    storage.mode(countsFiltered) <- "integer"

    .logMessages(
      date(),
      "...... Collapsing features into",
      L,
      "modules",
      logfile = logfile,
      append = TRUE,
      verbose = verbose
    )
    ## Celda clustering using recursive module splitting
    L <- min(L, nrow(countsFiltered))
    if (is.null(seed)) {
      initialModuleSplit <- recursiveSplitModule(countsFiltered,
        initialL = L, maxL = L, perplexity = FALSE, verbose = FALSE
      )
    } else {
      with_seed(seed,
                initialModuleSplit <- recursiveSplitModule(countsFiltered,
        initialL = L, maxL = L, perplexity = FALSE, verbose = FALSE
      ))
    }
    initialModel <- subsetCeldaList(initialModuleSplit, list(L = L))

    .logMessages(
      date(),
      "...... Reducing dimensionality with UMAP",
      logfile = logfile,
      append = TRUE,
      verbose = verbose
    )
    ## Louvan graph-based method to reduce dimension into 2 cluster
    nNeighbors <- min(15, ncol(countsFiltered))
    # resUmap <- uwot::umap(t(sqrt(fm)), n_neighbors = nNeighbors,
    #    min_dist = 0.01, spread = 1)
    # rm(fm)
    resUmap <- celdaUmap(countsFiltered, initialModel,
      minDist = 0.01, spread = 1, nNeighbors = nNeighbors, seed = seed
    )

    .logMessages(
      date(),
      " ...... Determining cell clusters with DBSCAN (Eps=",
      dbscanEps,
      ")",
      sep = "",
      logfile = logfile,
      append = TRUE,
      verbose = verbose
    )
    # Use dbSCAN on the UMAP to identify broad cell types
    totalClusters <- 1
    while (totalClusters <= 1 & dbscanEps > 0) {
      resDbscan <- dbscan::dbscan(resUmap, dbscanEps)
      dbscanEps <- dbscanEps - (0.25 * dbscanEps)
      totalClusters <- length(unique(resDbscan$cluster))
    }

    return(list(
      "z" = resDbscan$cluster,
      "umap" = resUmap
    ))
  }


## process varGenes
.processvarGenes <- function(varGenes) {
  if (is.null(varGenes)) {
    varGenes <- 5000
  } else {
    if (varGenes < 2 | length(varGenes) > 1) {
      stop("Parameter 'varGenes' must be an integer larger than 1.")
    }
  }
  return(varGenes)
}

## process dbscanEps for resolusion threshold using DBSCAN
.processdbscanEps <- function(dbscanEps) {
  if (is.null(dbscanEps)) {
    dbscanEps <- 1
  } else {
    if (dbscanEps < 0) {
      stop("Parameter 'dbscanEps' needs to be non-negative.")
    }
  }
  return(dbscanEps)
}

.checkDelta <- function(delta) {
  if (!is.numeric(delta) | length(delta) != 2 | any(delta < 0)) {
    stop("'delta' needs to be a numeric vector of length 2",
         " containing positive values.")
  }
  return(delta)
}




#########################
# Simulating Data       #
#########################

#' @title Simulate contaminated count matrix
#' @description This function generates a list containing two count matrices --
#'  one for real expression, the other one for contamination, as well as other
#'  parameters used in the simulation which can be useful for running
#'  decontamination.
#' @param C Integer. Number of cells to be simulated. Default \code{300}.
#' @param G Integer. Number of genes to be simulated. Default \code{100}.
#' @param K Integer. Number of cell populations to be simulated.
#' Default \code{3}.
#' @param NRange Integer vector. A vector of length 2 that specifies the lower
#'  and upper bounds of the number of counts generated for each cell. Default
#'  \code{c(500, 1000)}.
#' @param beta Numeric. Concentration parameter for Phi. Default \code{0.1}.
#' @param delta Numeric or Numeric vector. Concentration parameter for Theta.
#'  If input as a single numeric value, symmetric values for beta
#'  distribution are specified; if input as a vector of lenght 2, the two
#'  values will be the shape1 and shape2 paramters of the beta distribution
#'  respectively. Default \code{c(1, 5)}.
#' @param numMarkers Integer. Number of markers for each cell population.
#' Default \code{3}.
#' @param seed Integer. Passed to \code{\link[withr]{with_seed}}.
#' For reproducibility, a default value of 12345 is used. If NULL, no calls to
#'  \code{\link[withr]{with_seed}} are made.
#' @return A list containing the \code{nativeMatirx} (real expression),
#' \code{observedMatrix} (real expression + contamination), as well as other
#' parameters used in the simulation.
#' @author Shiyi Yang, Joshua Campbell
#' @examples
#' contaminationSim <- simulateContamination(K = 3, delta = c(1, 10))
#' @export
simulateContamination <- function(C = 300,
                                  G = 100,
                                  K = 3,
                                  NRange = c(500, 1000),
                                  beta = 0.1,
                                  delta = c(1, 10),
                                  numMarkers = 3,
                                  seed = 12345) {
  if (is.null(seed)) {
    res <- .simulateContaminatedMatrix(
      C = C,
      G = G,
      K = K,
      NRange = NRange,
      beta = beta,
      delta = delta,
      numMarkers = numMarkers
    )
  } else {
    with_seed(
      seed,
      res <- .simulateContaminatedMatrix(
        C = C,
        G = G,
        K = K,
        NRange = NRange,
        beta = beta,
        delta = delta,
        numMarkers = numMarkers
      )
    )
  }

  return(res)
}


.simulateContaminatedMatrix <- function(C = 300,
                                        G = 100,
                                        K = 3,
                                        NRange = c(500, 1000),
                                        beta = 0.5,
                                        delta = c(1, 2),
                                        numMarkers = 3) {
  if (length(delta) == 1) {
    cpByC <- stats::rbeta(
      n = C,
      shape1 = delta,
      shape2 = delta
    )
  } else {
    cpByC <- stats::rbeta(
      n = C,
      shape1 = delta[1],
      shape2 = delta[2]
    )
  }

  z <- sample(seq(K), size = C, replace = TRUE)
  if (length(unique(z)) < K) {
    warning(
      "Only ",
      length(unique(z)),
      " clusters are simulated. Try to increase numebr of cells 'C' if",
      " more clusters are needed"
    )
    K <- length(unique(z))
    z <- plyr::mapvalues(z, unique(z), seq(length(unique(z))))
  }

  NbyC <- sample(seq(min(NRange), max(NRange)),
    size = C,
    replace = TRUE
  )
  cNbyC <- vapply(seq(C), function(i) {
    stats::rbinom(
      n = 1,
      size = NbyC[i],
      p = cpByC[i]
    )
  }, integer(1))
  rNbyC <- NbyC - cNbyC

  phi <- .rdirichlet(K, rep(beta, G))

  ## Select random genes to be markers in each cell population
  ## by setting their values to zero.
  if (K * numMarkers > G) {
    stop("The number of markers ('numMarkers') times the number of cell",
         " populations ('K') cannot be greater than the number of",
         " genes ('G').")
  }
  markerKIndex <- rep(seq(K), each = numMarkers)
  markerRowIndex <- sample(seq(G), numMarkers * K)
  for (i in seq(K)) {
    ix <- markerRowIndex[markerKIndex == i]
    phi[i, ix] <- max(phi[i, ])
    for (j in setdiff(seq(K), i)) {
      phi[j, ix] <- 0
    }
  }
  phi <- prop.table(phi, margin = 1)

  ## sample real expressed count matrix
  cellRmat <- vapply(seq(C), function(i) {
    stats::rmultinom(1, size = rNbyC[i], prob = phi[z[i], ])
  }, integer(G))

  rownames(cellRmat) <- paste0("Gene_", seq(G))
  colnames(cellRmat) <- paste0("Cell_", seq(C))

  ## Get list of marker names
  markerNames <- list()
  for (i in seq(K)) {
    markerNames[[i]] <- rownames(cellRmat)[markerRowIndex[markerKIndex == i]]
  }
  names(markerNames) <- paste0("CellType_", seq(K), "_Markers")

  ## sample contamination count matrix
  nGByK <-
    rowSums(cellRmat) - .colSumByGroup(cellRmat, group = z, K = K)
  eta <- normalizeCounts(counts = nGByK, normalize = "proportion")

  cellCmat <- vapply(seq(C), function(i) {
    stats::rmultinom(1, size = cNbyC[i], prob = eta[, z[i]])
  }, integer(G))
  cellOmat <- cellRmat + cellCmat
  contamination <- colSums(cellCmat) / colSums(cellOmat)

  rownames(cellOmat) <- paste0("Gene_", seq(G))
  colnames(cellOmat) <- paste0("Cell_", seq(C))

  return(
    list(
      "nativeCounts" = cellRmat,
      "observedCounts" = cellOmat,
      "NByC" = NbyC,
      "z" = z,
      "eta" = eta,
      "phi" = t(phi),
      "markers" = markerNames,
      "numMarkers" = numMarkers,
      "contamination" = contamination
    )
  )
}
