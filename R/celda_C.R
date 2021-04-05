#' @title Cell clustering with Celda
#' @description Clusters the columns of a count matrix containing single-cell
#'  data into K subpopulations. The
#'  \code{useAssay} \link{assay} slot in
#'  \code{altExpName} \link{altExp} slot will be used if
#'  it exists. Otherwise, the \code{useAssay}
#'  \link{assay} slot in \code{x} will be used if
#'  \code{x} is a \linkS4class{SingleCellExperiment} object.
#' @param x A numeric \link{matrix} of counts or a
#'  \linkS4class{SingleCellExperiment}
#'  with the matrix located in the assay slot under \code{useAssay} in
#'  \code{altExp(x, altExpName)}.
#'  Rows represent features and columns represent cells.
#' @param useAssay A string specifying the name of the
#'  \link{assay} slot to use. Default "counts".
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default "featureSubset".
#' @param sampleLabel Vector or factor. Denotes the sample label for each cell
#'  (column) in the count matrix.
#' @param K Integer. Number of cell populations.
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount
#'  to each cell population in each sample. Default 1.
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to
#'  each feature in each cell population. Default 1.
#' @param algorithm String. Algorithm to use for clustering cell subpopulations.
#'  One of 'EM' or 'Gibbs'. The EM algorithm is faster, especially for larger
#'  numbers of cells. However, more chains may be required to ensure a good
#'  solution is found. If 'EM' is selected, then 'stopIter' will be
#'  automatically set to 1. Default 'EM'.
#' @param stopIter Integer. Number of iterations without improvement in the
#'  log likelihood to stop inference. Default 10.
#' @param maxIter Integer. Maximum number of iterations of Gibbs sampling or
#'  EM to perform. Default 200.
#' @param splitOnIter Integer. On every `splitOnIter` iteration, a heuristic
#'  will be applied to determine if a cell population should be reassigned and
#'  another cell population should be split into two clusters. To disable
#'  splitting, set to -1. Default 10.
#' @param splitOnLast Integer. After `stopIter` iterations have been
#'  performed without improvement, a heuristic will be applied to determine if
#'  a cell population should be reassigned and another cell population should be
#'  split into two clusters. If a split occurs, then `stopIter` will be reset.
#'  Default TRUE.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @param nchains Integer. Number of random cluster initializations. Default 3.
#' @param zInitialize Character. One of 'random', 'split', or 'predefined'.
#'  With 'random', cells are randomly assigned to a populations. With 'split',
#'  cells will be split into sqrt(K) populations and then each population will
#'  be subsequently split into another sqrt(K) populations. With 'predefined',
#'  values in `zInit` will be used to initialize `z`. Default 'split'.
#' @param zInit Integer vector. Sets initial starting values of z. 'zInit'
#'  is only used when `zInitialize = 'predfined'`. Default NULL.
#' @param countChecksum Character. An MD5 checksum for the `counts` matrix.
#'  Default NULL.
#' @param logfile Character. Messages will be redirected to a file named
#'  `logfile`. If NULL, messages will be printed to stdout.  Default NULL.
#' @param verbose Logical. Whether to print log messages. Default TRUE.
#' @param ... Ignored. Placeholder to prevent check warning.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object. Function
#'  parameter settings are stored in the \link{metadata}
#'  \code{"celda_parameters"} slot.
#'  Columns \code{celda_sample_label} and \code{celda_cell_cluster} in
#'  \link{colData} contain sample labels and celda cell
#'  population clusters.
#' @seealso \link{celda_G} for feature clustering and \link{celda_CG} for
#'  simultaneous clustering of features and cells. \link{celdaGridSearch} can
#'  be used to run multiple values of K and multiple chains in parallel.
#' @examples
#' data(celdaCSim)
#' sce <- celda_C(celdaCSim$counts,
#'     K = celdaCSim$K,
#'     sampleLabel = celdaCSim$sampleLabel,
#'     nchains = 1)
#' @import Rcpp RcppEigen
#' @importFrom withr with_seed
#' @export
setGeneric("celda_C", function(x, ...) {
    standardGeneric("celda_C")})


#' @rdname celda_C
#' @export
setMethod("celda_C",
    signature(x = "SingleCellExperiment"),
    function(x,
        useAssay = "counts",
        altExpName = "featureSubset",
        sampleLabel = NULL,
        K,
        alpha = 1,
        beta = 1,
        algorithm = c("EM", "Gibbs"),
        stopIter = 10,
        maxIter = 200,
        splitOnIter = 10,
        splitOnLast = TRUE,
        seed = 12345,
        nchains = 3,
        zInitialize = c("split", "random", "predefined"),
        countChecksum = NULL,
        zInit = NULL,
        logfile = NULL,
        verbose = TRUE) {

        xClass <- "SingleCellExperiment"

        if (!altExpName %in% SingleCellExperiment::altExpNames(x)) {
            stop(altExpName, " not in 'altExpNames(x)'. Run ",
                "selectFeatures(x) first!")
        }

        altExp <- SingleCellExperiment::altExp(x, altExpName)

        if (!useAssay %in% SummarizedExperiment::assayNames(altExp)) {
            stop(useAssay, " not in assayNames(altExp(x, altExpName))")
        }

        counts <- SummarizedExperiment::assay(altExp, i = useAssay)

        altExp <- .celdaCWithSeed(counts = counts,
            xClass = xClass,
            useAssay = useAssay,
            sce = altExp,
            sampleLabel = sampleLabel,
            K = K,
            alpha = alpha,
            beta = beta,
            algorithm = match.arg(algorithm),
            stopIter = stopIter,
            maxIter = maxIter,
            splitOnIter = splitOnIter,
            splitOnLast = splitOnLast,
            seed = seed,
            nchains = nchains,
            zInitialize = match.arg(zInitialize),
            countChecksum = countChecksum,
            zInit = zInit,
            logfile = logfile,
            verbose = verbose)
        SingleCellExperiment::altExp(x, altExpName) <- altExp
        return(x)
    }
)


#' @rdname celda_C
#' @export
setMethod("celda_C",
    signature(x = "matrix"),
    function(x,
        useAssay = "counts",
        altExpName = "featureSubset",
        sampleLabel = NULL,
        K,
        alpha = 1,
        beta = 1,
        algorithm = c("EM", "Gibbs"),
        stopIter = 10,
        maxIter = 200,
        splitOnIter = 10,
        splitOnLast = TRUE,
        seed = 12345,
        nchains = 3,
        zInitialize = c("split", "random", "predefined"),
        countChecksum = NULL,
        zInit = NULL,
        logfile = NULL,
        verbose = TRUE) {

        ls <- list()
        ls[[useAssay]] <- x
        sce <- SingleCellExperiment::SingleCellExperiment(assays = ls)
        SingleCellExperiment::altExp(sce, altExpName) <- sce
        xClass <- "matrix"

        altExp <- .celdaCWithSeed(counts = x,
            xClass = xClass,
            useAssay = useAssay,
            sce = SingleCellExperiment::altExp(sce, altExpName),
            sampleLabel = sampleLabel,
            K = K,
            alpha = alpha,
            beta = beta,
            algorithm = match.arg(algorithm),
            stopIter = stopIter,
            maxIter = maxIter,
            splitOnIter = splitOnIter,
            splitOnLast = splitOnLast,
            seed = seed,
            nchains = nchains,
            zInitialize = match.arg(zInitialize),
            countChecksum = countChecksum,
            zInit = zInit,
            logfile = logfile,
            verbose = verbose)
        SingleCellExperiment::altExp(sce, altExpName) <- altExp
        return(sce)
    }
)


.celdaCWithSeed <- function(counts,
    xClass,
    useAssay,
    sce,
    sampleLabel,
    K,
    alpha,
    beta,
    algorithm,
    stopIter,
    maxIter,
    splitOnIter,
    splitOnLast,
    seed,
    nchains,
    zInitialize,
    countChecksum,
    zInit,
    logfile,
    verbose) {

    .validateCounts(counts)

    if (is.null(seed)) {
        celdaCMod <- .celda_C(counts = counts,
            sampleLabel = sampleLabel,
            K = K,
            alpha = alpha,
            beta = beta,
            algorithm = algorithm,
            stopIter = stopIter,
            maxIter = maxIter,
            splitOnIter = splitOnIter,
            splitOnLast = splitOnLast,
            nchains = nchains,
            zInitialize = zInitialize,
            countChecksum = countChecksum,
            zInit = zInit,
            logfile = logfile,
            verbose = verbose,
            reorder = TRUE)
    } else {
        with_seed(seed,
            celdaCMod <- .celda_C(counts = counts,
                sampleLabel = sampleLabel,
                K = K,
                alpha = alpha,
                beta = beta,
                algorithm = algorithm,
                stopIter = stopIter,
                maxIter = maxIter,
                splitOnIter = splitOnIter,
                splitOnLast = splitOnLast,
                nchains = nchains,
                zInitialize = zInitialize,
                countChecksum = countChecksum,
                zInit = zInit,
                logfile = logfile,
                verbose = verbose,
                reorder = TRUE))
    }

    sce <- .createSCEceldaC(celdaCMod = celdaCMod,
        sce = sce,
        xClass = xClass,
        useAssay = useAssay,
        algorithm = algorithm,
        stopIter = stopIter,
        maxIter = maxIter,
        splitOnIter = splitOnIter,
        splitOnLast = splitOnLast,
        nchains = nchains,
        zInitialize = zInitialize,
        zInit = zInit,
        logfile = logfile,
        verbose = verbose)
    return(sce)
}


# celda_C main function
.celda_C <- function(counts,
    sampleLabel = NULL,
    K,
    alpha = 1,
    beta = 1,
    algorithm = c("EM", "Gibbs"),
    stopIter = 10,
    maxIter = 200,
    splitOnIter = 10,
    splitOnLast = TRUE,
    nchains = 3,
    zInitialize = c("split", "random", "predefined"),
    countChecksum = NULL,
    zInit = NULL,
    logfile = NULL,
    verbose = TRUE,
    reorder = TRUE) {

    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = FALSE,
        verbose = verbose)

    .logMessages("Starting Celda_C: Clustering cells.",
        logfile = logfile,
        append = TRUE,
        verbose = verbose)

    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)

    startTime <- Sys.time()

    ## Error checking and variable processing
    #counts <- .processCounts(counts)
    if (is.null(countChecksum)) {
        countChecksum <- .createCountChecksum(counts)
    }

    sampleLabel <- .processSampleLabels(sampleLabel, ncol(counts))
    s <- as.integer(sampleLabel)

    algorithm <- match.arg(algorithm)
    if (algorithm == "EM") {
      stopIter <- 1
    }

    algorithmFun <- ifelse(algorithm == "Gibbs",
      ".cCCalcGibbsProbZ",
      ".cCCalcEMProbZ"
    )
    zInitialize <- match.arg(zInitialize)

    allChains <- seq(nchains)

    bestResult <- NULL
    for (i in allChains) {
      ## Initialize cluster labels
      .logMessages(date(),
        ".. Initializing 'z' in chain",
        i,
        "with",
        paste0("'", zInitialize, "' "),
        logfile = logfile,
        append = TRUE,
        verbose = verbose
      )

      if (zInitialize == "predefined") {
        if (is.null(zInit)) {
          stop("'zInit' needs to specified when initilize.z == 'given'.")
        }

      z <- .initializeCluster(K,
        ncol(counts),
        initial = zInit,
        fixed = NULL
      )
    } else if (zInitialize == "split") {
      z <- .initializeSplitZ(counts,
        K = K,
        alpha = alpha,
        beta = beta
      )
    } else {
      z <- .initializeCluster(K,
        ncol(counts),
        initial = NULL,
        fixed = NULL
      )
    }

    zBest <- z

    ## Calculate counts one time up front
    p <- .cCDecomposeCounts(counts, s, z, K)
    nS <- p$nS
    nG <- p$nG
    nM <- p$nM
    mCPByS <- p$mCPByS
    nGByCP <- p$nGByCP
    nCP <- p$nCP
    nByC <- p$nByC

    ll <- .cCCalcLL(
      mCPByS = mCPByS,
      nGByCP = nGByCP,
      s = s,
      K = K,
      nS = nS,
      nG = nG,
      alpha = alpha,
      beta = beta
    )

    iter <- 1L
    numIterWithoutImprovement <- 0L
    doCellSplit <- TRUE
    while (iter <= maxIter & numIterWithoutImprovement <= stopIter) {
      nextZ <- do.call(algorithmFun, list(
        counts = counts,
        mCPByS = mCPByS,
        nGByCP = nGByCP,
        nByC = nByC,
        nCP = nCP,
        z = z,
        s = s,
        K = K,
        nG = nG,
        nM = nM,
        alpha = alpha,
        beta = beta
      ))

      mCPByS <- nextZ$mCPByS
      nGByCP <- nextZ$nGByCP
      nCP <- nextZ$nCP
      z <- nextZ$z

      ## Perform split on i-th iteration of no improvement in log
      ## likelihood
      tempLl <- .cCCalcLL(
        mCPByS = mCPByS,
        nGByCP = nGByCP,
        s = s,
        K = K,
        nS = nS,
        nG = nG,
        alpha = alpha,
        beta = beta
      )

      if (K > 2 & iter != maxIter &
        ((((numIterWithoutImprovement == stopIter &
          !all(tempLl >= ll))) & isTRUE(splitOnLast)) |
          (splitOnIter > 0 & iter %% splitOnIter == 0 &
            isTRUE(doCellSplit)))) {
        .logMessages(date(),
          " .... Determining if any cell clusters should be split.",
          logfile = logfile,
          append = TRUE,
          sep = "",
          verbose = verbose
        )

        res <- .cCSplitZ(
          counts,
          mCPByS,
          nGByCP,
          nCP,
          s,
          z,
          K,
          nS,
          nG,
          alpha,
          beta,
          zProb = t(as.matrix(nextZ$probs)),
          maxClustersToTry = K,
          minCell = 3
        )

        .logMessages(res$message,
          logfile = logfile,
          append = TRUE,
          verbose = verbose
        )

        # Reset convergence counter if a split occured
        if (!isTRUE(all.equal(z, res$z))) {
          numIterWithoutImprovement <- 0L
          doCellSplit <- TRUE
        } else {
          doCellSplit <- FALSE
        }

        ## Re-calculate variables
        z <- res$z
        mCPByS <- res$mCPByS
        nGByCP <- res$nGByCP
        nCP <- res$nCP
      }

      ## Calculate complete likelihood
      tempLl <- .cCCalcLL(
        mCPByS = mCPByS,
        nGByCP = nGByCP,
        s = s,
        K = K,
        nS = nS,
        nG = nG,
        alpha = alpha,
        beta = beta
      )

      if ((all(tempLl > ll)) | iter == 1) {
        zBest <- z
        llBest <- tempLl
        numIterWithoutImprovement <- 1L
      } else {
        numIterWithoutImprovement <- numIterWithoutImprovement + 1L
      }

      ll <- c(ll, tempLl)

      .logMessages(date(),
        ".... Completed iteration:",
        iter,
        "| logLik:",
        tempLl,
        logfile = logfile,
        append = TRUE,
        verbose = verbose
      )
      iter <- iter + 1
    }

    names <- list(
      row = rownames(counts),
      column = colnames(counts),
      sample = levels(sampleLabel)
    )

    result <- list(
      z = zBest,
      completeLogLik = ll,
      finalLogLik = llBest,
      K = K,
      sampleLabel = sampleLabel,
      alpha = alpha,
      beta = beta,
      countChecksum = countChecksum,
      names = names
    )

    if (is.null(bestResult) ||
      result$finalLogLik > bestResult$finalLogLik) {
      bestResult <- result
    }

    .logMessages(date(),
      ".. Finished chain",
      i,
      logfile = logfile,
      append = TRUE,
      verbose = verbose
    )
  }

  bestResult <- methods::new("celda_C",
    clusters = list(z = bestResult$z),
    params = list(
      K = as.integer(bestResult$K),
      alpha = bestResult$alpha,
      beta = bestResult$beta,
      countChecksum = bestResult$countChecksum
    ),
    sampleLabel = bestResult$sampleLabel,
    completeLogLik = bestResult$completeLogLik,
    finalLogLik = bestResult$finalLogLik,
    names = bestResult$names
  )

  if (isTRUE(reorder)) {
    bestResult <- .reorderCeldaC(counts = counts, res = bestResult)
  }

  endTime <- Sys.time()
  .logMessages(paste(rep("-", 50), collapse = ""),
    logfile = logfile,
    append = TRUE,
    verbose = verbose
  )

  .logMessages("Completed Celda_C. Total time:",
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

  return(bestResult)
}


# Gibbs sampling for the celda_C Model
.cCCalcGibbsProbZ <- function(counts,
                              mCPByS,
                              nGByCP,
                              nByC,
                              nCP,
                              z,
                              s,
                              K,
                              nG,
                              nM,
                              alpha,
                              beta,
                              doSample = TRUE) {

  ## Set variables up front outside of loop
  probs <- matrix(NA, ncol = nM, nrow = K)

  ix <- sample(seq(nM))
  for (i in ix) {
    ## Subtract cell counts from current population assignment
    # nGByCP1 <- nGByCP
    # nGByCP1[, z[i]] <- nGByCP[, z[i]] - counts[, i]
    # nGByCP1 <- .colSums(lgamma(nGByCP1 + beta), nrow(nGByCP), ncol(nGByCP))

    # nCP1 <- nCP
    # nCP1[z[i]] <- nCP1[z[i]] - nByC[i]
    # nCP1 <- lgamma(nCP1 + (nG * beta))

    ## Add cell counts to all other populations
    # nGByCP2 <- nGByCP
    # otherIx <- seq(K)[-z[i]]
    # nGByCP2[, otherIx] <- nGByCP2[, otherIx] + counts[, i]
    # nGByCP2 <- .colSums(lgamma(nGByCP2 + beta), nrow(nGByCP), ncol(nGByCP))

    # nCP2 <- nCP
    # nCP2[otherIx] <- nCP2[otherIx] + nByC[i]
    # nCP2 <- lgamma(nCP2 + (nG * beta))


    mCPByS[z[i], s[i]] <- mCPByS[z[i], s[i]] - 1L

    ## Calculate probabilities for each state
    ## when consider a specific cluster fo this cell,
    ##   no need to calculate cells in other cluster
    for (j in seq_len(K)) {
      # otherIx <- seq(K)[-j]
      if (j != z[i]) { # when j is not current population assignment
        ## Theta simplified
        probs[j, i] <- log(mCPByS[j, s[i]] + alpha) +
          # if adding this cell -- Phi Numerator
          sum(lgamma(nGByCP[, j] + counts[, i] + beta)) -
          # if adding this cell -- Phi Denominator
          lgamma(nCP[j] + nByC[i] + nG * beta) -
          # if without this cell -- Phi Numerator
          sum(lgamma(nGByCP[, j] + beta)) +
          # if without this cell -- Phi Denominator
          lgamma(nCP[j] + nG * beta)
        # sum(nGByCP1[otherIx]) + ## Phi Numerator (other cells)
        # nGByCP2[j] - ## Phi Numerator (current cell)
        # sum(nCP1[otherIx]) - ## Phi Denominator (other cells)
        # nCP2[j] - ## Phi Denominator (current cell)
      } else { # when j is current population assignment
        ## Theta simplified
        probs[j, i] <- log(mCPByS[j, s[i]] + alpha) +
          sum(lgamma(nGByCP[, j] + beta)) -
          lgamma(nCP[j] + nG * beta) -
          sum(lgamma(nGByCP[, j] - counts[, i] + beta)) +
          lgamma(nCP[j] - nByC[i] + nG * beta)
      }
    }

    ## Sample next state and add back counts
    prevZ <- z[i]
    if (isTRUE(doSample)) {
      z[i] <- .sampleLl(probs[, i])
    }

    if (prevZ != z[i]) {
      nGByCP[, prevZ] <- nGByCP[, prevZ] - counts[, i]
      nGByCP[, z[i]] <- nGByCP[, z[i]] + counts[, i]

      nCP[prevZ] <- nCP[prevZ] - nByC[i]
      nCP[z[i]] <- nCP[z[i]] + nByC[i]
    }
    mCPByS[z[i], s[i]] <- mCPByS[z[i], s[i]] + 1L
  }

  return(list(
    mCPByS = mCPByS,
    nGByCP = nGByCP,
    nCP = nCP,
    z = z,
    probs = probs
  ))
}


.cCCalcEMProbZ <- function(counts,
                           mCPByS,
                           nGByCP,
                           nByC,
                           nCP,
                           z,
                           s,
                           K,
                           nG,
                           nM,
                           alpha,
                           beta,
                           doSample = TRUE) {

  ## Expectation given current cell population labels
  theta <- fastNormPropLog(mCPByS, alpha)
  phi <- fastNormPropLog(nGByCP, beta)

  ## Maximization to find best label for each cell
  if (inherits(counts, "matrix") & is.integer(counts)) {
    probs <- eigenMatMultInt(phi, counts) + theta[, s]
  } else if (inherits(counts, "matrix") & is.numeric(counts)) {
    probs <- eigenMatMultNumeric(phi, counts) + theta[, s] 
  } else if (inherits(counts, "dgCMatrix")) {
    probs <- (t(phi) %*% counts) + theta[, s] 
  } else {
    stop("'counts' must be an integer, numeric, or dgCMatrix matrix.")
  }    
  
  if (isTRUE(doSample)) {
    zPrevious <- z
    z <- apply(probs, 2, which.max)

    ## Recalculate counts based on new label
    p <- .cCReDecomposeCounts(counts, s, z, zPrevious, nGByCP, K)
    mCPByS <- p$mCPByS
    nGByCP <- p$nGByCP
    nCP <- p$nCP
  }

  return(list(
    mCPByS = mCPByS,
    nGByCP = nGByCP,
    nCP = nCP,
    z = z,
    probs = probs
  ))
}


# Calculate log-likelihood for celda_C model
.cCCalcLL <- function(mCPByS,
                      nGByCP,
                      s,
                      z,
                      K,
                      nS,
                      nG,
                      alpha,
                      beta) {

  ## Calculate for "Theta" component
  a <- nS * lgamma(K * alpha)
  b <- sum(lgamma(mCPByS + alpha))
  c <- -nS * K * lgamma(alpha)
  d <- -sum(lgamma(colSums(mCPByS + alpha)))

  thetaLl <- a + b + c + d

  ## Calculate for "Phi" component
  a <- K * lgamma(nG * beta)
  b <- sum(lgamma(nGByCP + beta))
  c <- -K * nG * lgamma(beta)
  d <- -sum(lgamma(colSums(nGByCP + beta)))

  phiLl <- a + b + c + d

  final <- thetaLl + phiLl
  return(final)
}


# Takes raw counts matrix and converts it to a series of matrices needed for
# log likelihood calculation
# @param counts Integer matrix. Rows represent features and columns represent
# cells.
# @param s Integer vector. Contains the sample label for each cell (column) in
# the count matrix.
# @param z Numeric vector. Denotes cell population labels.
# @param K Integer. Number of cell populations.
.cCDecomposeCounts <- function(counts, s, z, K) {
  nS <- length(unique(s))
  nG <- nrow(counts)
  nM <- ncol(counts)

  mCPByS <- matrix(as.integer(table(factor(z, levels = seq(K)), s)),
    ncol = nS
  )

  if (inherits(counts, "matrix") & is.integer(counts)) {
    nGByCP <- .colSumByGroup(counts, group = z, K = K)
    nByC <- as.integer(colSums(counts))
  } else if (inherits(counts, "matrix") & is.numeric(counts)) {
    nGByCP <- .colSumByGroupNumeric(counts, group = z, K = K)
    nByC <- as.integer(colSums(counts))
  } else if (inherits(counts, "dgCMatrix")) {
    nGByCP <- colSumByGroupSparse(counts, group = z)
    nByC <- as.integer(Matrix::colSums(counts))
  } else {
    stop("'counts' must be an integer, numeric, or dgCMatrix matrix.")
  }
  
  nCP <- as.integer(colSums(nGByCP))

  return(list(
    mCPByS = mCPByS,
    nGByCP = nGByCP,
    nCP = nCP,
    nByC = nByC,
    nS = nS,
    nG = nG,
    nM = nM
  ))
}


.cCReDecomposeCounts <- function(counts, s, z, previousZ, nGByCP, K) {
  ## Recalculate counts based on new label
  if (inherits(counts, "matrix") & is.integer(counts)) {
    nGByCP <- .colSumByGroupChange(counts, nGByCP, z, previousZ, K)
    nCP <- as.integer(colSums(nGByCP))
  } else if (inherits(counts, "matrix") & is.numeric(counts)) {
    nGByCP <- .colSumByGroupChangeNumeric(counts, nGByCP, z, previousZ, K)
    nCP <- as.integer(colSums(nGByCP))
  } else if (inherits(counts, "dgCMatrix")) {
    nGByCP <- colSumByGroupChangeSparse(counts, nGByCP, group = z,
                                        pgroup = previousZ)
    nCP <- as.integer(Matrix::colSums(nGByCP))
  }
  
  nS <- length(unique(s))
  mCPByS <- matrix(as.integer(table(factor(z, levels = seq(K)), s)),
    ncol = nS
  )

  return(list(
    mCPByS = mCPByS,
    nGByCP = nGByCP,
    nCP = nCP
  ))
}


.prepareCountsForDimReductionCeldaC <- function(sce,
    useAssay,
    maxCells,
    minClusterSize,
    normalize,
    scaleFactor,
    transformationFun) {

    counts <- SummarizedExperiment::assay(sce, i = useAssay)
    counts <- .processCounts(counts)

    ## Checking if maxCells and minClusterSize will work
    if (!is.null(maxCells)) {
        if ((maxCells < ncol(counts)) &
                (maxCells / minClusterSize <
                        S4Vectors::metadata(sce)$celda_parameters$K)) {

            stop("Cannot distribute ",
                maxCells,
                " cells among ",
                S4Vectors::metadata(sce)$celda_parameters$K,
                " clusters while maintaining a minumum of ",
                minClusterSize,
                " cells per cluster. Try increasing 'maxCells' or decreasing",
                " 'minClusterSize'.")
        }
    } else {
        maxCells <- ncol(counts)
    }

    ## Select a subset of cells to sample if greater than 'maxCells'
    totalCellsToRemove <- ncol(counts) - maxCells
    zInclude <- rep(TRUE, ncol(counts))

    if (totalCellsToRemove > 0) {
        zTa <- tabulate(SummarizedExperiment::colData(sce)$celda_cell_cluster,
            S4Vectors::metadata(sce)$celda_parameters$K)

        ## Number of cells that can be sampled from each cluster without
        ## going below the minimum threshold
        clusterCellsToSample <- zTa - minClusterSize
        clusterCellsToSample[clusterCellsToSample < 0] <- 0

        ## Number of cells to sample after exluding smaller clusters
        ## Rounding can cause number to be off by a few, so ceiling is
        ## used with a second round of subtraction
        clusterNToSample <- ceiling((clusterCellsToSample /
                sum(clusterCellsToSample)) * totalCellsToRemove)
        diff <- sum(clusterNToSample) - totalCellsToRemove
        clusterNToSample[which.max(clusterNToSample)] <-
            clusterNToSample[which.max(clusterNToSample)] - diff

        ## Perform sampling for each cluster
        for (i in which(clusterNToSample > 0)) {
            zInclude[sample(which(
                SummarizedExperiment::colData(sce)$celda_cell_cluster == i),
                clusterNToSample[i])] <- FALSE
        }
    }

    cellIx <- which(zInclude)
    norm <- t(normalizeCounts(counts[, cellIx],
        normalize = normalize,
        scaleFactor = scaleFactor,
        transformationFun = transformationFun))
    return(list(norm = norm, cellIx = cellIx))
}


.createSCEceldaC <- function(celdaCMod,
    sce,
    xClass,
    useAssay,
    algorithm,
    stopIter,
    maxIter,
    splitOnIter,
    splitOnLast,
    nchains,
    zInitialize,
    zInit,
    logfile,
    verbose) {

    # add metadata
    S4Vectors::metadata(sce)[["celda_parameters"]] <- list(
        model = "celda_C",
        xClass = xClass,
        useAssay = useAssay,
        sampleLevels = celdaCMod@names$sample,
        K = celdaCMod@params$K,
        alpha = celdaCMod@params$alpha,
        beta = celdaCMod@params$beta,
        algorithm = algorithm,
        stopIter = stopIter,
        maxIter = maxIter,
        splitOnIter = splitOnIter,
        splitOnLast = splitOnLast,
        seed = celdaCMod@params$seed,
        nchains = nchains,
        zInitialize = zInitialize,
        countChecksum = celdaCMod@params$countChecksum,
        zInit = zInit,
        logfile = logfile,
        verbose = verbose,
        completeLogLik = celdaCMod@completeLogLik,
        finalLogLik = celdaCMod@finalLogLik,
        cellClusterLevels = sort(unique(celdaClusters(celdaCMod)$z)))

    SummarizedExperiment::rowData(sce)["rownames"] <- celdaCMod@names$row
    SummarizedExperiment::colData(sce)["colnames"] <-
        celdaCMod@names$column
    SummarizedExperiment::colData(sce)["celda_sample_label"] <-
        celdaCMod@sampleLabel
    SummarizedExperiment::colData(sce)["celda_cell_cluster"] <-
        celdaClusters(celdaCMod)$z

    return(sce)
}
