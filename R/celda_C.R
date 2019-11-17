#' @title Cell clustering with Celda
#' @description Clusters the columns of a count matrix containing single-cell
#'  data into K subpopulations.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells.
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
#' @param zInitialize Chararacter. One of 'random', 'split', or 'predefined'.
#'  With 'random', cells are randomly assigned to a populations. With 'split',
#'  cells will be split into sqrt(K) populations and then each popluation will
#'  be subsequently split into another sqrt(K) populations. With 'predefined',
#'  values in `zInit` will be used to initialize `z`. Default 'split'.
#' @param zInit Integer vector. Sets initial starting values of z. If NULL,
#'  starting values for each cell will be randomly sampled from `1:K`. 'zInit'
#'  can only be used when `initialize = 'random'`. Default NULL.
#' @param countChecksum "Character. An MD5 checksum for the `counts` matrix.
#'  Default NULL.
#' @param logfile Character. Messages will be redirected to a file named
#'  `logfile`. If NULL, messages will be printed to stdout.  Default NULL.
#' @param verbose Logical. Whether to print log messages. Default TRUE.
#' @return An object of class `celda_C` with the cell population clusters
#'  stored in `z`.
#' @seealso `celda_G()` for feature clustering and `celda_CG()` for simultaneous
#'  clustering of features and cells. `celdaGridSearch()` can be used to run
#'  multiple values of K and multiple chains in parallel.
#' @examples
#' data(celdaCSim)
#' celdaCMod <- celda_C(celdaCSim$counts,
#'     K = celdaCSim$K,
#'     sampleLabel = celdaCSim$sampleLabel)
#' @import Rcpp RcppEigen
#' @rawNamespace import(gridExtra, except = c(combine))
#' @importFrom withr with_seed
#' @export
celda_C <- function(counts,
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

    .validateCounts(counts)
    if (is.null(seed)) {
        res <- .celda_C(counts,
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
            verbose,
            reorder = TRUE)
    } else {
        with_seed(seed,
            res <- .celda_C(counts,
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
                verbose,
                reorder = TRUE))
    }

    return(res)
}


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
    seed = 12345,
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
    counts <- .processCounts(counts)
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
        ".cCCalcEMProbZ")
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
            verbose = verbose)

        if (zInitialize == "predefined") {
            if (is.null(zInit)) {
                stop("'zInit' needs to specified when initilize.z == 'given'.")
            }

            z <- .initializeCluster(K,
                ncol(counts),
                initial = zInit,
                fixed = NULL)
        } else if (zInitialize == "split") {
            z <- .initializeSplitZ(counts,
                    K = K,
                    alpha = alpha,
                    beta = beta)
        } else {
            z <- .initializeCluster(K,
                    ncol(counts),
                    initial = NULL,
                    fixed = NULL)
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

        ll <- .cCCalcLL(mCPByS = mCPByS,
                nGByCP = nGByCP,
                s = s,
                K = K,
                nS = nS,
                nG = nG,
                alpha = alpha,
                beta = beta)

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
                beta = beta))

            mCPByS <- nextZ$mCPByS
            nGByCP <- nextZ$nGByCP
            nCP <- nextZ$nCP
            z <- nextZ$z

            ## Perform split on i-th iteration of no improvement in log
            ## likelihood
            tempLl <- .cCCalcLL(mCPByS = mCPByS,
                    nGByCP = nGByCP,
                    s = s,
                    K = K,
                    nS = nS,
                    nG = nG,
                    alpha = alpha,
                    beta = beta)

            if (K > 2 & iter != maxIter &
                ((((numIterWithoutImprovement == stopIter &
                    !all(tempLl > ll))) & isTRUE(splitOnLast)) |
                        (splitOnIter > 0 & iter %% splitOnIter == 0 &
                            isTRUE(doCellSplit)))) {

                .logMessages(date(),
                    " .... Determining if any cell clusters should be split.",
                    logfile = logfile,
                    append = TRUE,
                    sep = "",
                    verbose = verbose)

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
                    zProb = t(nextZ$probs),
                    maxClustersToTry = K,
                    minCell = 3)

                .logMessages(res$message,
                    logfile = logfile,
                    append = TRUE,
                    verbose = verbose)

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
            tempLl <- .cCCalcLL(mCPByS = mCPByS,
                    nGByCP = nGByCP,
                    s = s,
                    K = K,
                    nS = nS,
                    nG = nG,
                    alpha = alpha,
                    beta = beta)

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
                verbose = verbose)
            iter <- iter + 1
        }

        names <- list(row = rownames(counts),
            column = colnames(counts),
            sample = levels(sampleLabel))

        result <- list(z = zBest,
            completeLogLik = ll,
            finalLogLik = llBest,
            K = K,
            sampleLabel = sampleLabel,
            alpha = alpha,
            beta = beta,
            countChecksum = countChecksum,
            names = names)

        if (is.null(bestResult) ||
            result$finalLogLik > bestResult$finalLogLik) {

            bestResult <- result
        }

        .logMessages(date(),
            ".. Finished chain",
            i,
            logfile = logfile,
            append = TRUE,
            verbose = verbose)
    }

    bestResult <- methods::new("celda_C",
        clusters = list(z = bestResult$z),
        params = list(K = as.integer(bestResult$K),
            alpha = bestResult$alpha,
            beta = bestResult$beta,
            countChecksum = bestResult$countChecksum),
        sampleLabel = bestResult$sampleLabel,
        completeLogLik = bestResult$completeLogLik,
        finalLogLik = bestResult$finalLogLik,
        names = bestResult$names)

    if (isTRUE(reorder)) {
        bestResult <- .reorderCelda_C(counts = counts, res = bestResult)
    }

    endTime <- Sys.time()
    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)

    .logMessages("Completed Celda_C. Total time:",
        format(difftime(endTime, startTime)),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)

    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)

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
        #nGByCP1 <- nGByCP
        #nGByCP1[, z[i]] <- nGByCP[, z[i]] - counts[, i]
        #nGByCP1 <- .colSums(lgamma(nGByCP1 + beta), nrow(nGByCP), ncol(nGByCP))

        #nCP1 <- nCP
        #nCP1[z[i]] <- nCP1[z[i]] - nByC[i]
        #nCP1 <- lgamma(nCP1 + (nG * beta))

        ## Add cell counts to all other populations
        #nGByCP2 <- nGByCP
        #otherIx <- seq(K)[-z[i]]
        #nGByCP2[, otherIx] <- nGByCP2[, otherIx] + counts[, i]
        #nGByCP2 <- .colSums(lgamma(nGByCP2 + beta), nrow(nGByCP), ncol(nGByCP))

        #nCP2 <- nCP
        #nCP2[otherIx] <- nCP2[otherIx] + nByC[i]
        #nCP2 <- lgamma(nCP2 + (nG * beta))


        mCPByS[z[i], s[i]] <- mCPByS[z[i], s[i]] - 1L

        ## Calculate probabilities for each state
		## when consider a specific cluster fo this cell, no need to calculate cells in other cluster	
        for (j in seq_len(K)) {
            #otherIx <- seq(K)[-j]
		if (j != z[i]) { # when j is not current population assignment
				probs[j, i] <- log(mCPByS[j, s[i]] + alpha) + ## Theta simplified
				sum(lgamma(nGByCP[, j] + counts[, i] + beta)) -   # if adding this cell -- Phi Numerator
				lgamma(nCP[j] + nByC[i] + nG * beta) -    # if adding this cell -- Phi Denominator
				sum(lgamma(nGByCP[, j] + beta)) +   # if without this cell -- Phi Numerator
				lgamma(nCP[j] + nG * beta  )    # if without this cell -- Phi Denominator
				#sum(nGByCP1[otherIx]) + ## Phi Numerator (other cells)
				#nGByCP2[j] - ## Phi Numerator (current cell)
				#sum(nCP1[otherIx]) - ## Phi Denominator (other cells)
				#nCP2[j] - ## Phi Denominator (current cell)
		}
		else {  # when j is current population assignment 
				probs[j, i] <- log(mCPByS[j, s[i]] + alpha) + ## Theta simplified
				sum(lgamma(nGByCP[, j] + beta)) - 
				lgamma(nCP[j] + nG * beta) -
				sum(lgamma(nGByCP[, j] - counts[, i] + beta)) +
				lgamma(nCP[j] - nByC[i] + nG * beta) 

		}

        }

        ## Sample next state and add back counts
        prevZ <- z[i]
        if (isTRUE(doSample))
            z[i] <- .sampleLl(probs[, i])

        if (prevZ != z[i]) {
            nGByCP[, prevZ] <- nGByCP[, prevZ] - counts[, i]
            nGByCP[, z[i]] <- nGByCP[, z[i]] + counts[, i]

            nCP[prevZ] <- nCP[prevZ] - nByC[i]
            nCP[z[i]] <- nCP[z[i]] + nByC[i]
        }
        mCPByS[z[i], s[i]] <- mCPByS[z[i], s[i]] + 1L
    }

    return(list(mCPByS = mCPByS,
        nGByCP = nGByCP,
        nCP = nCP,
        z = z,
        probs = probs))
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
    probs <- eigenMatMultInt(phi, counts) + theta[, s]

    if (isTRUE(doSample)) {
        zPrevious <- z
        z <- apply(probs, 2, which.max)

        ## Recalculate counts based on new label
        p <- .cCReDecomposeCounts(counts, s, z, zPrevious, nGByCP, K)
        mCPByS <- p$mCPByS
        nGByCP <- p$nGByCP
        nCP <- p$nCP
    }

    return(list(mCPByS = mCPByS,
        nGByCP = nGByCP,
        nCP = nCP,
        z = z,
        probs = probs))
}


#' @title Simulate cells from the celda_C model
#' @description Generates a simulated counts matrix, cell subpopulation
#'  clusters, and sample labels according to the generative process of the
#'  celda_C model.
#' @param model Character. Options available in `celda::availableModels`.
#' @param S Integer. Number of samples to simulate. Default 5.
#' @param CRange Vector of length 2 given the range (min, max) of number of
#'  cells for each sample to be randomly generated from the uniform
#'  distribution. Default c(50, 100).
#' @param NRange Integer vector. A vector of length 2 that specifies the lower
#'  and upper bounds of the number of counts generated for each cell. Default
#'  c(500, 1000).
#' @param G Integer. The total number of features to be simulated. Default 100.
#' @param K Integer. Number of cell populations. Default 5.
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount
#'  to each cell population in each sample. Default 1.
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to
#'  each feature in each cell population. Default 1.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @param ... Additional parameters.
#' @return List. Contains the simulated matrix `counts`, cell population
#'  clusters `z`, sample assignments `sampleLabel`, and input parameters.
#' @seealso `celda_G()` for simulating feature modules and `celda_CG()` for
#'  simulating feature modules and cell populations.
#' @examples
#' celdaCSim <- simulateCells(model = "celda_C", K = 10)
#' simCounts <- celdaCSim$counts
#' @rawNamespace import(stats, except = c(start, end))
#' @export
simulateCellscelda_C <- function(model,
    S = 5,
    CRange = c(50, 100),
    NRange = c(500, 1000),
    G = 100,
    K = 5,
    alpha = 1,
    beta = 1,
    seed = 12345,
    ...) {

    if (is.null(seed)) {
        res <- .simulateCellscelda_C(model = model,
            S = S,
            CRange = CRange,
            NRange = NRange,
            G = G,
            K = K,
            alpha = alpha,
            beta = beta,
            ...)
    } else {
        res <- with_seed(seed,
            .simulateCellscelda_C(model = model,
                S = S,
                CRange = CRange,
                NRange = NRange,
                G = G,
                K = K,
                alpha = alpha,
                beta = beta,
                ...))
    }

    return(res)
}


.simulateCellscelda_C <- function(model,
    S = 5,
    CRange = c(50, 100),
    NRange = c(500, 1000),
    G = 100,
    K = 5,
    alpha = 1,
    beta = 1,
    ...) {

    phi <- .rdirichlet(K, rep(beta, G))
    theta <- .rdirichlet(S, rep(alpha, K))

    ## Select the number of cells per sample
    nC <- sample(seq(CRange[1], CRange[2]), size = S, replace = TRUE)
    cellSampleLabel <- rep(seq(S), nC)

    ## Select state of the cells
    z <- unlist(lapply(seq(S), function(i) {
        sample(seq(K),
            size = nC[i],
            prob = theta[i, ],
            replace = TRUE)
    }))

    ## Select number of transcripts per cell
    nN <- sample(seq(NRange[1], NRange[2]),
        size = length(cellSampleLabel),
        replace = TRUE)

    ## Select transcript distribution for each cell
    cellCounts <- vapply(seq(length(cellSampleLabel)), function(i) {
        stats::rmultinom(1, size = nN[i], prob = phi[z[i], ])
    }, integer(G))

    rownames(cellCounts) <- paste0("Gene_", seq(nrow(cellCounts)))
    colnames(cellCounts) <- paste0("Cell_", seq(ncol(cellCounts)))
    cellSampleLabel <- paste0("Sample_", seq(S))[cellSampleLabel]
    cellSampleLabel <- factor(cellSampleLabel,
        levels = paste0("Sample_", seq(S)))

    ## Peform reordering on final Z and Y assigments:
    cellCounts <- .processCounts(cellCounts)
    names <- list(row = rownames(cellCounts),
        column = colnames(cellCounts),
        sample = unique(cellSampleLabel))
    countChecksum <- .createCountChecksum(cellCounts)
    result <- methods::new("celda_C",
        clusters = list(z = z),
        params = list(K = as.integer(K),
            alpha = alpha,
            beta = beta,
            countChecksum = countChecksum),
        sampleLabel = cellSampleLabel,
        names = names)
    class(result) <- "celda_C"
    result <- .reorderCelda_C(counts = cellCounts, res = result)

    return(list(z = clusters(result)$z,
        counts = .processCounts(cellCounts),
        sampleLabel = cellSampleLabel,
        K = K,
        alpha = alpha,
        beta = beta,
        CRange = CRange,
        NRange = NRange,
        S = S))
}


#' @title Matrix factorization for results from celda_C()
#' @description Generates factorized matrices showing the contribution of each
#'  feature in each cell population or each cell population in each sample.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda object of class "celda_C".
#' @param type Character vector. A vector containing one or more of "counts",
#'  "proportion", or "posterior". "counts" returns the raw number of counts for
#'  each factorized matrix. "proportions" returns the normalized probabilities
#'  for each factorized matrix, which are calculated by dividing the raw counts
#'  in each factorized matrix by the total counts in each column. "posterior"
#'  returns the posterior estimates. Default
#'  `c("counts", "proportion", "posterior")`.
#' @examples
#' data(celdaCSim, celdaCMod)
#' factorizedMatrices <- factorizeMatrix(celdaCSim$counts,
#'     celdaCMod, "posterior")
#' @return A list with elements for `counts`, `proportions`, or `posterior`
#'  probabilities. Each element will be a list containing factorized matrices
#'  for `module` and `sample`.
#' @seealso `celda_C()` for clustering cells
#' @export
setMethod("factorizeMatrix", signature(celdaMod = "celda_C"),
    function(counts,
        celdaMod,
        type = c("counts", "proportion", "posterior")) {

        counts <- .processCounts(counts)
        compareCountMatrix(counts, celdaMod)

        K <- params(celdaMod)$K
        z <- clusters(celdaMod)$z
        alpha <- params(celdaMod)$alpha
        beta <- params(celdaMod)$beta
        sampleLabel <- sampleLabel(celdaMod)
        s <- as.integer(sampleLabel)

        p <- .cCDecomposeCounts(counts, s, z, K)
        mCPByS <- p$mCPByS
        nGByCP <- p$nGByCP

        KNames <- paste0("K", seq(K))
        rownames(nGByCP) <- matrixNames(celdaMod)$row
        colnames(nGByCP) <- KNames
        rownames(mCPByS) <- KNames
        colnames(mCPByS) <- matrixNames(celdaMod)$sample

        countsList <- c()
        propList <- c()
        postList <- c()
        res <- list()

        if (any("counts" %in% type)) {
            countsList <- list(sample = mCPByS, module = nGByCP)
            res <- c(res, list(counts = countsList))
        }

        if (any("proportion" %in% type)) {
            ## Need to avoid normalizing cell/gene states with zero cells/genes
            uniqueZ <- sort(unique(z))
            tempNGByCP <- nGByCP
            tempNGByCP[, uniqueZ] <- normalizeCounts(tempNGByCP[, uniqueZ],
                normalize = "proportion")

            propList <- list(sample = normalizeCounts(mCPByS,
                normalize = "proportion"),
                module = tempNGByCP)
            res <- c(res, list(proportions = propList))
        }

        if (any("posterior" %in% type)) {
            postList <- list(sample = normalizeCounts(mCPByS + alpha,
                normalize = "proportion"),
                module = normalizeCounts(nGByCP + beta,
                    normalize = "proportion"))

            res <- c(res, posterior = list(postList))
        }

        return(res)
    })


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


#' @title Calculate Celda_C log likelihood
#' @description Calculates the log likelihood for user-provided cell population
#'  clusters using the `celda_C()` model.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells.
#' @param sampleLabel Vector or factor. Denotes the sample label for each cell
#'  (column) in the count matrix.
#' @param z Numeric vector. Denotes cell population labels.
#' @param K Integer. Number of cell populations.
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount
#'  to each cell population in each sample. Default 1.
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to
#'  each feature in each cell population. Default 1.
#' @return Numeric. The log likelihood for the given cluster assignments
#' @seealso `celda_C()` for clustering cells
#' @examples
#' data(celdaCSim)
#' loglik <- logLikelihoodcelda_C(celdaCSim$counts,
#'     sampleLabel = celdaCSim$sampleLabel,
#'     z = celdaCSim$z,
#'     K = celdaCSim$K,
#'     alpha = celdaCSim$alpha,
#'     beta = celdaCSim$beta)
#'
#' loglik <- logLikelihood(celdaCSim$counts,
#'     model = "celda_C",
#'     sampleLabel = celdaCSim$sampleLabel,
#'     z = celdaCSim$z,
#'     K = celdaCSim$K,
#'     alpha = celdaCSim$alpha,
#'     beta = celdaCSim$beta)
#' @export
logLikelihoodcelda_C <- function(counts, sampleLabel, z, K, alpha, beta) {

    if (sum(z > K) > 0) {
        stop("An entry in z contains a value greater than the provided K.")
    }
    sampleLabel <- .processSampleLabels(sampleLabel, ncol(counts))
    s <- as.integer(sampleLabel)
    p <- .cCDecomposeCounts(counts, s, z, K)
    final <- .cCCalcLL(mCPByS = p$mCPByS,
        nGByCP = p$nGByCP,
        s = s,
        z = z,
        K = K,
        nS = p$nS,
        nG = p$nG,
        alpha = alpha,
        beta = beta)
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
        ncol = nS)
    nGByCP <- .colSumByGroup(counts, group = z, K = K)
    nCP <- as.integer(colSums(nGByCP))
    nByC <- as.integer(colSums(counts))

    return(list(mCPByS = mCPByS,
            nGByCP = nGByCP,
            nCP = nCP,
            nByC = nByC,
            nS = nS,
            nG = nG,
            nM = nM))
}


.cCReDecomposeCounts <- function(counts, s, z, previousZ, nGByCP, K) {
    ## Recalculate counts based on new label
    nGByCP <- .colSumByGroupChange(counts, nGByCP, z, previousZ, K)
    nS <- length(unique(s))
    mCPByS <- matrix(as.integer(table(factor(z, levels = seq(K)), s)),
        ncol = nS)
    nCP <- as.integer(colSums(nGByCP))

    return(list(mCPByS = mCPByS,
        nGByCP = nGByCP,
        nCP = nCP))
}


#' @title Conditional probabilities for cells in subpopulations from a Celda_C
#'  model
#' @description Calculates the conditional probability of each cell belonging to
#'  each subpopulation given all other cell cluster assignments in a `celda_C()`
#'  result.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda object of class `celda_C`.
#' @param log Logical. If FALSE, then the normalized conditional probabilities
#'  will be returned. If TRUE, then the unnormalized log probabilities will be
#'  returned. Default FALSE.
#' @param ... Additional parameters.
#' @return A list containging a matrix for the conditional cell subpopulation
#'  cluster probabilities.
#' @seealso `celda_C()` for clustering cells
#' @examples
#' data(celdaCSim, celdaCMod)
#' clusterProb <- clusterProbability(celdaCSim$counts, celdaCMod)
#' @export
setMethod("clusterProbability", signature(celdaMod = "celda_C"),
    function(counts, celdaMod, log = FALSE, ...) {
        z <- clusters(celdaMod)$z
        sampleLabel <- sampleLabel(celdaMod)
        s <- as.integer(sampleLabel)

        K <- params(celdaMod)$K
        alpha <- params(celdaMod)$alpha
        beta <- params(celdaMod)$beta

        p <- .cCDecomposeCounts(counts, s, z, K)

        nextZ <- .cCCalcGibbsProbZ(counts = counts,
            mCPByS = p$mCPByS,
            nGByCP = p$nGByCP,
            nByC = p$nByC,
            nCP = p$nCP,
            z = z,
            s = s,
            K = K,
            nG = p$nG,
            nM = p$nM,
            alpha = alpha,
            beta = beta,
            doSample = FALSE)
        zProb <- t(nextZ$probs)

        if (!isTRUE(log)) {
            zProb <- .normalizeLogProbs(zProb)
        }

        return(list(zProbability = zProb))
    })


#' @title Calculate the perplexity on new data with a celda_C model
#' @description Perplexity is a statistical measure of how well a probability
#'  model can predict new data. Lower perplexity indicates a better model.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda object of class "celda_C"
#' @param newCounts A new counts matrix used to calculate perplexity. If NULL,
#'  perplexity will be calculated for the 'counts' matrix. Default NULL.
#' @return Numeric. The perplexity for the provided count data and model.
#' @seealso `celda_C()` for clustering cells
#' @examples
#' data(celdaCSim, celdaCMod)
#' perplexity <- perplexity(celdaCSim$counts, celdaCMod)
#' @importFrom matrixStats logSumExp
#' @export
setMethod("perplexity", signature(celdaMod = "celda_C"),
    function(counts, celdaMod, newCounts = NULL) {
        if (!("celda_C" %in% class(celdaMod))) {
            stop("The celdaMod provided was not of class celda_C.")
        }

        counts <- .processCounts(counts)
        compareCountMatrix(counts, celdaMod)

        if (is.null(newCounts)) {
            newCounts <- counts
        } else {
            newCounts <- .processCounts(newCounts)
        }

        if (nrow(newCounts) != nrow(counts)) {
            stop("newCounts should have the same number of rows as counts.")
        }

        factorized <- factorizeMatrix(counts = counts,
            celdaMod = celdaMod,
            type = "posterior")
        theta <- log(factorized$posterior$sample)
        phi <- log(factorized$posterior$module)
        s <- as.integer(sampleLabel(celdaMod))

        # inner.log.prob = (t(phi) %*% newCounts) + theta[, s]
        inner.log.prob <- eigenMatMultInt(phi, newCounts) + theta[, s]
        logPx <- sum(apply(inner.log.prob, 2, matrixStats::logSumExp))

        perplexity <- exp(- (logPx / sum(newCounts)))
        return(perplexity)
    })


.reorderCelda_C <- function(counts, res) {
    if (params(res)$K > 2 & isTRUE(length(unique(clusters(res)$z)) > 1)) {
        res@clusters$z <- as.integer(as.factor(clusters(res)$z))
        fm <- factorizeMatrix(counts = counts, celdaMod = res)
        uniqueZ <- sort(unique(clusters(res)$z))
        d <- .cosineDist(fm$posterior$module[, uniqueZ])
        h <- stats::hclust(d, method = "complete")
        res <- recodeClusterZ(res,
            from = h$order,
            to = seq(length(h$order)))
    }
    return(res)
}


#' @title Heatmap for celda_C
#' @description Renders an expression heatmap to visualize `celda_C()` results.
#'  Features to include in the heatmap must be supplied.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda object of class `celda_C`.
#' @param featureIx Integer vector. Indices of features to plot, such the top
#'  features from a differential expression analysis.
#' @param ... Additional parameters.
#' @seealso `celda_C()` for clustering cells and `celdaTsne()` for generating
#'  2-dimensional coordinates
#' @examples
#' data(celdaCSim, celdaCMod)
#' celdaHeatmap(celdaCSim$counts, celdaCMod)
#' @return list A list containing dendrograms and the heatmap grob
#' @export
setMethod("celdaHeatmap", signature(celdaMod = "celda_C"),
    function(counts, celdaMod, featureIx, ...) {
        norm <- normalizeCounts(counts,
            normalize = "proportion",
            transformationFun = sqrt)
        plotHeatmap(norm[featureIx, ], z = clusters(celdaMod)$z, ...)
    })


#' @title tSNE for celda_C
#' @description Embeds cells in two dimensions using tSNE based on a `celda_C`
#'  model. PCA on the normalized counts is used to reduce the number of
#'  features before applying tSNE.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda object of class `celda_C`.
#' @param maxCells Integer. Maximum number of cells to plot. Cells will be
#'  randomly subsampled if ncol(counts) > maxCells. Larger numbers of cells
#'  requires more memory. If NULL, no subsampling will be performed.
#'  Default NULL.
#' @param minClusterSize Integer. Do not subsample cell clusters below this
#'  threshold. Default 100.
#' @param initialDims Integer. PCA will be used to reduce the dimentionality
#'  of the dataset. The top 'initialDims' principal components will be used
#'  for tSNE. Default 20.
#' @param perplexity Numeric. Perplexity parameter for tSNE. Default 20.
#' @param maxIter Integer. Maximum number of iterations in tSNE generation.
#'  Default 2500.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @seealso `celda_C()` for clustering cells and `celdaHeatmap()` for displaying
#'  expression
#' @examples
#' data(celdaCSim, celdaCMod)
#' tsneRes <- celdaTsne(celdaCSim$counts, celdaCMod)
#' @return A two column matrix of t-SNE coordinates
#' @export
setMethod("celdaTsne", signature(celdaMod = "celda_C"),
    function(counts,
        celdaMod,
        maxCells = NULL,
        minClusterSize = 100,
        initialDims = 20,
        perplexity = 20,
        maxIter = 2500,
        seed = 12345) {

        if (is.null(seed)) {
            res <- .celdaTsneC(counts = counts,
                celdaMod = celdaMod,
                maxCells = maxCells,
                minClusterSize = minClusterSize,
                initialDims = initialDims,
                perplexity = perplexity,
                maxIter = maxIter)
        } else {
            with_seed(seed,
                res <- .celdaTsneC(counts = counts,
                    celdaMod = celdaMod,
                    maxCells = maxCells,
                    minClusterSize = minClusterSize,
                    initialDims = initialDims,
                    perplexity = perplexity,
                    maxIter = maxIter))
        }

        return(res)
    })


.celdaTsneC <- function(counts,
    celdaMod,
    maxCells = NULL,
    minClusterSize = 100,
    initialDims = 20,
    perplexity = 20,
    maxIter = 2500) {

    preparedCountInfo <- .prepareCountsForDimReductionCeldaC(counts,
        celdaMod,
        maxCells,
        minClusterSize)

    res <- .calculateTsne(preparedCountInfo$norm,
        perplexity = perplexity,
        maxIter = maxIter,
        doPca = TRUE,
        initialDims = initialDims)

    final <- matrix(NA, nrow = ncol(counts), ncol = 2)
    final[preparedCountInfo$cellIx, ] <- res
    rownames(final) <- colnames(counts)
    colnames(final) <- c("tSNE1", "tSNE_2")
    return(final)
}


#' @title umap for celda_C
#' @description Embeds cells in two dimensions using umap based on a `celda_C`
#'  model. PCA on the normalized counts is used to reduce the number of features
#'  before applying umap.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda object of class `celda_C`.
#' @param maxCells Integer. Maximum number of cells to plot. Cells will be
#'  randomly subsampled if ncol(counts) > maxCells. Larger numbers of cells
#'  requires more memory. If NULL, no subsampling will be performed.
#'  Default NULL.
#' @param minClusterSize Integer. Do not subsample cell clusters below this
#'  threshold. Default 100.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @param nNeighbors The size of local neighborhood used for
#'   manifold approximation. Larger values result in more global
#'   views of the manifold, while smaller values result in more
#'   local data being preserved. Default 30.
#'   See `?uwot::umap` for more information.
#' @param minDist The effective minimum distance between embedded points.
#'   Smaller values will result in a more clustered/clumped
#'   embedding where nearby points on the manifold are drawn
#'   closer together, while larger values will result on a more
#'   even dispersal of points. Default 0.2.
#'   See `?uwot::umap` for more information.
#' @param spread The effective scale of embedded points. In combination with
#'   ‘min_dist’, this determines how clustered/clumped the
#'   embedded points are. Default 1. See `?uwot::umap` for more information.
#' @param pca Logical. Whether to perform
#' dimensionality reduction with PCA before UMAP.
#' @param initialDims Integer. Number of dimensions from PCA to use as
#' input in UMAP. Default 50.
#' @param cores Number of threads to use. Default 1.
#' @param ... Other parameters to pass to `uwot::umap`.
#' @seealso `celda_C()` for clustering cells and `celdaHeatmap()` for displaying
#'  expression.
#' @examples
#' data(celdaCSim, celdaCMod)
#' umapRes <- celdaUmap(celdaCSim$counts, celdaCMod)
#' @return A two column matrix of UMAP coordinates
#' @export
setMethod("celdaUmap", signature(celdaMod = "celda_C"),
    function(counts,
        celdaMod,
        maxCells = NULL,
        minClusterSize = 100,
        seed = 12345,
        nNeighbors = 30,
        minDist = 0.75,
        spread = 1,
        pca = TRUE,
        initialDims = 50,
        cores = 1,
        ...) {

        if (is.null(seed)) {
            res <- .celdaUmapC(counts = counts,
                celdaMod = celdaMod,
                maxCells = maxCells,
                minClusterSize = minClusterSize,
                nNeighbors = nNeighbors,
                minDist = minDist,
                spread = spread,
                pca = pca,
                initialDims = initialDims,
                cores = cores,
                ...)
        } else {
            with_seed(seed,
                res <- .celdaUmapC(counts = counts,
                    celdaMod = celdaMod,
                    maxCells = maxCells,
                    minClusterSize = minClusterSize,
                    nNeighbors = nNeighbors,
                    minDist = minDist,
                    spread = spread,
                    pca = pca,
                    initialDims = initialDims,
                    cores = cores,
                    ...))
        }

        return(res)
    })


.celdaUmapC <- function(counts,
    celdaMod,
    maxCells = NULL,
    minClusterSize = 100,
    nNeighbors = 30,
    minDist = 0.2,
    spread = 1,
    pca = TRUE,
    initialDims = 50,
    cores = 1,
    ...) {

    preparedCountInfo <- .prepareCountsForDimReductionCeldaC(counts,
        celdaMod,
        maxCells,
        minClusterSize)
    umapRes <- .calculateUmap(preparedCountInfo$norm,
        nNeighbors = nNeighbors,
        minDist = minDist,
        spread = spread,
        pca = pca,
        initialDims = initialDims,
        cores = cores,
        ...)

    final <- matrix(NA, nrow = ncol(counts), ncol = 2)
    final[preparedCountInfo$cellIx, ] <- umapRes
    rownames(final) <- colnames(counts)
    colnames(final) <- c("UMAP_1", "UMAP_2")
    return(final)
}


.prepareCountsForDimReductionCeldaC <- function(counts,
    celdaMod,
    maxCells = NULL,
    minClusterSize = 100) {

    counts <- .processCounts(counts)
    compareCountMatrix(counts, celdaMod)

    ## Checking if maxCells and minClusterSize will work
    if (!is.null(maxCells)) {
      if ((maxCells < ncol(counts)) &
            (maxCells / minClusterSize < params(celdaMod)$K)) {

        stop("Cannot distribute ",
            maxCells,
            " cells among ",
            params(celdaMod)$K,
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
        zTa <- tabulate(clusters(celdaMod)$z, params(celdaMod)$K)

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
            zInclude[sample(which(clusters(celdaMod)$z == i),
                clusterNToSample[i])] <- FALSE
        }
    }

    cellIx <- which(zInclude)
    norm <- t(normalizeCounts(counts[, cellIx],
        normalize = "proportion",
        transformationFun = sqrt))
    return(list(norm = norm, cellIx = cellIx))
}


#' @title Probability map for a celda_C model
#' @description Renders probability and relative expression heatmaps to
#'  visualize the relationship between cell populations and samples.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda object of class `celda_C`.
#' @param level Character. 'sample' will display the absolute probabilities and
#'  relative normalized abundance of each cell population in each sample.
#'  Default 'sample'.
#' @param ... Additional parameters.
#' @seealso `celda_C()` for clustering cells
#' @examples
#' data(celdaCSim, celdaCMod)
#' celdaProbabilityMap(celdaCSim$counts, celdaCMod)
#' @return A grob containing the specified plots
#' @importFrom gridExtra grid.arrange
#' @export
setMethod("celdaProbabilityMap", signature(celdaMod = "celda_C"),
    function(counts, celdaMod, level = c("sample"), ...) {
        counts <- .processCounts(counts)
        compareCountMatrix(counts, celdaMod)

        zInclude <- which(tabulate(clusters(celdaMod)$z,
            params(celdaMod)$K) > 0)

        level <- match.arg(level)
        factorized <- factorizeMatrix(celdaMod = celdaMod, counts = counts)

        samp <- factorized$proportions$sample[zInclude, , drop = FALSE]
        col <- colorRampPalette(c("white",
            "blue",
            "#08306B",
            "#006D2C",
            "yellowgreen",
            "yellow",
            "orange",
            "red"))(100)
        breaks <- seq(0, 1, length.out = length(col))
        g1 <- plotHeatmap(samp,
            colorScheme = "sequential",
            scaleRow = NULL,
            clusterCell = FALSE,
            clusterFeature = FALSE,
            showNamesCell = TRUE,
            showNamesFeature = TRUE,
            breaks = breaks,
            col = col,
            main = "Absolute Probability",
            silent = TRUE)

        if (ncol(samp) > 1) {
            sampNorm <- normalizeCounts(samp,
                normalize = "proportion",
                transformationFun = sqrt,
                scaleFun = base::scale)
            g2 <- plotHeatmap(sampNorm,
                colorScheme = "divergent",
                clusterCell = FALSE,
                clusterFeature = FALSE,
                showNamesCell = TRUE,
                showNamesFeature = TRUE,
                main = "Relative Abundance",
                silent = TRUE)
            gridExtra::grid.arrange(g1$gtable, g2$gtable, ncol = 2)
        } else {
            gridExtra::grid.arrange(g1$gtable)
        }
    })


#' @title Lookup the module of a feature
#' @description Finds the module assignments of given features in a `celda_C()`
#'  model.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Model of class `celda_C`.
#' @param feature Character vector. The module assignemnts will be found for
#'  feature names in this vector.
#' @param exactMatch Logical. Whether an exact match or a partial match using
#'  `grep()` is required to look up the feature in the rownames of the counts
#'  matrix. Default TRUE.
#' @return List. Each element contains the module of the provided feature.
#' @export
setMethod("featureModuleLookup", signature(celdaMod = "celda_C"),
    function(counts, celdaMod, feature, exactMatch) {
        stop("Celda_C models do not contain feature modules.")
    })
