.celda_C_mm <- function(counts,
    counts2,
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
    zInitialize = c("random", "predefined"),
    countChecksum = NULL,
    zInit = NULL,
    logfile = NULL,
    verbose = TRUE) {

    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = FALSE,
        verbose = verbose)

    .logMessages("Starting multi-modal Celda_C: Clustering cells.",
        logfile = logfile,
        append = TRUE,
        verbose = verbose)

    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)

    startTime <- Sys.time()

    if (ncol(counts) != ncol(counts2)) {
        stop("'counts2' should have the same number of columns as 'counts'!")
    }

    if (isTRUE(any(counts != counts2))) {
        stop("'counts2' should have the same column names as 'counts'!")
    }

    ## Error checking and variable processing
    counts <- .processCounts(counts)
    counts2 <- .processCounts(counts2)
    if (is.null(countChecksum)) {
        countChecksum <- c(.createCountChecksum(counts),
            .createCountChecksum(counts2))
    }

    sampleLabel <- .processSampleLabels(sampleLabel, ncol(counts))
    s <- as.integer(sampleLabel)

    algorithm <- match.arg(algorithm)
    if (algorithm == "EM") {
        stopIter <- 1
    }

    if (algorithm == "Gibbs") {
        stop("Gibbs sampling has not been implemented for multi-modal",
        " clustering yet. Use 'EM' instead.")
    } else if (algorithm == "EM") {
        algorithmFun <- ".cCCalcEMProbZMM"
    } else {
        stop("Invalid input '", algorithm, "' for 'algorithm'")
    }

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
                stop("'zInit' needs to specified when zInitialize ==",
                    " 'predefined'.")
            }
            z <- .initializeCluster(K,
                ncol(counts),
                initial = zInit,
                fixed = NULL)
            } else {
            z <- .initializeCluster(K,
                ncol(counts),
                initial = NULL,
                fixed = NULL)
        }

        zBest <- z

        ## Calculate counts one time up front
        p <- .cCDecomposeCountsMM(counts, counts2, s, z, K)
        nS <- p$nS
        nG <- p$nG
        nM <- p$nM
        mCPByS <- p$mCPByS
        nGByCP <- p$nGByCP
        nCP <- p$nCP
        nByC <- p$nByC
        nGByCP2 <- p$nGByCP2
        nCP2 <- p$nCP2
        nByC2 <- p$nByC2
        nG2 <- p$nG2

        ll <- .cCCalcLLMM(mCPByS = mCPByS,
            nGByCP = nGByCP,
            nGByCP2 = nGByCP2,
            s = s,
            K = K,
            nS = nS,
            nG = nG,
            nG2 = nG2,
            alpha = alpha,
            beta = beta)

        iter <- 1L
        numIterWithoutImprovement <- 0L
        doCellSplit <- TRUE
        while (iter <= maxIter & numIterWithoutImprovement <= stopIter) {
            #.cCCalcEMProbZMM
            nextZ <- do.call(algorithmFun, list(
                counts = counts,
                counts2 = counts2,
                mCPByS = mCPByS,
                nGByCP = nGByCP,
                nByC = nByC,
                nCP = nCP,
                nGByCP2 = nGByCP2,
                nByC2 = nByC2,
                nCP2 = nCP2,
                z = z,
                s = s,
                K = K,
                nG = nG,
                nG2 = nG2,
                nM = nM,
                alpha = alpha,
                beta = beta))

            mCPByS <- nextZ$mCPByS
            nGByCP <- nextZ$nGByCP
            nCP <- nextZ$nCP
            nGByCP2 <- nextZ$nGByCP2
            nCP2 <- nextZ$nCP2

            z <- nextZ$z

            ## Perform split on i-th iteration of no improvement in log
            ## likelihood
            tempLl <- .cCCalcLLMM(mCPByS = mCPByS,
                nGByCP = nGByCP,
                nGByCP2 = nGByCP2,
                s = s,
                K = K,
                nS = nS,
                nG = nG,
                nG2 = nG2,
                alpha = alpha,
                beta = beta)

            if (K > 2 & iter != maxIter &
                    (((numIterWithoutImprovement == stopIter &
                            !all(tempLl > ll)) & isTRUE(splitOnLast)) |
                            (splitOnIter > 0 & iter %% splitOnIter == 0 &
                                    isTRUE(doCellSplit)))) {

                .logMessages(date(),
                    " .... Determining if any cell clusters should be split.",
                    logfile = logfile,
                    append = TRUE,
                    sep = "",
                    verbose = verbose)

                # the splitting of clusters do not consider protein counts
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
            tempLl <- .cCCalcLLMM(mCPByS = mCPByS,
                nGByCP = nGByCP,
                nGByCP2 = nGByCP2,
                s = s,
                K = K,
                nS = nS,
                nG = nG,
                nG2 = nG2,
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

    # if (isTRUE(reorder)) {
    #     bestResult <- .reorderCelda_C(counts = counts, res = bestResult)
    # }

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


.cCCalcEMProbZMM <- function(counts,
    counts2,
    mCPByS,
    nGByCP,
    nByC,
    nCP,
    nGByCP2,
    nByC2,
    nCP2,
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
    phi2 <- fastNormPropLog(nGByCP2, beta)

    ## Maximization to find best label for each cell
    probs <- eigenMatMultInt(phi, counts) +
        eigenMatMultInt(phi2, counts2) + theta[, s]

    if (isTRUE(doSample)) {
        zPrevious <- z
        z <- apply(probs, 2, which.max)

        ## Recalculate counts based on new label
        p <- .cCReDecomposeCountsMM(counts,
            counts2, s, z, zPrevious, nGByCP, nGByCP2, K)
        mCPByS <- p$mCPByS
        nGByCP <- p$nGByCP
        nCP <- p$nCP
        nGByCP2 <- p$nGByCP2
        nCP2 <- nCP2
    }

    return(list(mCPByS = mCPByS,
        nGByCP = nGByCP,
        nCP = nCP,
        z = z,
        probs = probs))
}


# Calculate log-likelihood for celda_C model
.cCCalcLLMM <- function(mCPByS,
    nGByCP,
    nGByCP2,
    s,
    z,
    K,
    nS,
    nG,
    nG2,
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

    ## Calculate for "Phi2" component
    a <- K * lgamma(nG2 * beta)
    b <- sum(lgamma(nGByCP2 + beta))
    c <- -K * nG2 * lgamma(beta)
    d <- -sum(lgamma(colSums(nGByCP2 + beta)))

    phi2Ll <- a + b + c + d

    final <- thetaLl + phiLl + phi2Ll
    return(final)
}


# Takes raw counts matrix and converts it to a series of matrices needed for
# log likelihood calculation
# @param counts Integer matrix. Rows represent features and columns represent
# cells.
# @param counts2 Integer matrix. It's needed only for the multi-modal celda
#  clustering algorithm. This count matrix represents the additional
#  quantification of the other layer for multi-modality data (i.e. CITE-seq).
#  For example, for CITE-seq data, the protein counts should be provided here.
# @param s Integer vector. Contains the sample label for each cell (column) in
# the count matrix.
# @param z Numeric vector. Denotes cell population labels.
# @param K Integer. Number of cell populations.
.cCDecomposeCountsMM <- function(counts, counts2, s, z, K) {
    nS <- length(unique(s))
    nG <- nrow(counts)
    nM <- ncol(counts)
    nG2 <- nrow(counts2)

    mCPByS <- matrix(as.integer(table(factor(z, levels = seq(K)), s)),
        ncol = nS)
    nGByCP <- .colSumByGroup(counts, group = z, K = K)
    nCP <- as.integer(colSums(nGByCP))
    nByC <- as.integer(colSums(counts))

    nGByCP2 <- .colSumByGroup(counts2, group = z, K = K)
    nCP2 <- as.integer(colSums(nGByCP2))
    nbyC2 <- as.integer(colSums(counts2))

    return(list(mCPByS = mCPByS,
        nGByCP = nGByCP,
        nCP = nCP,
        nByC = nByC,
        nGByCP2 = nGByCP2,
        nCP2 = nCP2,
        nbyC2 = nbyC2,
        nS = nS,
        nG = nG,
        nG2 = nG2,
        nM = nM))
}


.cCReDecomposeCountsMM <- function(counts,
    counts2, s, z, previousZ, nGByCP, nGByCP2, K) {

    ## Recalculate counts based on new label
    nGByCP <- .colSumByGroupChange(counts, nGByCP, z, previousZ, K)
    nS <- length(unique(s))
    mCPByS <- matrix(as.integer(table(factor(z, levels = seq(K)), s)),
        ncol = nS)
    nCP <- as.integer(colSums(nGByCP))

    nGByCP2 <- .colSumByGroupChange(counts2, nGByCP2, z, previousZ, K)
    nCP2 <- as.integer(colSums(nGByCP2))

    return(list(mCPByS = mCPByS,
        nGByCP = nGByCP,
        nCP = nCP,
        nGByCP2 = nGByCP2,
        nCP2 = nCP2))
}
