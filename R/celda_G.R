#' @title Feature clustering with Celda
#' @description Clusters the rows of a count matrix containing single-cell data
#'  into L modules.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells.
#' @param L Integer. Number of feature modules.
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to
#'  each feature module in each cell. Default 1.
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to
#'  each feature in each module. Default 1.
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to
#'  the number of features in each module. Default 1.
#' @param stopIter Integer. Number of iterations without improvement in the
#'  log likelihood to stop inference. Default 10.
#' @param maxIter Integer. Maximum number of iterations of Gibbs sampling to
#'  perform. Default 200.
#' @param splitOnIter Integer. On every `splitOnIter` iteration, a heuristic
#'  will be applied to determine if a feature module should be reassigned and
#'  another feature module should be split into two clusters. To disable
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
#' @param yInitialize Chararacter. One of 'random', 'split', or 'predefined'.
#'  With 'random', features are randomly assigned to a modules. With 'split',
#'  features will be split into sqrt(L) modules and then each module will be
#'  subsequently split into another sqrt(L) modules. With 'predefined', values
#'  in `yInit` will be used to initialize `y`. Default 'split'.
#' @param yInit Integer vector. Sets initial starting values of y. If NULL,
#'  starting values for each feature will be randomly sampled from `1:L`.
#'  `yInit` can only be used when `initialize = 'random'`. Default NULL.
#' @param countChecksum Character. An MD5 checksum for the `counts` matrix.
#'  Default NULL.
#' @param logfile Character. Messages will be redirected to a file named
#'  `logfile`. If NULL, messages will be printed to stdout.  Default NULL.
#' @param verbose Logical. Whether to print log messages. Default TRUE.
#' @return An object of class `celda_G` with the feature module clusters stored
#'  in `y`.
#' @seealso `celda_C()` for cell clustering and `celda_CG()` for simultaneous
#'  clustering of features and cells. `celdaGridSearch()` can be used to run
#'  multiple values of L and multiple chains in parallel.
#' @examples
#' data(celdaGSim)
#' celdaMod <- celda_G(celdaGSim$counts, L = celdaGSim$L)
#' @export
celda_G <- function(counts,
    L,
    beta = 1,
    delta = 1,
    gamma = 1,
    stopIter = 10,
    maxIter = 200,
    splitOnIter = 10,
    splitOnLast = TRUE,
    seed = 12345,
    nchains = 3,
    yInitialize = c("split", "random", "predefined"),
    countChecksum = NULL,
    yInit = NULL,
    logfile = NULL,
    verbose = TRUE) {

    .validateCounts(counts)
    if (is.null(seed)) {
        res <- .celda_G(counts,
            L,
            beta,
            delta,
            gamma,
            stopIter,
            maxIter,
            splitOnIter,
            splitOnLast,
            nchains,
            yInitialize,
            countChecksum,
            yInit,
            logfile,
            verbose,
            reorder = TRUE)
    } else {
        with_seed(seed,
            res <- .celda_G(counts,
                L,
                beta,
                delta,
                gamma,
                stopIter,
                maxIter,
                splitOnIter,
                splitOnLast,
                nchains,
                yInitialize,
                countChecksum,
                yInit,
                logfile,
                verbose,
                reorder = TRUE))
    }

    return(res)
}


.celda_G <- function(counts,
    L,
    beta = 1,
    delta = 1,
    gamma = 1,
    stopIter = 10,
    maxIter = 200,
    splitOnIter = 10,
    splitOnLast = TRUE,
    nchains = 3,
    yInitialize = c("split", "random", "predefined"),
    countChecksum = NULL,
    yInit = NULL,
    logfile = NULL,
    verbose = TRUE,
    reorder = TRUE) {

    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = FALSE,
        verbose = verbose)
    .logMessages("Starting Celda_G: Clustering genes.",
        logfile = logfile,
        append = TRUE,
        verbose = verbose)
    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)
    start.time <- Sys.time()

    ## Error checking and variable processing
    counts <- .processCounts(counts)
    if (is.null(countChecksum)) {
        countChecksum <- .createCountChecksum(counts)
    }
    yInitialize <- match.arg(yInitialize)

    allChains <- seq(nchains)

    # Pre-compute lgamma values
    lgbeta <- lgamma(seq(0, max(.colSums(counts,
        nrow(counts), ncol(counts)))) + beta)
    lggamma <- lgamma(seq(0, nrow(counts) + L) + gamma)
    lgdelta <- c(NA, lgamma((seq(nrow(counts) + L) * delta)))

    bestResult <- NULL
    for (i in allChains) {
        ## Randomly select y or y to supplied initial values
        ## Initialize cluster labels
        .logMessages(date(),
            ".. Initializing 'y' in chain",
            i,
            "with",
            paste0("'", yInitialize, "' "),
            logfile = logfile,
            append = TRUE,
            verbose = verbose)

        if (yInitialize == "predefined") {
            if (is.null(yInit)) {
                stop("'yInit' needs to specified when initilize.y == 'given'.")
            }
            y <- .initializeCluster(L,
                nrow(counts),
                initial = yInit,
                fixed = NULL)
        } else if (yInitialize == "split") {
            y <- .initializeSplitY(counts,
                L,
                beta = beta,
                delta = delta,
                gamma = gamma)
        } else {
            y <- .initializeCluster(L,
                    nrow(counts),
                    initial = NULL,
                    fixed = NULL)
        }
        yBest <- y

        ## Calculate counts one time up front
        p <- .cGDecomposeCounts(counts = counts, y = y, L = L)
        nTSByC <- p$nTSByC
        nByG <- p$nByG
        nByTS <- p$nByTS
        nGByTS <- p$nGByTS
        nM <- p$nM
        nG <- p$nG
        rm(p)

        ## Calculate initial log likelihood
        ll <- .cGCalcLL(nTSByC = nTSByC,
                nByTS = nByTS,
                nByG = nByG,
                nGByTS = nGByTS,
                nM = nM,
                nG = nG,
                L = L,
                beta = beta,
                delta = delta,
                gamma = gamma)

        iter <- 1L
        numIterWithoutImprovement <- 0L
        doGeneSplit <- TRUE
        while (iter <= maxIter & numIterWithoutImprovement <= stopIter) {
            nextY <- .cGCalcGibbsProbY(counts = counts,
                nTSByC = nTSByC,
                nByTS = nByTS,
                nGByTS = nGByTS,
                nByG = nByG,
                y = y,
                nG = nG,
                L = L,
                beta = beta,
                delta = delta,
                gamma = gamma,
                lgbeta = lgbeta,
                lggamma = lggamma,
                lgdelta = lgdelta)
            nTSByC <- nextY$nTSByC
            nGByTS <- nextY$nGByTS
            nByTS <- nextY$nByTS
            y <- nextY$y

            ## Perform split on i-th iteration of no improvement in log
            ## likelihood
            tempLl <- .cGCalcLL(nTSByC = nTSByC,
                    nByTS = nByTS,
                    nByG = nByG,
                    nGByTS = nGByTS,
                    nM = nM,
                    nG = nG,
                    L = L,
                    beta = beta,
                    delta = delta,
                    gamma = gamma)
            if (L > 2 & iter != maxIter &
                ((((numIterWithoutImprovement == stopIter &
                    !all(tempLl > ll))) & isTRUE(splitOnLast)) |
                        (splitOnIter > 0 & iter %% splitOnIter == 0 &
                            isTRUE(doGeneSplit)))) {
                .logMessages(date(),
                    " .... Determining if any gene clusters should be split.",
                    logfile = logfile,
                    append = TRUE,
                    sep = "",
                    verbose = verbose)
                res <- .cGSplitY(counts,
                    y,
                    nTSByC,
                    nByTS,
                    nByG,
                    nGByTS,
                    nM,
                    nG,
                    L,
                    beta,
                    delta,
                    gamma,
                    yProb = t(nextY$probs),
                    minFeature = 3,
                    maxClustersToTry = max(L / 2, 10))
                .logMessages(res$message,
                    logfile = logfile,
                    append = TRUE,
                    verbose = verbose)

                # Reset convergence counter if a split occured
                if (!isTRUE(all.equal(y, res$y))) {
                    numIterWithoutImprovement <- 1L
                    doGeneSplit <- TRUE
                } else {
                    doGeneSplit <- FALSE
                }

                ## Re-calculate variables
                y <- res$y
                nTSByC <- res$nTSByC
                nByTS <- res$nByTS
                nGByTS <- res$nGByTS
            }

            ## Calculate complete likelihood
            tempLl <- .cGCalcLL(nTSByC = nTSByC,
                    nByTS = nByTS,
                    nByG = nByG,
                    nGByTS = nGByTS,
                    nM = nM,
                    nG = nG,
                    L = L,
                    beta = beta,
                    delta = delta,
                    gamma = gamma)
            if ((all(tempLl > ll)) | iter == 1) {
                yBest <- y
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

        names <- list(row = rownames(counts), column = colnames(counts))

        result <- list(y = yBest,
            completeLogLik = ll,
            finalLogLik = llBest,
            L = L,
            beta = beta,
            delta = delta,
            gamma = gamma,
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

    bestResult <- methods::new("celda_G",
        clusters = list(y = yBest),
        params = list(L = as.integer(L),
            beta = beta,
            delta = delta,
            gamma = gamma,
            countChecksum = countChecksum),
        completeLogLik = ll,
        finalLogLik = llBest,
        names = names)
    if (isTRUE(reorder)) {
        bestResult <- .reorderCeldaG(counts = counts, res = bestResult)
    }

    endTime <- Sys.time()
    .logMessages(paste0(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)
    .logMessages("Completed Celda_G. Total time:",
        format(difftime(endTime, start.time)),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)
    .logMessages(paste0(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)

    return(bestResult)
}


# Calculate Log Likelihood For Single Set of Cluster Assignments
# (Gene Clustering)
# This function calculates the log-likelihood of a given set of cluster
# assigments for the samples
# represented in the provided count matrix.
# @param nTSByC Number of counts in each Transcriptional State per Cell.
# @param nByTS Number of counts per Transcriptional State.
# @param nGByTS Number of genes in each Transcriptional State.
# @param nG.in.Y  Number of genes in each of the cell cluster.
# @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to
# the number of features in each module. Default 1.
# @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to
# each feature in each module. Default 1.
# @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to
# each feature module in each cell. Default 1.
# @keywords log likelihood
.cGCalcGibbsProbY <- function(counts,
    nTSByC,
    nByTS,
    nGByTS,
    nByG,
    y,
    L,
    nG,
    beta,
    delta,
    gamma,
    lgbeta,
    lggamma,
    lgdelta,
    doSample = TRUE) {

    ## Set variables up front outside of loop
    probs <- matrix(NA, ncol = nG, nrow = L)
    ix <- sample(seq(nG))
    for (i in ix) {
        probs[, i] <- cG_CalcGibbsProbY_fast(
            index = i,
            counts = counts,
            nTSbyC = nTSByC,
            nbyTS = nByTS,
            nGbyTS = nGByTS,
            nbyG = nByG,
            y = y,
            L = L,
            nG = nG,
            lg_beta = lgbeta,
            lg_gamma = lggamma,
            lg_delta = lgdelta,
            delta = delta
        )
        ## Sample next state and add back counts
        if (isTRUE(doSample)) {
            prevY <- y[i]
            y[i] <- .sampleLl(probs[, i])

            if (prevY != y[i]) {
                nTSByC[prevY, ] <- nTSByC[prevY, ] - counts[i, ]
                nGByTS[prevY] <- nGByTS[prevY] - 1L
                nByTS[prevY] <- nByTS[prevY] - nByG[i]

                nTSByC[y[i], ] <- nTSByC[y[i], ] + counts[i, ]
                nGByTS[y[i]] <- nGByTS[y[i]] + 1L
                nByTS[y[i]] <- nByTS[y[i]] + nByG[i]
            }
        }
    }

    return(list(nTSByC = nTSByC,
        nGByTS = nGByTS,
        nByTS = nByTS,
        y = y,
        probs = probs))
}


#' @title Simulate cells from the celda_G model
#' @description Generates a simulated counts matrix and feature module clusters
#'  according to the generative process of the celda_G model.
#' @param model Character. Options available in `celda::availableModels`.
#' @param C Integer. Number of cells to simulate. Default 100.
#' @param L Integer. Number of feature modules. Default 10.
#' @param NRange Integer vector. A vector of length 2 that specifies the lower
#'  and upper bounds of the number of counts generated for each cell. Default
#'  c(500, 5000).
#' @param G Integer. The total number of features to be simulated. Default 100.
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to
#'  each feature module in each cell. Default 1.
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to
#'  each feature in each module. Default 1.
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to
#'  the number of features in each module. Default 5.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @param ... Additional parameters.
#' @return List. Contains the simulated matrix `counts`, feature module clusters
#'  `y`, and input parameters.
#' @seealso `celda_C()` for simulating cell subpopulations and `celda_CG()` for
#'  simulating feature modules and cell populations.
#' @examples
#' celdaGSim <- simulateCells(model = "celda_G")
#' @export
simulateCellscelda_G <- function(model,
    C = 100,
    NRange = c(500, 1000),
    G = 100,
    L = 10,
    beta = 1,
    gamma = 5,
    delta = 1,
    seed = 12345,
    ...) {

    if (is.null(seed)) {
        res <- .simulateCellscelda_G(model = model,
            C = C,
            NRange = NRange,
            G = G,
            L = L,
            beta = beta,
            gamma = gamma,
            delta = delta,
            ...)
    } else {
        with_seed(seed,
            res <- .simulateCellscelda_G(model = model,
                C = C,
                NRange = NRange,
                G = G,
                L = L,
                beta = beta,
                gamma = gamma,
                delta = delta,
                ...))
    }

    return(res)
}


.simulateCellscelda_G <- function(model,
    C = 100,
    NRange = c(500, 1000),
    G = 100,
    L = 10,
    beta = 1,
    gamma = 5,
    delta = 1,
    ...) {

    eta <- .rdirichlet(1, rep(gamma, L))

    y <- sample(seq(L),
        size = G,
        prob = eta,
        replace = TRUE)
    if (length(table(y)) < L) {
        stop("Some states did not receive any genes after sampling. Try",
            " increasing G and/or setting gamma > 1.")
    }

    psi <- matrix(0, nrow = G, ncol = L)
    for (i in seq(L)) {
        ind <- y == i
        psi[ind, i] <- .rdirichlet(1, rep(delta, sum(ind)))
    }

    phi <- .rdirichlet(C, rep(beta, L))

    ## Select number of transcripts per cell
    nN <- sample(seq(NRange[1], NRange[2]), size = C, replace = TRUE)

    ## Select transcript distribution for each cell
    cellCounts <- matrix(0, nrow = G, ncol = C)
    for (i in seq(C)) {
        cellDist <- stats::rmultinom(1, size = nN[i], prob = phi[i, ])
        for (j in seq(L)) {
            cellCounts[, i] <- cellCounts[, i] + stats::rmultinom(1,
                size = cellDist[j], prob = psi[, j])
        }
    }

    ## Ensure that there are no all-0 rows in the counts matrix, which violates
    ## a celda modeling
    ## constraint (columns are guarnteed at least one count):
    zeroRowIdx <- which(rowSums(cellCounts) == 0)
    if (length(zeroRowIdx > 0)) {
        cellCounts <- cellCounts[-zeroRowIdx, ]
        y <- y[-zeroRowIdx]
    }

    rownames(cellCounts) <- paste0("Gene_", seq(nrow(cellCounts)))
    colnames(cellCounts) <- paste0("Cell_", seq(ncol(cellCounts)))

    ## Peform reordering on final Z and Y assigments:
    cellCounts <- .processCounts(cellCounts)
    names <- list(row = rownames(cellCounts),
        column = colnames(cellCounts))
    countChecksum <- .createCountChecksum(cellCounts)
    result <- methods::new("celda_G",
        clusters = list(y = y),
        params = list(L = as.integer(L),
            beta = beta,
            delta = delta,
            gamma = gamma,
            countChecksum = countChecksum),
        names = names
    )
    result <- .reorderCeldaG(counts = cellCounts, res = result)

    return(list(y = clusters(result)$y,
        counts = .processCounts(cellCounts),
        L = L,
        beta = beta,
        delta = delta,
        gamma = gamma))
}


#' @title Matrix factorization for results from celda_G
#' @description Generates factorized matrices showing the contribution of each
#'  feature in each module and each module in each cell.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda object of class "celda_G".
#' @param type Character vector. A vector containing one or more of "counts",
#'  "proportion", or "posterior". "counts" returns the raw number of counts for
#'  each factorized matrix. "proportions" returns the normalized probabilities
#'  for each factorized matrix, which are calculated by dividing the raw counts
#'  in each factorized matrix by the total counts in each column. "posterior"
#'  returns the posterior estimates. Default
#'  `c("counts", "proportion", "posterior")`.
#' @return A list with elements for `counts`, `proportions`, or `posterior`
#'  probabilities. Each element will be a list containing factorized matrices
#'  for `module` and `cell`.
#' @seealso `celda_G()` for clustering features
#' @examples
#' data(celdaGSim, celdaGMod)
#' factorizedMatrices <- factorizeMatrix(celdaGSim$counts,
#'     celdaGMod, "posterior")
#' @export
setMethod("factorizeMatrix", signature(celdaMod = "celda_G"),
    function(counts,
        celdaMod,
        type = c("counts", "proportion", "posterior")) {

        counts <- .processCounts(counts)
        # compareCountMatrix(counts, celdaMod)

        L <- params(celdaMod)$L
        y <- clusters(celdaMod)$y
        beta <- params(celdaMod)$beta
        delta <- params(celdaMod)$delta
        gamma <- params(celdaMod)$gamma

        p <- .cGDecomposeCounts(counts = counts, y = y, L = L)
        nTSByC <- p$nTSByC
        nByG <- p$nByG
        nByTS <- p$nByTS
        nGByTS <- p$nGByTS
        nM <- p$nM
        nG <- p$nG
        rm(p)

        nGByTS[nGByTS == 0] <- 1
        nGByTS <- matrix(0, nrow = length(y), ncol = L)
        nGByTS[cbind(seq(nG), y)] <- nByG

        LNames <- paste0("L", seq(L))
        colnames(nTSByC) <- matrixNames(celdaMod)$column
        rownames(nTSByC) <- LNames
        colnames(nGByTS) <- LNames
        rownames(nGByTS) <- matrixNames(celdaMod)$row
        names(nGByTS) <- LNames

        countsList <- c()
        propList <- c()
        postList <- c()
        res <- list()

        if (any("counts" %in% type)) {
            countsList <- list(cell = nTSByC,
                module = nGByTS,
                geneDistribution = nGByTS)
            res <- c(res, list(counts = countsList))
        }

        if (any("proportion" %in% type)) {
            ## Need to avoid normalizing cell/gene states with zero cells/genes
            uniqueY <- sort(unique(y))
            tempNGByTS <- nGByTS
            tempNGByTS[, uniqueY] <- normalizeCounts(tempNGByTS[, uniqueY],
                normalize = "proportion")
            tempNGByTS <- nGByTS / sum(nGByTS)

            propList <- list(cell = normalizeCounts(nTSByC,
                normalize = "proportion"),
                module = tempNGByTS,
                geneDistribution = tempNGByTS)
            res <- c(res, list(proportions = propList))
        }

        if (any("posterior" %in% type)) {
            gs <- nGByTS
            gs[cbind(seq(nG), y)] <- gs[cbind(seq(nG), y)] + delta
            gs <- normalizeCounts(gs, normalize = "proportion")
            tempNGByTS <- (nGByTS + gamma) / sum(nGByTS + gamma)

            postList <- list(cell = normalizeCounts(nTSByC + beta,
                normalize = "proportion"),
                module = gs,
                geneDistribution = tempNGByTS)
            res <- c(res, posterior = list(postList))
        }

        return(res)
    })


# Calculate log-likelihood of celda_CG model
.cGCalcLL <- function(nTSByC,
    nByTS,
    nByG,
    nGByTS,
    nM,
    nG,
    L,
    beta,
    delta,
    gamma) {

    nG <- sum(nGByTS)

    ## Calculate for "Phi" component
    a <- nM * lgamma(L * beta)
    b <- sum(lgamma(nTSByC + beta))
    c <- -nM * L * lgamma(beta)
    d <- -sum(lgamma(colSums(nTSByC + beta)))

    phiLl <- a + b + c + d

    ## Calculate for "Psi" component
    a <- sum(lgamma(nGByTS * delta))
    b <- sum(lgamma(nByG + delta))
    c <- -nG * lgamma(delta)
    d <- -sum(lgamma(nByTS + (nGByTS * delta)))

    psiLl <- a + b + c + d

    ## Calculate for "Eta" component
    a <- lgamma(L * gamma)
    b <- sum(lgamma(nGByTS + gamma))
    c <- -L * lgamma(gamma)
    d <- -sum(lgamma(sum(nGByTS + gamma)))

    etaLl <- a + b + c + d

    final <- phiLl + psiLl + etaLl
    return(final)
}


#' @title Calculate Celda_G log likelihood
#' @description Calculates the log likelihood for user-provided feature module
#'  clusters using the `celda_G()` model.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells.
#' @param y Numeric vector. Denotes feature module labels.
#' @param L Integer. Number of feature modules.
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to
#'  each feature module in each cell. Default 1.
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to
#'  each feature in each module. Default 1.
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to
#'  the number of features in each module. Default 1.
#' @keywords log likelihood
#' @return The log-likelihood for the given cluster assignments.
#' @seealso `celda_G()` for clustering features
#' @examples
#' data(celdaGSim)
#' loglik <- logLikelihoodcelda_G(celdaGSim$counts,
#'     y = celdaGSim$y,
#'     L = celdaGSim$L,
#'     beta = celdaGSim$beta,
#'     delta = celdaGSim$delta,
#'     gamma = celdaGSim$gamma)
#'
#' loglik <- logLikelihood(celdaGSim$counts,
#'     model = "celda_G",
#'     y = celdaGSim$y,
#'     L = celdaGSim$L,
#'     beta = celdaGSim$beta,
#'     delta = celdaGSim$delta,
#'     gamma = celdaGSim$gamma)
#' @export
logLikelihoodcelda_G <- function(counts, y, L, beta, delta, gamma) {
    if (sum(y > L) > 0) {
        stop("An entry in y contains a value greater than the provided L.")
    }
    p <- .cGDecomposeCounts(counts = counts, y = y, L = L)
    final <- .cGCalcLL(nTSByC = p$nTSByC,
        nByTS = p$nByTS,
        nByG = p$nByG,
        nGByTS = p$nGByTS,
        nM = p$nM,
        nG = p$nG,
        L = L,
        beta = beta,
        delta = delta,
        gamma = gamma)

    return(final)
}


# Takes raw counts matrix and converts it to a series of matrices needed for
# log likelihood calculation
# @param counts Integer matrix. Rows represent features and columns represent
# cells.
# @param y Numeric vector. Denotes feature module labels.
# @param L Integer. Number of feature modules.
.cGDecomposeCounts <- function(counts, y, L) {
    if (any(y > L)) {
        stop("Entries in the module clusters 'y' are greater than L.")
    }
    nTSByC <- .rowSumByGroup(counts, group = y, L = L)
    nByG <- as.integer(.rowSums(counts, nrow(counts), ncol(counts)))
    nByTS <- as.integer(.rowSumByGroup(matrix(nByG, ncol = 1),
        group = y, L = L))
    nGByTS <- tabulate(y, L) + 1 ## Add pseudogene to each state
    nM <- ncol(counts)
    nG <- nrow(counts)

    return(list(nTSByC = nTSByC,
        nByG = nByG,
        nByTS = nByTS,
        nGByTS = nGByTS,
        nM = nM,
        nG = nG))
}


.cGReDecomposeCounts <- function(counts, y, previousY, nTSByC, nByG, L) {
    ## Recalculate counts based on new label
    nTSByC <- .rowSumByGroupChange(counts, nTSByC, y, previousY, L)
    nByTS <- as.integer(.rowSumByGroup(matrix(nByG, ncol = 1),
        group = y, L = L))
    nGByTS <- tabulate(y, L) + 1

    return(list(nTSByC = nTSByC,
        nByTS = nByTS,
        nGByTS = nGByTS))
}


#' @title Conditional probabilities for features in modules from a Celda_G model
#' @description Calculates the conditional probability of each feature belonging
#'  to each module given all other feature cluster assignments in a `celda_G()`
#'  result.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda object of class `celda_G`.
#' @param log Logical. If FALSE, then the normalized conditional probabilities
#'  will be returned. If TRUE, then the unnormalized log probabilities will be
#'  returned. Default FALSE.
#' @param ... Additional parameters.
#' @return A list containging a matrix for the conditional cell cluster
#'  probabilities.
#' @seealso `celda_G()` for clustering features
#' @examples
#' data(celdaGSim, celdaGMod)
#' clusterProb <- clusterProbability(celdaGSim$counts, celdaGMod)
#' @export
setMethod("clusterProbability", signature(celdaMod = "celda_G"),
    function(counts, celdaMod, log = FALSE, ...) {
        y <- clusters(celdaMod)$y
        L <- params(celdaMod)$L
        delta <- params(celdaMod)$delta
        beta <- params(celdaMod)$beta
        gamma <- params(celdaMod)$gamma

        ## Calculate counts one time up front
        p <- .cGDecomposeCounts(counts = counts, y = y, L = L)
        lgbeta <- lgamma(seq(0, max(.colSums(counts,
            nrow(counts), ncol(counts)))) + beta)
        lggamma <- lgamma(seq(0, nrow(counts) + L) + gamma)
        lgdelta <- c(NA, lgamma(seq(nrow(counts) + L) * delta))

        nextY <- .cGCalcGibbsProbY(counts = counts,
            nTSByC = p$nTSByC,
            nByTS = p$nByTS,
            nGByTS = p$nGByTS,
            nByG = p$nByG,
            y = y,
            nG = p$nG,
            L = L,
            lgbeta = lgbeta,
            lgdelta = lgdelta,
            lggamma = lggamma,
            delta = delta,
            doSample = FALSE)
        yProb <- t(nextY$probs)

        if (!isTRUE(log)) {
            yProb <- .normalizeLogProbs(yProb)
        }

        return(list(yProbability = yProb))
    })


#' @title Calculate the perplexity on new data with a celda_G model
#' @description Perplexity is a statistical measure of how well a probability
#'  model can predict new data. Lower perplexity indicates a better model.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda object of class "celda_C"
#' @param newCounts A new counts matrix used to calculate perplexity. If NULL,
#'  perplexity will be calculated for the 'counts' matrix. Default NULL.
#' @return Numeric. The perplexity for the provided count data and model.
#' @seealso `celda_G()` for clustering features
#' @examples
#' data(celdaGSim, celdaGMod)
#' perplexity <- perplexity(celdaGSim$counts, celdaGMod)
#' @export
setMethod("perplexity", signature(celdaMod = "celda_G"),
    function(counts, celdaMod, newCounts = NULL) {
        counts <- .processCounts(counts)
        # compareCountMatrix(counts, celdaMod)

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
            type = c("posterior", "counts"))
        psi <- factorized$posterior$module
        phi <- factorized$posterior$cell
        eta <- factorized$posterior$geneDistribution
        nGByTS <- factorized$counts$geneDistribution

        etaProb <- log(eta) * nGByTS
        # gene.by.cell.prob = log(psi %*% phi)
        # logPx = sum(gene.by.cell.prob * newCounts) # + sum(etaProb)
        logPx <- .perplexityGLogPx(newCounts,
            phi,
            psi,
            clusters(celdaMod)$y,
            params(celdaMod)$L) # + sum(etaProb)
        perplexity <- exp(- (logPx / sum(newCounts)))
        return(perplexity)
    })


.reorderCeldaG <- function(counts, res) {
    if (params(res)$L > 2 & isTRUE(length(unique(clusters(res)$y)) > 1)) {
        res@clusters$y <- as.integer(as.factor(clusters(res)$y))
        fm <- factorizeMatrix(counts = counts, celdaMod = res)
        uniqueY <- sort(unique(clusters(res)$y))
        cs <- prop.table(t(fm$posterior$cell[uniqueY, ]), 2)
        d <- .cosineDist(cs)
        h <- stats::hclust(d, method = "complete")
        res <- recodeClusterY(res,
            from = h$order,
            to = seq(length(h$order)))
    }
    return(res)
}


#' @title Heatmap for celda_CG
#' @description Renders an expression heatmap to visualize `celda_CG()` results.
#'  The top `nfeatures` for each module will be included in the heatmap.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda object of class `celda_G`.
#' @param nfeatures Integer. Maximum number of features to select for each
#'  module. Default 25.
#' @param ... Additional parameters.
#' @seealso `celda_G()` for clustering features and `celdaTsne()` for
#'  generating 2-dimensional coordinates.
#' @examples
#' data(celdaGSim, celdaGMod)
#' celdaHeatmap(celdaGSim$counts, celdaGMod)
#' @return list A list containing the dendrograms and the heatmap grob.
#' @export
setMethod("celdaHeatmap", signature(celdaMod = "celda_G"),
    function(counts, celdaMod, nfeatures = 25, ...) {
        fm <- factorizeMatrix(counts, celdaMod, type = "proportion")
        top <- topRank(fm$proportions$module, n = nfeatures)
        ix <- unlist(top$index)
        norm <- normalizeCounts(counts,
            normalize = "proportion",
            transformationFun = sqrt)
        plotHeatmap(norm[ix, ], y = clusters(celdaMod)$y[ix], ...)
    })

#' @title tSNE for celda_G
#' @description Embeds cells in two dimensions using tSNE based on a `celda_G`
#'  model. tSNE is run on module probabilities to reduce the number of features
#'  instead of using PCA. Module probabilities square-root trasformed before
#'  applying tSNE.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda object of class `celda_G`.
#' @param maxCells Integer. Maximum number of cells to plot. Cells will be
#'  randomly subsampled if ncol(conts) > maxCells. Larger numbers of cells
#'  requires more memory. If NULL, no subsampling will be performed.
#'  Default NULL.
#' @param minClusterSize Integer. Do not subsample cell clusters below this
#'  threshold. Default 100.
#' @param initialDims Integer. PCA will be used to reduce the dimentionality of
#'  the dataset. The top 'initialDims' principal components will be used for
#'  tSNE. Default 20.
#' @param modules Integer vector. Determines which feature modules to use for
#'  tSNE. If NULL, all modules will be used. Default NULL.
#' @param perplexity Numeric. Perplexity parameter for tSNE. Default 20.
#' @param maxIter Integer. Maximum number of iterations in tSNE generation.
#'  Default 2500.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @seealso `celda_G()` for clustering features and `celdaHeatmap()` for
#'  displaying expression
#' @examples
#' data(celdaGSim, celdaGMod)
#' tsneRes <- celdaTsne(celdaGSim$counts, celdaGMod)
#' @return A two column matrix of t-SNE coordinates.
#' @export
setMethod("celdaTsne", signature(celdaMod = "celda_G"),
    function(counts,
        celdaMod,
        maxCells = NULL,
        minClusterSize = 100,
        initialDims = 20,
        modules = NULL,
        perplexity = 20,
        maxIter = 2500,
        seed = 12345) {

        if (is.null(seed)) {
            res <- .celdaTsneG(counts = counts,
                celdaMod = celdaMod,
                maxCells = maxCells,
                minClusterSize = minClusterSize,
                initialDims = initialDims,
                modules = modules,
                perplexity = perplexity,
                maxIter = maxIter)
        } else {
            with_seed(seed,
                res <- .celdaTsneG(counts = counts,
                    celdaMod = celdaMod,
                    maxCells = maxCells,
                    minClusterSize = minClusterSize,
                    initialDims = initialDims,
                    modules = modules,
                    perplexity = perplexity,
                    maxIter = maxIter))
        }

        return(res)

    })


.celdaTsneG <- function(counts,
    celdaMod,
    maxCells = NULL,
    minClusterSize = 100,
    initialDims = 20,
    modules = NULL,
    perplexity = 20,
    maxIter = 2500) {

    preparedCountInfo <- .prepareCountsForDimReductionCeldaCG(counts,
        celdaMod,
        maxCells,
        minClusterSize,
        modules)
    res <- .calculateTsne(preparedCountInfo$norm,
        doPca = FALSE,
        perplexity = perplexity,
        maxIter = maxIter)
    final <- matrix(NA, nrow = ncol(counts), ncol = 2)
    final[preparedCountInfo$cellIx, ] <- res
    rownames(final) <- colnames(counts)
    colnames(final) <- c("tSNE_1", "tSNE_2")
    return(final)
}


#' @title umap for celda_G
#' @description Embeds cells in two dimensions using umap based on a `celda_G`
#'  model. umap is run on module probabilities to reduce the number of features
#'  instead of using PCA. Module probabilities square-root trasformed before
#'  applying tSNE.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda object of class `celda_CG`.
#' @param maxCells Integer. Maximum number of cells to plot. Cells will be
#'  randomly subsampled if ncol(counts) > maxCells. Larger numbers of cells
#'  requires more memory. If NULL, no subsampling will be performed.
#'  Default NULL.
#' @param minClusterSize Integer. Do not subsample cell clusters below this
#'  threshold. Default 100.
#' @param modules Integer vector. Determines which features modules to use for
#'  UMAP. If NULL, all modules will be used. Default NULL.
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
#'   embedded points are. Default 1.
#'   See `?uwot::umap` for more information.
#' @param cores Number of threads to use. Default 1.
#' @param ... Other parameters to pass to `uwot::umap`.
#' @seealso `celda_G()` for clustering features and cells  and `celdaHeatmap()`
#'  for displaying expression
#' @examples
#' data(celdaGSim, celdaGMod)
#' umapRes <- celdaUmap(celdaGSim$counts, celdaGMod)
#' @return A two column matrix of umap coordinates
#' @export
setMethod("celdaUmap", signature(celdaMod = "celda_G"),
    function(counts,
        celdaMod,
        maxCells = NULL,
        minClusterSize = 100,
        modules = NULL,
        seed = 12345,
        nNeighbors = 30,
        minDist = 0.2,
        spread = 1,
        cores = 1,
        ...) {

        if (is.null(seed)) {
            res <- .celdaUmapG(counts = counts,
                celdaMod = celdaMod,
                maxCells = maxCells,
                minClusterSize = minClusterSize,
                modules = modules,
                nNeighbors = nNeighbors,
                minDist = minDist,
                spread = spread,
                cores = cores,
                ...)
        } else {
            with_seed(seed,
                res <- .celdaUmapG(counts = counts,
                    celdaMod = celdaMod,
                    maxCells = maxCells,
                    minClusterSize = minClusterSize,
                    modules = modules,
                    nNeighbors = nNeighbors,
                    minDist = minDist,
                    spread = spread,
                    cores = cores,
                    ...))
        }

        return(res)
    })


.celdaUmapG <- function(counts,
    celdaMod,
    maxCells = NULL,
    minClusterSize = 100,
    modules = NULL,
    nNeighbors = nNeighbors,
    minDist = minDist,
    spread = spread,
    cores = cores,
    ...) {

    preparedCountInfo <- .prepareCountsForDimReductionCeldaCG(counts,
        celdaMod,
        maxCells,
        minClusterSize,
        modules)
    umapRes <- .calculateUmap(preparedCountInfo$norm,
        nNeighbors = nNeighbors,
        minDist = minDist,
        spread = spread,
        cores = cores,
        ...)

    final <- matrix(NA, nrow = ncol(counts), ncol = 2)
    final[preparedCountInfo$cellIx, ] <- umapRes
    rownames(final) <- colnames(counts)
    colnames(final) <- c("UMAP_1", "UMAP_2")
    return(final)
}


.prepareCountsForDimReductionCeldaCG <- function(counts,
    celdaMod,
    maxCells = NULL,
    minClusterSize = 100,
    modules = NULL) {

    if (is.null(maxCells) || maxCells > ncol(counts)) {
      maxCells <- ncol(counts)
      cellIx <- seq_len(ncol(counts))
    } else {
      cellIx <- sample(seq(ncol(counts)), maxCells)
    }

    fm <- factorizeMatrix(counts = counts,
        celdaMod = celdaMod,
        type = "counts")

    modulesToUse <- seq(nrow(fm$counts$cell))
    if (!is.null(modules)) {
        if (!all(modules %in% modulesToUse)) {
            stop("'modules' must be a vector of numbers between 1 and ",
                modulesToUse,
                ".")
        }
        modulesToUse <- modules
    }

    norm <- t(normalizeCounts(fm$counts$cell[modulesToUse, cellIx],
        normalize = "proportion",
        transformationFun = sqrt))
    return(list(norm = norm, cellIx = cellIx))
}


#' @title Lookup the module of a feature
#' @description Finds the module assignments of given features in a `celda_G()`
#'  model.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Model of class `celda_G`.
#' @param feature Character vector. The module assignemnts will be found for
#'  feature names in this vector.
#' @param exactMatch Logical. Whether an exact match or a partial match using
#'  `grep()` is used to look up the feature in the rownames of the counts
#'  matrix. Default TRUE.
#' @return List. Each element contains the module of the provided feature.
#' @seealso `celda_G()` for clustering features
#' @examples
#' data(celdaGSim, celdaGMod)
#' module <- featureModuleLookup(celdaGSim$counts,
#'     celdaGMod,
#'     c("Gene_1", "Gene_XXX"))
#' @export
setMethod("featureModuleLookup", signature(celdaMod = "celda_G"),
    function(counts, celdaMod, feature, exactMatch = TRUE) {
        if (!isTRUE(exactMatch)) {
            feature <- unlist(lapply(seq(length(feature)),
                function(x) {
                    rownames(counts)[grep(feature[x], rownames(counts))]
                }))
        }

        featList <- lapply(seq(length(feature)),
            function(x) {
                if (feature[x] %in% rownames(counts)) {
                    return(clusters(celdaMod)$y[which(rownames(counts) ==
                            feature[x])])
                } else {
                    return(paste0(
                        "No feature was identified matching '",
                        feature[x],
                        "'."
                    ))
                }
            })
        names(featList) <- feature
        return(featList)
    })
