#' @title Cell and feature clustering with Celda
#' @description Clusters the rows and columns of a count matrix containing
#'  single-cell data into L modules and K subpopulations, respectively.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells.
#' @param sampleLabel Vector or factor. Denotes the sample label for each cell
#'  (column) in the count matrix.
#' @param K Integer. Number of cell populations.
#' @param L Integer. Number of feature modules.
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount
#'  to each cell population in each sample. Default 1.
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to
#'  each feature module in each cell population. Default 1.
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to
#'  each feature in each module. Default 1.
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to
#'  the number of features in each module. Default 1.
#' @param algorithm String. Algorithm to use for clustering cell subpopulations.
#'  One of 'EM' or 'Gibbs'. The EM algorithm for cell clustering is faster,
#'  especially for larger numbers of cells. However, more chains may be required
#'  to ensure a good solution is found. Default 'EM'.
#' @param stopIter Integer. Number of iterations without improvement in the log
#'  likelihood to stop inference. Default 10.
#' @param maxIter Integer. Maximum number of iterations of Gibbs sampling to
#'  perform. Default 200.
#' @param splitOnIter Integer. On every `splitOnIter` iteration, a heuristic
#'  will be applied to determine if a cell population or feature module should
#'  be reassigned and another cell population or feature module should be split
#'  into two clusters. To disable splitting, set to -1. Default 10.
#' @param splitOnLast Integer. After `stopIter` iterations have been
#'  performed without improvement, a heuristic will be applied to determine if
#'  a cell population or feature module should be reassigned and another cell
#'  population or feature module should be split into two clusters. If a split
#'  occurs, then 'stopIter' will be reset. Default TRUE.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @param nchains Integer. Number of random cluster initializations. Default 3.
#' @param zInitialize Chararacter. One of 'random', 'split', or 'predefined'.
#'  With 'random', cells are randomly assigned to a populations. With 'split',
#'  cells will be split into sqrt(K) populations and then each popluation will
#'  be subsequently split into another sqrt(K) populations. With 'predefined',
#'  values in `zInit` will be used to initialize `z`. Default 'split'.
#' @param yInitialize Chararacter. One of 'random', 'split', or 'predefined'.
#'  With 'random', features are randomly assigned to a modules. With 'split',
#'  features will be split into sqrt(L) modules and then each module will be
#'  subsequently split into another sqrt(L) modules. With 'predefined', values
#'  in `yInit` will be used to initialize `y`. Default 'split'.
#' @param zInit Integer vector. Sets initial starting values of z. If NULL,
#'  starting values for each cell will be randomly sampled from 1:K. 'zInit'
#'  can only be used when `initialize' = 'random'`. Default NULL.
#' @param yInit Integer vector. Sets initial starting values of y. If NULL,
#'  starting values for each feature will be randomly sampled from 1:L.
#'  'yInit' can only be used when `initialize = 'random'`. Default NULL.
#' @param countChecksum Character. An MD5 checksum for the `counts` matrix.
#'  Default NULL.
#' @param logfile Character. Messages will be redirected to a file named
#'  `logfile`. If NULL, messages will be printed to stdout.  Default NULL.
#' @param verbose Logical. Whether to print log messages. Default TRUE.
#' @return An object of class `celda_CG` with the cell populations clusters
#'  stored in `z` and feature module clusters stored in `y`.
#' @seealso `celda_G()` for feature clustering and `celda_C()` for clustering
#'  cells. `celdaGridSearch()` can be used to run multiple values of K/L and
#'  multiple chains in parallel.
#' @examples
#' data(celdaCGSim)
#' celdaMod <- celda_CG(celdaCGSim$counts,
#'     K = celdaCGSim$K,
#'     L = celdaCGSim$L,
#'     sampleLabel = celdaCGSim$sampleLabel,
#'     nchains = 1)
#' @import Rcpp RcppEigen
#' @rawNamespace import(gridExtra, except = c(combine))
#' @export
celda_CG <- function(counts,
    sampleLabel = NULL,
    K,
    L,
    alpha = 1,
    beta = 1,
    delta = 1,
    gamma = 1,
    algorithm = c("EM", "Gibbs"),
    stopIter = 10,
    maxIter = 200,
    splitOnIter = 10,
    splitOnLast = TRUE,
    seed = 12345,
    nchains = 3,
    zInitialize = c("split", "random", "predefined"),
    yInitialize = c("split", "random", "predefined"),
    countChecksum = NULL,
    zInit = NULL,
    yInit = NULL,
    logfile = NULL,
    verbose = TRUE) {

    .validateCounts(counts)
    if (is.null(seed)) {
        res <- .celda_CG(
            counts,
            sampleLabel,
            K,
            L,
            alpha,
            beta,
            delta,
            gamma,
            algorithm,
            stopIter,
            maxIter,
            splitOnIter,
            splitOnLast,
            nchains,
            zInitialize,
            yInitialize,
            countChecksum,
            zInit,
            yInit,
            logfile,
            verbose,
            reorder = TRUE)
    } else {
        with_seed(seed,
            res <- .celda_CG(
                counts,
                sampleLabel,
                K,
                L,
                alpha,
                beta,
                delta,
                gamma,
                algorithm,
                stopIter,
                maxIter,
                splitOnIter,
                splitOnLast,
                nchains,
                zInitialize,
                yInitialize,
                countChecksum,
                zInit,
                yInit,
                logfile,
                verbose,
                reorder = TRUE))
    }

    return(res)
}

.celda_CG <- function(counts,
    sampleLabel = NULL,
    K,
    L,
    alpha = 1,
    beta = 1,
    delta = 1,
    gamma = 1,
    algorithm = c("EM", "Gibbs"),
    stopIter = 10,
    maxIter = 200,
    splitOnIter = 10,
    splitOnLast = TRUE,
    nchains = 3,
    zInitialize = c("split", "random", "predefined"),
    yInitialize = c("split", "random", "predefined"),
    countChecksum = NULL,
    zInit = NULL,
    yInit = NULL,
    logfile = NULL,
    verbose = TRUE,
    reorder = TRUE) {

    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = FALSE,
        verbose = verbose)

    .logMessages("Starting Celda_CG: Clustering cells and genes.",
        logfile = logfile,
        append = TRUE,
        verbose = verbose)

    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)

    startTime <- Sys.time()

    counts <- .processCounts(counts)
    if (is.null(countChecksum)) {
        countChecksum <- .createCountChecksum(counts)
    }

    sampleLabel <- .processSampleLabels(sampleLabel, ncol(counts))
    s <- as.integer(sampleLabel)

    algorithm <- match.arg(algorithm)
    algorithmFun <- ifelse(algorithm == "Gibbs",
        ".cCCalcGibbsProbZ",
        ".cCCalcEMProbZ")
    zInitialize <- match.arg(zInitialize)
    yInitialize <- match.arg(yInitialize)

    allChains <- seq(nchains)

	# delta needs to be integer for computational speed
	delta <- as.integer(round(delta))
	if(delta == 0) delta <- 1L

    # Pre-compute lgamma values
    lggamma <- lgamma(seq(0, nrow(counts) + L) + gamma)
    #lgdelta <- c(NA, lgamma((seq(nrow(counts) + L) * delta)))
	lgdelta <- c(NA, lgamma(1:(sum(counts)+nrow(counts))))


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

        .logMessages(date(),
            ".. Initializing 'y' in chain",
            i,
            "with",
            paste0("'", yInitialize, "' "),
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
            z <- .initializeSplitZ(
                counts,
                K = K,
                alpha = alpha,
                beta = beta)
        } else {
            z <- .initializeCluster(K,
                ncol(counts),
                initial = NULL,
                fixed = NULL)
        }

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

        zBest <- z
        yBest <- y

        ## Calculate counts one time up front
        p <- .cCGDecomposeCounts(counts, s, z, y, K, L)
        mCPByS <- p$mCPByS
        nTSByC <- p$nTSByC
        nTSByCP <- p$nTSByCP
        nCP <- p$nCP
        nByG <- p$nByG
        nByC <- p$nByC
        nByTS <- p$nByTS
        nGByTS <- p$nGByTS
        nGByCP <- p$nGByCP
        nM <- p$nM
        nG <- p$nG
        nS <- p$nS
        rm(p)

        ll <- .cCGCalcLL(K = K,
            L = L,
            mCPByS = mCPByS,
            nTSByCP = nTSByCP,
            nByG = nByG,
            nByTS = nByTS,
            nGByTS = nGByTS,
            nS = nS,
            nG = nG,
            alpha = alpha,
            beta = beta,
            delta = delta,
            gamma = gamma)

        iter <- 1L
        numIterWithoutImprovement <- 0L
        doCellSplit <- TRUE
        doGeneSplit <- TRUE
        while (iter <= maxIter & numIterWithoutImprovement <= stopIter) {
            ## Gibbs sampling for each gene
            lgbeta <- lgamma(seq(0, max(nCP)) + beta)
            nextY <- .cGCalcGibbsProbY(counts = nGByCP,
                    nTSByC = nTSByCP,
                    nByTS = nByTS,
                    nGByTS = nGByTS,
                    nByG = nByG,
                    y = y,
                    L = L,
                    nG = nG,
                    beta = beta,
                    delta = delta,
                    gamma = gamma,
                    lgbeta = lgbeta,
                    lggamma = lggamma,
                    lgdelta = lgdelta)
            nTSByCP <- nextY$nTSByC
            nGByTS <- nextY$nGByTS
            nByTS <- nextY$nByTS
            nTSByC <- .rowSumByGroupChange(counts, nTSByC, nextY$y, y, L)
            y <- nextY$y

            ## Gibbs or EM sampling for each cell
            nextZ <- do.call(algorithmFun, list(counts = nTSByC,
                mCPByS = mCPByS,
                nGByCP = nTSByCP,
                nCP = nCP,
                nByC = nByC,
                z = z,
                s = s,
                K = K,
                nG = L,
                nM = nM,
                alpha = alpha,
                beta = beta))
            mCPByS <- nextZ$mCPByS
            nTSByCP <- nextZ$nGByCP
            nCP <- nextZ$nCP
            nGByCP <- .colSumByGroupChange(counts, nGByCP, nextZ$z, z, K)
            z <- nextZ$z

            ## Perform split on i-th iteration defined by splitOnIter
            tempLl <- .cCGCalcLL(K = K,
                L = L,
                mCPByS = mCPByS,
                nTSByCP = nTSByCP,
                nByG = nByG,
                nByTS = nByTS,
                nGByTS = nGByTS,
                nS = nS,
                nG = nG,
                alpha = alpha,
                beta = beta,
                delta = delta,
                gamma = gamma)

            if (L > 2 & iter != maxIter &
                (((numIterWithoutImprovement == stopIter &
                    !all(tempLl > ll)) & isTRUE(splitOnLast)) |
                        (splitOnIter > 0 & iter %% splitOnIter == 0 &
                            isTRUE(doGeneSplit)))) {
                .logMessages(date(),
                    " .... Determining if any gene clusters should be split.",
                    logfile = logfile,
                    append = TRUE,
                    sep = "",
                    verbose = verbose)
                res <- .cCGSplitY(counts,
                        y,
                        mCPByS,
                        nGByCP,
                        nTSByC,
                        nTSByCP,
                        nByG,
                        nByTS,
                        nGByTS,
                        nCP,
                        s,
                        z,
                        K,
                        L,
                        nS,
                        nG,
                        alpha,
                        beta,
                        delta,
                        gamma,
                        yProb = t(nextY$probs),
                        maxClustersToTry = max(L / 2, 10),
                        minCell = 3)
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
                nTSByCP <- res$nTSByCP
                nByTS <- res$nByTS
                nGByTS <- res$nGByTS
                nTSByC <- .rowSumByGroup(counts, group = y, L = L)
            }

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
                res <- .cCGSplitZ(counts,
                        mCPByS,
                        nTSByC,
                        nTSByCP,
                        nByG,
                        nByTS,
                        nGByTS,
                        nCP,
                        s,
                        z,
                        K,
                        L,
                        nS,
                        nG,
                        alpha,
                        beta,
                        delta,
                        gamma,
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
                nTSByCP <- res$nTSByCP
                nCP <- res$nCP
                nGByCP <- .colSumByGroup(counts, group = z, K = K)
            }

            ## Calculate complete likelihood
            tempLl <- .cCGCalcLL(K = K,
                    L = L,
                    mCPByS = mCPByS,
                    nTSByCP = nTSByCP,
                    nByG = nByG,
                    nByTS = nByTS,
                    nGByTS = nGByTS,
                    nS = nS,
                    nG = nG,
                    alpha = alpha,
                    beta = beta,
                    delta = delta,
                    gamma = gamma)
            if ((all(tempLl > ll)) | iter == 1) {
                zBest <- z
                yBest <- y
                llBest <- tempLl
                numIterWithoutImprovement <- 1L
            } else {
                numIterWithoutImprovement <- numIterWithoutImprovement + 1L
            }
            ll <- c(ll, tempLl)

            .logMessages(date(),
                " .... Completed iteration: ",
                iter,
                " | logLik: ",
                tempLl,
                logfile = logfile,
                append = TRUE,
                sep = "",
                verbose = verbose)
            iter <- iter + 1L
        }

        names <- list(row = rownames(counts),
            column = colnames(counts),
            sample = levels(sampleLabel))

        result <- list(z = zBest,
            y = yBest,
            completeLogLik = ll,
            finalLogLik = llBest,
            K = K,
            L = L,
            alpha = alpha,
            beta = beta,
            delta = delta,
            gamma = gamma,
            sampleLabel = sampleLabel,
            names = names,
            countChecksum = countChecksum)

        class(result) <- "celda_CG"

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

    ## Peform reordering on final Z and Y assigments:
    bestResult <- methods::new("celda_CG",
        clusters = list(z = zBest, y = yBest),
        params = list(K = as.integer(K),
            L = as.integer(L),
            alpha = alpha,
            beta = beta,
            delta = delta,
            gamma = gamma,
            countChecksum = countChecksum),
        completeLogLik = ll,
        finalLogLik = llBest,
        sampleLabel = sampleLabel,
        names = names)
    if (isTRUE(reorder)) {
        bestResult <- .reorderCeldaCG(counts = counts, res = bestResult)
    }

    endTime <- Sys.time()
    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)
    .logMessages("Completed Celda_CG. Total time:",
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


#' @title Simulate cells from the celda_CG model
#' @description Generates a simulated counts matrix, cell subpopulation
#'  clusters, sample labels, and feature module clusters according to the
#'  generative process of the celda_CG model.
#' @param model Character. Options available in `celda::availableModels`.
#' @param S Integer. Number of samples to simulate. Default 5.
#' @param CRange Integer vector. A vector of length 2 that specifies the lower
#'  and upper bounds of the number of cells to be generated in each sample.
#'  Default c(50, 100).
#' @param NRange Integer vector. A vector of length 2 that specifies the lower
#'  and upper bounds of the number of counts generated for each cell. Default
#'  c(500, 1000).
#' @param G Integer. The total number of features to be simulated. Default 100.
#' @param K Integer. Number of cell populations. Default 5.
#' @param L Integer. Number of feature modules. Default 10.
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount
#'  to each cell population in each sample. Default 1.
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to
#'  each feature module in each cell population. Default 1.
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to
#'  the number of features in each module. Default 5.
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to
#'  each feature in each module. Default 1.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @param ... Additional parameters.
#' @return List. Contains the simulated matrix `counts`, cell population
#'  clusters `z`, feature module clusters `y`, sample assignments `sampleLabel`,
#'  and input parameters.
#' @seealso `celda_C()` for simulating cell subpopulations and `celda_G()` for
#'  simulating feature modules.
#' @examples
#' celdaCGSim <- simulateCells(model = "celda_CG")
#' @export
simulateCellscelda_CG <- function(model,
    S = 5,
    CRange = c(50, 100),
    NRange = c(500, 1000),
    G = 100,
    K = 5,
    L = 10,
    alpha = 1,
    beta = 1,
    gamma = 5,
    delta = 1,
    seed = 12345,
    ...) {

    if (is.null(seed)) {
        res <- .simulateCellscelda_CG(model = model,
            S = S,
            CRange = CRange,
            NRange = NRange,
            G = G,
            K = K,
            L = L,
            alpha = alpha,
            beta = beta,
            gamma = gamma,
            delta = delta,
            ...)
    } else {
        with_seed(seed,
            res <- .simulateCellscelda_CG(model = model,
                S = S,
                CRange = CRange,
                NRange = NRange,
                G = G,
                K = K,
                L = L,
                alpha = alpha,
                beta = beta,
                gamma = gamma,
                delta = delta,
                ...))
    }
    return(res)
}


.simulateCellscelda_CG <- function(model = model,
    S = S,
    CRange = CRange,
    NRange = NRange,
    G = G,
    K = K,
    L = L,
    alpha = alpha,
    beta = beta,
    gamma = gamma,
    delta = delta,
    ...) {

    ## Number of cells per sample
    nC <- sample(seq(CRange[1], CRange[2]), size = S, replace = TRUE)
    nCSum <- sum(nC)
    cellSampleLabel <- rep(seq(S), nC)

    ## Select number of transcripts per cell
    nN <- sample(seq(NRange[1], NRange[2]),
        size = length(cellSampleLabel),
        replace = TRUE)

    ## Generate cell population distribution for each sample
    theta <- t(.rdirichlet(S, rep(alpha, K)))

    ## Assign cells to cellular subpopulations
    z <- unlist(lapply(seq(S), function(i) {
        sample(seq(K),
            size = nC[i],
            prob = theta[, i],
            replace = TRUE)
    }))

    ## Generate transcriptional state distribution for each cell subpopulation
    phi <- .rdirichlet(K, rep(beta, L))

    ## Assign genes to gene modules
    eta <- .rdirichlet(1, rep(gamma, L))
    y <- sample(seq(L),
        size = G,
        prob = eta,
        replace = TRUE)
    if (length(table(y)) < L) {
        warning("Some gene modules did not receive any genes after sampling.",
            " Try increasing G and/or making gamma larger.")
        L <- length(table(y))
        y <- as.integer(as.factor(y))
    }

    psi <- matrix(0, nrow = G, ncol = L)
    for (i in seq(L)) {
        ind <- y == i
        psi[ind, i] <- .rdirichlet(1, rep(delta, sum(ind)))
    }

    ## Select transcript distribution for each cell
    cellCounts <- matrix(0, nrow = G, ncol = nCSum)
    for (i in seq(nCSum)) {
        transcriptionalStateDist <- as.integer(stats::rmultinom(1,
            size = nN[i], prob = phi[z[i], ]))
        for (j in seq(L)) {
            if (transcriptionalStateDist[j] > 0) {
                cellCounts[, i] <- cellCounts[, i] + stats::rmultinom(1,
                    size = transcriptionalStateDist[j], prob = psi[, j])
            }
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

    ## Assign gene/cell/sample names
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
    result <- methods::new("celda_CG",
        clusters = list(z = z, y = y),
        params = list(K = as.integer(K),
            L = as.integer(L),
            alpha = alpha,
            beta = beta,
            delta = delta,
            gamma = gamma,
            countChecksum = countChecksum),
        sampleLabel = cellSampleLabel,
        names = names
    )

    result <- .reorderCeldaCG(counts = cellCounts, res = result)

    return(list(z = clusters(result)$z,
        y = clusters(result)$y,
        sampleLabel = cellSampleLabel,
        counts = cellCounts,
        K = K,
        L = L,
        CRange = CRange,
        NRange = NRange,
        S = S,
        alpha = alpha,
        beta = beta,
        gamma = gamma,
        delta = delta)
    )
}


#' @title Matrix factorization for results from celda_CG
#' @description Generates factorized matrices showing the contribution of each
#'  feature in each module, each module in each cell and/or cell population,
#'  and each cell population in each sample.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda model. Options are "celda_C" or "celda_CG". Celda
#'  object of class "celda_CG".
#' @param type Character vector. A vector containing one or more of "counts",
#'  "proportion", or "posterior". "counts" returns the raw number of counts for
#'  each factorized matrix. "proportions" returns the normalized probabilities
#'  for each factorized matrix, which are calculated by dividing the raw counts
#'  in each factorized matrix by the total counts in each column. "posterior"
#'  returns the posterior estimates. Default
#'  `c("counts", "proportion", "posterior")`.
#' @return A list with elements for `counts`, `proportions`, or `posterior`
#'  probabilities. Each element will be a list containing factorized matrices
#'  for `module`, `cellPopulation`, and `sample`. Additionally, the
#'  contribution of each module in each individual cell will be included in the
#'  `cell` element of `counts` and `proportions` elements.
#' @seealso `celda_CG()` for clustering features and cells
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' factorizedMatrices <- factorizeMatrix(celdaCGSim$counts,
#'     celdaCGMod, "posterior")
#' @export
setMethod("factorizeMatrix", signature(celdaMod = "celda_CG"),
    function(counts,
        celdaMod,
        type = c("counts", "proportion", "posterior")) {

        counts <- .processCounts(counts)
        compareCountMatrix(counts, celdaMod)

        K <- params(celdaMod)$K
        L <- params(celdaMod)$L
        z <- clusters(celdaMod)$z
        y <- clusters(celdaMod)$y
        alpha <- params(celdaMod)$alpha
        beta <- params(celdaMod)$beta
        delta <- params(celdaMod)$delta
        gamma <- params(celdaMod)$gamma
        sampleLabel <- sampleLabel(celdaMod)
        s <- as.integer(sampleLabel)

        ## Calculate counts one time up front
        p <- .cCGDecomposeCounts(counts, s, z, y, K, L)
        nS <- p$nS
        nG <- p$nG
        nM <- p$nM
        mCPByS <- p$mCPByS
        nTSByC <- p$nTSByC
        nTSByCP <- p$nTSByCP
        nByG <- p$nByG
        nByTS <- p$nByTS
        nGByTS <- p$nGByTS
        nGByTS[nGByTS == 0] <- 1

        nGByTS <- matrix(0, nrow = length(y), ncol = L)
        nGByTS[cbind(seq(nG), y)] <- p$nByG

        LNames <- paste0("L", seq(L))
        KNames <- paste0("K", seq(K))
        colnames(nTSByC) <- matrixNames(celdaMod)$column
        rownames(nTSByC) <- LNames
        colnames(nGByTS) <- LNames
        rownames(nGByTS) <- matrixNames(celdaMod)$row
        rownames(mCPByS) <- KNames
        colnames(mCPByS) <- matrixNames(celdaMod)$sample
        colnames(nTSByCP) <- KNames
        rownames(nTSByCP) <- LNames

        countsList <- c()
        propList <- c()
        postList <- c()
        res <- list()

        if (any("counts" %in% type)) {
            countsList <- list(sample = mCPByS,
                cellPopulation = nTSByCP,
                cell = nTSByC,
                module = nGByTS,
                geneDistribution = nGByTS)
            res <- c(res, list(counts = countsList))
        }

        if (any("proportion" %in% type)) {
            ## Need to avoid normalizing cell/gene states with zero cells/genes
            uniqueZ <- sort(unique(z))
            tempNTSByCP <- nTSByCP
            tempNTSByCP[, uniqueZ] <- normalizeCounts(tempNTSByCP[, uniqueZ],
                normalize = "proportion")

            uniqueY <- sort(unique(y))
            tempNGByTS <- nGByTS
            tempNGByTS[, uniqueY] <- normalizeCounts(tempNGByTS[, uniqueY],
                normalize = "proportion")
            tempNGByTS <- nGByTS / sum(nGByTS)

            propList <- list(sample = normalizeCounts(mCPByS,
                normalize = "proportion"),
                cellPopulation = tempNTSByCP,
                cell = normalizeCounts(nTSByC, normalize = "proportion"),
                module = tempNGByTS,
                geneDistribution = tempNGByTS)
            res <- c(res, list(proportions = propList))
        }

        if (any("posterior" %in% type)) {
            gs <- nGByTS
            gs[cbind(seq(nG), y)] <- gs[cbind(seq(nG), y)] + delta
            gs <- normalizeCounts(gs, normalize = "proportion")
            tempNGByTS <- (nGByTS + gamma) / sum(nGByTS + gamma)

            postList <- list(sample = normalizeCounts(mCPByS + alpha,
                normalize = "proportion"),
                cellPopulation = normalizeCounts(nTSByCP + beta,
                    normalize = "proportion"),
                module = gs,
                geneDistribution = tempNGByTS)
            res <- c(res, posterior = list(postList))
        }

        return(res)
    })

# Calculate the loglikelihood for the celda_CG model
.cCGCalcLL <- function(K,
    L,
    mCPByS,
    nTSByCP,
    nByG,
    nByTS,
    nGByTS,
    nS,
    nG,
    alpha,
    beta,
    delta,
    gamma) {

    nG <- sum(nGByTS)

    ## Calculate for "Theta" component
    a <- nS * lgamma(K * alpha)
    b <- sum(lgamma(mCPByS + alpha))
    c <- -nS * K * lgamma(alpha)
    d <- -sum(lgamma(colSums(mCPByS + alpha)))

    thetaLl <- a + b + c + d

    ## Calculate for "Phi" component
    a <- K * lgamma(L * beta)
    b <- sum(lgamma(nTSByCP + beta))
    c <- -K * L * lgamma(beta)
    d <- -sum(lgamma(colSums(nTSByCP + beta)))

    phiLl <- a + b + c + d

    ## Calculate for "Psi" component
    a <- sum(lgamma(nGByTS * delta))
    b <- sum(lgamma(nByG + delta))
    c <- -nG * lgamma(delta)
    d <- -sum(lgamma(nByTS + (nGByTS * delta)))

    psiLl <- a + b + c + d

    ## Calculate for "Eta" side
    a <- lgamma(L * gamma)
    b <- sum(lgamma(nGByTS + gamma))
    c <- -L * lgamma(gamma)
    d <- -lgamma(sum(nGByTS + gamma))

    etaLl <- a + b + c + d

    final <- thetaLl + phiLl + psiLl + etaLl
    return(final)
}


#' @title Calculate Celda_CG log likelihood
#' @description Calculates the log likelihood for user-provided cell population
#'  and feature module clusters using the `celda_CG()` model.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells.
#' @param sampleLabel Vector or factor. Denotes the sample label for each cell
#'  (column) in the count matrix.
#' @param z Numeric vector. Denotes cell population labels.
#' @param y Numeric vector. Denotes feature module labels.
#' @param K Integer. Number of cell populations.
#' @param L Integer. Number of feature modules.
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount
#'  to each cell population in each sample. Default 1.
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to
#'  each feature module in each cell population. Default 1.
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to
#'  each feature in each module. Default 1.
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to
#'  the number of features in each module. Default 1.
#' @return The log likelihood for the given cluster assignments
#' @seealso `celda_CG()` for clustering features and cells
#' @examples
#' data(celdaCGSim)
#' loglik <- logLikelihoodcelda_CG(celdaCGSim$counts,
#'     sampleLabel = celdaCGSim$sampleLabel,
#'     z = celdaCGSim$z,
#'     y = celdaCGSim$y,
#'     K = celdaCGSim$K,
#'     L = celdaCGSim$L,
#'     alpha = celdaCGSim$alpha,
#'     beta = celdaCGSim$beta,
#'     gamma = celdaCGSim$gamma,
#'     delta = celdaCGSim$delta)
#'
#' loglik <- logLikelihood(celdaCGSim$counts,
#'     model = "celda_CG",
#'     sampleLabel = celdaCGSim$sampleLabel,
#'     z = celdaCGSim$z,
#'     y = celdaCGSim$y,
#'     K = celdaCGSim$K,
#'     L = celdaCGSim$L,
#'     alpha = celdaCGSim$alpha,
#'     beta = celdaCGSim$beta,
#'     gamma = celdaCGSim$gamma,
#'     delta = celdaCGSim$delta)
#' @export
logLikelihoodcelda_CG <- function(counts,
    sampleLabel,
    z,
    y,
    K,
    L,
    alpha,
    beta,
    delta,
    gamma) {

    if (sum(z > K) > 0) {
        stop("An entry in z contains a value greater than the provided K.")
    }
    if (sum(y > L) > 0) {
        stop("An entry in y contains a value greater than the provided L.")
    }

    sampleLabel <- .processSampleLabels(sampleLabel, ncol(counts))
    s <- as.integer(sampleLabel)
    p <- .cCGDecomposeCounts(counts, s, z, y, K, L)
    final <- .cCGCalcLL(K = K,
            L = L,
            mCPByS = p$mCPByS,
            nTSByCP = p$nTSByCP,
            nByG = p$nByG,
            nByTS = p$nByTS,
            nGByTS = p$nGByTS,
            nS = p$nS,
            nG = p$nG,
            alpha = alpha,
            beta = beta,
            delta = delta,
            gamma = gamma)
    return(final)
}


# Takes raw counts matrix and converts it to a series of matrices needed for
# log likelihood calculation
# @param counts Integer matrix. Rows represent features and columns represent
# cells.
# @param s Integer vector. Contains the sample label for each cell (column) in
# the count matrix.
# @param z Numeric vector. Denotes cell population labels.
# @param y Numeric vector. Denotes feature module labels.
# @param K Integer. Number of cell populations.
# @param L Integer. Number of feature modules.
.cCGDecomposeCounts <- function(counts, s, z, y, K, L) {
    nS <- length(unique(s))
    mCPByS <- matrix(as.integer(table(factor(z, levels = seq(K)), s)),
        ncol = nS)
    nTSByC <- .rowSumByGroup(counts, group = y, L = L)
    nTSByCP <- .colSumByGroup(nTSByC, group = z, K = K)
    nCP <- as.integer(colSums(nTSByCP))
    nByG <- as.integer(rowSums(counts))
    nByC <- as.integer(colSums(counts))
    nByTS <- as.integer(.rowSumByGroup(matrix(nByG, ncol = 1),
        group = y, L = L))
    nGByTS <- tabulate(y, L) + 1 ## Add pseudogene to each module
    nGByCP <- .colSumByGroup(counts, group = z, K = K)

    nG <- nrow(counts)
    nM <- ncol(counts)

    return(list(mCPByS = mCPByS,
            nTSByC = nTSByC,
            nTSByCP = nTSByCP,
            nCP = nCP,
            nByG = nByG,
            nByC = nByC,
            nByTS = nByTS,
            nGByTS = nGByTS,
            nGByCP = nGByCP,
            nM = nM,
            nG = nG,
            nS = nS))
}


#' @title Conditional probabilities for cells and features from a Celda_CG
#'  model
#' @description Calculates the conditional probability of each cell belonging
#'  to each subpopulation given all other cell cluster assignments as well as
#'  each feature belonging to each module given all other feature cluster
#'  assignments in a `celda_CG()` result.
#' @param celdaMod Celda object of class `celda_CG`.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param log Logical. If FALSE, then the normalized conditional probabilities
#'  will be returned. If TRUE, then the unnormalized log probabilities will be
#'  returned. Default FALSE.
#' @param ... Additional parameters.
#' @return A list containging a matrix for the conditional cell and feature
#'  cluster probabilities.
#' @seealso `celda_CG()` for clustering features and cells
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' clusterProb <- clusterProbability(celdaCGSim$counts, celdaCGMod)
#' @export
setMethod("clusterProbability", signature(celdaMod = "celda_CG"),
    function(counts, celdaMod, log = FALSE, ...) {
        s <- as.integer(sampleLabel(celdaMod))
        z <- clusters(celdaMod)$z
        K <- params(celdaMod)$K
        y <- clusters(celdaMod)$y
        L <- params(celdaMod)$L
        alpha <- params(celdaMod)$alpha
        delta <- params(celdaMod)$delta
        beta <- params(celdaMod)$beta
        gamma <- params(celdaMod)$gamma

        p <- .cCGDecomposeCounts(counts, s, z, y, K, L)
        lgbeta <- lgamma(seq(0, max(p$nCP)) + beta)
        lggamma <- lgamma(seq(0, nrow(counts) + L) + gamma)
        lgdelta <- c(NA, lgamma((seq(nrow(counts) + L) * delta)))

        nextZ <- .cCCalcGibbsProbZ(counts = p$nTSByC,
            mCPByS = p$mCPByS,
            nGByCP = p$nTSByCP,
            nCP = p$nCP,
            nByC = p$nByC,
            z = z,
            s = s,
            K = K,
            nG = L,
            nM = p$nM,
            alpha = alpha,
            beta = beta,
            doSample = FALSE)
        zProb <- t(nextZ$probs)

        ## Gibbs sampling for each gene
        nextY <- .cGCalcGibbsProbY(counts = p$nGByCP,
            nTSByC = p$nTSByCP,
            nByTS = p$nByTS,
            nGByTS = p$nGByTS,
            nByG = p$nByG,
            y = y,
            L = L,
            nG = p$nG,
            lgbeta = lgbeta,
            lgdelta = lgdelta,
            lggamma = lggamma,
            delta = delta,
            doSample = FALSE)

        yProb <- t(nextY$probs)

        if (!isTRUE(log)) {
            zProb <- .normalizeLogProbs(zProb)
            yProb <- .normalizeLogProbs(yProb)
        }

        return(list(zProbability = zProb, yProbability = yProb))
    })


#' @title Calculate the perplexity on new data with a celda_CG model
#' @description Perplexity is a statistical measure of how well a probability
#'  model can predict new data. Lower perplexity indicates a better model.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda object of class "celda_C", "celda_G" or "celda_CG".
#' @param newCounts A new counts matrix used to calculate perplexity. If NULL,
#'  perplexity will be calculated for the 'counts' matrix. Default NULL.
#' @return Numeric. The perplexity for the provided count data and model.
#' @seealso `celda_CG()` for clustering features and cells
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' perplexity <- perplexity(celdaCGSim$counts, celdaCGMod)
#' @importFrom matrixStats logSumExp
#' @export
setMethod("perplexity", signature(celdaMod = "celda_CG"),
    function(counts, celdaMod, newCounts = NULL) {

        if (!("celda_CG" %in% class(celdaMod))) {
            stop("The celdaMod provided was not of class celda_CG.")
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
            type = c("posterior", "counts"))

        theta <- log(factorized$posterior$sample)
        phi <- factorized$posterior$cellPopulation
        psi <- factorized$posterior$module
        s <- as.integer(sampleLabel(celdaMod))
        eta <- factorized$posterior$geneDistribution
        nGByTS <- factorized$counts$geneDistribution

        etaProb <- log(eta) * nGByTS
        geneByPopProb <- log(psi %*% phi)
        innerLogProb <- eigenMatMultInt(geneByPopProb, newCounts) + theta[, s]
        # innerLogProb = (t(geneByPopProb) %*% newCounts) + theta[, s]

        log.px <- sum(apply(innerLogProb, 2, matrixStats::logSumExp))
        # + sum(etaProb)
        perplexity <- exp(- (log.px / sum(newCounts)))
        return(perplexity)
    })


.reorderCeldaCG <- function(counts, res) {
    # Reorder K
    if (params(res)$K > 2 & isTRUE(length(unique(clusters(res)$z)) > 1)) {

        res@clusters$z <- as.integer(as.factor(clusters(res)$z))
        fm <- factorizeMatrix(counts = counts,
            celdaMod = res,
            type = "posterior")
        uniqueZ <- sort(unique(clusters(res)$z))
        d <- .cosineDist(fm$posterior$cellPopulation[, uniqueZ])
        h <- stats::hclust(d, method = "complete")

        res <- recodeClusterZ(res,
            from = h$order,
            to = seq(length(h$order)))
    }

    # Reorder L
    if (params(res)$L > 2 & isTRUE(length(unique(clusters(res)$y)) > 1)) {

        res@clusters$y <- as.integer(as.factor(clusters(res)$y))
        fm <- factorizeMatrix(counts = counts,
                celdaMod = res,
                type = "posterior")
        uniqueY <- sort(unique(clusters(res)$y))
        cs <- prop.table(t(fm$posterior$cellPopulation[uniqueY, ]), 2)
        d <- .cosineDist(cs)
        h <- stats::hclust(d, method = "complete")

        res <- recodeClusterY(res, from = h$order, to = seq(length(h$order)))
    }
    return(res)
}


#' @title Heatmap for celda_CG
#' @description Renders an expression heatmap to visualize `celda_CG()` results.
#'  The top `nfeatures` for each module will be included in the heatmap.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate.
#'  `celdaMod`.
#' @param celdaMod Celda object of class `celda_CG`.
#' @param nfeatures Integer. Maximum number of features to select for each
#'  module. Default 25.
#' @param ... Additional parameters.
#' @seealso `celda_CG()` for clustering features and cells and `celdaTsne()`
#'  for generating 2-dimensional coordinates.
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' celdaHeatmap(celdaCGSim$counts, celdaCGMod)
#' @return A list containing dendrograms and the heatmap grob
#' @export
setMethod("celdaHeatmap", signature(celdaMod = "celda_CG"),
    function(counts, celdaMod, nfeatures = 25, ...) {
        fm <- factorizeMatrix(counts, celdaMod, type = "proportion")
        top <- celda::topRank(fm$proportions$module, n = nfeatures)
        ix <- unlist(top$index)
        norm <- normalizeCounts(counts,
                normalize = "proportion",
                transformationFun = sqrt)
        plotHeatmap(norm[ix, ],
            z = clusters(celdaMod)$z,
            y = clusters(celdaMod)$y[ix],
            ...)
    })


#' @title tSNE for celda_CG
#' @description Embeds cells in two dimensions using tSNE based on a `celda_CG`
#'  model. tSNE is run on module probabilities to reduce the number of features
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
#' @param initialDims Integer. PCA will be used to reduce the dimentionality
#'  of the dataset. The top 'initialDims' principal components will be used
#'  for tSNE. Default 20.
#' @param modules Integer vector. Determines which features modules to use for
#'  tSNE. If NULL, all modules will be used. Default NULL.
#' @param perplexity Numeric. Perplexity parameter for tSNE. Default 20.
#' @param maxIter Integer. Maximum number of iterations in tSNE generation.
#'  Default 2500.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @seealso `celda_CG()` for clustering features and cells  and `celdaHeatmap()`
#'  for displaying expression
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' tsneRes <- celdaTsne(celdaCGSim$counts, celdaCGMod)
#' @return A two column matrix of t-SNE coordinates
#' @export
setMethod("celdaTsne", signature(celdaMod = "celda_CG"),
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
            res <- .celdaTsneCG(counts = counts,
                celdaMod = celdaMod,
                maxCells = maxCells,
                minClusterSize = minClusterSize,
                initialDims = initialDims,
                modules = modules,
                perplexity = perplexity,
                maxIter = maxIter)
        } else {
            with_seed(seed,
                res <- .celdaTsneCG(counts = counts,
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


.celdaTsneCG <- function(counts,
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
    norm <- preparedCountInfo$norm
    res <- .calculateTsne(norm,
        doPca = FALSE,
        perplexity = perplexity,
        maxIter = maxIter,
        initialDims = initialDims)
    final <- matrix(NA, nrow = ncol(counts), ncol = 2)
    final[preparedCountInfo$cellIx, ] <- res
    rownames(final) <- colnames(counts)
    colnames(final) <- c("tSNE_1", "tSNE_2")
    return(final)
}


#' @title umap for celda_CG
#' @description Embeds cells in two dimensions using umap based on a `celda_CG`
#'  model. umap is run on module probabilities to reduce the number of features
#'  instead of using PCA. Module probabilities square-root trasformed before
#'  applying tSNE.
#'
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
#' @seealso `celda_CG()` for clustering features and cells and `celdaHeatmap()`
#'  for displaying expression.
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' umapRes <- celdaUmap(celdaCGSim$counts, celdaCGMod)
#' @return A two column matrix of umap coordinates
#' @export
setMethod("celdaUmap",
    signature(celdaMod = "celda_CG"), function(counts,
        celdaMod,
        maxCells = NULL,
        minClusterSize = 100,
        modules = NULL,
        seed = 12345,
        nNeighbors = 30,
        minDist = 0.75,
        spread = 1,
        cores = 1,
        ...) {

        if (is.null(seed)) {
            res <- .celdaUmapCG(counts = counts,
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
                res <- .celdaUmapCG(counts = counts,
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


.celdaUmapCG <- function(counts,
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
            " cells per cluster. Try increasing 'maxCells' or",
            " decreasing 'minClusterSize'.")
      }
    } else {
      maxCells <- ncol(counts)
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
        ## Rounding can cause number to be off by a few, so ceiling is used
        ## with a second round of subtraction
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

    norm <- t(normalizeCounts(fm$counts$cell[modulesToUse, cellIx],
        normalize = "proportion",
        transformationFun = sqrt))
    return(list(norm = norm, cellIx = cellIx))
}


#' @title Probability map for a celda_CG model
#' @description Renders probability and relative expression heatmaps to
#'  visualize the relationship between features and cell populations or cell
#'  populations and samples.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda object of class `celda_CG`.
#' @param level Character. One of 'cellPopulation' or 'sample'.
#'  'cellPopulation' will display the absolute probabilities and relative
#'  normalized expression of each module in each cell population. 'sample'
#'  will display the absolute probabilities and relative normalized abundance
#'  of each cell population in each sample. Default 'cellPopulation'.
#' @param ... Additional parameters.
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' celdaProbabilityMap(celdaCGSim$counts, celdaCGMod)
#' @return A grob containing the specified plots
#' @importFrom gridExtra grid.arrange
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @seealso `celda_CG()` for clustering features and cells
#' @export
setMethod("celdaProbabilityMap", signature(celdaMod = "celda_CG"),
    function(counts, celdaMod, level = c("cellPopulation", "sample"), ...) {

        counts <- .processCounts(counts)
        compareCountMatrix(counts, celdaMod)

        level <- match.arg(level)
        factorized <- factorizeMatrix(celdaMod = celdaMod, counts = counts)
        zInclude <- which(tabulate(clusters(celdaMod)$z,
            params(celdaMod)$K) > 0)
        yInclude <- which(tabulate(clusters(celdaMod)$y,
            params(celdaMod)$L) > 0)

        if (level == "cellPopulation") {
            pop <- factorized$proportions$cellPopulation[yInclude,
                zInclude, drop = FALSE]
            popNorm <- normalizeCounts(pop,
                normalize = "proportion",
                transformationFun = sqrt,
                scaleFun = base::scale)

            percentile9 <- round(stats::quantile(pop, .9), digits = 2) * 100
            col1 <- grDevices::colorRampPalette(c("#FFFFFF",
                RColorBrewer::brewer.pal(n = 9, name = "Blues")))(percentile9)
            col2 <- grDevices::colorRampPalette(c("#08306B",
                c("#006D2C", "Yellowgreen", "Yellow", "Orange",
                    "Red")))(100 - percentile9)
            col <- c(col1, col2)
            breaks <- seq(0, 1, length.out = length(col))

            g1 <- plotHeatmap(pop,
                colorScheme = "sequential",
                scaleRow = NULL,
                clusterCell = FALSE,
                clusterFeature = FALSE,
                showNamesCell = TRUE,
                showNamesFeature = TRUE,
                breaks = breaks,
                col = col,
                main = "Absolute Probability",
                silent = TRUE
            )
            g2 <- plotHeatmap(popNorm,
                colorScheme = "divergent",
                clusterCell = FALSE,
                clusterFeature = FALSE,
                showNamesCell = TRUE,
                showNamesFeature = TRUE,
                main = "Relative Expression",
                silent = TRUE)
            gridExtra::grid.arrange(g1$gtable, g2$gtable, ncol = 2)
        } else {
            samp <- factorized$proportions$sample
            col <- grDevices::colorRampPalette(c("white",
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
                sampNorm <- normalizeCounts(factorized$counts$sample,
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
        }
    })


#' @title Lookup the module of a feature
#' @description Finds the module assignments of given features in a `celda_G()`
#'  model
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Model of class `celda_CG`.
#' @param feature Character vector. The module assignemnts will be found for
#'  feature names in this vector.
#' @param exactMatch Logical. Whether an exact match or a partial match using
#'  `grep()` is used to look up the feature in the rownames of the counts
#'  matrix. Default TRUE.
#' @return List. Each element contains the module of the provided feature.
#' @seealso `celda_CG()` for clustering features and cells
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' module <- featureModuleLookup(celdaCGSim$counts,
#'     celdaCGMod,
#'     c("Gene_1", "Gene_XXX"))
#' @export
setMethod("featureModuleLookup", signature(celdaMod = "celda_CG"),
    function(counts, celdaMod, feature, exactMatch = TRUE) {
        list <- list()
        if (!isTRUE(exactMatch)) {
            featureGrep <- c()
            for (x in seq(length(feature))) {
                featureGrep <- c(featureGrep, rownames(counts)[grep(feature[x],
                    rownames(counts))])
            }
            feature <- featureGrep
        }
        for (x in seq(length(feature))) {
            if (feature[x] %in% rownames(counts)) {
                list[x] <- clusters(celdaMod)$y[which(rownames(counts) ==
                    feature[x])]
            } else {
                list[x] <- paste0("No feature was identified matching '",
                        feature[x],
                        "'.")
            }
        }
        names(list) <- feature
        return(list)
    })
