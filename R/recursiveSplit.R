.singleSplitZ <- function(counts,
    z,
    s,
    K,
    minCell = 3,
    alpha = 1,
    beta = 1) {

    zTa <- tabulate(z, K)
    zToSplit <- which(zTa > minCell)
    bestZ <- z
    bestLl <- -Inf
    for (i in zToSplit) {
        clustLabel <- .celda_C(
            counts[, z == i, drop = FALSE],
            K = 2,
            zInitialize = "random",
            splitOnIter = -1,
            splitOnLast = FALSE,
            verbose = FALSE)

        if (length(unique(clusters(clustLabel)$z)) == 2) {
            ix <- z == i
            newZ <- z
            newZ[ix] <- ifelse(clusters(clustLabel)$z == 2, i, K)
            ll <- logLikelihood.celda_C(counts, s, newZ, K, alpha, beta)

            if (ll > bestLl) {
                bestZ <- newZ
                bestLl <- ll
            }
        }
    }
    return(list(ll = bestLl, z = bestZ))
}

.singleSplitY <- function(counts,
    y,
    L,
    minFeature = 3,
    beta = 1,
    delta = 1,
    gamma = 1) {

    yTa <- tabulate(y, L)
    yToSplit <- which(yTa > minFeature)

    bestY <- y
    bestLl <- -Inf
    # previousY <- y
    for (i in yToSplit) {
        clustLabel <- .celda_G(counts[y == i, , drop = FALSE],
            L = 2,
            yInitialize = "random",
            splitOnIter = -1,
            splitOnLast = FALSE,
            nchains = 1,
            verbose = FALSE)

        if (length(unique(clusters(clustLabel)$y)) == 2) {
            ix <- y == i
            newY <- y
            newY[ix] <- ifelse(clusters(clustLabel)$y == 2, i, L)
            ll <- logLikelihood.celda_G(counts, newY, L, beta, delta, gamma)

            if (ll > bestLl) {
                bestY <- newY
                bestLl <- ll
            }
        }
    }
    return(list(ll = bestLl, y = bestY))
}


#' @title Recursive cell splitting
#' @description Uses the `celda_C` model to cluster cells into population for
#'  range of possible K's. The cell population labels of the previous "K-1"
#'  model are used as the initial values in the current model with K cell
#'  populations. The best split of an existing cell population is found to
#'  create the K-th cluster. This procedure is much faster than randomly
#'  initializing each model with a different K. If module labels for each
#'  feature are given in 'yInit', the `celda_CG` model will be used to split
#'  cell populations based on those modules instead of individual features.
#'  Module labels will also be updated during sampling and thus may end up
#'  slightly different than `yInit`.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells.
#' @param sampleLabel Vector or factor. Denotes the sample label for each cell
#'  (column) in the count matrix.
#' @param initialK Integer. Minimum number of cell populations to try.
#' @param maxK Integer. Maximum number of cell populations to try.
#' @param tempL Integer. Number of temporary modules to identify and use in cell
#'  splitting. Only used if `yInit = NULL`. Collapsing features to a relatively
#'  smaller number of modules will increase the speed of clustering and tend to
#'  produce better cell populations. This number should be larger than the
#'  number of true modules expected in the dataset. Default NULL.
#' @param yInit Integer vector. Module labels for features. Cells will be
#'  clustered using the `celda_CG` model based on the modules specified in
#'  `yInit` rather than the counts of individual features. While the features
#'  will be initialized to the module labels in `yInit`, the labels will be
#'  allowed to move within each new model with a different K.
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount
#'  to each cell population in each sample. Default 1.
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to
#'  each feature in each cell (if `yInit` is NULL) or to each module in each
#'  cell population (if `yInit` is set). Default 1.
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount
#'  to each feature in each module. Only used if `yInit` is set. Default 1.
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount
#'  to the number of features in each module. Only used if `yInit` is set.
#'  Default 1.
#' @param minCell Integer. Only attempt to split cell populations with at
#'  least this many cells.
#' @param reorder Logical. Whether to reorder cell populations using
#'  hierarchical clustering after each model has been created. If FALSE, cell
#'  populations numbers will correspond to the split which created the cell
#'  populations (i.e. 'K15' was created at split 15, 'K16' was created at split
#'  16, etc.). Default TRUE.
#' @param perplexity Logical. Whether to calculate perplexity for each model.
#'  If FALSE, then perplexity can be calculated later with
#'  `resamplePerplexity()`. Default TRUE.
#' @param verbose Logical. Whether to print log messages. Default TRUE.
#' @param logfile Character. Messages will be redirected to a file named
#'  `logfile`. If NULL, messages will be printed to stdout.  Default NULL.
#' @return Object of class `celda_list`, which contains results for all model
#'  parameter combinations and summaries of the run parameters. The models in
#'  the list will be of class `celda_C` if `yInit = NULL` or `celda_CG` if
#'  `zInit` is set.
#' @seealso `recursiveSplitModule()` for recursive splitting of cell
#'  populations.
#' @examples
#' data(celdaCGSim, celdaCSim)
#' ## Create models that range from K = 3 to K = 7 by recursively splitting
#' ## cell populations into two to produce `celda_C` cell clustering models
#' testZ <- recursiveSplitCell(celdaCSim$counts, initialK = 3, maxK = 7)
#'
#' ## Alternatively, first identify features modules using
#' ## `recursiveSplitModule()`
#' moduleSplit <- recursiveSplitModule(celdaCGSim$counts,
#'     initialL = 3, maxL = 15)
#' plotGridSearchPerplexity(moduleSplit)
#' moduleSplitSelect <- subsetCeldaList(moduleSplit, list(L = 10))
#'
#' ## Then use module labels for initialization in `recursiveSplitCell()` to
#' ## produce `celda_CG` bi-clustering models
#' cellSplit <- recursiveSplitCell(celdaCGSim$counts,
#'     initialK = 3, maxK = 7, yInit = clusters(moduleSplitSelect)$y)
#' plotGridSearchPerplexity(cellSplit)
#' celdaMod <- subsetCeldaList(cellSplit, list(K = 5, L = 10))
#' @export
recursiveSplitCell <- function(counts,
    sampleLabel = NULL,
    initialK = 5,
    maxK = 25,
    tempL = NULL,
    yInit = NULL,
    alpha = 1,
    beta = 1,
    delta = 1,
    gamma = 1,
    minCell = 3,
    reorder = TRUE,
    perplexity = TRUE,
    logfile = NULL,
    verbose = TRUE) {

    .logMessages(paste(rep("=", 50), collapse = ""),
        logfile = logfile,
        append = FALSE,
        verbose = verbose)
    .logMessages("Starting recursive cell population splitting.",
        logfile = logfile,
        append = TRUE,
        verbose = verbose)
    .logMessages(paste(rep("=", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)

    startTime <- Sys.time()
    counts <- .processCounts(counts)
    countChecksum <- .createCountChecksum(counts)

    sampleLabel <- .processSampleLabels(sampleLabel, numCells = ncol(counts))
    s <- as.integer(sampleLabel)
    names <- list(row = rownames(counts),
        column = colnames(counts),
        sample = levels(sampleLabel))

    if (!is.null(yInit)) {
        # Create collapsed module matrix
        L <- length(unique(yInit))
        .logMessages(date(),
            ".. Collapsing to",
            L,
            "modules",
            append = TRUE,
            verbose = verbose,
            logfile = logfile)
        overallY <- .initializeCluster(L, nrow(counts), initial = yInit)
        countsY <- .rowSumByGroup(counts, overallY, L)

        # Create initial model with initialK and predifined y labels
        .logMessages(date(),
            ".. Initializing with",
            initialK,
            "populations",
            append = TRUE,
            verbose = verbose,
            logfile = logfile)
        modelInitial <- .celda_CG(counts,
            sampleLabel = s,
            K = as.integer(initialK),
            L = as.integer(L),
            zInitialize = "split",
            yInitialize = "predefined",
            nchains = 1,
            yInit = overallY,
            alpha = alpha,
            beta = beta,
            gamma = gamma,
            delta = delta,
            verbose = FALSE,
            reorder = reorder)
        currentK <- length(unique(clusters(modelInitial)$z)) + 1
        overallZ <- clusters(modelInitial)$z
        resList <- list(modelInitial)
        while (currentK <= maxK) {
            # previousY <- overallY
            tempSplit <- .singleSplitZ(countsY,
                overallZ,
                s,
                currentK,
                minCell = 3,
                alpha = alpha,
                beta = beta)
            tempModel <- .celda_CG(counts,
                sampleLabel = s,
                K = as.integer(currentK),
                L = as.integer(L),
                yInit = overallY,
                zInit = tempSplit$z,
                nchains = 1,
                zInitialize = "predefined",
                yInitialize = "predefined",
                splitOnLast = FALSE,
                stopIter = 5,
                alpha = alpha,
                beta = beta,
                gamma = gamma,
                delta = delta,
                verbose = FALSE,
                reorder = reorder)

            # Calculate new decomposed counts matrix with new module labels
            # overallY = clusters(tempModel)$y
            # p = .cGReDecomposeCounts(counts, overallY, previousY, countsY,
            # nByG, L = as.integer(L))
            # countsY = p$nTSByC

            # If the number of clusters is still "currentK", then keep the
            # reordering, otherwise keep the previous configuration
            if (length(unique(clusters(tempModel)$z)) == currentK) {
                overallZ <- clusters(tempModel)$z
            } else {
                overallZ <- tempSplit$z
                ll <- logLikelihood.celda_CG(counts,
                    s,
                    overallZ,
                    clusters(tempModel)$y,
                    currentK,
                    L,
                    alpha,
                    beta,
                    delta,
                    gamma
                )
                tempModel <- methods::new("celda_CG",
                    clusters = list(z = overallZ, y = clusters(tempModel)$y),
                    params = list(K = as.integer(currentK),
                        L = as.integer(L),
                        alpha = alpha,
                        beta = beta,
                        delta = delta,
                        gamma = gamma,
                        countChecksum = countChecksum
                    ),
                    finalLogLik = ll,
                    sampleLabel = sampleLabel,
                    names = names
                )
            }

            resList <- c(resList, list(tempModel))
            .logMessages(date(),
                ".. Current cell population",
                currentK,
                "| logLik:",
                bestLogLikelihood(tempModel),
                append = TRUE,
                verbose = verbose,
                logfile = logfile
            )
            currentK <- length(unique(overallZ)) + 1
        }

        runK <- vapply(resList, function(mod) {
            params(mod)$K
        }, integer(1))
        runL <- vapply(resList, function(mod) {
            params(mod)$L
        }, integer(1))
        runParams <- data.frame(index = seq.int(1, length(resList)),
            L = as.integer(runL),
            K = as.integer(runK),
            stringsAsFactors = FALSE)
    } else if (!is.null(tempL)) {
        L <- tempL
        .logMessages(date(),
            ".. Collapsing to",
            L,
            "temporary modules",
            append = TRUE,
            verbose = verbose,
            logfile = logfile
        )
        tempY <- .initializeSplitY(counts,
            L = as.integer(L),
            tempK = max(100, maxK),
            minFeature = 3)
        tempY <- as.integer(as.factor(tempY))
        L <- length(unique(tempY)) # Recalculate in case some modules are empty
        countsY <- .rowSumByGroup(counts, tempY, L)

        # Create initial model with initialK
        .logMessages(date(),
            ".. Initializing with",
            initialK,
            "populations",
            append = TRUE,
            verbose = verbose,
            logfile = logfile
        )
        modelInitial <- .celda_C(countsY,
            sampleLabel = s,
            K = as.integer(initialK),
            zInitialize = "split",
            nchains = 1,
            alpha = alpha,
            beta = beta,
            verbose = FALSE,
            reorder = reorder)
        currentK <- length(unique(clusters(modelInitial)$z)) + 1
        overallZ <- clusters(modelInitial)$z
        ll <- logLikelihood.celda_C(counts, s, overallZ, currentK,
            alpha, beta)
        params(modelInitial)$countChecksum <- countChecksum
        logLikelihoodHistory(modelInitial) <- ll
        bestLogLikelihood(modelInitial) <- ll

        resList <- list(modelInitial)
        while (currentK <= maxK) {
            # Find next best split, then do a new celda_C run with that split
            tempSplit <- .singleSplitZ(countsY,
                overallZ,
                s,
                currentK,
                minCell = 3,
                alpha = alpha,
                beta = beta)
            tempModel <- .celda_C(countsY,
                sampleLabel = s,
                K = as.integer(currentK),
                nchains = 1,
                zInitialize = "random",
                alpha = alpha,
                beta = beta,
                stopIter = 5,
                splitOnLast = FALSE,
                verbose = FALSE,
                zInit = tempSplit$z,
                reorder = reorder)

            # Handle rare cases where a population has no cells after running
            # the model
            if (length(unique(clusters(tempModel)$z)) == currentK) {
                overallZ <- clusters(tempModel)$z
            } else {
                overallZ <- tempSplit$z
            }

            # Need to change below line to use decompose counts to save time
            ll <- logLikelihood.celda_C(counts, s, overallZ, currentK,
                alpha, beta)
            tempModel <- methods::new("celda_C",
                clusters = list(z = overallZ),
                params = list(K = as.integer(currentK),
                    alpha = alpha,
                    beta = beta,
                    countChecksum = countChecksum),
                finalLogLik = ll,
                sampleLabel = sampleLabel,
                names = names
            )

            resList <- c(resList, list(tempModel))
            .logMessages(date(),
                ".. Current cell population",
                currentK,
                "| logLik:",
                bestLogLikelihood(tempModel),
                append = TRUE,
                verbose = verbose,
                logfile = logfile
            )
            currentK <- length(unique(overallZ)) + 1
        }

        runK <- vapply(resList, function(mod) {
            params(mod)$K
        }, integer(1))
        runParams <- data.frame(index = seq.int(1, length(resList)),
            K = as.integer(runK),
            stringsAsFactors = FALSE
        )
    } else {
        # Create initial model with initialK
        .logMessages(date(),
            ".. Initializing with",
            initialK,
            "populations",
            append = TRUE,
            verbose = verbose,
            logfile = logfile)
        modelInitial <- .celda_C(counts,
            sampleLabel = s,
            K = as.integer(initialK),
            zInitialize = "split",
            nchains = 1,
            alpha = alpha,
            beta = beta,
            verbose = FALSE,
            reorder = reorder)
        currentK <- length(unique(clusters(modelInitial)$z)) + 1
        overallZ <- clusters(modelInitial)$z
        resList <- list(modelInitial)
        while (currentK <= maxK) {
            tempSplit <- .singleSplitZ(counts,
                overallZ,
                s,
                currentK,
                minCell = 3,
                alpha = alpha,
                beta = beta)
            tempModel <- .celda_C(counts,
                sampleLabel = s,
                K = as.integer(currentK),
                nchains = 1,
                zInitialize = "random",
                alpha = alpha,
                beta = beta,
                stopIter = 5,
                splitOnLast = FALSE,
                verbose = FALSE,
                zInit = tempSplit$z,
                reorder = reorder)

            if (length(unique(clusters(tempModel)$z)) == currentK) {
                overallZ <- clusters(tempModel)$z
            } else {
                overallZ <- tempSplit$z
                ll <-
                    logLikelihood.celda_C(counts, s, overallZ,
                        currentK, alpha, beta)
                tempModel <- methods::new("celda_C",
                    clusters = list(z = overallZ),
                    params = list(K = as.integer(currentK),
                        alpha = alpha,
                        beta = beta,
                        countChecksum = countChecksum),
                    finalLogLik = ll,
                    sampleLabel = sampleLabel,
                    names = names
                )
            }

            resList <- c(resList, list(tempModel))
            .logMessages(date(),
                ".. Current cell population",
                currentK,
                "| logLik:",
                bestLogLikelihood(tempModel),
                append = TRUE,
                verbose = verbose,
                logfile = logfile
            )
            currentK <- length(unique(overallZ)) + 1
        }

        runK <- vapply(resList, function(mod) {
            params(mod)$K
        }, integer(1))
        runParams <- data.frame(index = seq.int(1, length(resList)),
            K = as.integer(runK),
            stringsAsFactors = FALSE
        )
    }

    # Summarize paramters of different models
    logliks <- vapply(resList, function(mod) {
        bestLogLikelihood(mod)
    }, double(1))
    runParams <- data.frame(runParams,
        log_likelihood = logliks,
        stringsAsFactors = FALSE)

    celdaRes <- methods::new("celdaList",
        runParams = runParams,
        resList = resList,
        countChecksum = countChecksum
    )

    if (isTRUE(perplexity)) {
        .logMessages(date(),
            ".. Calculating perplexity",
            append = TRUE,
            verbose = verbose,
            logfile = NULL
        )
        celdaRes <- resamplePerplexity(counts, celdaRes)
    }
    endTime <- Sys.time()
    .logMessages(
        paste(rep("=", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose
    )
    .logMessages(
        "Completed recursive cell population splitting. Total time:",
        format(difftime(endTime, startTime)),
        logfile = logfile,
        append = TRUE,
        verbose = verbose
    )
    .logMessages(
        paste(rep("=", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose
    )

    return(celdaRes)
}


#' @title Recursive module splitting
#'
#' @description Uses the `celda_G` model to cluster features into modules for
#'  a range of possible L's. The module labels of the previous "L-1" model are
#'  used as the initial values in the current model with L modules. The best
#'  split of an existing module is found to create the L-th module. This
#'  procedure is much faster than randomly initializing each model with a
#'  different L.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells.
#' @param initialL Integer. Minimum number of modules to try.
#' @param maxL Integer. Maximum number of modules to try.
#' @param tempK Integer. Number of temporary cell populations to identify and
#'  use in module splitting. Only used if `zInit=NULL` Collapsing cells to a
#'  relatively smaller number of cell popluations will increase the speed of
#'  module clustering and tend to produce better modules. This number should be
#'  larger than the number of true cell populations expected in the dataset.
#'  Default 100.
#' @param zInit Integer vector. Collapse cells to cell populations based on
#'  labels in `zInit` and then perform module splitting. If NULL, no
#'  collapasing will be performed unless `tempK` is specified. Default NULL.
#' @param sampleLabel Vector or factor. Denotes the sample label for each cell
#'  (column) in the count matrix. Only used if `zInit` is set.
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount
#'  to each cell population in each sample. Only used if `zInit` is set.
#'  Default 1.
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount
#'  to each feature module in each cell. Default 1.
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount
#'  to each feature in each module. Default 1.
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount
#'  to the number of features in each module. Default 1.
#' @param minFeature Integer. Only attempt to split modules with at least this
#'  many features.
#' @param reorder Logical. Whether to reorder modules using hierarchical
#'  clustering after each model has been created. If FALSE, module numbers will
#'  correspond to the split which created the module (i.e. 'L15' was created at
#'  split 15, 'L16' was created at split 16, etc.). Default TRUE.
#' @param perplexity Logical. Whether to calculate perplexity for each model.
#'  If FALSE, then perplexity can be calculated later with
#'  `resamplePerplexity()`. Default TRUE.
#' @param verbose Logical. Whether to print log messages. Default TRUE.
#' @param logfile Character. Messages will be redirected to a file named
#'  `logfile`. If NULL, messages will be printed to stdout.  Default NULL.
#' @return Object of class `celda_list`, which contains results for all model
#'  parameter combinations and summaries of the run parameters. The models in
#'  the list will be of class `celda_G` if `zInit=NULL` or `celda_CG` if
#'  `zInit` is set.
#' @seealso `recursiveSplitCell()` for recursive splitting of cell populations.
#' @examples
#' data(celdaCGSim)
#' ## Create models that range from L=3 to L=20 by recursively splitting modules
#' ## into two
#' moduleSplit <- recursiveSplitModule(celdaCGSim$counts,
#'     initialL = 3, maxL = 20)
#'
#' ## Example results with perplexity
#' plotGridSearchPerplexity(moduleSplit)
#'
#' ## Select model for downstream analysis
#' celdaMod <- subsetCeldaList(moduleSplit, list(L = 10))
#' @export
recursiveSplitModule <- function(counts,
    initialL = 10,
    maxL = 100,
    tempK = 100,
    zInit = NULL,
    sampleLabel = NULL,
    alpha = 1,
    beta = 1,
    delta = 1,
    gamma = 1,
    minFeature = 3,
    reorder = TRUE,
    perplexity = TRUE,
    verbose = TRUE,
    logfile = NULL) {

    .logMessages(paste(rep("=", 50), collapse = ""),
        logfile = logfile,
        append = FALSE,
        verbose = verbose)
    .logMessages("Starting recursive module splitting.",
        logfile = logfile,
        append = TRUE,
        verbose = verbose)
    .logMessages(paste(rep("=", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)
    startTime <- Sys.time()

    counts <- .processCounts(counts)
    countChecksum <- .createCountChecksum(counts)
    names <-
        list(row = rownames(counts), column = colnames(counts))
    sampleLabel <- .processSampleLabels(sampleLabel,
        numCells = ncol(counts))
    s <- as.integer(sampleLabel)

    if (!is.null(zInit)) {
        # Create collapsed module matrix
        K <- length(unique(zInit))
        .logMessages(
            date(),
            ".. Collapsing to",
            K,
            "cell populations",
            append = TRUE,
            verbose = verbose,
            logfile = logfile
        )
        overallZ <- .initializeCluster(N = K,
            len = ncol(counts),
            initial = zInit)
        countsZ <- .colSumByGroup(counts, overallZ, K)

        # Create initial model with initialL and predifined z labels
        .logMessages(date(),
            ".. Initializing with",
            initialL,
            "modules",
            append = TRUE,
            verbose = verbose,
            logfile = logfile)
        modelInitial <- .celda_CG(counts,
            sampleLabel = s,
            L = initialL,
            K = K,
            zInitialize = "predefined",
            yInitialize = "split",
            nchains = 1,
            zInit = overallZ,
            alpha = alpha,
            beta = beta,
            gamma = gamma,
            delta = delta,
            verbose = FALSE,
            reorder = reorder)
        currentL <- length(unique(clusters(modelInitial)$y)) + 1
        overallY <- clusters(modelInitial)$y

        resList <- list(modelInitial)
        while (currentL <= maxL) {
            # Allow features to cluster further with celda_CG
            tempSplit <- .singleSplitY(
                countsZ,
                overallY,
                currentL,
                minFeature = 3,
                beta = beta,
                delta = delta,
                gamma = gamma)
            tempModel <- .celda_CG(
                counts,
                L = currentL,
                K = K,
                stopIter = 5,
                splitOnIter = -1,
                splitOnLast = FALSE,
                nchains = 1,
                verbose = FALSE,
                yInitialize = "predefined",
                zInitialize = "predefined",
                yInit = tempSplit$y,
                zInit = overallZ,
                reorder = reorder
            )
            overallY <- clusters(tempModel)$y

            ## Add new model to results list and increment L
            .logMessages(
                date(),
                ".. Created module",
                currentL,
                "| logLik:",
                bestLogLikelihood(tempModel),
                append = TRUE,
                verbose = verbose,
                logfile = NULL
            )
            resList <- c(resList, list(tempModel))
            currentL <- currentL + 1
        }

        runK <- vapply(resList, function(mod) {
            params(mod)$K
        }, integer(1))
        runL <- vapply(resList, function(mod) {
            params(mod)$L
        }, integer(1))
        runParams <- data.frame(
            index = seq.int(1, length(resList)),
            L = runL,
            K = runK,
            stringsAsFactors = FALSE
        )
    } else if (!is.null(tempK)) {
        K <- tempK
        .logMessages(
            date(),
            ".. Collapsing to",
            K,
            "temporary cell populations",
            append = TRUE,
            verbose = verbose,
            logfile = logfile
        )
        z <- .initializeSplitZ(counts,
            K = K,
            minCell = 3)
        countsZ <- .colSumByGroup(counts, z, length(unique(z)))

        .logMessages(
            date(),
            ".. Initializing with",
            initialL,
            "modules",
            append = TRUE,
            verbose = verbose,
            logfile = logfile
        )
        modelInitial <- .celda_G(
            countsZ,
            L = initialL,
            nchains = 1,
            verbose = FALSE)

        currentL <- length(unique(clusters(modelInitial)$y)) + 1
        overallY <- clusters(modelInitial)$y

        ## Decomposed counts for full count matrix
        p <- .cGDecomposeCounts(counts, overallY, currentL)
        nTSByC <- p$nTSByC
        nByTS <- p$nByTS
        nGByTS <- p$nGByTS
        nByG <- p$nByG
        nG <- p$nG
        nM <- p$nM

        resList <- list(modelInitial)
        while (currentL <= maxL) {
            # Allow features to cluster further
            previousY <- overallY
            tempSplit <- .singleSplitY(
                countsZ,
                overallY,
                currentL,
                minFeature = 3,
                beta = beta,
                delta = delta,
                gamma = gamma)
            tempModel <- .celda_G(
                countsZ,
                L = currentL,
                stopIter = 5,
                splitOnIter = -1,
                splitOnLast = FALSE,
                nchains = 1,
                verbose = FALSE,
                yInitialize = "predefined",
                yInit = tempSplit$y,
                reorder = reorder
            )
            overallY <- clusters(tempModel)$y

            # Adjust decomposed count matrices
            p <- .cGReDecomposeCounts(counts,
                overallY,
                previousY,
                nTSByC,
                nByG,
                L = currentL)
            nTSByC <- p$nTSByC
            nByTS <- p$nByTS
            nGByTS <- p$nGByTS
            previousY <- overallY

            ## Create the final model object with correct info on full counts
            ## matrix
            bestLogLikelihood(tempModel) <- .cGCalcLL(
                nTSByC = nTSByC,
                nByTS = nByTS,
                nByG = nByG,
                nGByTS = nGByTS,
                nM = nM,
                nG = nG,
                L = currentL,
                beta = beta,
                delta = delta,
                gamma = gamma
            )
            logLikelihoodHistory(tempModel) <- bestLogLikelihood(tempModel)
            params(tempModel)$countChecksum <- countChecksum
            matrixNames(tempModel) <- names

            ## Add extra row/column for next round of L
            nTSByC <- rbind(nTSByC, rep(0L, ncol(nTSByC)))
            nByTS <- c(nByTS, 0L)
            nGByTS <- c(nGByTS, 0L)

            ## Add new model to results list and increment L
            .logMessages(
                date(),
                ".. Created module",
                currentL,
                "| logLik:",
                bestLogLikelihood(tempModel),
                append = TRUE,
                verbose = verbose,
                logfile = NULL
            )
            resList <- c(resList, list(tempModel))
            currentL <- currentL + 1
        }

        runL <- vapply(resList, function(mod) {
            params(mod)$L
        }, integer(1))
        runParams <- data.frame(
            index = seq.int(1, length(resList)),
            L = runL,
            stringsAsFactors = FALSE
        )
    } else {
        .logMessages(
            date(),
            ".. Initializing with",
            initialL,
            "modules",
            append = TRUE,
            verbose = verbose,
            logfile = logfile
        )
        modelInitial <- .celda_G(
            counts,
            L = initialL,
            maxIter = 20,
            nchains = 1,
            verbose = FALSE)
        overallY <- clusters(modelInitial)$y
        currentL <- length(unique(overallY)) + 1

        ## Perform splitting for y labels
        resList <- list(modelInitial)
        while (currentL <= maxL) {
            # Allow features to cluster further
            previousY <- overallY
            tempSplit <- .singleSplitY(
                counts,
                overallY,
                currentL,
                minFeature = 3,
                beta = beta,
                delta = delta,
                gamma = gamma)
            tempModel <- .celda_G(
                counts,
                L = currentL,
                stopIter = 5,
                splitOnIter = -1,
                splitOnLast = FALSE,
                nchains = 1,
                verbose = FALSE,
                yInitialize = "predefined",
                yInit = tempSplit$y,
                reorder = reorder
            )
            overallY <- clusters(tempModel)$y

            ## Add new model to results list and increment L
            .logMessages(
                date(),
                ".. Created module",
                currentL,
                "| logLik:",
                bestLogLikelihood(tempModel),
                append = TRUE,
                verbose = verbose,
                logfile = NULL
            )
            resList <- c(resList, list(tempModel))
            currentL <- currentL + 1
        }

        runL <- vapply(resList, function(mod) {
            params(mod)$L
        }, integer(1))
        runParams <- data.frame(
            index = seq.int(1, length(resList)),
            L = runL,
            stringsAsFactors = FALSE
        )
    }

    ## Summarize paramters of different models
    logliks <- vapply(resList, function(mod) {
        bestLogLikelihood(mod)
    }, double(1))
    runParams <- data.frame(runParams,
        log_likelihood = logliks,
        stringsAsFactors = FALSE)

    celdaRes <- methods::new(
        "celdaList",
        runParams = runParams,
        resList = resList,
        countChecksum = countChecksum
    )

    if (isTRUE(perplexity)) {
        .logMessages(
            date(),
            ".. Calculating perplexity",
            append = TRUE,
            verbose = verbose,
            logfile = NULL
        )
        celdaRes <- resamplePerplexity(counts, celdaRes)
    }

    endTime <- Sys.time()
    .logMessages(paste(rep("=", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)
    .logMessages("Completed recursive module splitting. Total time:",
        format(difftime(endTime, startTime)),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)
    .logMessages(paste(rep("=", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)

    return(celdaRes)
}
