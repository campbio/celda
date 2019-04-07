.singleSplitZ <- function(counts,
    z,
    s,
    K,
    min.cell = 3,
    alpha = 1,
    beta = 1,
    seed) {

    z.ta <- tabulate(z, K)
    z.to.split <- which(z.ta > min.cell)
    best.z <- z
    best.ll <- -Inf
    for (i in z.to.split) {
        clustLabel <- .celda_C(
            counts[, z == i, drop = FALSE],
            K = 2,
            z.initialize = "random",
            split.on.iter = -1,
            split.on.last = FALSE,
            verbose = FALSE,
            seed = seed
        )
        if (length(unique(clustLabel@clusters$z)) == 2) {
            ix <- z == i
            new.z <- z
            new.z[ix] <- ifelse(clustLabel@clusters$z == 2, i, K)
            ll <-
                logLikelihood.celda_C(counts, s, new.z, K, alpha, beta)

            if (ll > best.ll) {
                best.z <- new.z
                best.ll <- ll
            }
        }
    }
    return(list(ll = best.ll, z = best.z))
}

.singleSplitY <- function(counts,
    y,
    L,
    min.feature = 3,
    beta = 1,
    delta = 1,
    gamma = 1,
    seed) {
    y.ta <- tabulate(y, L)
    y.to.split <- which(y.ta > min.feature)

    best.y <- y
    best.ll <- -Inf
    previous.y <- y
    for (i in y.to.split) {
        clustLabel <- .celda_G(counts[y == i, , drop = FALSE],
            L = 2,
            y.initialize = "random",
            split.on.iter = -1,
            split.on.last = FALSE,
            nchains = 1,
            verbose = FALSE,
            seed = seed
        )

        if (length(unique(clustLabel@clusters$y)) == 2) {
            ix <- y == i
            new.y <- y
            new.y[ix] <- ifelse(clustLabel@clusters$y == 2, i, L)
            ll <- logLikelihood.celda_G(counts, new.y, L, beta, delta, gamma)

            if (ll > best.ll) {
                best.y <- new.y
                best.ll <- ll
            }
        }
    }
    return(list(ll = best.ll, y = best.y))
}

#' @title Recursive cell splitting
#' @description Uses the `celda_C` model to cluster cells into population for range of possible K's. The cell population labels of the previous "K-1" model are used as the initial values in the current model with K cell populations. The best split of an existing cell population is found to create the K-th cluster. This procedure is much faster than randomly initializing each model with a different K. If module labels for each feature are given in 'y.init', the `celda_CG` model will be used to split cell populations based on those modules instead of individual features. Module labels will also be updated during sampling and thus may end up slightly different than `y.init`.
#' @param counts Integer matrix. Rows represent features and columns represent cells.
#' @param sample.label Vector or factor. Denotes the sample label for each cell (column) in the count matrix.
#' @param initial.K Integer. Minimum number of cell populations to try.
#' @param max.K Integer. Maximum number of cell populations to try.
#' @param temp.L Integer. Number of temporary modules to identify and use in cell splitting. Only used if `y.init=NULL`. Collapsing features to a relatively smaller number of modules will increase the speed of clustering and tend to produce better cell populations. This number should be larger than the number of true modules expected in the dataset. Default NULL.
#' @param y.init Integer vector. Module labels for features. Cells will be clustered using the `celda_CG` model based on the modules specified in `y.init` rather than the counts of individual features. While the features will be initialized to the module labels in `y.init`, the labels will be allowed to move within each new model with a different K.
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount to each cell population in each sample. Default 1.
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to each feature in each cell (if `y.init` is NULL) or to each module in each cell population (if `y.init` is set). Default 1.
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to each feature in each module. Only used if `y.init` is set. Default 1.
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to the number of features in each module. Only used if `y.init` is set. Default 1.
#' @param min.cell Integer. Only attempt to split cell populations with at least this many cells.
#' @param reorder Logical. Whether to reorder cell populations using hierarchical clustering after each model has been created. If FALSE, cell populations numbers will correspond to the split which created the cell populations (i.e. 'K15' was created at split 15, 'K16' was created at split 16, etc.). Default TRUE.
#' @param perplexity Logical. Whether to calculate perplexity for each model. If FALSE, then perplexity can be calculated later with `resamplePerplexity()`. Default TRUE.
#' @param seed Integer. Passed to `set.seed()`. Default 12345. If NULL, no calls to `set.seed()` are made.
#' @param verbose Logical. Whether to print log messages. Default TRUE.
#' @param logfile Character. Messages will be redirected to a file named `logfile`. If NULL, messages will be printed to stdout.  Default NULL.
#' @return Object of class `celda_list`, which contains results for all model parameter combinations and summaries of the run parameters. The models in the list will be of class `celda_C` if `y.init=NULL` or `celda_CG` if `z.init` is set.
#' @seealso `recursiveSplitModule()` for recursive splitting of cell populations.
#' @examples
#' ## Create models that range from K=3 to K=10 by recursively splitting cell populations into two to produce `celda_C` cell clustering models
#' testZ <- recursiveSplitCell(celda.C.sim$counts, initial.K = 3, max.K = 10)
#'
#' ## Alternatively, first identify features modules usinge `recursiveSplitModule()`
#' module.split <- recursiveSplitModule(celda.CG.sim$counts, initial.L = 3, max.L = 20)
#' plotGridSearchPerplexity(module.split)
#' module.split.select <- subsetCeldaList(module.split, list(L = 10))
#'
#' ## Then use module labels for initialization in `recursiveSplitCell()` to produce `celda_CG` bi-clustering models
#' cell.split <- recursiveSplitCell(celda.CG.sim$counts, initial.K = 3, max.K = 20, y.init = clusters(module.split.select)$y)
#' plotGridSearchPerplexity(cell.split)
#' celda.mod <- subsetCeldaList(cell.split, list(K = 5, L = 10))
#' @export
recursiveSplitCell <- function(counts,
    sample.label = NULL,
    initial.K = 5,
    max.K = 25,
    temp.L = NULL,
    y.init = NULL,
    alpha = 1,
    beta = 1,
    delta = 1,
    gamma = 1,
    min.cell = 3,
    reorder = TRUE,
    perplexity = TRUE,
    seed = 12345,
    logfile = NULL,
    verbose = TRUE) {

    logMessages(
        paste(rep("=", 50), collapse = "" ),
        logfile = logfile,
        append = FALSE,
        verbose = verbose
    )
    logMessages(
        "Starting recursive cell population splitting.",
        logfile = logfile,
        append = TRUE,
        verbose = verbose
    )
    logMessages(
        paste(rep("=", 50), collapse = "" ),
        logfile = logfile,
        append = TRUE,
        verbose = verbose
    )

    start.time <- Sys.time()
    counts <- processCounts(counts)
    count.checksum <- createCountChecksum(counts)

    sample.label <- processSampleLabels(sample.label,
        num.cells = ncol(counts))
    s <- as.integer(sample.label)
    names <- list(
        row = rownames(counts),
        column = colnames(counts),
        sample = levels(sample.label)
    )

    if (!is.null(y.init)) {
        # Create collapsed module matrix
        L <- length(unique(y.init))
        logMessages(date(),
            ".. Collapsing to",
            L,
            "modules",
            append = TRUE,
            verbose = verbose,
            logfile = logfile
        )
        overall.y <- initializeCluster(L, nrow(counts), initial = y.init,
                seed = seed)
        counts.y <- .rowSumByGroup(counts, overall.y, L)

        # Create initial model with initial.K and predifined y labels
        logMessages(date(),
            ".. Initializing with",
            initial.K,
            "populations",
            append = TRUE,
            verbose = verbose,
            logfile = logfile
        )
        model.initial <- .celda_CG(counts,
            sample.label = s,
            K = initial.K,
            L = L,
            z.initialize = "split",
            y.initialize = "predefined",
            nchains = 1,
            y.init = overall.y,
            alpha = alpha,
            beta = beta,
            gamma = gamma,
            delta = delta,
            verbose = FALSE,
            seed = seed,
            reorder = reorder
        )
        current.K <- length(unique(model.initial@clusters$z)) + 1
        overall.z <- model.initial@clusters$z
        res.list <- list(model.initial)
        while (current.K <= max.K) {
            previous.y <- overall.y
            temp.split <- singleSplitZ(counts.y,
                overall.z,
                s,
                current.K,
                min.cell = 3,
                alpha = alpha,
                beta = beta,
                seed = seed
            )
            temp.model <- .celda_CG(counts,
                sample.label = s,
                K = current.K,
                L = L,
                y.init = overall.y,
                z.init = temp.split$z,
                nchains = 1,
                z.initialize = "predefined",
                y.initialize = "predefined",
                split.on.last = FALSE,
                stop.iter = 5,
                alpha = alpha,
                beta = beta,
                gamma = gamma,
                delta = delta,
                verbose = FALSE,
                seed = seed,
                reorder = reorder
            )

            # Calculate new decomposed counts matrix with new module labels
            # overall.y = temp.model@clusters$y
            # p = cG.reDecomposeCounts(counts, overall.y, previous.y, counts.y, n.by.G, L = L)
            # counts.y = p$n.TS.by.C

            # If the number of clusters is still "current.K", then keep the reordering, otherwise keep the previous configuration
            if (length(unique(temp.model@clusters$z)) == current.K) {
                overall.z <- temp.model@clusters$z
            } else {
                overall.z <- temp.split$z
                ll <- logLikelihood.celda_CG(counts,
                    s,
                    overall.z,
                    temp.model@clusters$y,
                    current.K,
                    L,
                    alpha,
                    beta,
                    delta,
                    gamma
                )
                temp.model <- methods::new("celda_CG",
                    clusters = list(z = overall.z, y = temp.model@clusters$y),
                    params = list(K = current.K,
                        L = L,
                        alpha = alpha,
                        beta = beta,
                        delta = delta,
                        gamma = gamma,
                        seed = seed,
                        count.checksum = count.checksum
                    ),
                    finalLogLik = ll,
                    sample.label = sample.label,
                    names = names
                )
            }

            res.list <- c(res.list, list(temp.model))
            logMessages(date(),
                ".. Current cell population",
                current.K,
                "| logLik:",
                temp.model@finalLogLik,
                append = TRUE,
                verbose = verbose,
                logfile = logfile
            )
            current.K <- length(unique(overall.z)) + 1
        }

        runK <- sapply(res.list, function(mod) {
            mod@params$K
        })
        runL <- sapply(res.list, function(mod) {
            mod@params$L
        })
        run.params <-data.frame(index = seq.int(1, length(res.list)),
            L = runL,
            K = runK,
            stringsAsFactors = FALSE
            )
    } else if (!is.null(temp.L)) {
        L <- temp.L
        logMessages(date(),
            ".. Collapsing to",
            L,
            "temporary modules",
            append = TRUE,
            verbose = verbose,
            logfile = logfile
        )
        temp.y <- initialize.splitY(counts,
            L = L,
            temp.K = max(100, max.K),
            min.feature = 3,
            seed = seed
        )
        temp.y <- as.integer(as.factor(temp.y))
        L <- length(unique(temp.y)) # Recalculate in case some modules are empty
        counts.y <- .rowSumByGroup(counts, temp.y, L)

        # Create initial model with initial.K
        logMessages(date(),
            ".. Initializing with",
            initial.K,
            "populations",
            append = TRUE,
            verbose = verbose,
            logfile = logfile
        )
        model.initial <- .celda_C(counts.y,
            sample.label = s,
            K = initial.K,
            z.initialize = "split",
            nchains = 1,
            alpha = alpha,
            beta = beta,
            verbose = FALSE,
            seed = seed,
            reorder = reorder
        )
        current.K <- length(unique(model.initial@clusters$z)) + 1
        overall.z <- model.initial@clusters$z
        ll <- logLikelihood.celda_C(counts, s, overall.z, current.K,
                alpha, beta)
        model.initial@params$count.checksum <- count.checksum
        model.initial@completeLogLik <- ll
        model.initial@finalLogLik <- ll

        res.list <- list(model.initial)
        while (current.K <= max.K) {
            # Find next best split, then seed a new celda_C run with that split
            temp.split <- singleSplitZ(counts.y,
                overall.z,
                s,
                current.K,
                min.cell = 3,
                alpha = alpha,
                beta = beta,
                seed = seed
            )
            temp.model <- .celda_C(counts.y,
                sample.label = s,
                K = current.K,
                nchains = 1,
                z.initialize = "random",
                alpha = alpha,
                beta = beta,
                stop.iter = 5,
                split.on.last = FALSE,
                verbose = FALSE,
                z.init = temp.split$z,
                seed = seed,
                reorder = reorder
            )

            # Handle rare cases where a population has no cells after running the model
            if (length(unique(temp.model@clusters$z)) == current.K) {
                overall.z <- temp.model@clusters$z
            } else {
                overall.z <- temp.split$z
            }

            # Need to change below line to use decompose counts to save time
            ll <- logLikelihood.celda_C(counts, s, overall.z, current.K,
                    alpha, beta)
            temp.model <- methods::new("celda_C",
                clusters = list(z = overall.z),
                params = list(K = current.K,
                    alpha = alpha,
                    beta = beta,
                    seed = seed,
                    count.checksum = count.checksum
                ),
                finalLogLik = ll,
                sample.label = sample.label,
                names = names
            )

            res.list <- c(res.list, list(temp.model))
            logMessages(date(),
                ".. Current cell population",
                current.K,
                "| logLik:",
                temp.model@finalLogLik,
                append = TRUE,
                verbose = verbose,
                logfile = logfile
            )
            current.K <- length(unique(overall.z)) + 1
        }

        runK <- sapply(res.list, function(mod) {
            mod@params$K
        })
        run.params <- data.frame(index = seq.int(1, length(res.list)),
                K = runK,
                stringsAsFactors = FALSE
            )
    } else {
        # Create initial model with initial.K
        logMessages(date(),
            ".. Initializing with",
            initial.K,
            "populations",
            append = TRUE,
            verbose = verbose,
            logfile = logfile
        )
        model.initial <- .celda_C(counts,
            sample.label = s,
            K = initial.K,
            z.initialize = "split",
            nchains = 1,
            alpha = alpha,
            beta = beta,
            verbose = FALSE,
            seed = seed,
            reorder = reorder
        )
        current.K <- length(unique(model.initial@clusters$z)) + 1
        overall.z <- model.initial@clusters$z
        res.list <- list(model.initial)
        while (current.K <= max.K) {
            temp.split <- singleSplitZ(counts,
                overall.z,
                s,
                current.K,
                min.cell = 3,
                alpha = alpha,
                beta = beta,
                seed = seed
            )
            temp.model <- .celda_C(counts,
                sample.label = s,
                K = current.K,
                nchains = 1,
                z.initialize = "random",
                alpha = alpha,
                beta = beta,
                stop.iter = 5,
                split.on.last = FALSE,
                verbose = FALSE,
                z.init = temp.split$z,
                seed = seed,
                reorder = reorder
            )

            if (length(unique(temp.model@clusters$z)) == current.K) {
                overall.z <- temp.model@clusters$z
            } else {
                overall.z <- temp.split$z
                ll <-
                    logLikelihood.celda_C(counts, s, overall.z,
                        current.K, alpha, beta)
                temp.model <- methods::new("celda_C",
                    clusters = list(z = overall.z),
                    params = list(K = current.K,
                        alpha = alpha,
                        beta = beta,
                        seed = seed,
                        count.checksum = count.checksum
                    ),
                    finalLogLik = ll,
                    sample.label = sample.label,
                    names = names
                )
            }

            res.list <- c(res.list, list(temp.model))
            logMessages(date(),
                ".. Current cell population",
                current.K,
                "| logLik:",
                temp.model@finalLogLik,
                append = TRUE,
                verbose = verbose,
                logfile = logfile
            )
            current.K <- length(unique(overall.z)) + 1
        }

        runK <- sapply(res.list, function(mod) {
            mod@params$K
        })
        run.params <- data.frame(index = seq.int(1, length(res.list)),
                K = runK,
                stringsAsFactors = FALSE
            )
    }

    # Summarize paramters of different models
    logliks <- sapply(res.list, function(mod) {
        mod@finalLogLik
    })
    run.params <- data.frame(run.params,
        log_likelihood = logliks,
        stringsAsFactors = FALSE)

    celda.res <- methods::new("celdaList",
        run.params = run.params,
        res.list = res.list,
        count.checksum = count.checksum
    )

    if (isTRUE(perplexity)) {
        logMessages(date(),
            ".. Calculating perplexity",
            append = TRUE,
            verbose = verbose,
            logfile = NULL
        )
        celda.res <- resamplePerplexity(counts, celda.res)
    }
    end.time <- Sys.time()
    logMessages(
        paste(rep("=", 50), collapse = "" ),
        logfile = logfile,
        append = TRUE,
        verbose = verbose
    )
    logMessages(
        "Completed recursive cell population splitting. Total time:",
        format(difftime(end.time, start.time)),
        logfile = logfile,
        append = TRUE,
        verbose = verbose
    )
    logMessages(
        paste(rep("=", 50), collapse = "" ),
        logfile = logfile,
        append = TRUE,
        verbose = verbose
    )

    return(celda.res)
}


#' @title Recursive module splitting
#'
#' @description Uses the `celda_G` model to cluster features into modules for a range of possible L's. The module labels of the previous "L-1" model are used as the initial values in the current model with L modules. The best split of an existing module is found to create the L-th module. This procedure is much faster than randomly initializing each model with a different L.
#'
#' @param counts Integer matrix. Rows represent features and columns represent cells.
#' @param initial.L Integer. Minimum number of modules to try.
#' @param max.L Integer. Maximum number of modules to try.
#' @param temp.K Integer. Number of temporary cell populations to identify and use in module splitting. Only used if `z.init=NULL` Collapsing cells to a relatively smaller number of cell popluations will increase the speed of module clustering and tend to produce better modules. This number should be larger than the number of true cell populations expected in the dataset. Default 100.
#' @param z.init Integer vector. Collapse cells to cell populations based on labels in `z.init` and then perform module splitting. If NULL, no collapasing will be performed unless `temp.K` is specified. Default NULL.
#' @param sample.label Vector or factor. Denotes the sample label for each cell (column) in the count matrix. Only used if `z.init` is set.
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount to each cell population in each sample. Only used if `z.init` is set. Default 1.
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to each feature module in each cell. Default 1.
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to each feature in each module. Default 1.
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to the number of features in each module. Default 1.
#' @param min.feature Integer. Only attempt to split modules with at least this many features.
#' @param reorder Logical. Whether to reorder modules using hierarchical clustering after each model has been created. If FALSE, module numbers will correspond to the split which created the module (i.e. 'L15' was created at split 15, 'L16' was created at split 16, etc.). Default TRUE.
#' @param perplexity Logical. Whether to calculate perplexity for each model. If FALSE, then perplexity can be calculated later with `resamplePerplexity()`. Default TRUE.
#' @param seed Integer. Passed to `set.seed()`. Default 12345. If NULL, no calls to `set.seed()` are made.
#' @param verbose Logical. Whether to print log messages. Default TRUE.
#' @param logfile Character. Messages will be redirected to a file named `logfile`. If NULL, messages will be printed to stdout.  Default NULL.
#' @return Object of class `celda_list`, which contains results for all model parameter combinations and summaries of the run parameters. The models in the list will be of class `celda_G` if `z.init=NULL` or `celda_CG` if `z.init` is set.
#' @seealso `recursiveSplitCell()` for recursive splitting of cell populations.
#' @examples
#' ## Create models that range from L=3 to L=20 by recursively splitting modules into two
#' module.split <- recursiveSplitModule(celda.CG.sim$counts, initial.L = 3, max.L = 20)
#'
#' ## Example results with perplexity
#' plotGridSearchPerplexity(module.split)
#'
#' ## Select model for downstream analysis
#' celda.mod <- subsetCeldaList(module.split, list(L = 10))
#' @export
recursiveSplitModule <- function(counts,
    initial.L = 10,
    max.L = 100,
    temp.K = 100,
    z.init = NULL,
    sample.label = NULL,
    alpha = 1,
    beta = 1,
    delta = 1,
    gamma = 1,
    min.feature = 3,
    reorder = TRUE,
    perplexity = TRUE,
    seed = 12345,
    verbose = TRUE,
    logfile = NULL) {
    logMessages(
        paste(rep("=", 50), collapse = "" ),
        logfile = logfile,
        append = FALSE,
        verbose = verbose
    )
    logMessages(
        "Starting recursive module splitting.",
        logfile = logfile,
        append = TRUE,
        verbose = verbose
    )
    logMessages(
        paste(rep("=", 50), collapse = "" ),
        logfile = logfile,
        append = TRUE,
        verbose = verbose
    )
    start.time <- Sys.time()

    counts <- processCounts(counts)
    count.checksum <- createCountChecksum(counts)
    names <-
        list(row = rownames(counts), column = colnames(counts))
    sample.label <- processSampleLabels(sample.label,
        num.cells = ncol(counts))
    s <- as.integer(sample.label)

    if (!is.null(z.init)) {
        # Create collapsed module matrix
        K <- length(unique(z.init))
        logMessages(
            date(),
            ".. Collapsing to",
            K,
            "cell populations",
            append = TRUE,
            verbose = verbose,
            logfile = logfile
        )
        overall.z <- initialize.cluster(
            N = K,
            len = ncol(counts),
            initial = z.init,
            seed = seed
        )
        counts.z <- .colSumByGroup(counts, overall.z, K)

        # Create initial model with initial.L and predifined z labels
        logMessages(
            date(),
            ".. Initializing with",
            initial.L,
            "modules",
            append = TRUE,
            verbose = verbose,
            logfile = logfile
        )
        model.initial <- .celda_CG(
            counts,
            sample.label = s,
            L = initial.L,
            K = K,
            z.initialize = "predefined",
            y.initialize = "split",
            nchains = 1,
            z.init = overall.z,
            alpha = alpha,
            beta = beta,
            gamma = gamma,
            delta = delta,
            verbose = FALSE,
            seed = seed,
            reorder = reorder
        )
        current.L <- length(unique(model.initial@clusters$y)) + 1
        overall.y <- model.initial@clusters$y

        res.list <- list(model.initial)
        while (current.L <= max.L) {
            # Allow features to cluster further with celda_CG
            temp.split <- singleSplitY(
                counts.z,
                overall.y,
                current.L,
                min.feature = 3,
                beta = beta,
                delta = delta,
                gamma = gamma,
                seed = seed
            )
            temp.model <- .celda_CG(
                counts,
                L = current.L,
                K = K,
                stop.iter = 5,
                split.on.iter = -1,
                split.on.last = FALSE,
                nchains = 1,
                verbose = FALSE,
                y.initialize = "predefined",
                z.initialize = "predefined",
                y.init = temp.split$y,
                z.init = overall.z,
                reorder = reorder
            )
            overall.y <- temp.model@clusters$y

            ## Add new model to results list and increment L
            logMessages(
                date(),
                ".. Created module",
                current.L,
                "| logLik:",
                temp.model@finalLogLik,
                append = TRUE,
                verbose = verbose,
                logfile = NULL
            )
            res.list <- c(res.list, list(temp.model))
            current.L <- current.L + 1
        }

        runK <- sapply(res.list, function(mod) {
            mod@params$K
        })
        runL <- sapply(res.list, function(mod) {
            mod@params$L
        })
        run.params <- data.frame(
            index = seq.int(1, length(res.list)),
            L = runL,
            K = runK,
            stringsAsFactors = FALSE
        )
    } else if (!is.null(temp.K)) {
        K <- temp.K
        logMessages(
            date(),
            ".. Collapsing to",
            K,
            "temporary cell populations",
            append = TRUE,
            verbose = verbose,
            logfile = logfile
        )
        z <- initialize.splitZ(counts,
            K = K,
            min.cell = 3,
            seed = seed)
        counts.z <- .colSumByGroup(counts, z, length(unique(z)))

        logMessages(
            date(),
            ".. Initializing with",
            initial.L,
            "modules",
            append = TRUE,
            verbose = verbose,
            logfile = logfile
        )
        model.initial <- .celda_G(
            counts.z,
            L = initial.L,
            nchains = 1,
            verbose = FALSE,
            seed = seed
        )

        current.L <- length(unique(model.initial@clusters$y)) + 1
        overall.y <- model.initial@clusters$y

        ## Decomposed counts for full count matrix
        p <- cG.decomposeCounts(counts, overall.y, current.L)
        n.TS.by.C <- p$n.TS.by.C
        n.by.TS <- p$n.by.TS
        nG.by.TS <- p$nG.by.TS
        n.by.G <- p$n.by.G
        nG <- p$nG
        nM <- p$nM

        res.list <- list(model.initial)
        while (current.L <= max.L) {
            # Allow features to cluster further
            previous.y <- overall.y
            temp.split <- .singleSplitY(
                counts.z,
                overall.y,
                current.L,
                min.feature = 3,
                beta = beta,
                delta = delta,
                gamma = gamma,
                seed = seed
            )
            temp.model <- .celda_G(
                counts.z,
                L = current.L,
                stop.iter = 5,
                split.on.iter = -1,
                split.on.last = FALSE,
                nchains = 1,
                verbose = FALSE,
                y.initialize = "predefined",
                y.init = temp.split$y,
                reorder = reorder
            )
            overall.y <- temp.model@clusters$y

            # Adjust decomposed count matrices
            p <- cG.reDecomposeCounts(counts,
                overall.y,
                previous.y,
                n.TS.by.C,
                n.by.G,
                L = current.L)
            n.TS.by.C <- p$n.TS.by.C
            n.by.TS <- p$n.by.TS
            nG.by.TS <- p$nG.by.TS
            previous.y <- overall.y

            ## Create the final model object with correct info on full counts matrix
            temp.model@finalLogLik <- cG.calcLL(
                n.TS.by.C = n.TS.by.C,
                n.by.TS = n.by.TS,
                n.by.G = n.by.G,
                nG.by.TS = nG.by.TS,
                nM = nM,
                nG = nG,
                L = current.L,
                beta = beta,
                delta = delta,
                gamma = gamma
            )
            temp.model@completeLogLik <- temp.model@finalLogLik
            temp.model@params$count.checksum <- count.checksum
            temp.model@names <- names

            ## Add extra row/column for next round of L
            n.TS.by.C <- rbind(n.TS.by.C, rep(0L, ncol(n.TS.by.C)))
            n.by.TS <- c(n.by.TS, 0L)
            nG.by.TS <- c(nG.by.TS, 0L)

            ## Add new model to results list and increment L
            logMessages(
                date(),
                ".. Created module",
                current.L,
                "| logLik:",
                temp.model@finalLogLik,
                append = TRUE,
                verbose = verbose,
                logfile = NULL
            )
            res.list <- c(res.list, list(temp.model))
            current.L <- current.L + 1
        }

        runL <- sapply(res.list, function(mod) {
            mod@params$L
        })
        run.params <- data.frame(
            index = seq.int(1, length(res.list)),
            L = runL,
            stringsAsFactors = FALSE
        )
    } else {
        logMessages(
            date(),
            ".. Initializing with",
            initial.L,
            "modules",
            append = TRUE,
            verbose = verbose,
            logfile = logfile
        )
        model.initial <- .celda_G(
            counts,
            L = initial.L,
            max.iter = 20,
            nchains = 1,
            verbose = FALSE,
            seed = seed
        )
        overall.y <- model.initial@clusters$y
        current.L <- length(unique(overall.y)) + 1

        ## Perform splitting for y labels
        res.list <- list(model.initial)
        while (current.L <= max.L) {
            # Allow features to cluster further
            previous.y <- overall.y
            temp.split <- singleSplitY(
                counts,
                overall.y,
                current.L,
                min.feature = 3,
                beta = beta,
                delta = delta,
                gamma = gamma,
                seed = seed
            )
            temp.model <- .celda_G(
                counts,
                L = current.L,
                stop.iter = 5,
                split.on.iter = -1,
                split.on.last = FALSE,
                nchains = 1,
                verbose = FALSE,
                y.initialize = "predefined",
                y.init = temp.split$y,
                reorder = reorder
            )
            overall.y <- temp.model@clusters$y

            ## Add new model to results list and increment L
            logMessages(
                date(),
                ".. Created module",
                current.L,
                "| logLik:",
                temp.model@finalLogLik,
                append = TRUE,
                verbose = verbose,
                logfile = NULL
            )
            res.list <- c(res.list, list(temp.model))
            current.L <- current.L + 1
        }

        runL <- sapply(res.list, function(mod) {
            mod@params$L
        })
        run.params <- data.frame(
            index = seq.int(1, length(res.list)),
            L = runL,
            stringsAsFactors = FALSE
        )
    }

    ## Summarize paramters of different models
    logliks <- sapply(res.list, function(mod) {
        mod@finalLogLik
    })
    run.params <- data.frame(run.params,
        log_likelihood = logliks,
        stringsAsFactors = FALSE)

    celda.res <- methods::new(
        "celdaList",
        run.params = run.params,
        res.list = res.list,
        count.checksum = count.checksum
    )

    if (isTRUE(perplexity)) {
        logMessages(
            date(),
            ".. Calculating perplexity",
            append = TRUE,
            verbose = verbose,
            logfile = NULL
        )
        celda.res <- resamplePerplexity(counts, celda.res)
    }

    end.time <- Sys.time()
    logMessages(
        paste(rep("=", 50), collapse = "" ),
        logfile = logfile,
        append = TRUE,
        verbose = verbose
    )
    logMessages(
        "Completed recursive module splitting. Total time:",
        format(difftime(end.time, start.time)),
        logfile = logfile,
        append = TRUE,
        verbose = verbose
    )
    logMessages(
        paste(rep("=", 50), collapse = "" ),
        logfile = logfile,
        append = TRUE,
        verbose = verbose
    )

    return(celda.res)
}
