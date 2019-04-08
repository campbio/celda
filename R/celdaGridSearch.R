#' @title available models
#' @export
availableModels <- c("celda_C", "celda_G", "celda_CG")


#' @title Run Celda in parallel with multiple parameters
#' @description Run Celda with different combinations of parameters and
#'  multiple chains in parallel. The variable `availableModels` contains the
#'  potential models that can be utilized. Different parameters to be tested
#'  should be stored in a list and passed to the argument `paramsTest`. Fixed
#'  parameters to be used in all models, such as `sampleLabel`, can be passed
#'  as a list to the argument `paramsFixed`. When `verbose = TRUE`, output
#'  from each chain will be sent to a log file but not be displayed in stdout.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells.
#' @param model Celda model. Options available in `celda::availableModels`.
#' @param paramsTest List. A list denoting the combinations of parameters to
#'  run in a celda model. For example, `list(K = seq(5, 10), L = seq(15, 20))`
#'  will run all combinations of K from 5 to 10 and L from 15 to 20 in model
#'  `celda_CG()`.
#' @param paramsFixed List. A list denoting additional parameters to use in
#'  each celda model. Default NULL.
#' @param maxIter Integer. Maximum number of iterations of sampling to
#'  perform. Default 200.
#' @param nchains Integer. Number of random cluster initializations. Default 3.
#' @param cores Integer. The number of cores to use for parallel estimation of
#'  chains. Default 1.
#' @param bestOnly Logical. Whether to return only the chain with the highest
#'  log likelihood per combination of parameters or return all chains. Default
#'  TRUE.
#' @param seed Integer. Passed to `set.seed()`. Default 12345. If NULL, no
#'  calls to `set.seed()` are made.
#' @param perplexity Logical. Whether to calculate perplexity for each model.
#'  If FALSE, then perplexity can be calculated later with
#'  `resamplePerplexity()`. Default TRUE.
#' @param verbose Logical. Whether to print log messages during celda chain
#'  execution. Default TRUE.
#' @param logfilePrefix Character. Prefix for log files from worker threads
#'  and main process. Default "Celda".
#' @return Object of class `celdaList`, which contains results for all model
#'  parameter combinations and summaries of the run parameters
#' @seealso `celda_G()` for feature clustering, `celda_C()` for clustering of
#'  cells, and `celda_CG()` for simultaneous clustering of features and cells.
#'  `subsetCeldaList()` can subset the `celdaList` object. `selectBestModel()`
#'  can get the best model for each combination of parameters.
#' @examples
#' ## Run various combinations of parameters with 'celdaGridSearch'
#' cgs <- celdaGridSearch(celdaCGSim$counts,
#'     model = "celda_CG",
#'     paramsTest = list(K = seq(4, 6), L = seq(9, 11)),
#'     paramsFixed = list(sampleLabel = celdaCGSim$sampleLabel),
#'     bestOnly = TRUE,
#'     nchains = 1)
#' @import foreach
#' @export
celdaGridSearch <- function(counts,
    model,
    paramsTest,
    paramsFixed = NULL,
    maxIter = 200,
    nchains = 3,
    cores = 1,
    bestOnly = TRUE,
    seed = 12345,
    perplexity = TRUE,
    verbose = TRUE,
    logfilePrefix = "Celda") {

    ## Check parameters
    .validateCounts(counts)

    modelParams <- as.list(formals(model))
    if (!all(names(paramsTest) %in% names(modelParams))) {
        badParams <- setdiff(names(paramsTest), names(modelParams))
        stop("The following elements in 'paramsTest' are not arguments of '",
            model,
            "': ",
            paste(badParams, collapse = ","))
    }

    if (!is.null(paramsFixed) &&
            !all(names(paramsFixed) %in% names(modelParams))) {
        badParams <- setdiff(names(paramsFixed), names(modelParams))
        stop("The following elements in 'paramsFixed' are not arguments",
            " of '",
            model,
            "': ",
            paste(badParams, collapse = ","))
    }

    modelParamsRequired <- setdiff(names(modelParams[modelParams == ""]),
        "counts")

    if (!all(modelParamsRequired %in% c(names(paramsTest),
        names(paramsFixed)))) {
        missing.params <- setdiff(modelParamsRequired,
            c(names(paramsTest), names(paramsFixed)))
        stop("The following arguments are not in 'paramsTest' or 'paramsFixed'",
            " but are required for '",
            model,
            "': ",
            paste(missing.params, collapse = ","))
    }

    if (any(c("z.init", "y.init", "sampleLabel") %in% names(paramsTest))) {
        stop("Setting parameters such as 'z.init', 'y.init', and 'sampleLabel'",
            " in 'paramsTest' is not currently supported.")
    }

    if (any(c("nchains") %in% names(paramsTest))) {
        warning("Parameter 'nchains' should not be used within the paramsTest",
            " list")
        paramsTest[["nchains"]] <- NULL
    }

    # Set up parameter combinations for each individual chain
    runParams <- base::expand.grid(c(chain = list(seq_len(nchains)),
        paramsTest))
    runParams <- cbind(index = seq_len(nrow(runParams)), runParams)

    # Pre-generate a set of random seeds to be used for each chain
    allSeeds <- seq(seed, (seed + nchains - 1))

    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = NULL,
        append = FALSE,
        verbose = verbose)

    .logMessages("Starting celdaGridSearch with",
        model,
        logfile = NULL,
        append = TRUE,
        verbose = verbose)

    .logMessages("Number of cores:",
        cores,
        logfile = NULL,
        append = TRUE,
        verbose = verbose)

    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = NULL,
        append = TRUE,
        verbose = verbose)

    startTime <- Sys.time()

    # An MD5 checksum of the count matrix. Passed to models so
    # later on, we can check on celda_* model objects which
    # count matrix was used.
    counts <- .processCounts(counts)
    countChecksum <- .createCountChecksum(counts)

    ## Use DoParallel to loop through each combination of parameters
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    i <- NULL # Setting visible binding for R CMD CHECK
    resList <- foreach(i = seq_len(nrow(runParams)),
        .export = model,
        .combine = c,
        .multicombine = TRUE) %dopar% {

            ## Set up chain parameter list
            current.run <- c(runParams[i, ])
            chainParams <- list()
            for (j in names(paramsTest)) {
                chainParams[[j]] <- current.run[[j]]
            }
            chainParams$counts <- counts
            chainParams$seed <- allSeeds[ifelse(i %% nchains == 0,
                nchains, i %% nchains)]
            chainParams$maxIter <- maxIter
            chainParams$nchain <- 1
            chainParams$countChecksum <- countChecksum
            chainParams$verbose <- verbose
            chainParams$logfile <- paste0(logfilePrefix,
                "_",
                paste(paste(
                    colnames(runParams), runParams[i, ], sep = "-"
                ), collapse = "_"),
                "_Seed-",
                chainParams$seed,
                "_log.txt"
            )

            ## Run model
            res <- do.call(model, c(chainParams, paramsFixed))
            return(list(res))
        }
    parallel::stopCluster(cl)

    logliks <- vapply(resList, function(mod) {
        mod@finalLogLik
    }, double(1))
    runParams <- cbind(runParams, logLikelihood = logliks)

    celdaRes <- methods::new(
        "celdaList",
        runParams = runParams,
        resList = resList,
        countChecksum = countChecksum
    )

    if (isTRUE(bestOnly)) {
        celdaRes <- selectBestModel(celdaRes)
    }

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
    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = NULL,
        append = TRUE,
        verbose = verbose)
    .logMessages("Completed celdaGridSearch. Total time:",
        format(difftime(endTime, startTime)),
        logfile = NULL,
        append = TRUE,
        verbose = verbose)
    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = NULL,
        append = TRUE,
        verbose = verbose)

    return(celdaRes)
}


################################################################################
# Methods for manipulating celdaList objects                                  #
################################################################################
#' @title Subset celdaList object from celdaGridSearch
#' @description Select a subset of models from a `celdaList` object generated
#'  by `celdaGridSearch()` that match the criteria in the argument `params`.
#' @param celdaList celdaList Object of class `celdaList`. An object
#'  containing celda models returned from `celdaGridSearch`.
#' @param params List. List of parameters used to subset celdaList.
#' @return A new `celdaList` object containing all models matching the
#'  provided criteria in `params`. If only one item in the `celdaList` matches
#'  the given criteria, the matching model will be returned directly instead of
#'  a `celdaList` object.
#' @seealso `celdaGridSearch()` can run Celda with multiple parameters and
#'  chains in parallel. `selectBestModel()` can get the best model for each
#'  combination of parameters.
#' @examples
#' resK5L10 <- subsetCeldaList(celdaCGGridSearchRes, params = list(K = 5,
#'     L = 10))
#' @export
subsetCeldaList <- function(celdaList, params) {
    if (!methods::is(celdaList, "celdaList")) {
        stop("celdaList parameter was not of class celdaList.")
    }

    ## Check for bad parameter names
    if (!all(names(params) %in% colnames(celdaList@runParams))) {
        badParams <- setdiff(names(params), colnames(celdaList@runParams))
        stop("The following elements in 'params' are not columns in runParams",
            " (celdaList) ",
            paste(badParams, collapse = ","))
    }

    ## Subset 'runParams' based on items in 'params'
    new.runParams <- celdaList@runParams
    for (i in names(params)) {
        newRunParams <-
            subset(newRunParams, newRunParams[, i] %in% params[[i]])

        if (nrow(newRunParams) == 0) {
            stop("No runs matched the criteria given in 'params'. Check",
                " 'runParams(celdaList)' for complete list of parameters used to",
                " generate 'celdaList'.")
        }
    }

    ## Get index of selected models, subset celdaList, and return
    ix <- match(newRunParams$index, celdaList@runParams$index)
    if (length(ix) == 1) {
        return(celdaList@resList[[ix]])
    } else {
        celdaList@runParams <- as.data.frame(newRunParams)
        celdaList@resList <- celdaList@resList[ix]
        return(celdaList)
    }
}


#' @title Select best chain within each combination of parameters
#' @description Select the chain with the best log likelihood for each
#'  combination of tested parameters from a `celdaList` object gererated by
#'  `celdaGridSearch()`.
#' @param celdaList Object of class `celdaList`. An object containing celda
#'  models returned from `celdaGridSearch()`.
#' @return A new `celdaList` object containing one model with the best log
#'  likelihood for each set of parameters. If only one set of parameters is in
#'  the `celdaList`, the best model will be returned directly instead of a
#'  `celdaList` object.
#' @seealso `celdaGridSearch()` can run Celda with multiple parameters and
#'  chains in parallel. `subsetCeldaList()` can subset the `celdaList` object.
#' @examples
#' ## Returns same result as running celdaGridSearch with "bestOnly = TRUE"
#' cgsBest <- selectBestModel(celdaCGGridSearchRes)
#' @import data.table
#' @export
selectBestModel <- function(celdaList) {
    if (!methods::is(celdaList, "celdaList"))
        stop("celdaList parameter was not of class celdaList.")

    logLikelihood <- NULL
    group <- setdiff(colnames(celdaList@runParams),
        c("index", "chain", "logLikelihood", "mean_perplexity"))
    dt <- data.table::as.data.table(celdaList@runParams)
    newRunParams <- as.data.frame(dt[, .SD[which.max(logLikelihood)],
        by = group])
    newRunParams <- newRunParams[, colnames(celdaList@runParams)]

    ix <- match(newRunParams$index, celdaList@runParams$index)
    if (nrow(newRunParams) == 1) {
        return(celdaList@resList[[ix]])
    } else {
        celdaList@runParams <- as.data.frame(newRunParams)
        celdaList@resList <- celdaList@resList[ix]
        return(celdaList)
    }
}


#' @title Celda models
#' @description List of available Celda models with correpsonding descriptions.
#' @export
#' @return None
celda <- function() {
    message("celda_C: Clusters the columns of a count matrix containing",
        " single-cell data into K subpopulations.")
    message("celda_G: Clusters the rows of a count matrix containing",
        " single-cell data into L modules.")
    message("celda_CG: Clusters the rows and columns of a count matrix",
        " containing single-cell data into L modules and K subpopulations,",
        " respectively.")
    message("celdaGridSearch: Run Celda with different combinations of",
        " parameters and multiple chains in parallel.")
}
