#' @title Run Celda in parallel with multiple parameters
#' @description Run Celda with different combinations of parameters and
#'  multiple chains in parallel. The variable \link{availableModels} contains
#'  the potential models that can be utilized. Different parameters to be tested
#'  should be stored in a list and passed to the argument \code{paramsTest}.
#'  Fixed parameters to be used in all models, such as \code{sampleLabel}, can
#'  be passed as a list to the argument \code{paramsFixed}. When
#'  \code{verbose = TRUE}, output from each chain will be sent to a log file
#'  but not be displayed in stdout.
#' @param x A numeric \link{matrix} of counts or a
#'  \linkS4class{SingleCellExperiment}
#'  with the matrix located in the assay slot under \code{useAssay}.
#'  Rows represent features and columns represent cells.
#' @param useAssay A string specifying which \link[SummarizedExperiment]{assay}
#'  slot to use if \code{x} is a
#'  \link[SingleCellExperiment]{SingleCellExperiment} object. Default "counts".
#' @param model Celda model. Options available in \link{availableModels}.
#' @param paramsTest List. A list denoting the combinations of parameters to
#'  run in a celda model. For example,
#'  \code{list(K = seq(5, 10), L = seq(15, 20))}
#'  will run all combinations of K from 5 to 10 and L from 15 to 20 in model
#'  \link{celda_CG}.
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
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. Seed values
#'  \code{seq(seed, (seed + nchains - 1))} will be supplied to each chain in
#'  \code{nchains} If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @param perplexity Logical. Whether to calculate perplexity for each model.
#'  If FALSE, then perplexity can be calculated later with
#'  `resamplePerplexity()`. Default TRUE.
#' @param verbose Logical. Whether to print log messages during celda chain
#'  execution. Default TRUE.
#' @param logfilePrefix Character. Prefix for log files from worker threads
#'  and main process. Default "Celda".
#' @return A \linkS4class{SingleCellExperiment} object. Function
#'  parameter settings and celda model results are stored in the
#'  \link[S4Vectors]{metadata} \code{"celda_grid_search"} slot.
#' @seealso \link{celda_G} for feature clustering, \link{celda_C} for
#'  clustering of cells, and \link{celda_CG} for simultaneous clustering of
#'  features and cells. \link{subsetCeldaList} can subset the \link{celdaList}
#'  object. \link{selectBestModel} can get the best model for each combination
#'  of parameters.
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom methods is
#' @examples \dontrun{
#' data(celdaCGSim)
#' ## Run various combinations of parameters with 'celdaGridSearch'
#' celdaCGGridSearchRes <- celdaGridSearch(celdaCGSim$counts,
#'   model = "celda_CG",
#'   paramsTest = list(K = seq(4, 6), L = seq(9, 11)),
#'   paramsFixed = list(sampleLabel = celdaCGSim$sampleLabel),
#'   bestOnly = TRUE,
#'   nchains = 1,
#'   cores = 1)}
#' @export
setGeneric("celdaGridSearch", function(x, ...) {
    standardGeneric("celdaGridSearch")})


#' @rdname celdaGridSearch
#' @export
setMethod("celdaGridSearch",
    signature(x = "SingleCellExperiment"),
    function(x,
        useAssay = "counts",
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

        xClass <- "SingleCellExperiment"
        counts <- SummarizedExperiment::assay(x, i = useAssay)

        celdaList <- .celdaGridSearch(counts = counts,
            model = model,
            paramsTest = paramsTest,
            paramsFixed = paramsFixed,
            maxIter = maxIter,
            nchains = nchains,
            cores = cores,
            bestOnly = bestOnly,
            seed = seed,
            perplexity = perplexity,
            verbose = verbose,
            logfilePrefix = logfilePrefix)

        sce <- .createSCEceldaGridSearch(celdaList = celdaList,
            sce = x,
            xClass = xClass,
            useAssay = useAssay,
            model = model,
            paramsTest = paramsTest,
            paramsFixed = paramsFixed,
            maxIter = maxIter,
            seed = seed,
            nchains = nchains,
            cores = cores,
            bestOnly = bestOnly,
            perplexity = perplexity,
            verbose = verbose,
            logfilePrefix = logfilePrefix)
        return(sce)
    })


#' @rdname celdaGridSearch
#' @export
setMethod("celdaGridSearch",
    signature(x = "matrix"),
    function(x,
        model,
        paramsTest,
        paramsFixed = NULL,
        maxIter = 200,
        nchains = 3,
        cores = 1,
        bestOnly = TRUE,
        perplexity = TRUE,
        verbose = TRUE,
        logfilePrefix = "Celda") {

        xClass <- "matrix"
        useAssay <- NULL
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = x))
        celdaList <- .celdaGridSearch(counts = x,
            model = model,
            paramsTest = paramsTest,
            paramsFixed = paramsFixed,
            maxIter = maxIter,
            nchains = nchains,
            cores = cores,
            bestOnly = bestOnly,
            seed = seed,
            perplexity = perplexity,
            verbose = verbose,
            logfilePrefix = logfilePrefix)

        sce <- .createSCEceldaGridSearch(celdaList = celdaList,
            sce = sce,
            xClass = xClass,
            useAssay = useAssay,
            model = model,
            paramsTest = paramsTest,
            paramsFixed = paramsFixed,
            maxIter = maxIter,
            seed = seed,
            nchains = nchains,
            cores = cores,
            bestOnly = bestOnly,
            perplexity = perplexity,
            verbose = verbose,
            logfilePrefix = logfilePrefix)
        return(sce)
    })


.celdaGridSearch <- function(counts,
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
    stop(
      "The following elements in 'paramsTest' are not arguments of '",
      model,
      "': ",
      paste(badParams, collapse = ",")
    )
  }

  if (!is.null(paramsFixed) &&
    !all(names(paramsFixed) %in% names(modelParams))) {
    badParams <- setdiff(names(paramsFixed), names(modelParams))
    stop(
      "The following elements in 'paramsFixed' are not arguments",
      " of '",
      model,
      "': ",
      paste(badParams, collapse = ",")
    )
  }

  modelParamsRequired <- setdiff(
    names(modelParams[modelParams == ""]),
    "counts"
  )

  if (!all(modelParamsRequired %in% c(
    names(paramsTest),
    names(paramsFixed)
  ))) {
    missing.params <- setdiff(
      modelParamsRequired,
      c(names(paramsTest), names(paramsFixed))
    )
    stop(
      "The following arguments are not in 'paramsTest' or 'paramsFixed'",
      " but are required for '",
      model,
      "': ",
      paste(missing.params, collapse = ",")
    )
  }

  if (any(c("z.init", "y.init", "sampleLabel") %in% names(paramsTest))) {
    stop(
      "Setting parameters such as 'z.init', 'y.init', and 'sampleLabel'",
      " in 'paramsTest' is not currently supported."
    )
  }

  if (any(c("nchains") %in% names(paramsTest))) {
    warning(
      "Parameter 'nchains' should not be used within the paramsTest",
      " list"
    )
    paramsTest[["nchains"]] <- NULL
  }

  # Pre-generate a set of random seeds to be used for each chain
  if (is.null(seed)) {
    allSeeds <- NULL
  } else {
    allSeeds <- seq(seed, (seed + nchains - 1))
  }

  # Set up parameter combinations for each individual chain
  runParams <- base::expand.grid(c(
    chain = list(seq_len(nchains)),
    paramsTest
  ))
  runParams <- cbind(index = seq_len(nrow(runParams)), runParams)

  if (is.null(allSeeds)) {
    runParams <- cbind(runParams,
      seed = rep("NULL", nrow(runParams)))
  } else {
    runParams <- cbind(runParams,
      seed = rep(allSeeds, nrow(runParams) / nchains))
  }

  .logMessages(paste(rep("-", 50), collapse = ""),
    logfile = NULL,
    append = FALSE,
    verbose = verbose
  )

  .logMessages("Starting celdaGridSearch with",
    model,
    logfile = NULL,
    append = TRUE,
    verbose = verbose
  )

  .logMessages("Number of cores:",
    cores,
    logfile = NULL,
    append = TRUE,
    verbose = verbose
  )

  .logMessages(paste(rep("-", 50), collapse = ""),
    logfile = NULL,
    append = TRUE,
    verbose = verbose
  )

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
  resList <- foreach(
    i = seq_len(nrow(runParams)),
    .export = model,
    .combine = c,
    .multicombine = TRUE
  ) %dopar% {

    ## Set up chain parameter list
    current.run <- c(runParams[i, ])
    chainParams <- list()
    for (j in names(paramsTest)) {
      chainParams[[j]] <- current.run[[j]]
    }
    chainParams$counts <- counts
    # silently ignored if allSeeds is NULL!
    chainParams$seed <- allSeeds[ifelse(i %% nchains == 0,
      nchains, i %% nchains)]
    chainParams$maxIter <- maxIter
    chainParams$nchain <- 1
    chainParams$countChecksum <- countChecksum
    chainParams$verbose <- verbose
    chainParams$logfile <- paste0(
      logfilePrefix,
      "_",
      paste(paste(
        colnames(runParams), runParams[i, ],
        sep = "-"
      ), collapse = "_"),
      "_Seed-",
      ifelse(is.null(chainParams$seed), "NULL", chainParams$seed),
      "_log.txt"
    )

    ## Run model
    if (is.null(seed)) {
      res <- do.call(model, c(chainParams, paramsFixed, list(seed = NULL)))
    } else {
      res <- do.call(model, c(chainParams, paramsFixed))
    }
    return(list(res))
  }
  parallel::stopCluster(cl)

  logliks <- vapply(resList, function(mod) {
    bestLogLikelihood(mod)
  }, double(1))
  runParams <- cbind(runParams, logLikelihood = logliks)

  celdaRes <- methods::new(
    "celdaList",
    runParams = runParams,
    resList = resList,
    countChecksum = countChecksum
  )

  if (isTRUE(bestOnly)) {
    celdaRes <- selectBestModel(celdaRes, asList = TRUE)
  }

  if (isTRUE(perplexity)) {
    .logMessages(
      date(),
      ".. Calculating perplexity",
      append = TRUE,
      verbose = verbose,
      logfile = NULL
    )
    celdaRes <- resamplePerplexity(counts, celdaRes, seed = seed)
  }

  endTime <- Sys.time()
  .logMessages(paste(rep("-", 50), collapse = ""),
    logfile = NULL,
    append = TRUE,
    verbose = verbose
  )
  .logMessages("Completed celdaGridSearch. Total time:",
    format(difftime(endTime, startTime)),
    logfile = NULL,
    append = TRUE,
    verbose = verbose
  )
  .logMessages(paste(rep("-", 50), collapse = ""),
    logfile = NULL,
    append = TRUE,
    verbose = verbose
  )

  return(celdaRes)
}

#'
#' ################################################################################
#' # Methods for manipulating celdaList objects                                  #
#' ################################################################################
#' #' @title Subset celda model from SCE object returned from
#' #'  \code{celdaGridSearch}
#' #' @description Select a subset of models from a
#' #'  \linkS4class{SingleCellExperiment} object generated by
#' #'  \link{celdaGridSearch} that match the criteria in the argument
#' #'  \code{params}.
#' #' @param x A \linkS4class{SingleCellExperiment} object returned from
#' #'  \code{celdaGridSearch}. Must contain a list named
#' #'  \code{"celda_grid_search"} in \code{metadata(x)}.
#' #' @param params List. List of parameters used to subset the matching celda
#' #'  models in list \code{"celda_grid_search"} in \code{metadata(x)}.
#' #' @return A new \linkS4class{SingleCellExperiment} object containing
#' #'  all models matching the
#' #'  provided criteria in \code{params}. If only one celda model result in the
#' #'  \code{"celda_grid_search"} slot in \code{metadata(x)} matches
#' #'  the given criteria, a new \linkS4class{SingleCellExperiment} object
#' #'  with the matching model stored in the
#' #'  \link[S4Vectors]{metadata}
#' #'  \code{"celda_parameters"} slot will be returned. Otherwise, a new
#' #'  \linkS4class{SingleCellExperiment} object with the subset models stored
#' #'  in the \link[S4Vectors]{metadata}
#' #'  \code{"celda_grid_search"} slot will be returned.
#' #' @seealso \link{celdaGridSearch} can run Celda with multiple parameters and
#' #'  chains in parallel. \link{selectBestModel} can get the best model for each
#' #'  combination of parameters.
#' #' @examples
#' #' data(celdaCGGridSearchRes)
#' #' resK5L10 <- subsetCeldaList(celdaCGGridSearchRes,
#' #'     params = list(K = 5, L = 10))
#' #' @export
#' setGeneric("subsetCeldaList", function(x, ...) {
#'     standardGeneric("subsetCeldaList")})
#'
#'
#' #' @rdname subsetCeldaList
#' #' @export
#' setMethod("subsetCeldaList",
#'     signature(x = "SingleCellExperiment"),
#'     function(x, params) {
#'
#'     ## Check for bad parameter names
#'     if (!all(names(params) %in% colnames(runParams(celdaList)))) {
#'         badParams <- setdiff(names(params), colnames(runParams(celdaList)))
#'         stop(
#'             "The following elements in 'params' are not columns in runParams",
#'             " (celdaList) ",
#'             paste(badParams, collapse = ",")
#'         )
#'     }
#'
#'     ## Subset 'runParams' based on items in 'params'
#'     newRunParams <- runParams(celdaList)
#'     for (i in names(params)) {
#'         newRunParams <-
#'             subset(newRunParams, newRunParams[, i] %in% params[[i]])
#'
#'         if (nrow(newRunParams) == 0) {
#'             stop(
#'                 "No runs matched the criteria given in 'params'. Check",
#'                 " 'runParams(celdaList)' for complete list of parameters used",
#'                 " to generate 'celdaList'."
#'             )
#'         }
#'     }
#'
#'     ## Get index of selected models, subset celdaList, and return
#'     ix <- match(newRunParams$index, runParams(celdaList)$index)
#'     if (length(ix) == 1) {
#'         return(resList(celdaList)[[ix]])
#'     } else {
#'         celdaList@runParams <- as.data.frame(newRunParams)
#'         celdaList@resList <- resList(celdaList)[ix]
#'         return(celdaList)
#'     }
#' }


################################################################################
# Methods for manipulating celdaList objects                                  #
################################################################################
#' @title Subset celdaList object from celdaGridSearch
#' @description Select a subset of models from a `celdaList` object generated
#'  by `celdaGridSearch()` that match the criteria in the argument `params`.
#' @param x celdaList Object of class `celdaList`. An object
#'  containing celda models returned from `celdaGridSearch` in older versions.
#' @param params List. List of parameters used to subset celdaList.
#' @return A new `celdaList` object containing all models matching the
#'  provided criteria in `params`. If only one item in the `celdaList` matches
#'  the given criteria, the matching model will be returned directly instead of
#'  a `celdaList` object.
#' @seealso `celdaGridSearch()` can run Celda with multiple parameters and
#'  chains in parallel. `selectBestModel()` can get the best model for each
#'  combination of parameters.
#' @examples
#' data(celdaCGGridSearchRes)
#' resK5L10 <- .subsetCeldaList(celdaCGGridSearchRes,
#'   params = list(K = 5, L = 10))
subsetCeldaList <- function(celdaList, params) {
  if (!methods::is(celdaList, "celdaList")) {
    stop("celdaList parameter was not of class celdaList.")
  }

  ## Check for bad parameter names
  if (!all(names(params) %in% colnames(runParams(celdaList)))) {
    badParams <- setdiff(names(params), colnames(runParams(celdaList)))
    stop(
      "The following elements in 'params' are not columns in runParams",
      " (celdaList) ",
      paste(badParams, collapse = ",")
    )
  }

  ## Subset 'runParams' based on items in 'params'
  newRunParams <- runParams(celdaList)
  for (i in names(params)) {
    newRunParams <-
      subset(newRunParams, newRunParams[, i] %in% params[[i]])

    if (nrow(newRunParams) == 0) {
      stop(
        "No runs matched the criteria given in 'params'. Check",
        " 'runParams(celdaList)' for complete list of parameters used",
        " to generate 'celdaList'."
      )
    }
  }

  ## Get index of selected models, subset celdaList, and return
  ix <- match(newRunParams$index, runParams(celdaList)$index)
  if (length(ix) == 1) {
    return(resList(celdaList)[[ix]])
  } else {
    celdaList@runParams <- as.data.frame(newRunParams)
    celdaList@resList <- resList(celdaList)[ix]
    return(celdaList)
  }
}


#' @title Select best chain within each combination of parameters
#' @description Select the chain with the best log likelihood for each
#'  combination of tested parameters from a `celdaList` object gererated by
#'  `celdaGridSearch()`.
#' @param celdaList Object of class `celdaList`. An object containing celda
#'  models returned from `celdaGridSearch()`.
#' @param asList `TRUE` or `FALSE`. Whether to return the best model as a
#'  `celdaList` object or not. If `FALSE`, return the best model as a
#'  corresponding `celda_C`, `celda_G` or `celda_CG` object.
#' @return A new `celdaList` object containing one model with the best log
#'  likelihood for each set of parameters. If only one set of parameters is in
#'  the `celdaList`, the best model will be returned directly instead of a
#'  `celdaList` object.
#' @seealso `celdaGridSearch()` can run Celda with multiple parameters and
#'  chains in parallel. `subsetCeldaList()` can subset the `celdaList` object.
#' @examples
#' data(celdaCGGridSearchRes)
#' ## Returns same result as running celdaGridSearch with "bestOnly = TRUE"
#' cgsBest <- selectBestModel(celdaCGGridSearchRes)
#' @importFrom data.table as.data.table
#' @export
selectBestModel <- function(celdaList, asList = FALSE) {
  if (!methods::is(celdaList, "celdaList")) {
    stop("celdaList parameter was not of class celdaList.")
  }

  logLikelihood <- NULL
  group <- setdiff(
    colnames(runParams(celdaList)),
    c("index", "chain", "logLikelihood", "mean_perplexity", "seed")
  )
  dt <- data.table::as.data.table(runParams(celdaList))
  newRunParams <- as.data.frame(dt[, .SD[which.max(logLikelihood)],
    by = group
  ])
  newRunParams <- newRunParams[, colnames(runParams(celdaList))]

  ix <- match(newRunParams$index, runParams(celdaList)$index)
  if (nrow(newRunParams) == 1 & !asList) {
    return(resList(celdaList)[[ix]])
  } else {
    celdaList@runParams <- as.data.frame(newRunParams)
    celdaList@resList <- resList(celdaList)[ix]
    return(celdaList)
  }
}


#' @title Celda models
#' @description List of available Celda models with correpsonding descriptions.
#' @export
#' @examples
#' celda()
#' @return None
celda <- function() {
  message(
    "celda_C: Clusters the columns of a count matrix containing",
    " single-cell data into K subpopulations."
  )
  message(
    "celda_G: Clusters the rows of a count matrix containing",
    " single-cell data into L modules."
  )
  message(
    "celda_CG: Clusters the rows and columns of a count matrix",
    " containing single-cell data into L modules and K subpopulations,",
    " respectively."
  )
  message(
    "celdaGridSearch: Run Celda with different combinations of",
    " parameters and multiple chains in parallel."
  )
}


.createSCEceldaGridSearch <- function(celdaList,
    sce,
    xClass,
    useAssay,
    model,
    paramsTest,
    paramsFixed,
    maxIter,
    nchains,
    cores,
    bestOnly,
    perplexity,
    verbose,
    logfilePrefix) {

    S4Vectors::metadata(sce)[["celda_grid_search"]] <- celdaList

    S4Vectors::metadata(sce)$celda_grid_search@celdaGridSearchParameters <-
        list(xClass = xClass,
            useAssay = useAssay,
            model = model,
            paramsTest = paramsTest,
            paramsFixed = paramsFixed,
            maxIter = maxIter,
            seed = seed,
            nchains = nchains,
            cores = cores,
            bestOnly = bestOnly,
            perplexity = perplexity,
            verbose = verbose,
            logfilePrefix = logfilePrefix)
    return(sce)
}
