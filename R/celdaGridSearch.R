#' @title Run Celda in parallel with multiple parameters
#' @description Run Celda with different combinations of parameters and
#'  multiple chains in parallel. The variable \link{availableModels} contains
#'  the potential models that can be utilized. Different parameters to be tested
#'  should be stored in a list and passed to the argument \code{paramsTest}.
#'  Fixed parameters to be used in all models, such as \code{sampleLabel}, can
#'  be passed as a list to the argument \code{paramsFixed}. When
#'  \code{verbose = TRUE}, output from each chain will be sent to a log file
#'  but not be displayed in \code{stdout}.
#' @param x A numeric \link{matrix} of counts or a
#'  \linkS4class{SingleCellExperiment}
#'  with the matrix located in the assay slot under \code{useAssay}.
#'  Rows represent features and columns represent cells.
#' @param useAssay A string specifying the name of the
#'  \link[SummarizedExperiment]{assay} slot to use. Default "counts".
#' @param altExpName The name for the \link[SingleCellExperiment]{altExp} slot
#'  to use. Default "featureSubset".
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
#'  \code{nchains}. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @param perplexity Logical. Whether to calculate perplexity for each model.
#'  If FALSE, then perplexity can be calculated later with
#'  \link{resamplePerplexity}. Default TRUE.
#' @param verbose Logical. Whether to print log messages during celda chain
#'  execution. Default TRUE.
#' @param logfilePrefix Character. Prefix for log files from worker threads
#'  and main process. Default "Celda".
#' @return A \linkS4class{SingleCellExperiment} object. Function
#'  parameter settings and celda model results are stored in the
#'  \link[S4Vectors]{metadata} \code{"celda_grid_search"} slot.
#' @seealso \link{celda_G} for feature clustering, \link{celda_C} for
#'  clustering of cells, and \link{celda_CG} for simultaneous clustering of
#'  features and cells. \link{subsetCeldaList} can subset the \code{celdaList}
#'  object. \link{selectBestModel} can get the best model for each combination
#'  of parameters.
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom methods is
#' @examples
#' \dontrun{
#' data(celdaCGSim)
#' ## Run various combinations of parameters with 'celdaGridSearch'
#' celdaCGGridSearchRes <- celdaGridSearch(celdaCGSim$counts,
#'   model = "celda_CG",
#'   paramsTest = list(K = seq(4, 6), L = seq(9, 11)),
#'   paramsFixed = list(sampleLabel = celdaCGSim$sampleLabel),
#'   bestOnly = TRUE,
#'   nchains = 1,
#'   cores = 1)
#' }
#' @export
setGeneric("celdaGridSearch", function(x, ...) {
    standardGeneric("celdaGridSearch")})


#' @rdname celdaGridSearch
#' @export
setMethod("celdaGridSearch",
    signature(x = "SingleCellExperiment"),
    function(x,
        useAssay = "counts",
        altExpName = "featureSubset",
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

        if (!altExpName %in% SingleCellExperiment::altExpNames(x)) {
            stop(altExpName, " not in 'altExpNames(x)'. Run ",
                "selectFeatures(x) first!")
        }

        altExp <- SingleCellExperiment::altExp(x, altExpName)

        if (!useAssay %in% SummarizedExperiment::assayNames(altExp)) {
            stop(useAssay, " not in assayNames(altExp(x, altExpName))")
        }

        counts <- SummarizedExperiment::assay(altExp, i = useAssay)

        celdaList <- .celdaGridSearch(counts = counts,
            model = paste0(".", model),
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

        altExp <- .createSCEceldaGridSearch(celdaList = celdaList,
            sce = altExp,
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
        SingleCellExperiment::altExp(x, altExpName) <- altExp
        return(x)
    })


#' @rdname celdaGridSearch
#' @export
setMethod("celdaGridSearch",
    signature(x = "matrix"),
    function(x,
        useAssay = "counts",
        altExpName = "featureSubset",
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

        ls <- list()
        ls[[useAssay]] <- x
        sce <- SingleCellExperiment::SingleCellExperiment(assays = ls)
        SingleCellExperiment::altExp(sce, altExpName) <- sce
        xClass <- "matrix"

        celdaList <- .celdaGridSearch(counts = x,
            model = paste0(".", model),
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

        altExp <- .createSCEceldaGridSearch(celdaList = celdaList,
            sce = SingleCellExperiment::altExp(sce, altExpName),
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
        SingleCellExperiment::altExp(sce, altExpName) <- altExp
        return(sce)
    })


.celdaGridSearch <- function(counts,
                            model,
                            paramsTest,
                            paramsFixed,
                            maxIter,
                            nchains,
                            cores,
                            bestOnly,
                            seed,
                            perplexity,
                            verbose,
                            logfilePrefix) {

  ## Check parameters
  .validateCounts(counts)

  modelParams <- as.list(formals(model))
  if (!all(names(paramsTest) %in% names(modelParams))) {
    badParams <- setdiff(names(paramsTest), names(modelParams))
    stop(
      "The following elements in 'paramsTest' are not arguments of '",
      substring(model, 2),
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
      substring(model, 2),
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
      substring(model, 2),
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
    substring(model, 2),
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
      res <- do.call(model, c(chainParams, paramsFixed))
    } else {
      chainSeed <- allSeeds[ifelse(i %% nchains == 0,
          nchains, i %% nchains)]
      res <- with_seed(chainSeed,
          do.call(model, c(chainParams, paramsFixed)))
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


#' @title Subset celda model from SCE object returned from
#'  \code{celdaGridSearch}
#' @description Select a subset of models from a
#'  \linkS4class{SingleCellExperiment} object generated by
#'  \link{celdaGridSearch} that match the criteria in the argument
#'  \code{params}.
#' @param x Can be one of
#' \itemize{
#'  \item A \linkS4class{SingleCellExperiment} object returned from
#'  \code{celdaGridSearch}, \code{recursiveSplitModule},
#'  or \code{recursiveSplitCell}. Must contain a list named
#'  \code{"celda_grid_search"} in \code{metadata(x)}.
#'  \item celdaList object.}
#' @param params List. List of parameters used to subset the matching celda
#'  models in list \code{"celda_grid_search"} in \code{metadata(x)}.
#' @param useAssay A string specifying which \code{assay}
#'  slot to use if \code{x} is a
#'  \linkS4class{SingleCellExperiment} object. Default "counts".
#' @param altExpName The name for the \link[SingleCellExperiment]{altExp} slot
#'  to use. Default "featureSubset".
#' @return One of
#' \itemize{
#'  \item A new \linkS4class{SingleCellExperiment} object containing
#'  all models matching the
#'  provided criteria in \code{params}. If only one celda model result in the
#'  \code{"celda_grid_search"} slot in \code{metadata(x)} matches
#'  the given criteria, a new \linkS4class{SingleCellExperiment} object
#'  with the matching model stored in the
#'  \link[S4Vectors]{metadata}
#'  \code{"celda_parameters"} slot will be returned. Otherwise, a new
#'  \linkS4class{SingleCellExperiment} object with the subset models stored
#'  in the \link[S4Vectors]{metadata}
#'  \code{"celda_grid_search"} slot will be returned.
#'  \item A new \code{celdaList} object containing all models matching the
#'  provided criteria in \code{params}. If only one item in the
#'  \code{celdaList} matches the given criteria, the matching model will be
#'  returned directly instead of a \code{celdaList} object.}
#' @seealso \link{celdaGridSearch} can run Celda with multiple parameters and
#'  chains in parallel. \link{selectBestModel} can get the best model for each
#'  combination of parameters.
#' @export
setGeneric("subsetCeldaList", function(x, ...) {
    standardGeneric("subsetCeldaList")})


#' @rdname subsetCeldaList
#' @examples
#' data(sceCeldaCGGridSearch)
#' sceK5L10 <- subsetCeldaList(sceCeldaCGGridSearch,
#'     params = list(K = 5, L = 10))
#' @export
setMethod("subsetCeldaList",
    signature(x = "SingleCellExperiment"),
    function(x, params, useAssay = "counts", altExpName = "featureSubset") {

        ## Check for bad parameter names
        if (!all(names(params) %in% colnames(runParams(x)))) {
            badParams <- setdiff(names(params), colnames(runParams(x)))
            stop("The following elements in 'params' are not columns in",
                " runParams(x) ",
                paste(badParams, collapse = ",")
            )
        }

        ## Subset 'runParams' based on items in 'params'
        newRunParams <- runParams(x)
        for (i in names(params)) {
            newRunParams <-
                subset(newRunParams, newRunParams[, i] %in% params[[i]])

            if (nrow(newRunParams) == 0) {
                stop("No runs matched the criteria given in 'params'. Check",
                    " 'runParams(x)' for complete list of parameters used",
                    " to generate 'x'.")
            }
        }

        ## Get index of selected models, subset celdaList, and return
        ix <- match(newRunParams$index, runParams(x)$index)
        altExp <- SingleCellExperiment::altExp(x, altExpName)

        if (length(ix) == 1) {
            altExp <- .subsetCeldaListSCE(altExp, ix)
        } else {
            altExp@metadata$celda_grid_search@runParams <-
                as.data.frame(newRunParams)
            altExp@metadata$celda_grid_search@resList <-
                altExp@metadata$celda_grid_search@resList[ix]
        }
        SingleCellExperiment::altExp(x, altExpName) <- altExp
        return(x)
    }
)


#' @rdname subsetCeldaList
#' @examples
#' data(celdaCGGridSearchRes)
#' resK5L10 <- subsetCeldaList(celdaCGGridSearchRes,
#'     params = list(K = 5, L = 10))
#' @export
setMethod("subsetCeldaList",
    signature(x = "celdaList"),
    function(x, params) {
        ## Check for bad parameter names
        if (!all(names(params) %in% colnames(runParams(x)))) {
            badParams <- setdiff(names(params), colnames(runParams(x)))
            stop("The following elements in 'params' are not columns in",
                " runParams (x) ",
                paste(badParams, collapse = ",")
            )
        }

        ## Subset 'runParams' based on items in 'params'
        newRunParams <- runParams(x)
        for (i in names(params)) {
            newRunParams <-
                subset(newRunParams, newRunParams[, i] %in% params[[i]])

            if (nrow(newRunParams) == 0) {
                stop("No runs matched the criteria given in 'params'. Check",
                    " 'runParams(x)' for complete list of parameters used",
                    " to generate 'x'.")
            }
        }

        ## Get index of selected models, subset celdaList, and return
        ix <- match(newRunParams$index, runParams(x)$index)
        if (length(ix) == 1) {
            return(resList(x)[[ix]])
        } else {
            x@runParams <- as.data.frame(newRunParams)
            x@resList <- resList(x)[ix]
            return(x)
        }
    }
)


#' @title Select best chain within each combination of parameters
#' @description Select the chain with the best log likelihood for each
#'  combination of tested parameters from a \code{SCE} object gererated by
#'  \link{celdaGridSearch} or from a \code{celdaList} object.
#' @param x Can be one of
#' \itemize{
#'  \item A \linkS4class{SingleCellExperiment} object returned from
#'  \code{celdaGridSearch}, \code{recursiveSplitModule},
#'  or \code{recursiveSplitCell}. Must contain a list named
#'  \code{"celda_grid_search"} in \code{metadata(x)}.
#'  \item celdaList object.}
#' @param asList \code{TRUE} or \code{FALSE}. Whether to return the
#'  best model as a
#'  \code{celdaList} object or not. If \code{FALSE}, return the best model as a
#'  corresponding celda model object.
#' @param useAssay A string specifying which \code{assay}
#'  slot to use if \code{x} is a
#'  \linkS4class{SingleCellExperiment} object. Default "counts".
#' @param altExpName The name for the \link[SingleCellExperiment]{altExp} slot
#'  to use. Default "featureSubset".
#' @return One of
#' \itemize{
#'  \item A new \linkS4class{SingleCellExperiment} object containing
#'  one model with the best log-likelihood for each set of parameters in
#'  \code{metadata(x)}. If there is only one set of parameters,
#'  a new \linkS4class{SingleCellExperiment} object
#'  with the matching model stored in the
#'  \link[S4Vectors]{metadata}
#'  \code{"celda_parameters"} slot will be returned. Otherwise, a new
#'  \linkS4class{SingleCellExperiment} object with the subset models stored
#'  in the \link[S4Vectors]{metadata}
#'  \code{"celda_grid_search"} slot will be returned.
#'  \item A new \code{celdaList} object containing one model with the best
#'  log-likelihood for each set of parameters. If only one set of parameters
#'  is in the \code{celdaList}, the best model will be returned directly
#'  instead of a \code{celdaList} object.}
#' @seealso \link{celdaGridSearch} \link{subsetCeldaList}
#' @export
setGeneric("selectBestModel", function(x, ...) {
    standardGeneric("selectBestModel")})


#' @rdname selectBestModel
#' @examples
#' data(sceCeldaCGGridSearch)
#' ## Returns same result as running celdaGridSearch with "bestOnly = TRUE"
#' sce <- selectBestModel(sceCeldaCGGridSearch)
#' @importFrom data.table as.data.table
#' @export
setMethod("selectBestModel", signature(x = "SingleCellExperiment"),
    function(x, asList = FALSE, useAssay = "counts",
        altExpName = "featureSubset") {

        altExp <- SingleCellExperiment::altExp(x, altExpName)
        logLikelihood <- NULL
        group <- setdiff(colnames(runParams(x)),
            c("index", "chain", "logLikelihood", "mean_perplexity", "seed"))
        runParams <- S4Vectors::metadata(altExp)$celda_grid_search@runParams
        dt <- data.table::as.data.table(runParams)
        newRunParams <- as.data.frame(dt[, .SD[which.max(logLikelihood)],
            by = group])
        newRunParams <- newRunParams[, colnames(runParams)]

        ix <- match(newRunParams$index, runParams$index)
        if (nrow(newRunParams) == 1 & !asList) {
            altExp <- .subsetCeldaListSCE(altExp, ix)
        } else {
            altExp@metadata$celda_grid_search@runParams <-
                as.data.frame(newRunParams)
            altExp@metadata$celda_grid_search@resList <-
                altExp@metadata$celda_grid_search@resList[ix]
        }
        SingleCellExperiment::altExp(x, altExpName) <- altExp
        return(x)
    }
)


#' @rdname selectBestModel
#' @examples
#' data(celdaCGGridSearchRes)
#' ## Returns same result as running celdaGridSearch with "bestOnly = TRUE"
#' cgsBest <- selectBestModel(celdaCGGridSearchRes)
#' @importFrom data.table as.data.table
#' @export
setMethod("selectBestModel", signature(x = "celdaList"),
    function(x, asList = FALSE) {
        logLikelihood <- NULL
        group <- setdiff(colnames(runParams(x)),
            c("index", "chain", "logLikelihood", "mean_perplexity", "seed"))
        dt <- data.table::as.data.table(runParams(x))
        newRunParams <- as.data.frame(dt[, .SD[which.max(logLikelihood)],
            by = group])
        newRunParams <- newRunParams[, colnames(runParams(x))]

        ix <- match(newRunParams$index, runParams(x)$index)
        if (nrow(newRunParams) == 1 & !asList) {
            return(resList(x)[[ix]])
        } else {
            x@runParams <- as.data.frame(newRunParams)
            x@resList <- resList(x)[ix]
            return(x)
        }
    }
)


.createSCEceldaGridSearch <- function(celdaList,
    sce,
    xClass,
    useAssay,
    model,
    paramsTest,
    paramsFixed,
    maxIter,
    seed,
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


.subsetCeldaListSCE <- function(x, ix) {
    cgsparam <- x@metadata$celda_grid_search@celdaGridSearchParameters
    if (cgsparam$model == "celda_C") {
        x <- .createSCEceldaC(celdaCMod =
                x@metadata$celda_grid_search@resList[[ix]],
            sce = x,
            xClass = cgsparam$xClass,
            useAssay = cgsparam$useAssay,
            algorithm = cgsparam$algorithm,
            stopIter = cgsparam$stopIter,
            maxIter = cgsparam$maxIter,
            splitOnIter = cgsparam$splitOnIter,
            splitOnLast = cgsparam$splitOnLast,
            nchains = cgsparam$nchains,
            zInitialize = cgsparam[["zInitialize"]],
            zInit = cgsparam[["zInit"]],
            logfile = cgsparam$logfile,
            verbose = cgsparam$verbose)
    } else if (cgsparam$model == "celda_G") {
        x <- .createSCEceldaG(celdaGMod =
                x@metadata$celda_grid_search@resList[[ix]],
            sce = x,
            xClass = cgsparam$xClass,
            useAssay = cgsparam$useAssay,
            stopIter = cgsparam$stopIter,
            maxIter = cgsparam$maxIter,
            splitOnIter = cgsparam$splitOnIter,
            splitOnLast = cgsparam$splitOnLast,
            nchains = cgsparam$nchains,
            yInitialize = cgsparam[["yInitialize"]],
            yInit = cgsparam[["yInit"]],
            logfile = cgsparam$logfile,
            verbose = cgsparam$verbose)
    } else if (cgsparam$model == "celda_CG") {
        x <- .createSCEceldaCG(celdaCGMod =
                x@metadata$celda_grid_search@resList[[ix]],
            sce = x,
            xClass = cgsparam$xClass,
            useAssay = cgsparam$useAssay,
            algorithm = cgsparam$algorithm,
            stopIter = cgsparam$stopIter,
            maxIter = cgsparam$maxIter,
            splitOnIter = cgsparam$splitOnIter,
            splitOnLast = cgsparam$splitOnLast,
            nchains = cgsparam$nchains,
            zInitialize = cgsparam[["zInitialize"]],
            yInitialize = cgsparam[["yInitialize"]],
            zInit = cgsparam[["zInit"]],
            yInit = cgsparam[["yInit"]],
            logfile = cgsparam$logfile,
            verbose = cgsparam$verbose)
    } else {
        stop("S4Vectors::metadata(altExp(x, altExpName))$celda_grid_search@",
            "celdaGridSearchParameters$model must be",
            " one of 'celda_C', 'celda_G', or 'celda_CG'")
    }
    return(x)
}
