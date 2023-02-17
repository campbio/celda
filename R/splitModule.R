#' @title Split celda feature module
#' @description Manually select a celda feature module to split into 2 or
#'  more modules. Useful for splitting up modules that show divergent
#'  expression of features in multiple cell clusters.
#' @param x A \linkS4class{SingleCellExperiment} object
#'  with the matrix located in the assay slot under \code{useAssay}.
#'  Rows represent features and columns represent cells.
#' @param module Integer. The module to be split.
#' @param useAssay A string specifying which \link{assay}
#'  slot to use for \code{x}. Default "counts".
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default \code{"featureSubset"}.
#' @param n Integer. How many modules should \code{module} be split into.
#'  Default \code{2}.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @return A updated \linkS4class{SingleCellExperiment} object with new
#'  feature modules stored in column \code{celda_feature_module} in
#'  \code{\link{rowData}(x)}.
#' @export
setGeneric("splitModule",
    function(x,
        module,             
        useAssay = "counts",
        altExpName = "featureSubset",
        n = 2,
        seed = 12345) {

        standardGeneric("splitModule")
    })


#' @rdname splitModule
#' @examples
#' data(sceCeldaCG)
#' # Split module 5 into 2 new modules.
#' sce <- splitModule(sceCeldaCG, module = 5)
#' @export
setMethod("splitModule", signature(x = "SingleCellExperiment"),
    function(x,
        module,             
        useAssay = "counts",
        altExpName = "featureSubset",
        n = 2,
        seed = 12345) {

        altExp <- SingleCellExperiment::altExp(x, altExpName)

        if (!module %in% celdaModules(x, altExpName = altExpName)) {
            stop("Module ", module, " is not found in celdaModules(x,",
                " altExpName = altExpName).",
                " Please specify a valid module.")
        }

        celdaGMod <- .splitModuleWithSeed(x = altExp,
            useAssay = useAssay,
            module = module,
            n = n,
            seed = seed)

        S4Vectors::metadata(altExp)[["celda_parameters"]]$L <-
            params(celdaGMod)$L
        S4Vectors::metadata(altExp)[["celda_parameters"]]$finalLogLik <-
            celdaGMod@finalLogLik
        S4Vectors::metadata(altExp)[["celda_parameters"]]$featureModuleLevels <-
            sort(unique(celdaClusters(celdaGMod)$y))
        SummarizedExperiment::rowData(altExp)["celda_feature_module"] <-
            as.factor(celdaClusters(celdaGMod)$y)
        SingleCellExperiment::altExp(x, altExpName) <- altExp
        return(x)
    }
)


.splitModuleWithSeed <- function(x,
    useAssay,
    module,
    n,
    seed) {

    if (is.null(seed)) {
        celdaGMod <- .splitModule(x, useAssay, module, n)
    } else {
        with_seed(seed, celdaGMod <- .splitModule(x, useAssay, module, n))
    }
    return(celdaGMod)
}


.splitModule <- function(x, useAssay, module, n) {
    counts <- SummarizedExperiment::assay(x, i = useAssay)
    counts <- .processCounts(counts)
    .validateCounts(counts)
    
    L <- S4Vectors::metadata(x)$celda_parameters$L
    y <- as.numeric(SummarizedExperiment::rowData(x)$celda_feature_module)
    ix <- y == module
    
    if (sum(ix) < n) {
      stop("Module ", module, " contains less than ", n, " features. ",
            "Module splitting was not performed.")
    }      
      
    tempModel <- .celda_G(
      counts = counts[ix, , drop = FALSE],
      L = n,
      yInitialize = "random",
      splitOnIter = -1,
      splitOnLast = FALSE,
      nchains = 1,
      verbose = FALSE
    )
    
    # Need to set some of the features to the original module number.
    # The remaining features need to have "L + something" as they represent
    # a new module. Note that there may be more than 1 new module. 
    splitY <-
      as.numeric(as.character(celdaClusters(tempModel)$y))
    splitIx <- splitY > 1
    splitY[splitIx] <- L + splitY[splitIx] - 1
    splitY[!splitIx] <- module
    
    # Set up new y and L
    newY <- y
    newY[ix] <- splitY
    newL <- max(newY)
    
    newLl <- .logLikelihoodcelda_G(
      counts = counts,
      y = newY,
      L = newL,
      beta = S4Vectors::metadata(x)$celda_parameters$beta,
      delta = S4Vectors::metadata(x)$celda_parameters$delta,
      gamma = S4Vectors::metadata(x)$celda_parameters$gamma
    )
    
    model <- methods::new(
      "celda_G",
      clusters = list(y = factor(newY, seq(newL))),
      params = list(
        L = newL,
        beta = S4Vectors::metadata(x)$celda_parameters$beta,
        delta = S4Vectors::metadata(x)$celda_parameters$delta,
        gamma = S4Vectors::metadata(x)$celda_parameters$gamma,
        countChecksum = .createCountChecksum(counts)
      ),
      names = list(
        row = rownames(x),
        column = colnames(x),
        sample = x@metadata$celda_parameters$sampleLevels
      ),
      finalLogLik = newLl
    )
    
    return(model)
}
