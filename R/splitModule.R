#' @title Split celda feature module
#' @description Manually select a celda feature module to split into 2 or
#'  more modules. Useful for splitting up modules that show divergent
#'  expression of features in multiple cell clusters.
#' @param x A \linkS4class{SingleCellExperiment} object
#'  with the matrix located in the assay slot under \code{useAssay}.
#'  Rows represent features and columns represent cells.
#' @param useAssay A string specifying which \link{assay}
#'  slot to use for \code{x}. Default "counts".
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default "featureSubset".
#' @param module Integer. The module to be split.
#' @param n Integer. How many modules should \code{module} be split into.
#'  Default 2.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @return A updated \linkS4class{SingleCellExperiment} object with new
#'  feature modules stored in column \code{celda_feature_module} in
#'  \code{\link{rowData}(x)}.
#' @export
setGeneric("splitModule",
    function(x,
        useAssay = "counts",
        altExpName = "featureSubset",
        module,
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
        useAssay = "counts",
        altExpName = "featureSubset",
        module,
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
    ix <- SummarizedExperiment::rowData(x)$celda_feature_module == module

    if (sum(ix) > 1) {
        tempModel <- .celda_G(
            counts = counts[ix, , drop = FALSE],
            L = n,
            yInitialize = "random",
            splitOnIter = -1,
            splitOnLast = FALSE,
            nchains = 1,
            verbose = FALSE)

        splitY <- celdaClusters(tempModel)$y
        splitIx <- celdaClusters(tempModel)$y > 1
        splitY[splitIx] <- S4Vectors::metadata(x)$celda_parameters$L +
            splitY[splitIx] - 1
        splitY[!splitIx] <- module

        newY <- as.integer(
            SummarizedExperiment::rowData(x)$celda_feature_module)
        newY[ix] <- splitY
        newL <- max(newY)

        newLl <- .logLikelihoodcelda_G(
            counts = counts,
            y = newY,
            L = newL,
            beta = S4Vectors::metadata(x)$celda_parameters$beta,
            delta = S4Vectors::metadata(x)$celda_parameters$delta,
            gamma = S4Vectors::metadata(x)$celda_parameters$gamma)
        model <- methods::new(
            "celda_G",
            clusters = list(y = newY),
            params = list(
                L = newL,
                beta = S4Vectors::metadata(x)$celda_parameters$beta,
                delta = S4Vectors::metadata(x)$celda_parameters$delta,
                gamma = S4Vectors::metadata(x)$celda_parameters$gamma,
                countChecksum = .createCountChecksum(counts)
            ),
            names = list(row = rownames(x),
                column = colnames(x),
                sample = x@metadata$celda_parameters$sampleLevels),
            finalLogLik = newLl
        )
    } else {
        stop("Module ", module, "contains <= 1 feature. No additional",
            " splitting was able to be performed.")
    }
    return(model)
}
