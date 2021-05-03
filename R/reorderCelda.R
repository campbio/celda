#' @title Reorder cells populations and/or features modules using
#'  hierarchical clustering
#' @description Apply hierarchical clustering to reorder the cell populations
#'  and/or feature modules and group similar ones together based
#'  on the cosine distance of the factorized matrix
#'  from \link{factorizeMatrix}.
#' @param x Can be one of
#'  \itemize{
#'  \item A \linkS4class{SingleCellExperiment} object returned by
#'  \link{celda_C}, \link{celda_G} or \link{celda_CG}, with the matrix
#'  located in the \code{useAssay} assay slot in \code{altExp(x, altExpName)}.
#'  Rows represent features and columns represent cells.
#'  \item Integer count matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  \code{celdaMod}.}
#' @param useAssay A string specifying which \link{assay}
#'  slot to use if \code{x} is a \linkS4class{SingleCellExperiment} object.
#'  Default "counts".
#' @param altExpName The name for the \link{altExp} slot.
#'  Default "featureSubset".
#' @param method Passed to \link{hclust}. The agglomeration method
#'  to be used to be used. Default "complete".
#' @param celdaMod Celda model object. Only works if \code{x} is an integer
#'  counts matrix. Ignored if \code{x} is a
#'  \linkS4class{SingleCellExperiment} object.
#' @return A \linkS4class{SingleCellExperiment} object (or Celda model object)
#'  with updated cell cluster and/or feature module labels.
#' @export
setGeneric("reorderCelda",
    function(x,
        celdaMod,
        useAssay = "counts",
        altExpName = "featureSubset",
        method = "complete") {

        standardGeneric("reorderCelda")})


#' @examples
#' data(sceCeldaCG)
#' reordersce <- reorderCelda(sceCeldaCG)
#' @rdname reorderCelda
#' @export
setMethod("reorderCelda", signature(x = "SingleCellExperiment"),
    function(x,
        useAssay = "counts",
        altExpName = "featureSubset",
        method = "complete") {

        altExp <- SingleCellExperiment::altExp(x, e = altExpName)

        if (!useAssay %in% SummarizedExperiment::assayNames(altExp)) {
            stop(useAssay, " not in assayNames(altExp(x, altExpName))")
        }

        if (celdaModel(x, altExpName = altExpName) == "celda_C") {
            sce <- .reorderCeldaCsce(sce = x,
                useAssay = useAssay,
                altExpName = altExpName,
                method = method)
        } else if (celdaModel(x, altExpName = altExpName) == "celda_CG") {
            sce <- .reorderCeldaCGsce(sce = x,
                useAssay = useAssay,
                altExpName = altExpName,
                method = method)
        } else if (celdaModel(x, altExpName = altExpName) == "celda_G") {
            sce <- .reorderCeldaGsce(sce = x,
                useAssay = useAssay,
                altExpName = altExpName,
                method = method)
        } else {
            stop("S4Vectors::metadata(altExp(x, altExpName))$",
                "celda_parameters$model must be",
                " one of 'celda_C', 'celda_G', or 'celda_CG'")
        }
        return(sce)
    })


#' @examples
#' data(celdaCGSim, celdaCGMod)
#' reorderCeldaCG <- reorderCelda(celdaCGSim$counts, celdaCGMod)
#' @rdname reorderCelda
#' @export
setMethod("reorderCelda", signature(x = "matrix", celdaMod = "celda_CG"),
    function(x, celdaMod, method = "complete") {
        res <- .reorderCeldaCG(x, celdaMod, method = method)
        return(res)
    })


#' @examples
#' data(celdaCSim, celdaCMod)
#' reorderCeldaC <- reorderCelda(celdaCSim$counts, celdaCMod)
#' @rdname reorderCelda
#' @export
setMethod("reorderCelda", signature(x = "matrix", celdaMod = "celda_C"),
    function(x, celdaMod, method = "complete") {
        res <- .reorderCeldaC(x, celdaMod, method = method)
        return(res)
    })


#' @examples
#' data(celdaGSim, celdaGMod)
#' reorderCeldaG <- reorderCelda(celdaGSim$counts, celdaGMod)
#' @rdname reorderCelda
#' @export
setMethod("reorderCelda", signature(x = "matrix", celdaMod = "celda_G"),
    function(x, celdaMod, method = "complete") {
        res <- .reorderCeldaG(x, celdaMod, method = method)
        return(res)
    })


.reorderCeldaC <- function(counts, res, method = "complete") {
    if (params(res)$K > 2 & isTRUE(length(unique(celdaClusters(res)$z)) > 1)) {
        res@clusters$z <- as.integer(as.factor(celdaClusters(res)$z))
        fm <- factorizeMatrix(counts, res, type = "posterior")
        uniqueZ <- sort(unique(celdaClusters(res)$z))
        d <- .cosineDist(fm$posterior$module[, uniqueZ])
        h <- stats::hclust(d, method = method)
        res <- .recodeClusterZ(res,
            from = h$order,
            to = seq(length(h$order))
        )
    }
    return(res)
}


.reorderCeldaCsce <- function(sce, useAssay, altExpName, method = "complete") {
    if (S4Vectors::metadata(SingleCellExperiment::altExp(sce,
        altExpName))$celda_parameters$K > 2 &
            isTRUE(length(unique(celdaClusters(sce, altExpName))) > 1)) {
        celdaClusters(sce, altExpName) <-
            as.integer(as.factor(celdaClusters(sce, altExpName)))
        fm <- factorizeMatrix(sce,
            useAssay = useAssay,
            altExpName = altExpName,
            type = "posterior")
        uniqueZ <- sort(unique(celdaClusters(sce)))
        d <- .cosineDist(fm$posterior$module[, uniqueZ])
        h <- stats::hclust(d, method = method)
        sce <- recodeClusterZ(sce,
            from = h$order,
            to = seq(length(h$order)),
            altExpName = altExpName)
    }
    return(sce)
}


.reorderCeldaG <- function(counts, res, method = "complete") {
    if (params(res)$L > 2 & isTRUE(length(unique(celdaClusters(res)$y)) > 1)) {
        res@clusters$y <- as.integer(as.factor(celdaClusters(res)$y))
        fm <- factorizeMatrix(counts, res, type = "posterior")
        uniqueY <- sort(unique(celdaClusters(res)$y))
        cs <- prop.table(t(fm$posterior$cell[uniqueY, ]), 2)
        d <- .cosineDist(cs)
        h <- stats::hclust(d, method = method)
        res <- .recodeClusterY(res, from = h$order, to = seq(length(h$order)))
    }
    return(res)
}


.reorderCeldaGsce <- function(sce, useAssay, altExpName, method = "complete") {
    if (S4Vectors::metadata(SingleCellExperiment::altExp(sce,
        altExpName))$celda_parameters$L > 2 &
            isTRUE(length(unique(celdaModules(sce, altExpName))) > 1)) {
        celdaModules(sce, altExpName) <-
            as.integer(as.factor(celdaModules(sce, altExpName)))
        fm <- factorizeMatrix(sce,
            useAssay = useAssay,
            altExpName = altExpName,
            type = "posterior")
        uniqueY <- sort(unique(celdaModules(sce)))
        cs <- prop.table(t(fm$posterior$cell[uniqueY, ]), 2)
        d <- .cosineDist(cs)
        h <- stats::hclust(d, method = method)
        sce <- recodeClusterY(sce,
            from = h$order,
            to = seq(length(h$order)),
            altExpName = altExpName)
    }
    return(sce)
}


.reorderCeldaCG <- function(counts, res, method = "complete") {
    # Reorder K
    if (params(res)$K > 2 & isTRUE(length(unique(celdaClusters(res)$z)) > 1)) {
        res@clusters$z <- as.integer(as.factor(celdaClusters(res)$z))
        fm <- factorizeMatrix(counts, res, type = "posterior")
        uniqueZ <- sort(unique(celdaClusters(res)$z))
        d <- .cosineDist(fm$posterior$cellPopulation[, uniqueZ])
        h <- stats::hclust(d, method = method)

        res <- .recodeClusterZ(res, from = h$order, to = seq(length(h$order)))
    }

    # Reorder L
    if (params(res)$L > 2 & isTRUE(length(unique(celdaClusters(res)$y)) > 1)) {
        res@clusters$y <- as.integer(as.factor(celdaClusters(res)$y))
        fm <- factorizeMatrix(counts, res, type = "posterior")
        uniqueY <- sort(unique(celdaClusters(res)$y))
        cs <- prop.table(t(fm$posterior$cellPopulation[uniqueY, ]), 2)
        d <- .cosineDist(cs)
        h <- stats::hclust(d, method = method)

        res <- .recodeClusterY(res, from = h$order, to = seq(length(h$order)))
    }
    return(res)
}


.reorderCeldaCGsce <- function(sce, useAssay, altExpName, method = "complete") {
    # Reorder K
    if (S4Vectors::metadata(SingleCellExperiment::altExp(sce,
        altExpName))$celda_parameters$K > 2 &
            isTRUE(length(unique(celdaClusters(sce, altExpName))) > 1)) {
        celdaClusters(sce, altExpName) <-
            as.integer(as.factor(celdaClusters(sce, altExpName)))
        fm <- factorizeMatrix(sce,
            useAssay = useAssay,
            altExpName = altExpName,
            type = "posterior")
        uniqueZ <- sort(unique(celdaClusters(sce)))
        d <- .cosineDist(fm$posterior$cellPopulation[, uniqueZ])
        h <- stats::hclust(d, method = method)
        sce <- recodeClusterZ(sce,
            from = h$order,
            to = seq(length(h$order)),
            altExpName = altExpName)
    }

    # Reorder L
    if (S4Vectors::metadata(SingleCellExperiment::altExp(sce,
        altExpName))$celda_parameters$L > 2 &
            isTRUE(length(unique(celdaModules(sce, altExpName))) > 1)) {
        celdaModules(sce, altExpName) <-
            as.integer(as.factor(celdaModules(sce, altExpName)))
        fm <- factorizeMatrix(sce,
            useAssay = useAssay,
            altExpName = altExpName,
            type = "posterior")
        uniqueY <- sort(unique(celdaModules(sce)))
        cs <- prop.table(t(fm$posterior$cellPopulation[uniqueY, ]), 2)
        d <- .cosineDist(cs)
        h <- stats::hclust(d, method = method)
        sce <- recodeClusterY(sce,
            from = h$order,
            to = seq(length(h$order)),
            altExpName = altExpName)
    }
    return(sce)
}
