#' @title Calculate the perplexity of a celda model
#' @description Perplexity is a statistical measure of how well a probability
#'  model can predict new data. Lower perplexity indicates a better model.
#' @param x Can be one of
#'  \itemize{
#'  \item A \linkS4class{SingleCellExperiment} object returned by
#'  \link{celda_C}, \link{celda_G} or \link{celda_CG}, with the matrix
#'  located in the \code{useAssay} assay slot.
#'  Rows represent features and columns represent cells.
#'  \item Integer counts matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  \code{celdaMod}.}
#' @param useAssay A string specifying which \link[SummarizedExperiment]{assay}
#'  slot to use if \code{x} is a \linkS4class{SingleCellExperiment} object.
#'  Default "counts".
#' @param celdaMod Celda model object. Only works if \code{x} is an integer
#'  counts matrix.
#' @param newCounts A new counts matrix used to calculate perplexity. If NULL,
#'  perplexity will be calculated for the matrix in \code{useAssay} slot in
#'  \code{x}. Default NULL.
#' @return Numeric. The perplexity for the provided \code{x} (and
#'  \code{celdaModel}).
#' @export
setGeneric("perplexity",
    function(x, celdaMod, ...) {standardGeneric("perplexity")})


#' @importFrom matrixStats logSumExp
#' @examples
#' data(sceCeldaCG)
#' perplexity <- perplexity(sceCeldaCG)
#' @rdname perplexity
#' @export
setMethod("perplexity", signature(x = "SingleCellExperiment"),
    function(x, useAssay = "counts", newCounts = NULL) {

        if (celdaModel(x) == "celda_C") {
            p <- .perplexityCelda_C(sce = x, useAssay = useAssay,
                newCounts = newCounts)
            return(p)
        } else if (celdaModel(x) == "celda_CG") {
            p <- .perplexityCelda_CG(sce = x, useAssay = useAssay,
                newCounts = newCounts)
            return(p)
        } else if (celdaModel(x) == "celda_G") {
            p <- .perplexityCelda_G(sce = x, useAssay = useAssay,
                newCounts = newCounts)
            return(p)
        } else {
            stop("S4Vectors::metadata(x)$celda_parameters$model must be",
                " one of 'celda_C', 'celda_G', or 'celda_CG'")
        }
    })


#' @importFrom matrixStats logSumExp
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' perplexity <- perplexity(celdaCGSim$counts, celdaCGMod)
#' @rdname perplexity
#' @export
setMethod("perplexity", signature(x = "matrix", celdaMod = "celda_CG"),
    function(x, celdaMod, newCounts = NULL) {
        if (!("celda_CG" %in% class(celdaMod))) {
            stop("The celdaMod provided was not of class celda_CG.")
        }

        counts <- .processCounts(x)
        compareCountMatrix(counts, celdaMod)

        if (is.null(newCounts)) {
            newCounts <- counts
        } else {
            newCounts <- .processCounts(newCounts)
        }
        if (nrow(newCounts) != nrow(counts)) {
            stop("newCounts should have the same number of rows as counts.")
        }

        factorized <- factorizeMatrix(
            x = counts,
            celdaMod = celdaMod,
            type = c("posterior", "counts")
        )

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
    }
)


#' @examples
#' data(celdaCSim, celdaCMod)
#' perplexity <- perplexity(celdaCSim$counts, celdaCMod)
#' @importFrom matrixStats logSumExp
#' @rdname perplexity
#' @export
setMethod(
    "perplexity", signature(x = "matrix", celdaMod = "celda_C"),
    function(x, celdaMod, newCounts = NULL) {
        if (!("celda_C" %in% class(celdaMod))) {
            stop("The celdaMod provided was not of class celda_C.")
        }

        counts <- .processCounts(x)
        compareCountMatrix(counts, celdaMod)

        if (is.null(newCounts)) {
            newCounts <- counts
        } else {
            newCounts <- .processCounts(newCounts)
        }

        if (nrow(newCounts) != nrow(counts)) {
            stop("newCounts should have the same number of rows as counts.")
        }

        factorized <- factorizeMatrix(
            x = counts,
            celdaMod = celdaMod,
            type = "posterior"
        )
        theta <- log(factorized$posterior$sample)
        phi <- log(factorized$posterior$module)
        s <- as.integer(sampleLabel(celdaMod))

        # inner.log.prob = (t(phi) %*% newCounts) + theta[, s]
        inner.log.prob <- eigenMatMultInt(phi, newCounts) + theta[, s]
        logPx <- sum(apply(inner.log.prob, 2, matrixStats::logSumExp))

        perplexity <- exp(- (logPx / sum(newCounts)))
        return(perplexity)
    }
)


#' @examples
#' data(celdaGSim, celdaGMod)
#' perplexity <- perplexity(celdaGSim$counts, celdaGMod)
#' @rdname perplexity
#' @export
setMethod(
    "perplexity", signature(x = "matrix", celdaMod = "celda_G"),
    function(x, celdaMod, newCounts = NULL) {
        counts <- .processCounts(x)
        # compareCountMatrix(counts, celdaMod)

        if (is.null(newCounts)) {
            newCounts <- counts
        } else {
            newCounts <- .processCounts(newCounts)
        }
        if (nrow(newCounts) != nrow(counts)) {
            stop("newCounts should have the same number of rows as counts.")
        }

        factorized <- factorizeMatrix(
            x = counts,
            celdaMod = celdaMod,
            type = c("posterior", "counts")
        )
        psi <- factorized$posterior$module
        phi <- factorized$posterior$cell
        eta <- factorized$posterior$geneDistribution
        nGByTS <- factorized$counts$geneDistribution

        etaProb <- log(eta) * nGByTS
        # gene.by.cell.prob = log(psi %*% phi)
        # logPx = sum(gene.by.cell.prob * newCounts) # + sum(etaProb)
        logPx <- .perplexityGLogPx(
            newCounts,
            phi,
            psi,
            celdaClusters(celdaMod)$y,
            params(celdaMod)$L
        ) # + sum(etaProb)
        perplexity <- exp(- (logPx / sum(newCounts)))
        return(perplexity)
    }
)
