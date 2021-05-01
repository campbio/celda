#' @title Get the conditional probabilities of cell in subpopulations from celda
#'  model
#' @description Calculate the conditional probability of each cell belonging to
#'  each subpopulation given all other cell cluster assignments and/or
#'  each feature belonging to each module given all other feature cluster
#'  assignments in a celda model.
#' @param sce A \linkS4class{SingleCellExperiment} object returned by
#'  \link{celda_C}, \link{celda_G}, or \link{celda_CG}, with the matrix
#'  located in the \code{useAssay} assay slot.
#'  Rows represent features and columns represent cells.
#' @param useAssay A string specifying which \link{assay}
#'  slot to use. Default "counts".
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default "featureSubset".
#' @param log Logical. If \code{FALSE}, then the normalized conditional
#'  probabilities will be returned. If \code{TRUE}, then the unnormalized log
#'  probabilities will be returned. Default \code{FALSE}.
#' @param ... Ignored. Placeholder to prevent check warning.
#' @examples
#' data(sceCeldaCG)
#' clusterProb <- clusterProbability(sceCeldaCG, log = TRUE)
#' @return A list containging a matrix for the conditional cell subpopulation
#'  cluster and/or feature module probabilities.
#' @export
setGeneric("clusterProbability",
    function(sce, ...) {
        standardGeneric("clusterProbability")
    })


#' @seealso `celda_C()` for clustering cells
#' @examples
#' data(sceCeldaC)
#' clusterProb <- clusterProbability(sceCeldaC)
#' @rdname clusterProbability
#' @export
setMethod("clusterProbability", signature(sce = "SingleCellExperiment"),
    function(sce,
        useAssay = "counts",
        altExpName = "featureSubset",
        log = FALSE) {

        model <- celdaModel(sce, altExpName = altExpName)
        altExp <- SingleCellExperiment::altExp(sce, altExpName)
        counts <- SummarizedExperiment::assay(altExp, i = useAssay)
        beta <- S4Vectors::metadata(altExp)$celda_parameters$beta

        if (model == "celda_C") {
            s <- as.integer(
                SummarizedExperiment::colData(altExp)$celda_sample_label)
            z <- SummarizedExperiment::colData(altExp)$celda_cell_cluster
            K <- S4Vectors::metadata(altExp)$celda_parameters$K
            alpha <- S4Vectors::metadata(altExp)$celda_parameters$alpha

            cp <- .clusterProbabilityCeldaC(
                counts = counts,
                z = z,
                s = s,
                K = K,
                alpha = alpha,
                beta = beta,
                log = log)
        } else if (model == "celda_CG") {
            s <- as.integer(
                SummarizedExperiment::colData(altExp)$celda_sample_label)
            z <- SummarizedExperiment::colData(altExp)$celda_cell_cluster
            K <- S4Vectors::metadata(altExp)$celda_parameters$K
            y <- SummarizedExperiment::rowData(altExp)$celda_feature_module
            L <- S4Vectors::metadata(altExp)$celda_parameters$L
            alpha <- S4Vectors::metadata(altExp)$celda_parameters$alpha
            delta <- S4Vectors::metadata(altExp)$celda_parameters$delta
            gamma <- S4Vectors::metadata(altExp)$celda_parameters$gamma

            cp <- .clusterProbabilityCeldaCG(
                counts = counts,
                s = s,
                z = z,
                y = y,
                K = K,
                L = L,
                alpha = alpha,
                delta = delta,
                beta = beta,
                gamma = gamma,
                log = log)
        } else if (model == "celda_G") {
            y <- SummarizedExperiment::rowData(altExp)$celda_feature_module
            L <- S4Vectors::metadata(altExp)$celda_parameters$L
            delta <- S4Vectors::metadata(altExp)$celda_parameters$delta
            gamma <- S4Vectors::metadata(altExp)$celda_parameters$gamma

            cp <- .clusterProbabilityCeldaG(
                counts = counts,
                y = y,
                L = L,
                delta = delta,
                beta = beta,
                gamma = gamma,
                log = log)
        } else {
            stop("S4Vectors::metadata(altExp(sce, altExpName))$",
                "celda_parameters$model must be",
                " one of 'celda_C', 'celda_G', or 'celda_CG'!")
        }
        return(cp)
    }
)


.clusterProbabilityCeldaC <- function(
    counts,
    z,
    s,
    K,
    alpha,
    beta,
    log) {

    p <- .cCDecomposeCounts(counts, s, z, K)

    nextZ <- .cCCalcGibbsProbZ(counts = counts,
        mCPByS = p$mCPByS,
        nGByCP = p$nGByCP,
        nByC = p$nByC,
        nCP = p$nCP,
        z = z,
        s = s,
        K = K,
        nG = p$nG,
        nM = p$nM,
        alpha = alpha,
        beta = beta,
        doSample = FALSE)
    zProb <- t(nextZ$probs)

    if (!isTRUE(log)) {
        zProb <- .normalizeLogProbs(zProb)
    }
    return(list(zProbability = zProb))
}


.clusterProbabilityCeldaCG <- function(
    counts,
    s,
    z,
    y,
    K,
    L,
    alpha,
    delta,
    beta,
    gamma,
    log) {

    p <- .cCGDecomposeCounts(counts, s, z, y, K, L)
    lgbeta <- lgamma(seq(0, max(p$nCP)) + beta)
    lggamma <- lgamma(seq(0, nrow(counts) + L) + gamma)
    lgdelta <- c(NA, lgamma((seq(nrow(counts) + L) * delta)))

    nextZ <- .cCCalcGibbsProbZ(
        counts = p$nTSByC,
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
        doSample = FALSE
    )
    zProb <- t(nextZ$probs)

    ## Gibbs sampling for each gene
    nextY <- .cGCalcGibbsProbY(
        counts = p$nGByCP,
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
        doSample = FALSE
    )

    yProb <- t(nextY$probs)

    if (!isTRUE(log)) {
        zProb <- .normalizeLogProbs(zProb)
        yProb <- .normalizeLogProbs(yProb)
    }

    return(list(zProbability = zProb, yProbability = yProb))
}


.clusterProbabilityCeldaG <- function(
    counts,
    y,
    L,
    delta,
    beta,
    gamma,
    log) {

    ## Calculate counts one time up front
    p <- .cGDecomposeCounts(counts = counts, y = y, L = L)
    lgbeta <- lgamma(seq(0, max(.colSums(
        counts,
        nrow(counts), ncol(counts)
    ))) + beta)
    lggamma <- lgamma(seq(0, nrow(counts) + L) + gamma)
    lgdelta <- c(NA, lgamma(seq(nrow(counts) + L) * delta))

    nextY <- .cGCalcGibbsProbY(
        counts = counts,
        nTSByC = p$nTSByC,
        nByTS = p$nByTS,
        nGByTS = p$nGByTS,
        nByG = p$nByG,
        y = y,
        nG = p$nG,
        L = L,
        lgbeta = lgbeta,
        lgdelta = lgdelta,
        lggamma = lggamma,
        delta = delta,
        doSample = FALSE
    )
    yProb <- t(nextY$probs)

    if (!isTRUE(log)) {
        yProb <- .normalizeLogProbs(yProb)
    }
    return(list(yProbability = yProb))
}
