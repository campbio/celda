#' @title Generate factorized matrices showing each feature's influence on cell
#'  / gene clustering
#' @description Generates factorized matrices showing the contribution of each
#'  feature in each cell population or each cell population in each sample.
#' @param x Can be one of
#'  \itemize{
#'  \item A \linkS4class{SingleCellExperiment} object returned by
#'  \link{celda_C}, \link{celda_G} or \link{celda_CG}, with the matrix
#'  located in the \code{useAssay} assay slot in \code{altExp(x, altExpName)}.
#'  Rows represent features and columns represent cells.
#'  \item Integer counts matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  \code{celdaMod}.}
#' @param useAssay A string specifying which \link{assay}
#'  slot to use if \code{x} is a \linkS4class{SingleCellExperiment} object.
#'  Default "counts".
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default "featureSubset".
#' @param celdaMod Celda model object. Only works if \code{x} is an integer
#'  counts matrix.
#' @param type Character vector. A vector containing one or more of "counts",
#'  "proportion", or "posterior". "counts" returns the raw number of counts for
#'  each factorized matrix. "proportions" returns the normalized probabilities
#'  for each factorized matrix, which are calculated by dividing the raw counts
#'  in each factorized matrix by the total counts in each column. "posterior"
#'  returns the posterior estimates which include the addition of the Dirichlet
#'  concentration parameter (essentially as a pseudocount). Default
#'  \code{"counts"}.
#' @param ... Ignored. Placeholder to prevent check warning.
#' @export
setGeneric("factorizeMatrix",
    function(x, celdaMod, ...) {
        standardGeneric("factorizeMatrix")})


#' @examples
#' data(sceCeldaCG)
#' factorizedMatrices <- factorizeMatrix(sceCeldaCG, type = "posterior")
#' @rdname factorizeMatrix
#' @export
setMethod("factorizeMatrix", signature(x = "SingleCellExperiment"),
    function(x,
        useAssay = "counts",
        altExpName = "featureSubset",
        type = c("counts", "proportion", "posterior")) {

        altExp <- SingleCellExperiment::altExp(x, e = altExpName)
        counts <- SummarizedExperiment::assay(altExp, i = useAssay)
        counts <- .processCounts(counts)

        beta <- S4Vectors::metadata(altExp)$celda_parameters$beta
        rNames <- rownames(altExp)

        if (celdaModel(x, altExpName = altExpName) == "celda_C") {
            z <- SummarizedExperiment::colData(altExp)$celda_cell_cluster
            K <- S4Vectors::metadata(altExp)$celda_parameters$K
            alpha <- S4Vectors::metadata(altExp)$celda_parameters$alpha
            sampleLabel <-
                SummarizedExperiment::colData(altExp)$celda_sample_label
            sNames <- S4Vectors::metadata(altExp)$celda_parameters$sampleLevels

            res <- .factorizeMatrixC(
                counts = counts,
                z = z,
                K = K,
                alpha = alpha,
                beta = beta,
                sampleLabel = sampleLabel,
                rNames = rNames,
                sNames = sNames,
                type = type)
        } else if (celdaModel(x, altExpName = altExpName) == "celda_CG") {
            K <- S4Vectors::metadata(altExp)$celda_parameters$K
            z <- SummarizedExperiment::colData(altExp)$celda_cell_cluster
            y <- SummarizedExperiment::rowData(altExp)$celda_feature_module
            L <- S4Vectors::metadata(altExp)$celda_parameters$L
            alpha <- S4Vectors::metadata(altExp)$celda_parameters$alpha
            delta <- S4Vectors::metadata(altExp)$celda_parameters$delta
            gamma <- S4Vectors::metadata(altExp)$celda_parameters$gamma
            sampleLabel <-
                SummarizedExperiment::colData(altExp)$celda_sample_label
            cNames <- colnames(altExp)
            sNames <- S4Vectors::metadata(altExp)$celda_parameters$sampleLevels

            res <- .factorizeMatrixCG(
                counts = counts,
                K = K,
                z = z,
                y = y,
                L = L,
                alpha = alpha,
                beta = beta,
                delta = delta,
                gamma = gamma,
                sampleLabel = sampleLabel,
                cNames = cNames,
                rNames = rNames,
                sNames = sNames,
                type = type)
        } else if (celdaModel(x, altExpName = altExpName) == "celda_G") {
            y <- SummarizedExperiment::rowData(altExp)$celda_feature_module
            L <- S4Vectors::metadata(altExp)$celda_parameters$L
            delta <- S4Vectors::metadata(altExp)$celda_parameters$delta
            gamma <- S4Vectors::metadata(altExp)$celda_parameters$gamma
            cNames <- colnames(altExp)

            res <- .factorizeMatrixG(
                counts = counts,
                y = y,
                L = L,
                beta = beta,
                delta = delta,
                gamma = gamma,
                cNames = cNames,
                rNames = rNames,
                type = type)
        } else {
            stop("S4Vectors::metadata(altExp(x, altExpName))$",
                "celda_parameters$model must be",
                " one of 'celda_C', 'celda_G', or 'celda_CG'")
        }
        return(res)
    })


#' @return For celda_CG model, A list with elements for "counts", "proportions",
#'  or "posterior" probabilities. Each element will be a list containing
#'  factorized matrices for "module", "cellPopulation", and "sample".
#'  Additionally, the contribution of each module in each individual cell will
#'  be included in the "cell" element of "counts" and "proportions" elements.
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' factorizedMatrices <- factorizeMatrix(
#'   celdaCGSim$counts,
#'   celdaCGMod,
#'   "posterior")
#' @rdname factorizeMatrix
#' @export
setMethod("factorizeMatrix", signature(x = "ANY", celdaMod = "celda_CG"),
    function(x,
        celdaMod,
        type = c("counts", "proportion", "posterior")) {

        counts <- .processCounts(x)
        compareCountMatrix(counts, celdaMod)

        z <- celdaClusters(celdaMod)$z
        y <- celdaClusters(celdaMod)$y
        # Sometimes, fewer clusters get returned by celda_C/G
        # Taking the max(z)/max(y) rather than
        # the original K/L will prevent errors
        # K <- params(celdaMod)$K; L <- params(celdaMod)$L
        K <- max(z)
        L <- max(y)
        alpha <- params(celdaMod)$alpha
        beta <- params(celdaMod)$beta
        delta <- params(celdaMod)$delta
        gamma <- params(celdaMod)$gamma
        sampleLabel <- sampleLabel(celdaMod)

        cNames <- matrixNames(celdaMod)$column
        rNames <- matrixNames(celdaMod)$row
        sNames <- matrixNames(celdaMod)$sample

        res <- .factorizeMatrixCG(
            counts = counts,
            K = K,
            z = z,
            y = y,
            L = L,
            alpha = alpha,
            beta = beta,
            delta = delta,
            gamma = gamma,
            sampleLabel = sampleLabel,
            cNames = cNames,
            rNames = rNames,
            sNames = sNames,
            type = type)
        return(res)
    }
)


.factorizeMatrixCG <- function(counts,
    K,
    z,
    y,
    L,
    alpha,
    beta,
    delta,
    gamma,
    sampleLabel,
    cNames,
    rNames,
    sNames,
    type) {

    s <- as.integer(sampleLabel)

    ## Calculate counts one time up front
    p <- .cCGDecomposeCounts(counts, s, z, y, K, L)
    nS <- p$nS
    nG <- p$nG
    nM <- p$nM
    mCPByS <- p$mCPByS
    nTSByC <- p$nTSByC
    nTSByCP <- p$nTSByCP
    nByG <- p$nByG
    nByTS <- p$nByTS
    nGByTS <- p$nGByTS
    nGByTS[nGByTS == 0] <- 1

    GByTS <- matrix(0, nrow = length(y), ncol = L)
    GByTS[cbind(seq(nG), y)] <- p$nByG

    LNames <- paste0("L", seq(L))
    KNames <- paste0("K", seq(K))
    colnames(nTSByC) <- cNames
    rownames(nTSByC) <- LNames
    colnames(GByTS) <- LNames
    rownames(GByTS) <- rNames
    rownames(mCPByS) <- KNames
    colnames(mCPByS) <- sNames
    colnames(nTSByCP) <- KNames
    rownames(nTSByCP) <- LNames

    countsList <- c()
    propList <- c()
    postList <- c()
    res <- list()

    if (any("counts" %in% type)) {
        countsList <- list(
            sample = mCPByS,
            cellPopulation = nTSByCP,
            cell = nTSByC,
            module = GByTS,
            geneDistribution = nGByTS
        )
        res <- c(res, list(counts = countsList))
    }

    if (any("proportion" %in% type)) {
        ## Need to avoid normalizing cell/gene states with zero cells/genes
        uniqueZ <- sort(unique(z))
        tempNTSByCP <- nTSByCP
        tempNTSByCP[, uniqueZ] <- normalizeCounts(tempNTSByCP[, uniqueZ],
            normalize = "proportion"
        )

        uniqueY <- sort(unique(y))
        tempGByTS <- GByTS
        tempGByTS[, uniqueY] <- normalizeCounts(tempGByTS[, uniqueY],
            normalize = "proportion"
        )
        tempNGByTS <- nGByTS / sum(nGByTS)

        propList <- list(
            sample = normalizeCounts(mCPByS,
                normalize = "proportion"
            ),
            cellPopulation = tempNTSByCP,
            cell = normalizeCounts(nTSByC, normalize = "proportion"),
            module = tempGByTS,
            geneDistribution = tempNGByTS
        )
        res <- c(res, list(proportions = propList))
    }

    if (any("posterior" %in% type)) {
        gs <- GByTS
        gs[cbind(seq(nG), y)] <- gs[cbind(seq(nG), y)] + delta
        gs <- normalizeCounts(gs, normalize = "proportion")
        tempNGByTS <- (nGByTS + gamma) / sum(nGByTS + gamma)

        postList <- list(
            sample = normalizeCounts(mCPByS + alpha,
                normalize = "proportion"
            ),
            cellPopulation = normalizeCounts(nTSByCP + beta,
                normalize = "proportion"
            ),
            module = gs,
            geneDistribution = tempNGByTS
        )
        res <- c(res, posterior = list(postList))
    }
    return(res)
}


#' @examples
#' data(celdaCSim, celdaCMod)
#' factorizedMatrices <- factorizeMatrix(
#'   celdaCSim$counts,
#'   celdaCMod, "posterior"
#' )
#' @return For celda_C model, a list with elements for "counts", "proportions",
#'  or "posterior" probabilities. Each element will be a list containing
#'  factorized matrices for "module" and "sample".
#' @rdname factorizeMatrix
#' @export
setMethod("factorizeMatrix", signature(x = "ANY", celdaMod = "celda_C"),
    function(x,
        celdaMod,
        type = c("counts", "proportion", "posterior")) {

        counts <- .processCounts(x)
        compareCountMatrix(counts, celdaMod)

        z <- celdaClusters(celdaMod)$z
        # Sometimes, fewer clusters get returned by celda_C
        # Taking the max(z) rather than the
        # original K will prevent errors
        # K <- params(celdaMod)$K
        K <- max(z)
        alpha <- params(celdaMod)$alpha
        beta <- params(celdaMod)$beta
        sampleLabel <- sampleLabel(celdaMod)
        rNames <- matrixNames(celdaMod)$row
        sNames <- matrixNames(celdaMod)$sample

        res <- .factorizeMatrixC(
            counts = counts,
            z = z,
            K = K,
            alpha = alpha,
            beta = beta,
            sampleLabel = sampleLabel,
            rNames = rNames,
            sNames = sNames,
            type = type)
        return(res)
    }
)


.factorizeMatrixC <- function(
    counts,
    z,
    K,
    alpha,
    beta,
    sampleLabel,
    rNames,
    sNames,
    type) {

    s <- as.integer(sampleLabel)

    p <- .cCDecomposeCounts(counts, s, z, K)
    mCPByS <- p$mCPByS
    nGByCP <- p$nGByCP
    KNames <- paste0("K", seq(K))
    rownames(nGByCP) <- rNames
    colnames(nGByCP) <- KNames
    rownames(mCPByS) <- KNames
    colnames(mCPByS) <- sNames

    countsList <- c()
    propList <- c()
    postList <- c()
    res <- list()

    if (any("counts" %in% type)) {
        countsList <- list(sample = mCPByS, module = nGByCP)
        res <- c(res, list(counts = countsList))
    }

    if (any("proportion" %in% type)) {
        ## Need to avoid normalizing cell/gene states with zero cells/genes
        uniqueZ <- sort(unique(z))
        tempNGByCP <- nGByCP
        tempNGByCP[, uniqueZ] <- normalizeCounts(tempNGByCP[, uniqueZ],
            normalize = "proportion"
        )

        propList <- list(
            sample = normalizeCounts(mCPByS,
                normalize = "proportion"
            ),
            module = tempNGByCP
        )
        res <- c(res, list(proportions = propList))
    }

    if (any("posterior" %in% type)) {
        postList <- list(
            sample = normalizeCounts(mCPByS + alpha,
                normalize = "proportion"
            ),
            module = normalizeCounts(nGByCP + beta,
                normalize = "proportion"
            )
        )

        res <- c(res, posterior = list(postList))
    }

    return(res)
}


#' @return For celda_G model, a list with elements for "counts", "proportions",
#'  or "posterior" probabilities. Each element will be a list containing
#'  factorized matrices for "module" and "cell".
#' @examples
#' data(celdaGSim, celdaGMod)
#' factorizedMatrices <- factorizeMatrix(
#'   celdaGSim$counts,
#'   celdaGMod, "posterior"
#' )
#' @rdname factorizeMatrix
#' @export
setMethod("factorizeMatrix", signature(x = "ANY", celdaMod = "celda_G"),
    function(x,
        celdaMod,
        type = c("counts", "proportion", "posterior")) {

        counts <- .processCounts(x)
        compareCountMatrix(counts, celdaMod)

        y <- celdaClusters(celdaMod)$y
        # Sometimes, fewer clusters get returned by celda_G
        # Taking the max(y) rather than the original
        # L will prevent errors
        # L <- params(celdaMod)$L
        L <- max(y)
        beta <- params(celdaMod)$beta
        delta <- params(celdaMod)$delta
        gamma <- params(celdaMod)$gamma
        cNames <- matrixNames(celdaMod)$column
        rNames <- matrixNames(celdaMod)$row

        res <- .factorizeMatrixG(
            counts = counts,
            y = y,
            L = L,
            beta = beta,
            delta = delta,
            gamma = gamma,
            cNames = cNames,
            rNames = rNames,
            type = type)
        return(res)
    }
)


.factorizeMatrixG <- function(
    counts,
    y,
    L,
    beta,
    delta,
    gamma,
    cNames,
    rNames,
    type) {

    p <- .cGDecomposeCounts(counts = counts, y = y, L = L)
    nTSByC <- p$nTSByC
    nByG <- p$nByG
    nByTS <- p$nByTS
    nGByTS <- p$nGByTS
    nGByTS[nGByTS == 0] <- 1
    nM <- p$nM
    nG <- p$nG
    rm(p)

    GByTS <- matrix(0, nrow = length(y), ncol = L)
    GByTS[cbind(seq(nG), y)] <- nByG

    LNames <- paste0("L", seq(L))
    colnames(nTSByC) <- cNames
    rownames(nTSByC) <- LNames
    colnames(GByTS) <- LNames
    rownames(GByTS) <- rNames
    names(nGByTS) <- LNames

    countsList <- c()
    propList <- c()
    postList <- c()
    res <- list()

    if (any("counts" %in% type)) {
        countsList <- list(
            cell = nTSByC,
            module = GByTS,
            geneDistribution = nGByTS
        )
        res <- c(res, list(counts = countsList))
    }

    if (any("proportion" %in% type)) {
        ## Need to avoid normalizing cell/gene states with zero cells/genes
        uniqueY <- sort(unique(y))
        tempGByTS <- GByTS
        tempGByTS[, uniqueY] <- normalizeCounts(tempGByTS[, uniqueY],
            normalize = "proportion"
        )
        tempNGByTS <- nGByTS / sum(nGByTS)

        propList <- list(
            cell = normalizeCounts(nTSByC,
                normalize = "proportion"
            ),
            module = tempGByTS,
            geneDistribution = tempNGByTS
        )
        res <- c(res, list(proportions = propList))
    }

    if (any("posterior" %in% type)) {
        gs <- GByTS
        gs[cbind(seq(nG), y)] <- gs[cbind(seq(nG), y)] + delta
        gs <- normalizeCounts(gs, normalize = "proportion")
        tempNGByTS <- (nGByTS + gamma) / sum(nGByTS + gamma)

        postList <- list(
            cell = normalizeCounts(nTSByC + beta,
                normalize = "proportion"
            ),
            module = gs,
            geneDistribution = tempNGByTS
        )
        res <- c(res, posterior = list(postList))
    }

    return(res)
}
