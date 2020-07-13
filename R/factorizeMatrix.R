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
#' @param useAssay A string specifying which \link[SummarizedExperiment]{assay}
#'  slot to use if \code{x} is a \linkS4class{SingleCellExperiment} object.
#'  Default "counts".
#' @param altExpName The name for the \link[SingleCellExperiment]{altExp} slot
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

        if (celdaModel(x) == "celda_C") {
            res <- .factorizeMatrixCelda_C(sce = altExp, useAssay = useAssay,
                type = type)
        } else if (celdaModel(x) == "celda_CG") {
            res <- .factorizeMatrixCelda_CG(sce = altExp, useAssay = useAssay,
                type = type)
        } else if (celdaModel(x) == "celda_G") {
            res <- .factorizeMatrixCelda_G(sce = altExp, useAssay = useAssay,
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
setMethod("factorizeMatrix", signature(x = "matrix", celdaMod = "celda_CG"),
    function(x,
        celdaMod,
        type = c("counts", "proportion", "posterior")) {

        counts <- .processCounts(x)
        compareCountMatrix(counts, celdaMod)

        K <- params(celdaMod)$K
        L <- params(celdaMod)$L
        z <- celdaClusters(celdaMod)$z
        y <- celdaClusters(celdaMod)$y
        alpha <- params(celdaMod)$alpha
        beta <- params(celdaMod)$beta
        delta <- params(celdaMod)$delta
        gamma <- params(celdaMod)$gamma
        sampleLabel <- sampleLabel(celdaMod)
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

        nGByTS <- matrix(0, nrow = length(y), ncol = L)
        nGByTS[cbind(seq(nG), y)] <- p$nByG

        LNames <- paste0("L", seq(L))
        KNames <- paste0("K", seq(K))
        colnames(nTSByC) <- matrixNames(celdaMod)$column
        rownames(nTSByC) <- LNames
        colnames(nGByTS) <- LNames
        rownames(nGByTS) <- matrixNames(celdaMod)$row
        rownames(mCPByS) <- KNames
        colnames(mCPByS) <- matrixNames(celdaMod)$sample
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
                module = nGByTS,
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
            tempNGByTS <- nGByTS
            tempNGByTS[, uniqueY] <- normalizeCounts(tempNGByTS[, uniqueY],
                normalize = "proportion"
            )
            tempNGByTS <- nGByTS / sum(nGByTS)

            propList <- list(
                sample = normalizeCounts(mCPByS,
                    normalize = "proportion"
                ),
                cellPopulation = tempNTSByCP,
                cell = normalizeCounts(nTSByC, normalize = "proportion"),
                module = tempNGByTS,
                geneDistribution = tempNGByTS
            )
            res <- c(res, list(proportions = propList))
        }

        if (any("posterior" %in% type)) {
            gs <- nGByTS
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
)


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
setMethod("factorizeMatrix", signature(x = "matrix", celdaMod = "celda_C"),
    function(x,
        celdaMod,
        type = c("counts", "proportion", "posterior")) {

        counts <- .processCounts(x)
        compareCountMatrix(counts, celdaMod)

        K <- params(celdaMod)$K
        z <- celdaClusters(celdaMod)$z
        alpha <- params(celdaMod)$alpha
        beta <- params(celdaMod)$beta
        sampleLabel <- sampleLabel(celdaMod)
        s <- as.integer(sampleLabel)

        p <- .cCDecomposeCounts(counts, s, z, K)
        mCPByS <- p$mCPByS
        nGByCP <- p$nGByCP

        KNames <- paste0("K", seq(K))
        rownames(nGByCP) <- matrixNames(celdaMod)$row
        colnames(nGByCP) <- KNames
        rownames(mCPByS) <- KNames
        colnames(mCPByS) <- matrixNames(celdaMod)$sample

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
)


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
setMethod("factorizeMatrix", signature(x = "matrix", celdaMod = "celda_G"),
    function(x,
        celdaMod,
        type = c("counts", "proportion", "posterior")) {

        counts <- .processCounts(x)
        # compareCountMatrix(counts, celdaMod)

        L <- params(celdaMod)$L
        y <- celdaClusters(celdaMod)$y
        beta <- params(celdaMod)$beta
        delta <- params(celdaMod)$delta
        gamma <- params(celdaMod)$gamma

        p <- .cGDecomposeCounts(counts = counts, y = y, L = L)
        nTSByC <- p$nTSByC
        nByG <- p$nByG
        nByTS <- p$nByTS
        nGByTS <- p$nGByTS
        nM <- p$nM
        nG <- p$nG
        rm(p)

        nGByTS[nGByTS == 0] <- 1
        nGByTS <- matrix(0, nrow = length(y), ncol = L)
        nGByTS[cbind(seq(nG), y)] <- nByG

        LNames <- paste0("L", seq(L))
        colnames(nTSByC) <- matrixNames(celdaMod)$column
        rownames(nTSByC) <- LNames
        colnames(nGByTS) <- LNames
        rownames(nGByTS) <- matrixNames(celdaMod)$row
        names(nGByTS) <- LNames

        countsList <- c()
        propList <- c()
        postList <- c()
        res <- list()

        if (any("counts" %in% type)) {
            countsList <- list(
                cell = nTSByC,
                module = nGByTS,
                geneDistribution = nGByTS
            )
            res <- c(res, list(counts = countsList))
        }

        if (any("proportion" %in% type)) {
            ## Need to avoid normalizing cell/gene states with zero cells/genes
            uniqueY <- sort(unique(y))
            tempNGByTS <- nGByTS
            tempNGByTS[, uniqueY] <- normalizeCounts(tempNGByTS[, uniqueY],
                normalize = "proportion"
            )
            tempNGByTS <- nGByTS / sum(nGByTS)

            propList <- list(
                cell = normalizeCounts(nTSByC,
                    normalize = "proportion"
                ),
                module = tempNGByTS,
                geneDistribution = tempNGByTS
            )
            res <- c(res, list(proportions = propList))
        }

        if (any("posterior" %in% type)) {
            gs <- nGByTS
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
)


.factorizeMatrixCelda_C <- function(sce, useAssay, type) {
    counts <- SummarizedExperiment::assay(sce, i = useAssay)
    counts <- .processCounts(counts)

    K <- S4Vectors::metadata(sce)$celda_parameters$K
    z <- SummarizedExperiment::colData(sce)$celda_cell_cluster
    alpha <- S4Vectors::metadata(sce)$celda_parameters$alpha
    beta <- S4Vectors::metadata(sce)$celda_parameters$beta
    sampleLabel <- SummarizedExperiment::colData(sce)$celda_sample_label
    s <- as.integer(sampleLabel)

    p <- .cCDecomposeCounts(counts, s, z, K)
    mCPByS <- p$mCPByS
    nGByCP <- p$nGByCP

    KNames <- paste0("K", seq(K))
    rownames(nGByCP) <- rownames(sce)
    colnames(nGByCP) <- KNames
    rownames(mCPByS) <- KNames
    colnames(mCPByS) <-
        S4Vectors::metadata(sce)$celda_parameters$sampleLevels

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
            normalize = "proportion")

        propList <- list(sample = normalizeCounts(mCPByS,
            normalize = "proportion"),
            module = tempNGByCP)
        res <- c(res, list(proportions = propList))
    }

    if (any("posterior" %in% type)) {
        postList <- list(sample = normalizeCounts(mCPByS + alpha,
            normalize = "proportion"),
            module = normalizeCounts(nGByCP + beta,
                normalize = "proportion"))

        res <- c(res, posterior = list(postList))
    }

    return(res)
}


.factorizeMatrixCelda_CG <- function(sce, useAssay, type) {
    counts <- SummarizedExperiment::assay(sce, i = useAssay)
    counts <- .processCounts(counts)

    K <- S4Vectors::metadata(sce)$celda_parameters$K
    L <- S4Vectors::metadata(sce)$celda_parameters$L
    z <- SummarizedExperiment::colData(sce)$celda_cell_cluster
    y <- SummarizedExperiment::rowData(sce)$celda_feature_module
    alpha <- S4Vectors::metadata(sce)$celda_parameters$alpha
    beta <- S4Vectors::metadata(sce)$celda_parameters$beta
    delta <- S4Vectors::metadata(sce)$celda_parameters$delta
    gamma <- S4Vectors::metadata(sce)$celda_parameters$gamma
    sampleLabel <- SummarizedExperiment::colData(sce)$celda_sample_label
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

    nGByTS <- matrix(0, nrow = length(y), ncol = L)
    nGByTS[cbind(seq(nG), y)] <- p$nByG

    LNames <- paste0("L", seq(L))
    KNames <- paste0("K", seq(K))
    colnames(nTSByC) <- colnames(sce)
    rownames(nTSByC) <- LNames
    colnames(nGByTS) <- LNames
    rownames(nGByTS) <- rownames(sce)
    rownames(mCPByS) <- KNames
    colnames(mCPByS) <- S4Vectors::metadata(sce)$celda_parameters$sampleLevels
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
            module = nGByTS,
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
        tempNGByTS <- nGByTS
        tempNGByTS[, uniqueY] <- normalizeCounts(tempNGByTS[, uniqueY],
            normalize = "proportion"
        )
        tempNGByTS <- nGByTS / sum(nGByTS)

        propList <- list(
            sample = normalizeCounts(mCPByS, normalize = "proportion"),
            cellPopulation = tempNTSByCP,
            cell = normalizeCounts(nTSByC, normalize = "proportion"),
            module = tempNGByTS,
            geneDistribution = tempNGByTS
        )
        res <- c(res, list(proportions = propList))
    }

    if (any("posterior" %in% type)) {
        gs <- nGByTS
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


.factorizeMatrixCelda_G <- function(sce, useAssay, type) {
    counts <- SummarizedExperiment::assay(sce, i = useAssay)
    counts <- .processCounts(counts)
    # compareCountMatrix(counts, celdaMod)

    L <- S4Vectors::metadata(sce)$celda_parameters$L
    y <- SummarizedExperiment::rowData(sce)$celda_feature_module
    beta <- S4Vectors::metadata(sce)$celda_parameters$beta
    delta <- S4Vectors::metadata(sce)$celda_parameters$delta
    gamma <- S4Vectors::metadata(sce)$celda_parameters$gamma

    p <- .cGDecomposeCounts(counts = counts, y = y, L = L)
    nTSByC <- p$nTSByC
    nByG <- p$nByG
    nByTS <- p$nByTS
    nGByTS <- p$nGByTS
    nM <- p$nM
    nG <- p$nG
    rm(p)

    nGByTS[nGByTS == 0] <- 1
    nGByTS <- matrix(0, nrow = length(y), ncol = L)
    nGByTS[cbind(seq(nG), y)] <- nByG

    LNames <- paste0("L", seq(L))
    colnames(nTSByC) <- colnames(sce)
    rownames(nTSByC) <- LNames
    colnames(nGByTS) <- LNames
    rownames(nGByTS) <- rownames(sce)
    names(nGByTS) <- LNames

    countsList <- c()
    propList <- c()
    postList <- c()
    res <- list()

    if (any("counts" %in% type)) {
        countsList <- list(
            cell = nTSByC,
            module = nGByTS,
            geneDistribution = nGByTS
        )
        res <- c(res, list(counts = countsList))
    }

    if (any("proportion" %in% type)) {
        ## Need to avoid normalizing cell/gene states with zero cells/genes
        uniqueY <- sort(unique(y))
        tempNGByTS <- nGByTS
        tempNGByTS[, uniqueY] <- normalizeCounts(tempNGByTS[, uniqueY],
            normalize = "proportion"
        )
        tempNGByTS <- nGByTS / sum(nGByTS)

        propList <- list(
            cell = normalizeCounts(nTSByC,
                normalize = "proportion"
            ),
            module = tempNGByTS,
            geneDistribution = tempNGByTS
        )
        res <- c(res, list(proportions = propList))
    }

    if (any("posterior" %in% type)) {
        gs <- nGByTS
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
