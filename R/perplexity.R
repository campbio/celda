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
#' @param useAssay A string specifying which \link{assay}
#'  slot to use if \code{x} is a \linkS4class{SingleCellExperiment} object.
#'  Default "counts".
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default "featureSubset".
#' @param celdaMod Celda model object. Only works if \code{x} is an integer
#'  counts matrix.
#' @param newCounts A new counts matrix used to calculate perplexity. If NULL,
#'  perplexity will be calculated for the matrix in \code{useAssay} slot in
#'  \code{x}. Default NULL.
#' @return Numeric. The perplexity for the provided \code{x} (and
#'  \code{celdaModel}).
#' @export
setGeneric("perplexity",
    function(x,
        celdaMod,
        useAssay = "counts",
        altExpName = "featureSubset",
        newCounts = NULL) {

        standardGeneric("perplexity")})


#' @importFrom matrixStats logSumExp
#' @examples
#' data(sceCeldaCG)
#' perplexity <- perplexity(sceCeldaCG)
#' @rdname perplexity
#' @export
setMethod("perplexity", signature(x = "SingleCellExperiment"),
    function(x,
        useAssay = "counts",
        altExpName = "featureSubset",
        newCounts = NULL) {

        altExp <- SingleCellExperiment::altExp(x, altExpName)

        counts <- SummarizedExperiment::assay(altExp, i = useAssay)
        counts <- .processCounts(counts)

        if (celdaModel(x, altExpName = altExpName) == "celda_C") {
            factorized <- factorizeMatrix(x = x,
                useAssay = useAssay,
                altExpName = altExpName,
                type = "posterior")
            s <- as.integer(sampleLabel(x, altExpName = altExpName))
            p <- .perplexityCelda_C(
                counts = counts,
                factorized = factorized,
                s = s,
                newCounts = newCounts)
        } else if (celdaModel(x, altExpName = altExpName) == "celda_CG") {
            factorized <- factorizeMatrix(x = x,
                useAssay = useAssay,
                altExpName = altExpName,
                type = c("posterior", "counts"))
            s <- as.integer(sampleLabel(x, altExpName = altExpName))
            p <- .perplexityCelda_CG(counts = counts,
                factorized = factorized,
                s = s,
                newCounts = newCounts)
        } else if (celdaModel(x, altExpName = altExpName) == "celda_G") {
            factorized <- factorizeMatrix(x = x,
                useAssay = useAssay,
                altExpName = altExpName,
                type = c("posterior", "counts"))
            L <- S4Vectors::metadata(altExp)$celda_parameters$L
            y <- celdaModules(x, altExpName = altExpName)
            p <- .perplexityCelda_G(counts,
                factorized,
                L,
                y,
                newCounts = newCounts)
        } else {
            stop("S4Vectors::metadata(altExp(x, altExpName))$",
                "celda_parameters$model must be",
                " one of 'celda_C', 'celda_G', or 'celda_CG'")
        }
        return(p)
    })


#' @importFrom matrixStats logSumExp
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' perplexity <- perplexity(celdaCGSim$counts, celdaCGMod)
#' @rdname perplexity
#' @export
setMethod("perplexity", signature(x = "ANY", celdaMod = "celda_CG"),
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
        innerLogProb <- .countsTimesProbs(newCounts, geneByPopProb) + theta[, s]
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
    "perplexity", signature(x = "ANY", celdaMod = "celda_C"),
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
        inner.log.prob <- .countsTimesProbs(newCounts, phi) + theta[, s]
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
    "perplexity", signature(x = "ANY", celdaMod = "celda_G"),
    function(x, celdaMod, newCounts = NULL) {
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


.perplexityCelda_C <- function(
    counts,
    factorized,
    s,
    newCounts) {

    if (is.null(newCounts)) {
        newCounts <- counts
    } else {
        newCounts <- .processCounts(newCounts)
    }

    if (nrow(newCounts) != nrow(counts)) {
        stop("'newCounts' should have the same number of rows as",
            " 'assay(altExp(x, altExpName), i = useAssay)'.")
    }

    theta <- log(factorized$posterior$sample)
    phi <- log(factorized$posterior$module)

    # inner.log.prob = (t(phi) %*% newCounts) + theta[, s]
    inner.log.prob <- .countsTimesProbs(newCounts, phi) + theta[, s]
    logPx <- sum(apply(inner.log.prob, 2, matrixStats::logSumExp))

    perplexity <- exp(- (logPx / sum(newCounts)))
    return(perplexity)
}


.perplexityCelda_CG <- function(
    counts,
    factorized,
    s,
    newCounts) {

    if (is.null(newCounts)) {
        newCounts <- counts
    } else {
        newCounts <- .processCounts(newCounts)
    }
    if (nrow(newCounts) != nrow(counts)) {
        stop("newCounts should have the same number of rows as",
            " 'assay(altExp(x, altExpName), i = useAssay)'.")
    }

    theta <- log(factorized$posterior$sample)
    phi <- factorized$posterior$cellPopulation
    psi <- factorized$posterior$module
    eta <- factorized$posterior$geneDistribution
    nGByTS <- factorized$counts$geneDistribution

    etaProb <- log(eta) * nGByTS
    geneByPopProb <- log(psi %*% phi)
    innerLogProb <- .countsTimesProbs(newCounts, geneByPopProb) + theta[, s]
    # innerLogProb = (t(geneByPopProb) %*% newCounts) + theta[, s]

    log.px <- sum(apply(innerLogProb, 2, matrixStats::logSumExp))
    # + sum(etaProb)
    perplexity <- exp(- (log.px / sum(newCounts)))
    return(perplexity)
}


.perplexityCelda_G <- function(counts,
    factorized,
    L,
    y,
    newCounts) {

    if (is.null(newCounts)) {
        newCounts <- counts
    } else {
        newCounts <- .processCounts(newCounts)
    }
    if (nrow(newCounts) != nrow(counts)) {
        stop("newCounts should have the same number of rows as",
            " 'assay(altExp(x, altExpName), i = useAssay)'.")
    }

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
        y,
        L
    ) # + sum(etaProb)
    perplexity <- exp(- (logPx / sum(newCounts)))
    return(perplexity)
}


#' @title Calculate and visualize perplexity of all models in a celdaList, with
#'  count resampling
#' @description Calculates the perplexity of each model's cluster assignments
#'  given the provided countMatrix, as well as resamplings of that count
#'  matrix, providing a distribution of perplexities and a better sense of the
#'  quality of a given K/L choice.
#' @param x A numeric \link{matrix} of counts or a
#'  \linkS4class{SingleCellExperiment} returned from \link{celdaGridSearch}
#'  with the matrix located in the assay slot under \code{useAssay}.
#'  Rows represent features and columns represent cells. Must contain
#'  "celda_grid_search" slot in \code{metadata(x)} if \code{x} is a
#'  \linkS4class{SingleCellExperiment} object.
#' @param useAssay A string specifying which \link{assay}
#'  slot to use if \code{x} is a
#'  \linkS4class{SingleCellExperiment} object. Default "counts".
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default "featureSubset".
#' @param celdaList Object of class 'celdaList'. Used only if \code{x} is a
#'  matrix object.
#' @param resample Integer. The number of times to resample the counts matrix
#'  for evaluating perplexity. Default 5.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @return A \linkS4class{SingleCellExperiment} object or
#'  \code{celdaList} object with a \code{perplexity}
#'  property, detailing the perplexity of all K/L combinations that appeared in
#'  the celdaList's models.
#' @export
setGeneric("resamplePerplexity",
    function(x,
        celdaList,
        useAssay = "counts",
        altExpName = "featureSubset",
        resample = 5,
        seed = 12345) {
    standardGeneric("resamplePerplexity")})


#' @rdname resamplePerplexity
#' @examples
#' data(sceCeldaCGGridSearch)
#' sce <- resamplePerplexity(sceCeldaCGGridSearch)
#' plotGridSearchPerplexity(sce)
#' @export
setMethod("resamplePerplexity",
    signature(x = "SingleCellExperiment"),
    function(x,
        useAssay = "counts",
        altExpName = "featureSubset",
        resample = 5,
        seed = 12345) {

        altExp <- SingleCellExperiment::altExp(x, altExpName)

        counts <- SummarizedExperiment::assay(altExp, i = useAssay)
        celdaList <- S4Vectors::metadata(altExp)$celda_grid_search
        if (is.null(seed)) {
            res <- .resamplePerplexity(
                counts = counts,
                celdaList = celdaList,
                resample = resample)
        } else {
            with_seed(seed,
                res <- .resamplePerplexity(
                    counts = counts,
                    celdaList = celdaList,
                    resample = resample)
            )
        }

        S4Vectors::metadata(altExp)$celda_grid_search <- res
        SingleCellExperiment::altExp(x, altExpName) <- altExp
        return(x)
    }
)


#' @rdname resamplePerplexity
#' @examples
#' data(celdaCGSim, celdaCGGridSearchRes)
#' celdaCGGridSearchRes <- resamplePerplexity(
#'   celdaCGSim$counts,
#'   celdaCGGridSearchRes
#' )
#' plotGridSearchPerplexity(celdaCGGridSearchRes)
#' @export
setMethod("resamplePerplexity",
    signature(x = "ANY"),
    function(x,
        celdaList,
        resample = 5,
        seed = 12345) {

        if (is.null(seed)) {
            res <- .resamplePerplexity(
                counts = x,
                celdaList = celdaList,
                resample = resample)
        } else {
            with_seed(seed,
                res <- .resamplePerplexity(
                    counts = x,
                    celdaList = celdaList,
                    resample = resample))
        }

        return(res)
    }
)


.resamplePerplexity <- function(counts,
    celdaList,
    resample = 5) {

    if (!methods::is(celdaList, "celdaList")) {
        stop("celdaList parameter was not of class celdaList.")
    }
    if (!isTRUE(is.numeric(resample))) {
        stop("Provided resample parameter was not numeric.")
    }

    perpRes <- matrix(NA, nrow = length(resList(celdaList)), ncol = resample)
    for (j in seq(resample)) {
        newCounts <- .resampleCountMatrix(counts)
        for (i in seq(length(resList(celdaList)))) {
            perpRes[i, j] <- perplexity(x = counts,
                celdaMod = resList(celdaList)[[i]],
                newCounts = newCounts)
        }
    }
    celdaList@perplexity <- perpRes

    ## Add mean perplexity to runParams
    perpMean <- apply(perpRes, 1, mean)
    celdaList@runParams$mean_perplexity <- perpMean

    return(celdaList)
}


#' @title Visualize perplexity of a list of celda models
#' @description Visualize perplexity of every model in a celdaList, by unique
#'  K/L combinations
#' @param x Can be one of
#' \itemize{
#'  \item A \linkS4class{SingleCellExperiment} object returned from
#'  \code{celdaGridSearch}, \code{recursiveSplitModule},
#'  or \code{recursiveSplitCell}. Must contain a list named
#'  \code{"celda_grid_search"} in \code{metadata(x)}.
#'  \item celdaList object.}
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default "featureSubset". Only works if \code{x} is a
#'  \linkS4class{SingleCellExperiment} object.
#' @param sep Numeric. Breaks in the x axis of the resulting plot.
#' @param alpha Numeric. Passed to \link{geom_jitter}. Opacity of the points.
#'  Values of alpha range from 0 to 1, with lower values corresponding
#'  to more transparent colors.
#' @return A ggplot plot object showing perplexity as a function of clustering
#'  parameters.
#' @export
setGeneric("plotGridSearchPerplexity",
    function(x, altExpName = "featureSubset", sep = 1, alpha = 0.5) {
    standardGeneric("plotGridSearchPerplexity")})


#' @rdname plotGridSearchPerplexity
#' @examples
#' data(sceCeldaCGGridSearch)
#' sce <- resamplePerplexity(sceCeldaCGGridSearch)
#' plotGridSearchPerplexity(sce)
#' @export
setMethod("plotGridSearchPerplexity",
    signature(x = "SingleCellExperiment"),
    function(x, altExpName = "featureSubset", sep = 1, alpha = 0.5) {
        altExp <- SingleCellExperiment::altExp(x, altExpName)
        celdaList <- S4Vectors::metadata(altExp)$celda_grid_search
        g <- do.call(paste0(".plotGridSearchPerplexity",
            as.character(class(resList(x, altExpName = altExpName)[[1]]))),
            args = list(celdaList, sep, alpha))
        return(g)
    }
)


#' @rdname plotGridSearchPerplexity
#' @examples
#' data(celdaCGSim, celdaCGGridSearchRes)
#' ## Run various combinations of parameters with 'celdaGridSearch'
#' celdaCGGridSearchRes <- resamplePerplexity(
#'   celdaCGSim$counts,
#'   celdaCGGridSearchRes)
#' plotGridSearchPerplexity(celdaCGGridSearchRes)
#' @export
setMethod("plotGridSearchPerplexity",
    signature(x = "celdaList"),
    function(x, sep = 1, alpha = 0.5) {
        g <- do.call(paste0(".plotGridSearchPerplexity",
            as.character(class(resList(x)[[1]]))),
            args = list(x, sep, alpha))
        return(g)
    }
)


.plotGridSearchPerplexitycelda_CG <- function(celdaList, sep, alpha) {
    if (!all(c("K", "L") %in% colnames(runParams(celdaList)))) {
        stop("runParams(celdaList) needs K and L columns.")
    }
    if (is.null(celdaPerplexity(celdaList))) {
        stop("No perplexity measurements available. First run",
            " 'resamplePerplexity' with celdaList object.")
    }

    ix1 <- rep(seq(nrow(celdaPerplexity(celdaList))),
        each = ncol(celdaPerplexity(celdaList))
    )
    ix2 <- rep(
        seq(ncol(celdaPerplexity(celdaList))),
        nrow(celdaPerplexity(celdaList))
    )
    df <- data.frame(runParams(celdaList)[ix1, ],
        perplexity = celdaPerplexity(celdaList)[cbind(ix1, ix2)]
    )
    df$K <- as.factor(df$K)
    df$L <- as.factor(df$L)

    lMeansByK <- stats::aggregate(df$perplexity,
        by = list(df$K, df$L),
        FUN = mean
    )
    colnames(lMeansByK) <- c("K", "L", "mean_perplexity")
    lMeansByK$K <- as.factor(lMeansByK$K)
    lMeansByK$L <- as.factor(lMeansByK$L)

    if (nlevels(df$K) > 1) {
        if (nlevels(df$L) > 1) {
            plot <- ggplot2::ggplot(
                df,
                ggplot2::aes_string(x = "K", y = "perplexity")
            ) +
                ggplot2::geom_jitter(
                    height = 0, width = 0.1, alpha = alpha,
                    ggplot2::aes_string(color = "L")
                ) +
                ggplot2::scale_color_discrete(name = "L") +
                ggplot2::geom_path(data = lMeansByK, ggplot2::aes_string(
                    x = "K",
                    y = "mean_perplexity", group = "L", color = "L"
                )) +
                ggplot2::ylab("Perplexity") +
                ggplot2::xlab("K") +
                ggplot2::scale_x_discrete(breaks = seq(
                    min(runParams(celdaList)$K),
                    max(runParams(celdaList)$K), sep
                )) +
                ggplot2::theme_bw() +
                ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                    panel.grid.minor = ggplot2::element_blank())
        } else {
            plot <- ggplot2::ggplot(
                df,
                ggplot2::aes_string(x = "K", y = "perplexity")) +
                ggplot2::geom_jitter(height = 0, width = 0.1,
                    color = "grey", alpha = alpha) +
                ggplot2::scale_color_manual(name = "L", values = "black") +
                ggplot2::geom_path(data = lMeansByK, ggplot2::aes_string(
                    x = "K",
                    y = "mean_perplexity", group = "L", color = "L"
                )) +
                ggplot2::ylab("Perplexity") +
                ggplot2::xlab("K") +
                ggplot2::scale_x_discrete(breaks = seq(
                    min(runParams(celdaList)$K),
                    max(runParams(celdaList)$K), sep
                )) +
                ggplot2::theme_bw() +
                ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                    panel.grid.minor = ggplot2::element_blank())
        }
    } else {
        plot <- ggplot2::ggplot(
            df,
            ggplot2::aes_string(x = "L", y = "perplexity")) +
            ggplot2::geom_jitter(height = 0, width = 0.1,
                color = "grey", alpha = alpha) +
            ggplot2::geom_path(data = lMeansByK,
                ggplot2::aes_string(x = "L",
                    y = "mean_perplexity",
                    group = "K",
                    color = "K")) +
            ggplot2::scale_color_manual(name = "K", values = "black") +
            ggplot2::ylab("Perplexity") +
            ggplot2::xlab("L") +
            ggplot2::scale_x_discrete(breaks = seq(min(runParams(celdaList)$L),
                max(runParams(celdaList)$L), sep)) +
            ggplot2::theme_bw() +
            ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank())
    }
    return(plot)
}


.plotGridSearchPerplexitycelda_C <- function(celdaList, sep, alpha) {
    if (!all(c("K") %in% colnames(runParams(celdaList)))) {
        stop("runParams(celdaList) needs the column K.")
    }
    if (is.null(celdaPerplexity(celdaList))) {
        stop(
            "No perplexity measurements available. First run",
            " 'resamplePerplexity' with celdaList object."
        )
    }

    ix1 <- rep(seq(nrow(celdaPerplexity(celdaList))),
        each = ncol(celdaPerplexity(celdaList))
    )
    ix2 <- rep(
        seq(ncol(celdaPerplexity(celdaList))),
        nrow(celdaPerplexity(celdaList))
    )
    df <- data.frame(runParams(celdaList)[ix1, ],
        perplexity = celdaPerplexity(celdaList)[cbind(ix1, ix2)]
    )
    df$K <- as.factor(df$K)

    meansByK <- stats::aggregate(df$perplexity, by = list(df$K), FUN = mean)
    colnames(meansByK) <- c("K", "mean_perplexity")
    meansByK$K <- as.factor(meansByK$K)

    plot <-
        ggplot2::ggplot(df, ggplot2::aes_string(x = "K", y = "perplexity")) +
        ggplot2::geom_jitter(height = 0, width = 0.1,
            color = "grey", alpha = alpha) +
        ggplot2::geom_path(
            data = meansByK,
            ggplot2::aes_string(x = "K", y = "mean_perplexity", group = 1)
        ) +
        ggplot2::ylab("Perplexity") +
        ggplot2::xlab("K") +
        ggplot2::scale_x_discrete(breaks = seq(
            min(runParams(celdaList)$K),
            max(runParams(celdaList)$K), sep
        )) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank())

    return(plot)
}


.plotGridSearchPerplexitycelda_G <- function(celdaList, sep, alpha) {
    if (!all(c("L") %in% colnames(runParams(celdaList)))) {
        stop("runParams(celdaList) needs the column L.")
    }
    if (length(celdaPerplexity(celdaList)) == 0) {
        stop(
            "No perplexity measurements available. First run",
            " 'resamplePerplexity' with celdaList object."
        )
    }

    ix1 <- rep(seq(nrow(celdaPerplexity(celdaList))),
        each = ncol(celdaPerplexity(celdaList))
    )
    ix2 <- rep(
        seq(ncol(celdaPerplexity(celdaList))),
        nrow(celdaPerplexity(celdaList))
    )
    df <- data.frame(runParams(celdaList)[ix1, ],
        perplexity = celdaPerplexity(celdaList)[cbind(ix1, ix2)]
    )
    df$L <- as.factor(df$L)

    meansByL <- stats::aggregate(df$perplexity, by = list(df$L), FUN = mean)
    colnames(meansByL) <- c("L", "mean_perplexity")
    meansByL$L <- as.factor(meansByL$L)

    plot <-
        ggplot2::ggplot(df, ggplot2::aes_string(x = "L", y = "perplexity")) +
        ggplot2::geom_jitter(height = 0, width = 0.1,
            color = "grey", alpha = alpha) +
        ggplot2::geom_path(
            data = meansByL,
            ggplot2::aes_string(x = "L", y = "mean_perplexity", group = 1)
        ) +
        ggplot2::ylab("Perplexity") +
        ggplot2::xlab("L") +
        ggplot2::scale_x_discrete(breaks = seq(
            min(runParams(celdaList)$L),
            max(runParams(celdaList)$L), sep
        )) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank())

    return(plot)
}


# Resample a counts matrix for evaluating perplexity
# Normalizes each column (cell) of a countMatrix by the column sum to
# create a distribution of observing a given number of counts for a given
# gene in that cell,
# then samples across all cells.
# This is primarily used to evaluate the stability of the perplexity for
# a given K/L combination.
# @param celda.mod A single celda run (usually from the _resList_ property
# of a celdaList).
# @return The perplexity for the provided chain as an mpfr number.
#' @importFrom Matrix colSums t
.resampleCountMatrix <- function(countMatrix) {
    colsums <- colSums(countMatrix)
    prob <- t(t(countMatrix) / colsums)
    resample <- vapply(seq(ncol(countMatrix)), function(idx) {
        stats::rmultinom(
            n = 1,
            size = colsums[idx],
            prob = prob[, idx]
        )
    }, integer(nrow(countMatrix)))
    return(resample)
}


#' @title Visualize perplexity differences of a list of celda models
#' @description Visualize perplexity differences of every model in a celdaList,
#'  by unique K/L combinations. Line represents centered moving average with
#'  windows of length \code{n}.
#' @param x Can be one of
#' \itemize{
#'  \item A \linkS4class{SingleCellExperiment} object returned from
#'  \code{celdaGridSearch}, \code{recursiveSplitModule},
#'  or \code{recursiveSplitCell}. Must contain a list named
#'  \code{"celda_grid_search"} in \code{metadata(x)}.
#'  \item celdaList object.}
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default "featureSubset".
#' @param sep Numeric. Breaks in the x axis of the resulting plot.
#' @param n Integer. Width of the rolling window. Default 10.
#' @param alpha Numeric. Passed to \link{geom_jitter}. Opacity of the points.
#'  Values of alpha range from 0 to 1, with lower values corresponding
#'  to more transparent colors.
#' @return A ggplot plot object showing perplexity diferences as a function of
#'  clustering parameters.
#' @export
setGeneric("plotRPC",
    function(x, altExpName = "featureSubset", sep = 1, n = 10, alpha = 0.5) {
    standardGeneric("plotRPC")})


#' @rdname plotRPC
#' @examples
#' data(sceCeldaCGGridSearch)
#' sce <- resamplePerplexity(sceCeldaCGGridSearch)
#' plotRPC(sce, n = 1)
#' @export
setMethod("plotRPC",
    signature(x = "SingleCellExperiment"),
    function(x, altExpName = "featureSubset", sep = 1, n = 10, alpha = 0.5) {
        altExp <- SingleCellExperiment::altExp(x, altExpName)
        model <- altExp@metadata$celda_grid_search@celdaGridSearchParameters$
            model
        celdaList <- S4Vectors::metadata(altExp)$celda_grid_search

        if (model == "celda_C") {
            g <- .plotRPCC(celdaList, sep, n = n,
                alpha = alpha)
        } else if (model == "celda_G") {
            g <- .plotRPCG(celdaList, sep, n = n,
                alpha = alpha)
        } else if (model == "celda_CG") {
            g <- .plotRPCCG(celdaList, sep, n = n,
                alpha = alpha)
        } else {
            stop("S4Vectors::metadata(altExp(x, altExpName))$",
                "celda_grid_search@",
                "celdaGridSearchParameters$model must be",
                " one of 'celda_C', 'celda_G', or 'celda_CG'")
        }
        return(g)
    }
)


#' @rdname plotRPC
#' @examples
#' data(celdaCGSim, celdaCGGridSearchRes)
#' ## Run various combinations of parameters with 'celdaGridSearch'
#' celdaCGGridSearchRes <- resamplePerplexity(
#'   celdaCGSim$counts,
#'   celdaCGGridSearchRes)
#' plotRPC(celdaCGGridSearchRes, n = 1)
#' @export
setMethod("plotRPC",
    signature(x = "celdaList"),
    function(x, sep = 1, n = 10, alpha = 0.5) {
        g <- do.call(paste0(".plotRPC",
            unlist(strsplit(as.character(class(resList(x)[[1]])), "_"))[[2]]),
            args = list(x, sep, n, alpha))
        return(g)
    }
)


.plotRPCCG <- function(celdaList, sep, n, alpha) {
    # fix check note
    K <- L <- perpdiffK <- meanperpdiffK <- perpdiffL <- meanperpdiffL <- NULL

    if (!all(c("K", "L") %in% colnames(runParams(celdaList)))) {
        stop("runParams(celdaList) needs K and L columns.")
    }
    if (is.null(celdaPerplexity(celdaList))) {
        stop("No perplexity measurements available. First run",
            " 'resamplePerplexity' with celdaList object.")
    }

    ix1 <- rep(seq(nrow(celdaPerplexity(celdaList))),
        each = ncol(celdaPerplexity(celdaList)))
    ix2 <- rep(seq(ncol(celdaPerplexity(celdaList))),
        nrow(celdaPerplexity(celdaList)))
    dt <- data.table::data.table(runParams(celdaList)[ix1, ],
        perplexity = celdaPerplexity(celdaList)[cbind(ix1, ix2)])
    dt$K <- as.factor(dt$K)
    dt$L <- as.factor(dt$L)

    if (nlevels(dt$K) > 1) {
        for (i in seq(nlevels(dt$L))) {
            for (j in seq(2, nlevels(dt$K))) {
                p1 <- dt[K == levels(dt$K)[j - 1] & L == levels(dt$L)[i],
                    perplexity]
                p2 <- dt[K == levels(dt$K)[j] & L == levels(dt$L)[i],
                    perplexity]
                dt[K == levels(dt$K)[j] & L == levels(dt$L)[i],
                    "perpdiffK"] <- p2 - p1
            }
        }

        diffMeansByK <- data.table::data.table(stats::aggregate(
            dt$perpdiffK,
            by = list(dt$K, dt$L),
            FUN = mean))
        colnames(diffMeansByK) <- c("K", "L", "meanperpdiffK")
        diffMeansByK$K <- as.factor(diffMeansByK$K)
        diffMeansByK$L <- as.factor(diffMeansByK$L)
        diffMeansByK$rollmean <- data.table::frollmean(
            diffMeansByK$meanperpdiffK, n = n, align = "center")
        diffMeansByK <- diffMeansByK[stats::complete.cases(diffMeansByK), ]

        if (nlevels(dt$L) > 1) {
            plot <- ggplot2::ggplot(dt[!is.na(perpdiffK), ],
                ggplot2::aes_string(x = "K",
                    y = "perpdiffK")) +
                ggplot2::geom_jitter(height = 0, width = 0.1, alpha = alpha,
                    ggplot2::aes_string(color = "L")) +
                ggplot2::scale_color_discrete(name = "L") +
                ggplot2::geom_path(data = diffMeansByK,
                    ggplot2::aes_string(x = "K", y = "rollmean", group = "L",
                        color = "L"), size = 1) +
                ggplot2::ylab("Rate of perplexity change") +
                ggplot2::xlab("K") +
                ggplot2::scale_x_discrete(
                    breaks = seq(min(as.integer(levels(dt$K))),
                        max(as.integer(levels(dt$K))), sep)) +
                ggplot2::theme_bw() +
                ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                    panel.grid.minor = ggplot2::element_blank())
        } else {
            plot <- ggplot2::ggplot(dt[!is.na(perpdiffK), ],
                ggplot2::aes_string(x = "K",
                    y = "perpdiffK")) +
                ggplot2::geom_jitter(height = 0, width = 0.1,
                    color = "grey", alpha = alpha) +
                ggplot2::scale_color_manual(name = "L", values = "black") +
                ggplot2::geom_path(data = diffMeansByK,
                    ggplot2::aes_string(x = "K", y = "rollmean", group = "L",
                        color = "L"), size = 1) +
                ggplot2::ylab("Rate of perplexity change") +
                ggplot2::xlab("K") +
                ggplot2::scale_x_discrete(
                    breaks = seq(min(as.integer(levels(dt$K))),
                        max(as.integer(levels(dt$K))), sep)) +
                ggplot2::theme_bw() +
                ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                    panel.grid.minor = ggplot2::element_blank())
        }
    } else if (nlevels(dt$L) > 1) {
        for (j in seq(2, nlevels(dt$L))) {
            p1 <- dt[K == levels(dt$K) & L == levels(dt$L)[j - 1],
                perplexity]
            p2 <- dt[K == levels(dt$K) & L == levels(dt$L)[j],
                perplexity]
            dt[K == levels(dt$K) & L == levels(dt$L)[j],
                "perpdiffL"] <- p2 - p1
        }

        diffMeansByL <- data.table::data.table(stats::aggregate(dt$perpdiffL,
            by = list(dt$K, dt$L),
            FUN = mean))
        colnames(diffMeansByL) <- c("K", "L", "meanperpdiffL")
        diffMeansByL$K <- as.factor(diffMeansByL$K)
        diffMeansByL$L <- as.factor(diffMeansByL$L)
        diffMeansByL$rollmean <- data.table::frollmean(
            diffMeansByL$meanperpdiffL, n = n, align = "center")
        diffMeansByL <- diffMeansByL[stats::complete.cases(diffMeansByL), ]

        plot <- ggplot2::ggplot(dt[!is.na(perpdiffL), ],
            ggplot2::aes_string(x = "L", y = "perpdiffL")) +
            ggplot2::geom_jitter(height = 0, width = 0.1,
                color = "grey", alpha = alpha) +
            ggplot2::scale_color_manual(name = "K", values = "black") +
            ggplot2::geom_path(
                data = diffMeansByL,
                ggplot2::aes_string(
                    x = "L", y = "rollmean", group = "K", color = "K"),
                size = 1) +
            ggplot2::ylab("Rate of perplexity change") +
            ggplot2::xlab("L") +
            ggplot2::scale_x_discrete(
                breaks = seq(min(as.integer(levels(dt$L))),
                max(as.integer(levels(dt$L))), sep)) +
            ggplot2::theme_bw() +
            ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank())
    } else {
        stop("Only one combination of K and L available! Unable to calculate",
            " perplexity differences.")
    }

    return(plot)
}


.plotRPCC <- function(celdaList, sep, n, alpha) {
    K <- perpdiffK <- meanperpdiffK <- NULL # fix check note
    if (!all(c("K") %in% colnames(runParams(celdaList)))) {
        stop("runParams(celdaList) needs the column K.")
    }
    if (is.null(celdaPerplexity(celdaList))) {
        stop(
            "No perplexity measurements available. First run",
            " 'resamplePerplexity' with celdaList object."
        )
    }

    ix1 <- rep(seq(nrow(celdaPerplexity(celdaList))),
        each = ncol(celdaPerplexity(celdaList)))
    ix2 <- rep(seq(ncol(celdaPerplexity(celdaList))),
        nrow(celdaPerplexity(celdaList)))
    dt <- data.table::data.table(runParams(celdaList)[ix1, ],
        perplexity = celdaPerplexity(celdaList)[cbind(ix1, ix2)])
    dt$K <- as.factor(dt$K)

    if (nlevels(dt$K) > 1) {
        for (i in seq(2, nlevels(dt$K))) {
            p1 <- dt[K == levels(dt$K)[i - 1], perplexity]
            p2 <- dt[K == levels(dt$K)[i], perplexity]
            dt[K == levels(dt$K)[i], "perpdiffK"] <- p2 - p1
        }

        diffMeansByK <- data.table::data.table(stats::aggregate(dt$perpdiffK,
            by = list(dt$K),
            FUN = mean))
        colnames(diffMeansByK) <- c("K", "meanperpdiffK")
        diffMeansByK$K <- as.factor(diffMeansByK$K)
        diffMeansByK$rollmean <- data.table::frollmean(
            diffMeansByK$meanperpdiffK, n = n, align = "center")
        diffMeansByK <- diffMeansByK[stats::complete.cases(diffMeansByK), ]

        plot <- ggplot2::ggplot(dt[!is.na(perpdiffK), ],
            ggplot2::aes_string(x = "K",
                y = "perpdiffK")) +
            ggplot2::geom_jitter(height = 0, width = 0.1,
                color = "grey", alpha = alpha) +
            ggplot2::geom_path(data = diffMeansByK,
                ggplot2::aes_string(x = "K", y = "rollmean", group = 1),
                size = 1) +
            ggplot2::ylab("Perplexity difference compared to previous K") +
            ggplot2::xlab("K") +
            ggplot2::scale_x_discrete(
                breaks = seq(min(as.integer(levels(dt$K))),
                max(as.integer(levels(dt$K))), sep)) +
            ggplot2::theme_bw() +
            ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank())
    } else {
        stop("Only one unique K value available! Unable to calculate",
            " perplexity differences.")
    }
    return(plot)
}


.plotRPCG <- function(celdaList, sep, n, alpha) {
    L <- perpdiffL <- meanperpdiffL <- NULL # fix check note
    if (!all(c("L") %in% colnames(runParams(celdaList)))) {
        stop("runParams(celdaList) needs the column L.")
    }
    if (length(celdaPerplexity(celdaList)) == 0) {
        stop(
            "No perplexity measurements available. First run",
            " 'resamplePerplexity' with celdaList object."
        )
    }

    ix1 <- rep(seq(nrow(celdaPerplexity(celdaList))),
        each = ncol(celdaPerplexity(celdaList)))
    ix2 <- rep(seq(ncol(celdaPerplexity(celdaList))),
        nrow(celdaPerplexity(celdaList)))
    dt <- data.table::data.table(runParams(celdaList)[ix1, ],
        perplexity = celdaPerplexity(celdaList)[cbind(ix1, ix2)])
    dt$L <- as.factor(dt$L)

    if (nlevels(dt$L) > 1) {
        for (i in seq(2, nlevels(dt$L))) {
            p1 <- dt[L == levels(dt$L)[i - 1], perplexity]
            p2 <- dt[L == levels(dt$L)[i], perplexity]
            dt[L == levels(dt$L)[i], "perpdiffL"] <- p2 - p1
        }

        diffMeansByL <- data.table::data.table(stats::aggregate(dt$perpdiffL,
            by = list(dt$L),
            FUN = mean))
        colnames(diffMeansByL) <- c("L", "meanperpdiffL")
        diffMeansByL$L <- as.factor(diffMeansByL$L)
        diffMeansByL$rollmean <- data.table::frollmean(
            diffMeansByL$meanperpdiffL, n = n, align = "center")
        diffMeansByL <- diffMeansByL[stats::complete.cases(diffMeansByL), ]

        plot <- ggplot2::ggplot(dt[!is.na(perpdiffL), ],
            ggplot2::aes_string(x = "L",
                y = "perpdiffL")) +
            ggplot2::geom_jitter(height = 0, width = 0.1,
                color = "grey", alpha = alpha) +
            ggplot2::geom_path(data = diffMeansByL,
                ggplot2::aes_string(x = "L", y = "rollmean", group = 1),
                size = 1) +
            ggplot2::ylab("Perplexity difference compared to previous L") +
            ggplot2::xlab("L") +
            ggplot2::scale_x_discrete(
                breaks = seq(min(as.integer(levels(dt$L))),
                max(as.integer(levels(dt$L))), sep)) +
            ggplot2::theme_bw() +
            ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank())
    } else {
        stop("Only one unique L value available! Unable to calculate",
            " perplexity differences.")
    }
    return(plot)
}
