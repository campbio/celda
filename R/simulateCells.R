#' @title Simulate count data from the celda generative models.
#' @description This function generates a \linkS4class{SingleCellExperiment}
#'  containing a simulated counts matrix in the \code{"counts"} assay slot, as
#'  well as various parameters used in the simulation which can be
#'  useful for running celda and are stored in \code{metadata} slot. The user
#'  must provide the desired model (one of celda_C, celda_G, celda_CG) as well
#'  as any desired tuning parameters for those model's simulation functions
#'  as detailed below.
#' @param model Character. Options available in \code{celda::availableModels}.
#'  Can be one of \code{"celda_CG"}, \code{"celda_C"}, or \code{"celda_G"}.
#'  Default \code{"celda_CG"}.
#' @param S Integer. Number of samples to simulate. Default 5. Only used if
#'  \code{model} is one of \code{"celda_CG"} or \code{"celda_C"}.
#' @param CRange Integer vector. A vector of length 2 that specifies the lower
#'  and upper bounds of the number of cells to be generated in each sample.
#'  Default c(50, 100). Only used if
#'  \code{model} is one of \code{"celda_CG"} or \code{"celda_C"}.
#' @param NRange Integer vector. A vector of length 2 that specifies the lower
#'  and upper bounds of the number of counts generated for each cell. Default
#'  c(500, 1000).
#' @param C Integer. Number of cells to simulate. Default 100. Only used if
#'  \code{model} is \code{"celda_G"}.
#' @param G Integer. The total number of features to be simulated. Default 100.
#' @param K Integer. Number of cell populations. Default 5. Only used if
#'  \code{model} is one of \code{"celda_CG"} or \code{"celda_C"}.
#' @param L Integer. Number of feature modules. Default 10. Only used if
#'  \code{model} is one of \code{"celda_CG"} or \code{"celda_G"}.
#' @param alpha Numeric. Concentration parameter for Theta. Adds a pseudocount
#'  to each cell population in each sample. Default 1. Only used if
#'  \code{model} is one of \code{"celda_CG"} or \code{"celda_C"}.
#' @param beta Numeric. Concentration parameter for Phi. Adds a pseudocount to
#'  each feature module in each cell population. Default 1.
#' @param gamma Numeric. Concentration parameter for Eta. Adds a pseudocount to
#'  the number of features in each module. Default 5. Only used if
#'  \code{model} is one of \code{"celda_CG"} or \code{"celda_G"}.
#' @param delta Numeric. Concentration parameter for Psi. Adds a pseudocount to
#'  each feature in each module. Default 1. Only used if
#'  \code{model} is one of \code{"celda_CG"} or \code{"celda_G"}.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with
#'  simulated count matrix stored in the "counts" assay slot. Function
#'  parameter settings are stored in the \link{metadata} slot. For
#'  \code{"celda_CG"} and \code{"celda_C"} models,
#'  columns \code{celda_sample_label} and \code{celda_cell_cluster} in
#'  \link{colData} contain simulated sample labels and
#'  cell population clusters. For \code{"celda_CG"} and \code{"celda_G"}
#'  models, column \code{celda_feature_module} in
#'  \link{rowData} contains simulated gene modules.
#' @examples
#' sce <- simulateCells()
#' @export
simulateCells <- function(
    model = c("celda_CG", "celda_C", "celda_G"),
    S = 5,
    CRange = c(50, 100),
    NRange = c(500, 1000),
    C = 100,
    G = 100,
    K = 5,
    L = 10,
    alpha = 1,
    beta = 1,
    gamma = 5,
    delta = 1,
    seed = 12345) {

    model <- match.arg(model)

    if (model == "celda_C") {
        sce <- .simulateCellsMaincelda_C(model = model,
            S = S,
            CRange = CRange,
            NRange = NRange,
            G = G,
            K = K,
            alpha = alpha,
            beta = beta,
            seed = seed)
    } else if (model == "celda_CG") {
        sce <- .simulateCellsMaincelda_CG(
            model = model,
            S = S,
            CRange = CRange,
            NRange = NRange,
            G = G,
            K = K,
            L = L,
            alpha = alpha,
            beta = beta,
            gamma = gamma,
            delta = delta,
            seed = seed)
    } else if (model == "celda_G") {
        sce <- .simulateCellsMaincelda_G(
            model = model,
            C = C,
            L = L,
            NRange = NRange,
            G = G,
            beta = beta,
            delta = delta,
            gamma = gamma,
            seed = seed)
    } else {
        stop("'model' must be one of 'celda_C', 'celda_G', or 'celda_CG'")
    }

    return(sce)
}


.simulateCellsMaincelda_C <- function(model,
    S,
    CRange,
    NRange,
    G,
    K,
    alpha,
    beta,
    seed) {

    if (is.null(seed)) {
        res <- .simulateCellscelda_C(model = model,
            S = S,
            CRange = CRange,
            NRange = NRange,
            G = G,
            K = K,
            alpha = alpha,
            beta = beta)
    } else {
        res <- with_seed(seed,
            .simulateCellscelda_C(model = model,
                S = S,
                CRange = CRange,
                NRange = NRange,
                G = G,
                K = K,
                alpha = alpha,
                beta = beta))
    }

    sce <- .createSCEsimulateCellsCeldaC(res, seed)

    return(sce)
}


.simulateCellscelda_C <- function(model,
    S = 5,
    CRange = c(50, 100),
    NRange = c(500, 1000),
    G = 100,
    K = 5,
    alpha = 1,
    beta = 1) {

    phi <- .rdirichlet(K, rep(beta, G))
    theta <- .rdirichlet(S, rep(alpha, K))

    ## Select the number of cells per sample
    nC <- sample(seq(CRange[1], CRange[2]), size = S, replace = TRUE)
    cellSampleLabel <- rep(seq(S), nC)

    ## Select state of the cells
    z <- unlist(lapply(seq(S), function(i) {
        sample(seq(K),
            size = nC[i],
            prob = theta[i, ],
            replace = TRUE)
    }))

    ## Select number of transcripts per cell
    nN <- sample(seq(NRange[1], NRange[2]),
        size = length(cellSampleLabel),
        replace = TRUE)

    ## Select transcript distribution for each cell
    cellCounts <- vapply(seq(length(cellSampleLabel)), function(i) {
        stats::rmultinom(1, size = nN[i], prob = phi[z[i], ])
    }, integer(G))

    rownames(cellCounts) <- paste0("Gene_", seq(nrow(cellCounts)))
    colnames(cellCounts) <- paste0("Cell_", seq(ncol(cellCounts)))
    cellSampleLabel <- paste0("Sample_", seq(S))[cellSampleLabel]
    cellSampleLabel <- factor(cellSampleLabel,
        levels = paste0("Sample_", seq(S)))

    ## Peform reordering on final Z and Y assigments:
    cellCounts <- .processCounts(cellCounts)
    names <- list(row = rownames(cellCounts),
        column = colnames(cellCounts),
        sample = unique(cellSampleLabel))
    countChecksum <- .createCountChecksum(cellCounts)
    result <- methods::new("celda_C",
        clusters = list(z = z),
        params = list(K = as.integer(K),
            alpha = alpha,
            beta = beta,
            countChecksum = countChecksum),
        sampleLabel = cellSampleLabel,
        names = names)
    class(result) <- "celda_C"
    result <- .reorderCeldaC(counts = cellCounts, res = result)

    return(list(z = celdaClusters(result)$z,
        counts = cellCounts,
        sampleLabel = cellSampleLabel,
        G = G,
        K = K,
        alpha = alpha,
        beta = beta,
        CRange = CRange,
        NRange = NRange,
        S = S))
}


.createSCEsimulateCellsCeldaC <- function(simList, seed) {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = simList$counts))

    # add metadata
    S4Vectors::metadata(sce)[["celda_simulateCellscelda_C"]] <- list(
        model = "celda_C",
        sampleLevels = as.character(unique(simList$sampleLabel)),
        cellClusterLevels = sort(unique(simList$z)),
        S = simList$S,
        CRange = simList$CRange,
        NRange = simList$NRange,
        G = simList$G,
        K = simList$K,
        alpha = simList$alpha,
        beta = simList$beta,
        seed = seed)

    SummarizedExperiment::rowData(sce)["rownames"] <- rownames(simList$counts)
    SummarizedExperiment::colData(sce)["colnames"] <-
        colnames(simList$counts)
    SummarizedExperiment::colData(sce)["celda_sample_label"] <-
        as.factor(simList$sampleLabel)
    SummarizedExperiment::colData(sce)["celda_cell_cluster"] <-
        as.factor(simList$z)

    return(sce)
}


.simulateCellsMaincelda_CG <- function(model,
    S,
    CRange,
    NRange,
    G,
    K,
    L,
    alpha,
    beta,
    gamma,
    delta,
    seed) {

    if (is.null(seed)) {
        res <- .simulateCellscelda_CG(
            model = model,
            S = S,
            CRange = CRange,
            NRange = NRange,
            G = G,
            K = K,
            L = L,
            alpha = alpha,
            beta = beta,
            gamma = gamma,
            delta = delta)
    } else {
        with_seed(
            seed,
            res <- .simulateCellscelda_CG(
                model = model,
                S = S,
                CRange = CRange,
                NRange = NRange,
                G = G,
                K = K,
                L = L,
                alpha = alpha,
                beta = beta,
                gamma = gamma,
                delta = delta
            )
        )
    }

    sce <- .createSCEsimulateCellsCeldaCG(res, seed)

    return(sce)
}


.simulateCellscelda_CG <- function(model = model,
    S = S,
    CRange = CRange,
    NRange = NRange,
    G = G,
    K = K,
    L = L,
    alpha = alpha,
    beta = beta,
    gamma = gamma,
    delta = delta) {

    ## Number of cells per sample
    nC <- sample(seq(CRange[1], CRange[2]), size = S, replace = TRUE)
    nCSum <- sum(nC)
    cellSampleLabel <- rep(seq(S), nC)

    ## Select number of transcripts per cell
    nN <- sample(seq(NRange[1], NRange[2]),
        size = length(cellSampleLabel),
        replace = TRUE
    )

    ## Generate cell population distribution for each sample
    theta <- t(.rdirichlet(S, rep(alpha, K)))

    ## Assign cells to cellular subpopulations
    z <- unlist(lapply(seq(S), function(i) {
        sample(seq(K),
            size = nC[i],
            prob = theta[, i],
            replace = TRUE
        )
    }))

    ## Generate transcriptional state distribution for each cell subpopulation
    phi <- .rdirichlet(K, rep(beta, L))

    ## Assign genes to gene modules
    eta <- .rdirichlet(1, rep(gamma, L))
    y <- sample(seq(L),
        size = G,
        prob = eta,
        replace = TRUE
    )
    if (length(table(y)) < L) {
        warning(
            "Some gene modules did not receive any genes after sampling.",
            " Try increasing G and/or making gamma larger."
        )
        L <- length(table(y))
        y <- as.integer(as.factor(y))
    }

    psi <- matrix(0, nrow = G, ncol = L)
    for (i in seq(L)) {
        ind <- y == i
        psi[ind, i] <- .rdirichlet(1, rep(delta, sum(ind)))
    }

    ## Select transcript distribution for each cell
    cellCounts <- matrix(0, nrow = G, ncol = nCSum)
    for (i in seq(nCSum)) {
        transcriptionalStateDist <- as.integer(stats::rmultinom(1,
            size = nN[i], prob = phi[z[i], ]
        ))
        for (j in seq(L)) {
            if (transcriptionalStateDist[j] > 0) {
                cellCounts[, i] <- cellCounts[, i] + stats::rmultinom(1,
                    size = transcriptionalStateDist[j], prob = psi[, j]
                )
            }
        }
    }

    ## Ensure that there are no all-0 rows in the counts matrix, which violates
    ## a celda modeling
    ## constraint (columns are guarnteed at least one count):
    zeroRowIdx <- which(rowSums(cellCounts) == 0)
    if (length(zeroRowIdx > 0)) {
        cellCounts <- cellCounts[-zeroRowIdx, ]
        y <- y[-zeroRowIdx]
    }

    ## Assign gene/cell/sample names
    rownames(cellCounts) <- paste0("Gene_", seq(nrow(cellCounts)))
    colnames(cellCounts) <- paste0("Cell_", seq(ncol(cellCounts)))
    cellSampleLabel <- paste0("Sample_", seq(S))[cellSampleLabel]
    cellSampleLabel <- factor(cellSampleLabel,
        levels = paste0("Sample_", seq(S))
    )

    ## Peform reordering on final Z and Y assigments:
    cellCounts <- .processCounts(cellCounts)
    names <- list(
        row = rownames(cellCounts),
        column = colnames(cellCounts),
        sample = unique(cellSampleLabel)
    )
    countChecksum <- .createCountChecksum(cellCounts)
    result <- methods::new("celda_CG",
        clusters = list(z = z, y = y),
        params = list(
            K = as.integer(K),
            L = as.integer(L),
            alpha = alpha,
            beta = beta,
            delta = delta,
            gamma = gamma,
            countChecksum = countChecksum
        ),
        sampleLabel = cellSampleLabel,
        names = names
    )

    result <- .reorderCeldaCG(counts = cellCounts, res = result)

    return(list(
        z = celdaClusters(result)$z,
        y = celdaClusters(result)$y,
        counts = cellCounts,
        sampleLabel = cellSampleLabel,
        G = G,
        K = K,
        L = L,
        CRange = CRange,
        NRange = NRange,
        S = S,
        alpha = alpha,
        beta = beta,
        gamma = gamma,
        delta = delta
    ))
}


.createSCEsimulateCellsCeldaCG <- function(simList, seed) {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = simList$counts))

    # add metadata
    S4Vectors::metadata(sce)[["celda_simulateCellscelda_CG"]] <- list(
        model = "celda_CG",
        sampleLevels = as.character(unique(simList$sampleLabel)),
        cellClusterLevels = sort(unique(simList$z)),
        featureModuleLevels = sort(unique(simList$y)),
        S = simList$S,
        CRange = simList$CRange,
        NRange = simList$NRange,
        G = simList$G,
        K = simList$K,
        L = simList$L,
        alpha = simList$alpha,
        beta = simList$beta,
        gamma = simList$gamma,
        delta = simList$delta,
        seed = seed)

    SummarizedExperiment::rowData(sce)["rownames"] <- rownames(simList$counts)
    SummarizedExperiment::colData(sce)["colnames"] <-
        colnames(simList$counts)
    SummarizedExperiment::colData(sce)["celda_sample_label"] <-
        as.factor(simList$sampleLabel)
    SummarizedExperiment::colData(sce)["celda_cell_cluster"] <-
        as.factor(simList$z)
    SummarizedExperiment::rowData(sce)["celda_feature_module"] <-
        as.factor(simList$y)

    return(sce)
}


.simulateCellsMaincelda_G <- function(model,
    C,
    L,
    NRange,
    G,
    beta,
    delta,
    gamma,
    seed) {

    if (is.null(seed)) {
        res <- .simulateCellscelda_G(
            model = model,
            C = C,
            L = L,
            NRange = NRange,
            G = G,
            beta = beta,
            delta = delta,
            gamma = gamma)
    } else {
        with_seed(
            seed,
            res <- .simulateCellscelda_G(
                model = model,
                C = C,
                L = L,
                NRange = NRange,
                G = G,
                beta = beta,
                delta = delta,
                gamma = gamma)
        )
    }

    sce <- .createSCEsimulateCellsCeldaG(res, seed)

    return(sce)
}


.simulateCellscelda_G <- function(model,
    C = 100,
    L = 10,
    NRange = c(500, 1000),
    G = 100,
    beta = 1,
    gamma = 5,
    delta = 1,
    ...) {

    eta <- .rdirichlet(1, rep(gamma, L))

    y <- sample(seq(L),
        size = G,
        prob = eta,
        replace = TRUE
    )
    if (length(table(y)) < L) {
        stop(
            "Some states did not receive any features after sampling. Try",
            " increasing G and/or setting gamma > 1."
        )
    }

    psi <- matrix(0, nrow = G, ncol = L)
    for (i in seq(L)) {
        ind <- y == i
        psi[ind, i] <- .rdirichlet(1, rep(delta, sum(ind)))
    }

    phi <- .rdirichlet(C, rep(beta, L))

    ## Select number of transcripts per cell
    nN <- sample(seq(NRange[1], NRange[2]), size = C, replace = TRUE)

    ## Select transcript distribution for each cell
    cellCounts <- matrix(0, nrow = G, ncol = C)
    for (i in seq(C)) {
        cellDist <- stats::rmultinom(1, size = nN[i], prob = phi[i, ])
        for (j in seq(L)) {
            cellCounts[, i] <- cellCounts[, i] + stats::rmultinom(1,
                size = cellDist[j], prob = psi[, j]
            )
        }
    }

    ## Ensure that there are no all-0 rows in the counts matrix, which violates
    ## a celda modeling
    ## constraint (columns are guarnteed at least one count):
    zeroRowIdx <- which(rowSums(cellCounts) == 0)
    if (length(zeroRowIdx > 0)) {
        cellCounts <- cellCounts[-zeroRowIdx, ]
        y <- y[-zeroRowIdx]
    }

    rownames(cellCounts) <- paste0("Gene_", seq(nrow(cellCounts)))
    colnames(cellCounts) <- paste0("Cell_", seq(ncol(cellCounts)))

    ## Peform reordering on final Z and Y assigments:
    cellCounts <- .processCounts(cellCounts)
    names <- list(
        row = rownames(cellCounts),
        column = colnames(cellCounts)
    )
    countChecksum <- .createCountChecksum(cellCounts)
    result <- methods::new("celda_G",
        clusters = list(y = y),
        params = list(
            L = as.integer(L),
            beta = beta,
            delta = delta,
            gamma = gamma,
            countChecksum = countChecksum
        ),
        names = names
    )
    result <- .reorderCeldaG(counts = cellCounts, res = result)

    return(list(
        y = celdaClusters(result)$y,
        counts = cellCounts,
        C = C,
        G = G,
        L = L,
        NRange = NRange,
        beta = beta,
        delta = delta,
        gamma = gamma
    ))
}


.createSCEsimulateCellsCeldaG <- function(simList, seed) {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = simList$counts))

    # add metadata
    S4Vectors::metadata(sce)[["celda_simulateCellscelda_G"]] <- list(
        model = "celda_G",
        featureModuleLevels = sort(unique(simList$y)),
        NRange = simList$NRange,
        C = simList$C,
        G = simList$G,
        L = simList$L,
        beta = simList$beta,
        gamma = simList$gamma,
        delta = simList$delta,
        seed = seed)

    SummarizedExperiment::rowData(sce)["rownames"] <- rownames(simList$counts)
    SummarizedExperiment::colData(sce)["colnames"] <-
        colnames(simList$counts)
    SummarizedExperiment::rowData(sce)["celda_feature_module"] <-
        as.factor(simList$y)

    return(sce)
}
