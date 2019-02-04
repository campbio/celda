#' @title Simulate observed count matrix with contamination
#' @description This function generates a list of two simulated count matrices,
#'  one for real expression and the other one for contamination. Numerous
#'  parameters are specified in the simulation process and can be useful for
#'  running decontamination.
#' @param C Integer. Number of cells to be simulated. Default 300.
#' @param G Integer. Number of genes to be simulated. Default 100.
#' @param K Integer. Number of cell populations to be simulated. Default 3.
#' @param NRange Integer vector. A vector of length 2 that specifies the lower
#'  and upper bounds of the number of counts generated for each cell. Default
#'  \code{c(500, 1000)}.
#' @param beta Numeric. Concentration parameter for Phi. Default 0.5.
#' @param delta Numeric or numeric vector. Concentration parameter for Theta.
#'  If the input is a single numeric value, the symmetric values for beta
#'  distribution are specified. If the input is a vector of length 2, the two
#'  values will be the shape1 and shape2 paramters of the beta distribution
#'  respectively.
#' @param seed Integer. Passed to \code{set.seed()}. Default 12345. If NULL,
#'  no calls to `set.seed()` are made.
#' @examples
#' contaminationSim <- simulateObservedMatrix(K=3, delta = c(1, 9))
#' contaminationSim <- simulateObservedMatrix(K=3, delta = 1)
#' @export
simulateObservedMatrix <- function(
    C = 300,
    G = 100,
    K = 3,
    NRange = c(500, 1000),
    beta = 0.5,
    delta = c(1, 2),
    seed = 12345) {

    if (!is.null(seed)) {
        set.seed(seed)
    }

    if (length(delta) == 1) {
        cpByC <- rbeta(n = C, shape1 = delta, shape2 = delta)
    } else {
        cpByC <- rbeta(n = C, shape1 = delta[1], shape2 = delta[2])
    }

    z <- sample(seq_len(K), size = C, replace = TRUE)
    if (length(unique(z)) < K) {
        warning("Only ",
            length(unique(z)),
            " clusters are simulated. Try increasing the numeber of cells 'C'",
            " if more clusters are needed", "\n")
        K <- length(unique(z))
        z <- plyr::mapvalues(z, unique(z), seq_along(unique(z)))
    }

    NByC <- sample(seq(min(NRange), max(NRange)), size = C, replace = TRUE)
    cNByC <- vapply(seq_len(C), function(i) {
        rbinom(n = 1, size = NByC[i], p = cpByC[i])
    }, integer(1))
    rNByC <- NByC - cNByC

    phi <- rdirichlet(K, rep(beta, G))

    # sample real expressed count matrix
    cellRmat <- vapply(seq_len(C), function(i) {
        stats::rmultinom(1, size = rNByC[i], prob = phi[z[i], ])
    }, integer(G))

    rownames(cellRmat) <- paste0("Gene_", seq_len(G))
    colnames(cellRmat) <- paste0("Cell_", seq_len(C))


    # sample contamination count matrix
    nGByK <- rowSums(cellRmat) - colSumByGroup(cellRmat, group = z, K = K)
    eta <- normalizeCounts(counts = nGByK, normalize = "proportion")

    cellCmat <- vapply(seq_len(C), function(i) {
        stats::rmultinom(1, size = cNByC[i], prob = eta[, z[i]])
    }, integer(G))
    rownames(cellCmat) <- paste0("Gene_", seq_len(G))
    colnames(cellCmat) <- paste0("Cell_", seq_len(C))

    return(list("rmat" = cellRmat, "cmat" = cellCmat,
        "N.by.C" = NByC,
        "z" = z,
        "eta" = eta, "phi" = t(phi)))
}


# This function calculates the log-likelihood
#
# omat Numeric/Integer matrix. Observed count matrix, rows represent features
#  and columns represent cells
# z Integer vector. Cell population labels
# phi Numeric matrix. Rows represent features and columns represent cell
#  populations
# eta Numeric matrix. Rows represent features and columns represent cell
#  populations
# theta Numeric vector. Proportion of truely expressed transcripts
decon.calcLL <- function(omat, z,  phi, eta, theta) {
    #ll = sum( t(omat) * log( (1-conP )*geneDist[z,] + conP * conDist[z, ] +
    # 1e-20 ) )  # when dist_mat are K x G matrices
    ll <- sum(t(omat) * log(theta * t(phi)[z, ] + (1 - theta) * t(eta)[z, ] +
            1e-20))
    return(ll)
}


# This function calculates the log-likelihood of background distribution
#  decontamination
#
# bgDist Numeric matrix. Rows represent feature and columns are the times that
#  the background-distribution has been replicated.
bg.calcLL <- function(omat, cellDist, bgDist, theta) {
    ll <- sum(t(omat) * log(theta * t(cellDist) + (1 - theta) * t(bgDist) +
            1e-20))
    return(ll)
}


# This function updates decontamination
#
#  phi Numeric matrix. Rows represent features and columns represent cell
#   populations
#  eta Numeric matrix. Rows represent features and columns represent cell
#   populations
#  theta Numeric vector. Proportion of truely expressed transctripts
cD.calcEMDecontamination <- function(
    omat,
    phi,
    eta,
    theta,
    z,
    K,
    beta,
    delta) {

    logPr <- log(t(phi)[z, ] + 1e-20) + log(theta + 1e-20)
    logPc <- log(t(eta)[z, ] + 1e-20) + log(1 - theta + 1e-20)

    Pr <- exp(logPr) / (exp(logPr) + exp(logPc))

    estRmat <- t(Pr) * omat
    rnGByK <- colSumByGroup.numeric(estRmat, z, K)
    cnGByK <- rowSums(rnGByK) - rnGByK

    #update parameters
    theta <- (colSums(estRmat) + delta) / (colSums(omat) + 2 * delta)
    phi <- normalizeCounts(rnGByK, normalize = "proportion",
        pseudocount.normalize = beta)
    eta <- normalizeCounts(cnGByK, normalize = "proportion",
        pseudocount.normalize = beta)

    return(list("est.rmat" = estRmat, "theta" = theta, "phi" = phi,
        "eta" = eta))
}


# This function updates decontamination using background distribution
#
cD.calcEMbgDecontamination <- function(
    omat,
    cellDist,
    bgDist,
    theta,
    beta,
    delta) {

    #meanNByC <- apply(omat, 2, mean)

    logPr <- log(t(cellDist) + 1e-20) + log(theta + 1e-20)
    # + log( t(omat) / meanNByC )   # better when without panelty

    logPc <- log(t(bgDist) + 1e-20) + log(1 - theta + 2e-20)

    Pr <- exp(logPr) / (exp(logPr) + exp(logPc))

    estRmat <- t(Pr) * omat

    # Update paramters
    theta <- (colSums(estRmat) + delta) / (colSums(omat) + 2 * delta)
    cellDist <- normalizeCounts(estRmat,
        normalize = "proportion",
        pseudocount.normalize = beta)

    return(list("est.rmat" = estRmat, "theta" = theta, "cellDist" = cellDist))
}


#' @title Decontaminate the observed count matrix
#' @description This function calculates the decontaminated count matrix from
#'  the observed count matrix. The log-likelihood and the contamination
#'  proportion are also reported.
#' @param omat Numeric/Integer matrix. Observed count matrix, rows represent
#'  features and columns represent cells.
#' @param z Integer vector. Cell population labels. Default \code{NULL}.
#' @param max.iter Integer. Maximum iterations of EM algorithm. Default 200.
#' @param beta Numeric. Concentration parameter for Phi. Default 1e-6.
#' @param delta Numeric. Symmetric concentration parameter for Theta. Default
#'  10.
#' @param logfile Character. Messages will be redirected to a file named
#'  `logfile`. If NULL, messages will be printed to standard output. Default
#'  \code{NULL}.
#' @param verbose Logical. Whether to print log messages. Default \code{TRUE}.
#' @param seed Integer. Passed to \code{set.seed()}. Default 1234567. If NULL,
#'  no calls to \code{set.seed()} are made.
#' @examples
#' decon.c <- DecontX(omat = contamination.sim$rmat + contamination.sim$cmat,
#'     z = contamination.sim$z, max.iter = 3)
#' decon.bg <- DecontX(omat = contamination.sim$rmat + contamination.sim$cmat,
#'     max.iter = 3)
#' @export
DecontX <- function(
    omat,
    z = NULL,
    max.iter = 200,
    beta = 1e-6,
    delta = 10,
    logfile = NULL,
    verbose = TRUE,
    seed = 1234567) {

    checkCounts.decon(omat)
    checkParameters.decon(proportion.prior = delta, distribution.prior = beta)

    #nG <- nrow(omat)
    nC <- ncol(omat)
    K <- length(unique(z))

    if (is.null(z)) {
        decon.method <- "background"
    } else {
        decon.method <- "clustering"
        z <- processCellLabels(z, num.cells = nC)
    }

    logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)
    logMessages("Start DecontX. Decontamination",
        logfile = logfile, append = TRUE, verbose = verbose)
    logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)
    start.time <- Sys.time()

    if (decon.method == "clustering") {

        # initialization
        if (!is.null(seed)) {
            set.seed(seed)
        }
        theta <- runif(nC, min = 0.1, max = 0.5)
        estRmat <- t (t(omat) * theta)
        phi <- colSumByGroup.numeric(estRmat, z, K)
        eta <- rowSums(phi) - phi
        phi <- normalizeCounts(phi, normalize = "proportion",
            pseudocount.normalize = beta)
        eta <- normalizeCounts(eta, normalize = "proportion",
            pseudocount.normalize = beta)
        ll <- c()

        # EM updates
        for (iteration in seq_len(max.iter)) {

            next.decon <- cD.calcEMDecontamination(omat = omat,
                phi = phi,
                eta = eta,
                theta = theta,
                z = z,
                K = K,
                beta = beta,
                delta = delta)

            theta <- next.decon$theta
            phi <- next.decon$phi
            eta <- next.decon$eta

            # Calculate log-likelihood
            ll.temp <- decon.calcLL(omat = omat,
                z = z,
                phi = phi,
                eta = eta,
                theta = theta)
            ll <- c(ll, ll.temp)
        }
    }

    if (decon.method == "background") {

        # Initialization
        if (!is.null(seed)) {
            set.seed(seed)
        }
        theta <- runif(nC, min = 0.1, max = 0.5)
        estRmat <- t(t(omat) * theta)
        bgDist <- rowSums(omat) / sum(omat)
        bgDist <- matrix(rep(bgDist, nC), ncol = nC)
        cellDist <- normalizeCounts(estRmat, normalize = "proportion",
            pseudocount.normalize = beta)
        ll <- c()

        # EM updates
        for (iteration in seq_len(max.iter)) {

            next.decon <- cD.calcEMbgDecontamination(omat = omat,
                cellDist = cellDist,
                bgDist = bgDist,
                theta = theta,
                beta = beta,
                delta = delta)

            theta <- next.decon$theta
            cellDist <- next.decon$cellDist

            # Calculate log-likelihood
            ll.temp <- bg.calcLL(omat = omat,
                cellDist = cellDist,
                bgDist = bgDist,
                theta = theta)
            ll  <- c(ll, ll.temp)
        }
    }

    resConp <- 1 - colSums(next.decon$estRmat) / colSums(omat)

    end.time <- Sys.time()
    logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)
    logMessages("Completed DecontX. Total time:",
        format(difftime(end.time, start.time)),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)
    logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)

    run.params <- list("beta" = beta, "delta" = delta, "iteration" = max.iter)

    res.list <- list("logLikelihood" = ll,
        "est.rmat" = next.decon$estRmat,
        "est.conp" = resConp,
        "theta" = theta)
    if (decon.method == "clustering") {
        posterior.params <- list("est.GeneDist" = phi, "est.ConDist" = eta)
        res.list <- append(res.list, posterior.params)
    }

    return(list("run.params" = run.params,
        "res.list" = res.list,
        "method" = decon.method))
}


### checking
# Make sure provided parameters are the right type and value range
checkParameters.decon <- function(proportion.prior, distribution.prior) {
    if (length(proportion.prior) > 1 || proportion.prior <= 0) {
        stop("'delta' should be a single positive value.")
    }
    if (length( distribution.prior) > 1 || distribution.prior <= 0) {
        stop("'beta' should be a single positive value.")
    }
}

# Make sure provided rmat is the right type
checkCounts.decon <- function(omat) {
    if (sum(is.na(omat)) > 0) {
        stop("Missing value in 'omat' matrix.")
    }
}

# Make sure provided cell labels are the right type
processCellLabels <- function(z, num.cells) {
    if (length(z) != num.cells) {
        stop("'z' must be of the same length as the number of cells in the",
            " 'counts' matrix.")
    }
    if (length(unique(z)) < 2) {
        stop("'z' must have at least 2 different values.")
        # Even though everything runs smoothly when length(unique(z)) == 1,
        #  result is not trustful
    }
    if (!is.factor(z)) {
        z <- plyr::mapvalues(z, unique(z), seq_along(unique(z)))
        z <- as.factor(z)
    }
    return(z)
}
