#' @title Simulate contaminated count matrix
#' @description This function generates a list containing two count matrices --
#'  one for real expression, the other one for contamination, as well as other
#'  parameters used in the simulation which can be useful for running
#'  decontamination.
#' @param C Integer. Number of cells to be simulated. Default to be 300.
#' @param G Integer. Number of genes to be simulated. Default to be 100.
#' @param K Integer. Number of cell populations to be simulated. Default to be
#'  3.
#' @param NRange Integer vector. A vector of length 2 that specifies the lower
#'  and upper bounds of the number of counts generated for each cell. Default to
#'  be c(500, 1000).
#' @param beta Numeric. Concentration parameter for Phi. Default to be 0.5.
#' @param delta Numeric or Numeric vector. Concentration parameter for Theta. If
#'  input as a single numeric value, symmetric values for beta distribution are
#'  specified; if input as a vector of lenght 2, the two values will be the
#'  shape1 and shape2 paramters of the beta distribution respectively.
#' @return A list object containing the real expression matrix and contamination
#'  expression matrix as well as other parameters used in the simulation.
#' @examples
#' contaminationSim <- simulateContaminatedMatrix(K = 3, delta = c(1, 9))
#' contaminationSim <- simulateContaminatedMatrix(K = 3, delta = 1)
#' @export
simulateContaminatedMatrix <- function(C = 300,
    G = 100,
    K = 3,
    NRange = c(500, 1000),
    beta = 0.5,
    delta = c(1, 2), 
    seed = 428) {

    set.seed(seed) 

    if (length(delta) == 1) {
        cpByC <- stats::rbeta(n = C, shape1 = delta, shape2 = delta)
    } else {
        cpByC <- stats::rbeta(n = C, shape1 = delta[1], shape2 = delta[2])
    }

    z <- sample(seq(K), size = C, replace = TRUE)
    if (length(unique(z)) < K) {
        warning("Only ",
            length(unique(z)),
            " clusters are simulated. Try to increase numebr of cells 'C' if",
            " more clusters are needed")
        K <- length(unique(z))
        z <- plyr::mapvalues(z, unique(z), seq(length(unique(z))))
    }

    NbyC <- sample(seq(min(NRange), max(NRange)),
        size = C,
        replace = TRUE)
    cNbyC <- vapply(seq(C), function(i) {
        stats::rbinom(n = 1, size = NbyC[i], p = cpByC[i])
    }, integer(1))
    rNbyC <- NbyC - cNbyC

    phi <- .rdirichlet(K, rep(beta, G))

    ## sample real expressed count matrix
    cellRmat <- vapply(seq(C), function(i) {
        stats::rmultinom(1, size = rNbyC[i], prob = phi[z[i], ])
    }, integer(G))

    rownames(cellRmat) <- paste0("Gene_", seq(G))
    colnames(cellRmat) <- paste0("Cell_", seq(C))

    ## sample contamination count matrix
    nGByK <- rowSums(cellRmat) - .colSumByGroup(cellRmat, group = z, K = K)
    eta <- normalizeCounts(counts = nGByK, normalize = "proportion")

    cellCmat <- vapply(seq(C), function(i) {
        stats::rmultinom(1, size = cNbyC[i], prob = eta[, z[i]])
    }, integer(G))
    cellOmat <- cellRmat + cellCmat

    rownames(cellOmat) <- paste0("Gene_", seq(G))
    colnames(cellOmat) <- paste0("Cell_", seq(C))

    return(list("nativeCounts" = cellRmat,
        "observedCounts" = cellOmat,
        "NByC" = NbyC,
        "z" = z,
        "eta" = eta,
        "phi" = t(phi)))
}


# This function calculates the log-likelihood
#
# counts Numeric/Integer matrix. Observed count matrix, rows represent features
# and columns represent cells
# z Integer vector. Cell population labels
# phi Numeric matrix. Rows represent features and columns represent cell
# populations
# eta Numeric matrix. Rows represent features and columns represent cell
# populations
# theta Numeric vector. Proportion of truely expressed transcripts
.deconCalcLL <- function(counts, z, phi, eta, theta) {
    # ll = sum( t(counts) * log( (1-conP )*geneDist[z,] + conP * conDist[z, ] +
    # 1e-20 ) )  # when dist_mat are K x G matrices
    ll <- sum(t(counts) * log(theta * t(phi)[z, ] +
            (1 - theta) * t(eta)[z, ] + 1e-20))
    return(ll)
}

# This function calculates the log-likelihood of background distribution
# decontamination
# bgDist Numeric matrix. Rows represent feature and columns are the times that
# the background-distribution has been replicated.
.bgCalcLL <- function(counts, cellDist, bgDist, theta) {
    ll <- sum(t(counts) * log(theta * t(cellDist) +
            (1 - theta) * t(bgDist) + 1e-20))
    return(ll)
}


# This function updates decontamination
#  phi Numeric matrix. Rows represent features and columns represent cell
# populations
#  eta Numeric matrix. Rows represent features and columns represent cell
# populations
#  theta Numeric vector. Proportion of truely expressed transctripts
#' @importFrom MCMCprecision fit_dirichlet
.cDCalcEMDecontamination <- function(counts,
    phi,
    eta,
    theta,
    z,
    K,
    beta,
    delta) {

    ## Notes: use fix-point iteration to update prior for theta, no need
    ## to feed delta anymore
    logPr <- log(t(phi)[z, ] + 1e-20) + log(theta + 1e-20)
    logPc <- log(t(eta)[z, ] + 1e-20) + log(1 - theta + 1e-20)

    Pr <- exp(logPr) / (exp(logPr) + exp(logPc))
    Pc <- 1 - Pr
    deltaV2 <- MCMCprecision::fit_dirichlet(matrix(c(Pr, Pc), ncol = 2))$alpha

    estRmat <- t(Pr) * counts
    rnGByK <- .colSumByGroupNumeric(estRmat, z, K)
    cnGByK <- rowSums(rnGByK) - rnGByK

    ## Update parameters
    theta <- (colSums(estRmat) + deltaV2[1]) / (colSums(counts) + sum(deltaV2))
    phi <- normalizeCounts(rnGByK,
        normalize = "proportion",
        pseudocountNormalize = beta)
    eta <- normalizeCounts(cnGByK,
        normalize = "proportion",
        pseudocountNormalize = beta)

    return(list("estRmat" = estRmat,
        "theta" = theta,
        "phi" = phi,
        "eta" = eta,
        "delta" = deltaV2))
}


# This function updates decontamination using background distribution
.cDCalcEMbgDecontamination <- function(counts, cellDist, bgDist, theta, beta) {
    # meanNByC <- apply(counts, 2, mean)
    logPr <- log(t(cellDist) + 1e-20) + log(theta + 1e-20) # +
    # log( t(counts) / meanNByC )   # better when without panelty
    logPc <- log(t(bgDist) + 1e-20) + log(1 - theta + 2e-20)

    Pr <- exp(logPr) / (exp(logPr) + exp(logPc))
    Pc <- 1 - Pr
    deltaV2 <- MCMCprecision::fit_dirichlet(matrix(c(Pr, Pc), ncol = 2))$alpha

    estRmat <- t(Pr) * counts

    ## Update paramters
    theta <- (colSums(estRmat) + deltaV2[1]) / (colSums(counts) + sum(deltaV2))
    cellDist <- normalizeCounts(estRmat,
        normalize = "proportion",
        pseudocountNormalize = beta)

    return(list("estRmat" = estRmat,
        "theta" = theta,
        "cellDist" = cellDist,
        "delta" = deltaV2))
}

#' @title Decontaminate count matrix
#' @description This function updates decontamination on dataset with multiple
#'  batches.
#' @param counts Numeric/Integer matrix. Observed count matrix, rows represent
#'  features and columns represent cells.
#' @param z Integer vector. Cell population labels. Default NULL.
#' @param batch Integer vector. Cell batch labels. Default NULL.
#' @param maxIter Integer. Maximum iterations of EM algorithm. Default to be
#'  200.
#' @param beta Numeric. Concentration parameter for Phi. Default to be 1e-6.
#' @param delta Numeric. Symmetric concentration parameter for Theta. Default
#'  to be 10.
#' @param logfile Character. Messages will be redirected to a file named
#'  `logfile`. If NULL, messages will be printed to stdout.  Default NULL.
#' @param verbose Logical. Whether to print log messages. Default TRUE.
#' @return A list object which contains the decontaminated count matrix and
#'  related parameters.
#' @examples
#' data(contaminationSim)
#' deconC <- decontX(counts = contaminationSim$rmat + contaminationSim$cmat,
#'  z = contaminationSim$z, maxIter = 3)
#' deconBg <- decontX(counts = contaminationSim$rmat + contaminationSim$cmat,
#'  maxIter = 3)
#' @export
decontX <- function(counts,
    z = NULL,
    batch = NULL,
    maxIter = 200,
    beta = 1e-6,
    delta = 10,
    logfile = NULL,
    verbose = TRUE) {

    if (!is.null(batch)) {
        ## Set result lists upfront for all cells from different batches
        logLikelihood <- c()
        estRmat <- matrix(NA,
            ncol = ncol(counts),
            nrow = nrow(counts),
            dimnames = list(rownames(counts), colnames(counts)))
        theta <- rep(NA, ncol(counts))
        estConp <- rep(NA, ncol(counts))

        batchIndex <- unique(batch)

        for (bat in batchIndex) {
            countsBat <- counts[, batch == bat]
            if (!is.null(z)) {
                zBat <- z[batch == bat]
            } else {
                zBat <- z
            }
            resBat <- .decontXoneBatch(counts = countsBat,
                z = zBat,
                batch = bat,
                maxIter = maxIter,
                beta = beta,
                delta = delta,
                logfile = logfile,
                verbose = verbose)

            estRmat[, batch == bat] <- resBat$resList$estNativeCounts
            estConp[batch == bat] <- resBat$resList$estConp
            theta[batch == bat] <- resBat$resList$theta

            if (is.null(logLikelihood)) {
                logLikelihood <- resBat$resList$logLikelihood
            } else {
                logLikelihood <- addLogLikelihood(logLikelihood,
                    resBat$resList$logLikelihood)
            }
        }

        runParams <- resBat$runParams
        method <- resBat$method
        resList <- list("logLikelihood" = logLikelihood,
            "estNativeCounts" = estRmat,
            "estConp" = estConp,
            "theta" = theta)

        return(list("runParams" = runParams,
            "resList" = resList,
            "method" = method))
    }

    return(.decontXoneBatch(counts = counts,
        z = z,
        maxIter = maxIter,
        beta = beta,
        delta = delta,
        logfile = logfile,
        verbose = verbose))
}


# This function updates decontamination for one batch
.decontXoneBatch <- function(counts,
    z = NULL,
    batch = NULL,
    maxIter = 200,
    beta = 1e-6,
    delta = 10,
    logfile = NULL,
    verbose = TRUE) {

    .checkCountsDecon(counts)
    .checkParametersDecon(proportionPrior = delta, distributionPrior = beta)

    # nG <- nrow(counts)
    nC <- ncol(counts)
    K <- length(unique(z))

    if (is.null(z)) {
        deconMethod <- "background"
    } else {
        deconMethod <- "clustering"
        z <- .processCellLabels(z, numCells = nC)
    }

    iter <- 1L
    numIterWithoutImprovement <- 0L
    stopIter <- 3L

    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)
    .logMessages("Start DecontX. Decontamination",
        logfile = logfile,
        append = TRUE,
        verbose = verbose)

    if (!is.null(batch)) {
        .logMessages("batch: ",
            batch,
            logfile = logfile,
            append = TRUE,
            verbose = verbose)
    }

    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)
    startTime <- Sys.time()

    if (deconMethod == "clustering") {
        ## Initialization
        deltaInit <- delta
        # theta  = stats::runif(nC, min = 0.1, max = 0.5)
        theta <- stats::rbeta(n = nC,
            shape1 = deltaInit,
            shape2 = deltaInit)
        estRmat <- t(t(counts) * theta)
        phi <- .colSumByGroupNumeric(estRmat, z, K)
        eta <- rowSums(phi) - phi
        phi <- normalizeCounts(phi,
            normalize = "proportion",
            pseudocountNormalize = beta)
        eta <- normalizeCounts(eta,
            normalize = "proportion",
            pseudocountNormalize = beta)
        ll <- c()

        llRound <- .deconCalcLL(counts = counts,
            z = z,
            phi = phi,
            eta = eta,
            theta = theta)

        ## EM updates
        while (iter <= maxIter & numIterWithoutImprovement <= stopIter) {
            nextDecon <- .cDCalcEMDecontamination(counts = counts,
                phi = phi,
                eta = eta,
                theta = theta,
                z = z,
                K = K,
                beta = beta,
                delta = delta)

            theta <- nextDecon$theta
            phi <- nextDecon$phi
            eta <- nextDecon$eta
            delta <- nextDecon$delta

            ## Calculate log-likelihood
            llTemp <- .deconCalcLL(counts = counts,
                z = z,
                phi = phi,
                eta = eta,
                theta = theta)
            ll <- c(ll, llTemp)
            llRound <- c(llRound, round(llTemp, 2))

            if (round(llTemp, 2) > llRound[iter] | iter == 1) {
                numIterWithoutImprovement <- 1L
            } else {
                numIterWithoutImprovement <- numIterWithoutImprovement + 1L
            }
            iter <- iter + 1L
        }
    }

    if (deconMethod == "background") {
        ## Initialization
        deltaInit <- delta
        theta <- stats::rbeta(n = nC,
            shape1 = deltaInit,
            shape2 = deltaInit)
        estRmat <- t(t(counts) * theta)
        bgDist <- rowSums(counts) / sum(counts)
        bgDist <- matrix(rep(bgDist, nC), ncol = nC)
        cellDist <- normalizeCounts(estRmat,
            normalize = "proportion",
            pseudocountNormalize = beta)
        ll <- c()

        llRound <- .bgCalcLL(counts = counts,
            cellDist = cellDist,
            bgDist = bgDist,
            theta = theta)

        ## EM updates
        while (iter <= maxIter & numIterWithoutImprovement <= stopIter) {
            nextDecon <- .cDCalcEMbgDecontamination(counts = counts,
                cellDist = cellDist,
                bgDist = bgDist,
                theta = theta,
                beta = beta)

            theta <- nextDecon$theta
            cellDist <- nextDecon$cellDist
            delta <- nextDecon$delta

            ## Calculate log-likelihood
            llTemp <- .bgCalcLL(counts = counts,
                cellDist = cellDist,
                bgDist = bgDist,
                theta = theta)
            ll <- c(ll, llTemp)
            llRound <- c(llRound, round(llTemp, 2))

            if (round(llTemp, 2) > llRound[iter] | iter == 1) {
                numIterWithoutImprovement <- 1L
            } else {
                numIterWithoutImprovement <- numIterWithoutImprovement + 1L
            }
            iter <- iter + 1L
        }
    }

    resConp <- 1 - colSums(nextDecon$estRmat) / colSums(counts)

    endTime <- Sys.time()
    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)
    .logMessages("Completed DecontX. Total time:",
        format(difftime(endTime, startTime)),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)
    if (!is.null(batch)) {
        .logMessages("batch: ",
            batch,
            logfile = logfile,
            append = TRUE,
            verbose = verbose)
    }
    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)

    runParams <- list("betaInit" = beta,
        "deltaInit" = deltaInit,
        "iteration" = iter - 1L)

    resList <- list("logLikelihood" = ll,
        "estNativeCounts" = nextDecon$estRmat,
        "estConp" = resConp,
        "theta" = theta,
        "delta" = delta)
    # if( deconMethod=="clustering" ) {
    #    posterior.params = list( "est.GeneDist"=phi,  "est.ConDist"=eta  )
    #    resList = append( resList , posterior.params )
    # }

    return(list("runParams" = runParams,
        "resList" = resList,
        "method" = deconMethod))
}


## Make sure provided parameters are the right type and value range
.checkParametersDecon <- function(proportionPrior, distributionPrior) {
    if (length(proportionPrior) > 1 | any(proportionPrior <= 0)) {
        stop("'delta' should be a single positive value.")
    }
    if (length(distributionPrior) > 1 | any(distributionPrior <= 0)) {
        stop("'beta' should be a single positive value.")
    }
}


## Make sure provided count matrix is the right type
.checkCountsDecon <- function(counts) {
    if (sum(is.na(counts)) > 0) {
        stop("Missing value in 'counts' matrix.")
    }
}


## Make sure provided cell labels are the right type
#' @importFrom plyr mapvalues
.processCellLabels <- function(z, numCells) {
    if (length(z) != numCells) {
        stop("'z' must be of the same length as the number of cells in the",
            " 'counts' matrix.")
    }
    if (length(unique(z)) < 2) {
        stop("'z' must have at least 2 different values.") # Even though
        # everything runs smoothly when length(unique(z)) == 1, result is not
        # trustful
    }
    if (!is.factor(z)) {
        z <- plyr::mapvalues(z, unique(z), seq(length(unique(z))))
        z <- as.factor(z)
    }
    return(z)
}


## Add two (veried-length) vectors of logLikelihood
addLogLikelihood <- function(llA, llB) {
    lengthA <- length(llA)
    lengthB <- length(llB)

    if (lengthA >= lengthB) {
        llB <- c(llB, rep(llB[lengthB], lengthA - lengthB))
        ll <- llA + llB
    } else {
        llA <- c(llA, rep(llA[lengthA], lengthB - lengthA))
        ll <- llA + llB
    }

    return(ll)
}
