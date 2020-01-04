

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
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
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
    seed = 12345) {

    if (is.null(seed)) {
        res <- .simulateContaminatedMatrix(C = C,
            G = G,
            K = K,
            NRange = NRange,
            beta = beta,
            delta = delta)
    } else {
        with_seed(seed,
            res <- .simulateContaminatedMatrix(C = C,
                G = G,
                K = K,
                NRange = NRange,
                beta = beta,
                delta = delta))
    }

    return(res)
}


.simulateContaminatedMatrix <- function(C = 300,
    G = 100,
    K = 3,
    NRange = c(500, 1000),
    beta = 0.5,
    delta = c(1, 2)) {
    if (length(delta) == 1) {
        cpByC <- stats::rbeta(n = C,
            shape1 = delta,
            shape2 = delta)
    } else {
        cpByC <- stats::rbeta(n = C,
            shape1 = delta[1],
            shape2 = delta[2])
    }

    z <- sample(seq(K), size = C, replace = TRUE)
    if (length(unique(z)) < K) {
        warning(
            "Only ",
            length(unique(z)),
            " clusters are simulated. Try to increase numebr of cells 'C' if",
            " more clusters are needed"
        )
        K <- length(unique(z))
        z <- plyr::mapvalues(z, unique(z), seq(length(unique(z))))
    }

    NbyC <- sample(seq(min(NRange), max(NRange)),
        size = C,
        replace = TRUE)
    cNbyC <- vapply(seq(C), function(i) {
        stats::rbinom(n = 1,
            size = NbyC[i],
            p = cpByC[i])
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
    nGByK <-
        rowSums(cellRmat) - .colSumByGroup(cellRmat, group = z, K = K)
    eta <- normalizeCounts(counts = nGByK, normalize = "proportion")

    cellCmat <- vapply(seq(C), function(i) {
        stats::rmultinom(1, size = cNbyC[i], prob = eta[, z[i]])
    }, integer(G))
    cellOmat <- cellRmat + cellCmat

    rownames(cellOmat) <- paste0("Gene_", seq(G))
    colnames(cellOmat) <- paste0("Cell_", seq(C))

    return(
        list(
            "nativeCounts" = cellRmat,
            "observedCounts" = cellOmat,
            "NByC" = NbyC,
            "z" = z,
            "eta" = eta,
            "phi" = t(phi)
        )
    )
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
    ll <- sum(Matrix::t(counts) * log(theta * t(phi)[z, ] +
            (1 - theta) * t(eta)[z, ] + 1e-20))
    return(ll)
}

# DEPRECATED. This is not used, but is kept as it might be useful in the future.
# This function calculates the log-likelihood of background distribution
# decontamination
# bgDist Numeric matrix. Rows represent feature and columns are the times that
# the background-distribution has been replicated.
.bgCalcLL <- function(counts, globalZ, cbZ, phi, eta, theta) {
    # ll <- sum(t(counts) * log(theta * t(cellDist) +
    #        (1 - theta) * t(bgDist) + 1e-20))
    ll <- sum(t(counts) * log(theta * t(phi)[cbZ, ] +
            (1 - theta) * t(eta)[globalZ, ] + 1e-20))
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
    delta) {
    ## Notes: use fix-point iteration to update prior for theta, no need
    ## to feed delta anymore

    logPr <- log(t(phi)[z, ] + 1e-20) + log(theta + 1e-20)
    logPc <- log(t(eta)[z, ] + 1e-20) + log(1 - theta + 1e-20)
    Pr.e <- exp(logPr)
    Pc.e <- exp(logPc)
    Pr <- Pr.e / (Pr.e + Pc.e)

    estRmat <- t(Pr) * counts
    rnGByK <- .colSumByGroupNumeric(estRmat, z, K)
    cnGByK <- rowSums(rnGByK) - rnGByK

    counts.cs <- colSums(counts)
    estRmat.cs <- colSums(estRmat)
    estRmat.cs.n <- estRmat.cs / counts.cs
    estCmat.cs.n <- 1 - estRmat.cs.n
    temp <- cbind(estRmat.cs.n, estCmat.cs.n)
    deltaV2 <- MCMCprecision::fit_dirichlet(temp)$alpha

    ## Update parameters
    theta <-
        (estRmat.cs + deltaV2[1]) / (counts.cs + sum(deltaV2))
    phi <- normalizeCounts(rnGByK,
        normalize = "proportion",
        pseudocountNormalize = 1e-20)
    eta <- normalizeCounts(cnGByK,
        normalize = "proportion",
        pseudocountNormalize = 1e-20)

    return(list(
        "estRmat" = estRmat,
        "theta" = theta,
        "phi" = phi,
        "eta" = eta,
        "delta" = deltaV2
    ))
}

# DEPRECATED. This is not used, but is kept as it might be useful in the
# feature.
# This function updates decontamination using background distribution
.cDCalcEMbgDecontamination <-
    function(counts, globalZ, cbZ, trZ, phi, eta, theta) {
        logPr <- log(t(phi)[cbZ, ] + 1e-20) + log(theta + 1e-20)
        logPc <-
            log(t(eta)[globalZ, ] + 1e-20) + log(1 - theta + 1e-20)

        Pr <- exp(logPr) / (exp(logPr) + exp(logPc))
        Pc <- 1 - Pr
        deltaV2 <-
            MCMCprecision::fit_dirichlet(matrix(c(Pr, Pc), ncol = 2))$alpha

        estRmat <- t(Pr) * counts
        phiUnnormalized <-
            .colSumByGroupNumeric(estRmat, cbZ, max(cbZ))
        etaUnnormalized <-
            rowSums(phiUnnormalized) - .colSumByGroupNumeric(phiUnnormalized,
                trZ, max(trZ))

        ## Update paramters
        theta <-
            (colSums(estRmat) + deltaV2[1]) / (colSums(counts) + sum(deltaV2))
        phi <-
            normalizeCounts(phiUnnormalized,
                normalize = "proportion",
                pseudocountNormalize = 1e-20)
        eta <-
            normalizeCounts(etaUnnormalized,
                normalize = "proportion",
                pseudocountNormalize = 1e-20)

        return(list(
            "estRmat" = estRmat,
            "theta" = theta,
            "phi" = phi,
            "eta" = eta,
            "delta" = deltaV2
        ))
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
#' @param delta Numeric. Symmetric concentration parameter for Theta. Default
#'  to be 10.
#' @param logfile Character. Messages will be redirected to a file named
#'  `logfile`. If NULL, messages will be printed to stdout.  Default NULL.
#' @param verbose Logical. Whether to print log messages. Default TRUE.
#' @param varGenes Positive Integer. Used only when z is not provided.
#' Need to be larger than 1. Default value is 5000 if not provided.
#' varGenes, being the number of most variable genes, is used to filter genes
#' based on the variability of gene's expression cross cells. While the
#' variability is calcualted using scran::trendVar() and scran::decomposeVar().
#' @param L Positive Integer. Used only when z is not provided.
#' Need to be larger than 1. Default value is 50 if not provided.
#' L, being the number of gene modules, is used on celda_CG clustering
#' to collapse genes into gene modules.
#' @param dbscanEps Numeric. Used only when z is not provided.
#' Need to be non-negative. Default is 1.0 if not provided.
#' dbscanEps is the clustering resolution parameter that is used to feed into
#' dbscan::dbscan() to estimate broad cell clusters.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @return A list object which contains the decontaminated count matrix and
#'  related parameters.
#' @examples
#' data(contaminationSim)
#' deconC <- decontX(
#'   counts = contaminationSim$rmat + contaminationSim$cmat,
#'   z = contaminationSim$z, maxIter = 3
#' )
#' deconBg <- decontX(
#'   counts = contaminationSim$rmat + contaminationSim$cmat,
#'   maxIter = 3
#' )
#' @export
decontX <- function(counts,
    z = NULL,
    batch = NULL,
    maxIter = 200,
    delta = 10,
    convergence = 0.001,
    logfile = NULL,
    verbose = TRUE,
    varGenes = NULL,
    L = NULL,
    dbscanEps = NULL,
    seed = 12345) {

    if (is.null(seed)) {
        res <- .decontX(counts = counts,
            z = z,
            batch = batch,
            maxIter = maxIter,
            delta = delta,
            convergence = convergence,
            logfile = logfile,
            verbose = verbose,
            varGenes = varGenes,
            L = L,
            dbscanEps = dbscanEps)
    } else {
        with_seed(seed,
            res <- .decontX(counts = counts,
                z = z,
                batch = batch,
                maxIter = maxIter,
                delta = delta,
                convergence = convergence,
                logfile = logfile,
                verbose = verbose,
                varGenes = varGenes,
                L = L,
                dbscanEps = dbscanEps))
    }

    return(res)
}


.decontX <- function(counts,
    z = NULL,
    batch = NULL,
    maxIter = 200,
    delta = 10,
    convergence = 0.001,
    logfile = NULL,
    verbose = TRUE,
    varGenes = NULL,
    dbscanEps = NULL,
    L = NULL) {

    startTime <- Sys.time()
    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)
    .logMessages("Starting DecontX",
        logfile = logfile,
        append = TRUE,
        verbose = verbose)
    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)

    # Convert to sparse matrix
    # After Celda can run on sparse matrix,
    # then we can just have this be required
    # as input
    counts <- as(counts, "dgCMatrix")

    ## Empty expression genes won't be used for estimation
    haveEmptyGenes <- FALSE
    totalGenes <- nrow(counts)
    totalCells <- ncol(counts)
    noneEmptyGeneIndex <- Matrix::rowSums(counts) != 0
    geneNames <- rownames(counts)
    if (sum(noneEmptyGeneIndex) != totalGenes) {
        counts <- counts[noneEmptyGeneIndex, ]
        haveEmptyGenes <- TRUE
    }

    nC <- ncol(counts)
    allCellNames <- colnames(counts)

    estRmat <- Matrix::Matrix(
      data = 0,
      ncol = totalCells,
      nrow = totalGenes,
      sparse = TRUE,
      dimnames = list(geneNames, allCellNames)
    )

    if (!is.null(batch)) {
        ## Set result lists upfront for all cells from different batches
        logLikelihood <- c()

        theta <- rep(NA, nC)
        estConp <- rep(NA, nC)
        returnZ <- rep(NA, nC)

        batchIndex <- unique(batch)

        for (bat in batchIndex) {
            .logMessages(
                date(),
                ".. Analyzing cells in batch",
                bat,
                logfile = logfile,
                append = TRUE,
                verbose = verbose
            )

            zBat <- NULL
            countsBat <- counts[, batch == bat]
            if (!is.null(z)) {
                zBat <- z[batch == bat]
            }
            resBat <- .decontXoneBatch(
                counts = countsBat,
                z = zBat,
                batch = bat,
                maxIter = maxIter,
                delta = delta,
                convergence = convergence,
                logfile = logfile,
                verbose = verbose,
                varGenes = varGenes,
                dbscanEps = dbscanEps,
                L = L
            )

            estRmat <- calculateNativeMatrix(
              counts = countsBat,
              native_counts = estRmat,
              theta = resBat$resList$theta,
              eta = resBat$resList$eta,
              row_index = which(noneEmptyGeneIndex),
              col_index = which(batch == bat),
              phi = resBat$resList$phi,
              z = as.integer(resBat$runParams$z),
              pseudocount = 1e-20)

            estConp[batch == bat] <- resBat$resList$estConp
            theta[batch == bat] <- resBat$resList$theta
            returnZ[batch == bat] <- resBat$runParams$z

            if (is.null(logLikelihood)) {
                logLikelihood <- resBat$resList$logLikelihood
            } else {
                logLikelihood <- addLogLikelihood(logLikelihood,
                    resBat$resList$logLikelihood)
            }
        }

        runParams <- resBat$runParams
        ## All batches share the same other parameters except cluster label z
        ## So update z in the final returned result
        runParams$z <- returnZ
        method <- resBat$method
        resList <- list(
            "logLikelihood" = logLikelihood,
            "estNativeCounts" = estRmat,
            "estConp" = estConp,
            "theta" = theta
        )

        returnResult <- list(
            "runParams" = runParams,
            "resList" = resList,
            "method" = method
        )
    } else { ## When there is only one batch
        returnResult <- .decontXoneBatch(
            counts = counts,
            z = z,
            maxIter = maxIter,
            delta = delta,
            convergence = convergence,
            logfile = logfile,
            verbose = verbose,
            varGenes = varGenes,
            dbscanEps = dbscanEps,
            L = L
        )

        estRmat <- calculateNativeMatrix(
          counts = counts,
          native_counts = estRmat,
          theta = returnResult$resList$theta,
          eta = returnResult$resList$eta,
          row_index = which(noneEmptyGeneIndex),
          col_index = seq(totalCells),
          phi = returnResult$resList$phi,
          z = as.integer(returnResult$runParams$z),
          pseudocount = 1e-20)
          returnResult$resList$estNativeCounts <- estRmat
    }

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
    .logMessages(paste(rep("-", 50), collapse = ""),
        logfile = logfile,
        append = TRUE,
        verbose = verbose)

    return(returnResult)

}


# This function updates decontamination for one batch
.decontXoneBatch <- function(counts,
    z = NULL,
    batch = NULL,
    maxIter = 200,
    delta = 10,
    convergence = 0.01,
    logfile = NULL,
    verbose = TRUE,
    varGenes = NULL,
    dbscanEps = NULL,
    L = NULL) {
    .checkCountsDecon(counts)
    .checkParametersDecon(proportionPrior = delta)

    # nG <- nrow(counts)
    nC <- ncol(counts)
    deconMethod <- "clustering"

    umap <- NULL
    if (is.null(z)) {
        .logMessages(
            date(),
            ".. Estimating cell types with Celda",
            logfile = logfile,
            append = TRUE,
            verbose = verbose
        )
        ## Always uses clusters for DecontX estimation
        #deconMethod <- "background"

        varGenes <- .processvarGenes(varGenes)
        dbscanEps <- .processdbscanEps(dbscanEps)
        L <- .processL(L)

        celda.init <- .decontxInitializeZ(object = counts,
            varGenes = varGenes,
            dbscanEps = dbscanEps,
            verbose = verbose,
            logfile = logfile)
        z <- celda.init$z
        umap <- celda.init$umap
        colnames(umap) <- c("DecontX_UMAP_1",
                            "DecontX_UMAP_2") 
        rownames(umap) <- colnames(counts)
    }

    z <- .processCellLabels(z, numCells = nC)
    K <- length(unique(z))

    iter <- 1L
    numIterWithoutImprovement <- 0L
    stopIter <- 3L

    .logMessages(
        date(),
        ".. Estimating contamination",
        logfile = logfile,
        append = TRUE,
        verbose = verbose
    )

    if (deconMethod == "clustering") {
        ## Initialization
        deltaInit <- delta
        # theta  = stats::runif(nC, min = 0.1, max = 0.5)
        theta <- stats::rbeta(n = nC,
            shape1 = deltaInit,
            shape2 = deltaInit)

        nextDecon <- decontXInitialize(
                       counts = counts,
                       theta = theta,
                       z = z,
                       pseudocount = 1e-20)
        phi <- nextDecon$phi
        eta <- nextDecon$eta

#        estRmat <- Matrix::t(Matrix::t(counts) * theta)
#        phi <- .colSumByGroupNumeric(as.matrix(estRmat), z, K)
#        eta <- rowSums(phi) - phi
#        phi <- normalizeCounts(phi,
#            normalize = "proportion",
#            pseudocountNormalize = 1e-20)
#        eta <- normalizeCounts(eta,
#            normalize = "proportion",
#            pseudocountNormalize = 1e-20)
        ll <- c()


#        llRound <- .deconCalcLL(
#            counts = counts,
#            z = z,
#            phi = phi,
#            eta = eta,
#            theta = theta
#        )
        llRound <- decontXLogLik(
            counts = counts,
            z = z,
            phi = phi,
            eta = eta,
            theta = theta,
            pseudocount = 1e-20)

        ## EM updates
        theta.previous <- theta
        converged <- FALSE
        counts.colsums <- Matrix::colSums(counts)
        while (iter <= maxIter & !isTRUE(converged) &
                numIterWithoutImprovement <= stopIter) {
#            nextDecon <- .cDCalcEMDecontamination(
#                counts = counts,
#                phi = phi,
#                eta = eta,
#                theta = theta,
#                z = z,
#                K = K,
#                delta = delta
#            )
            nextDecon <- decontXEM(counts = counts,
              counts_colsums = counts.colsums,
              phi = phi,
              eta = eta,
              theta = theta,
              z = z,
              pseudocount = 1e-20)

            theta <- nextDecon$theta
            phi <- nextDecon$phi
            eta <- nextDecon$eta
            delta <- nextDecon$delta

            ## Calculate log-likelihood
#            llTemp <- .deconCalcLL(
#                counts = counts,
#                z = z,
#                phi = phi,
#                eta = eta,
#                theta = theta
#            )
            llTemp <- decontXLogLik(
                counts = counts,
                z = z,
                phi = phi,
                eta = eta,
                theta = theta,
                pseudocount = 1e-20)

            ll <- c(ll, llTemp)
            llRound <- c(llRound, round(llTemp, 2))

            if (round(llTemp, 2) > llRound[iter] | iter == 1) {
                numIterWithoutImprovement <- 1L
            } else {
                numIterWithoutImprovement <- numIterWithoutImprovement + 1L
            }

            max.divergence <- max(abs(theta.previous - theta))
            if (max.divergence < convergence) {
              converged <- TRUE
            }
            theta.previous <- theta

            .logMessages(date(),
                ".... Completed iteration:",
                iter,
                "| converge:",
                signif(max.divergence, 4),
                logfile = logfile,
                append = TRUE,
                verbose = verbose)

            iter <- iter + 1L
        }
    }


#    resConp <- 1 - colSums(nextDecon$estRmat) / colSums(counts)
    resConp <- nextDecon$contamination
    names(resConp) <- colnames(counts)

    if (!is.null(batch)) {
        batchMessage <- paste(" ", "in batch ", batch, ".", sep = "")
    } else {
        batchMessage <- "."
    }

    runParams <- list("deltaInit" = deltaInit,
        "iteration" = iter - 1L,
        "z" = z)
    if(!is.null(umap)) {
      runParams[["UMAP"]] <- umap
    }
    
    resList <- list(
        "logLikelihood" = ll,
        "estNativeCounts" = nextDecon$estRmat,
        "estConp" = resConp,
        "theta" = theta,
        "delta" = delta,
        "phi" = phi,
        "eta" = eta
    )

    return(list(
        "runParams" = runParams,
        "resList" = resList
    ))
}


## Make sure provided parameters are the right type and value range
.checkParametersDecon <- function(proportionPrior) {
    if (length(proportionPrior) > 1 | any(proportionPrior <= 0)) {
        stop("'delta' should be a single positive value.")
    }
}


## Make sure provided count matrix is the right type
.checkCountsDecon <- function(counts) {
    if (sum(is.na(counts)) > 0) {
        stop("Missing value in 'counts' matrix.")
    }
    if (is.null(dim(counts))) {
        stop("At least 2 genes need to have non-zero expressions.")
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
        stop("No need to decontaminate when only one cluster",
            " is in the dataset.") # Even though
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



## Initialization of cell labels for DecontX when they are not given
.decontxInitializeZ <-
    function(object, # object is either a sce object or a count matrix
        varGenes = 5000,
        L = 50,
        dbscanEps = 1.0,
        verbose = TRUE,
        logfile = NULL) {

        if (!is(object, "SingleCellExperiment")) {
            sce <- SingleCellExperiment::SingleCellExperiment(assays =
                    list(counts = object))
        }

        ## Add the log2 normalized counts into sce object
        ## The normalized counts is also centered using library size in the
        ## original count matrix in scater::normalizeSCE()
        #sce <- suppressWarnings(scater::normalizeSCE(sce))
        sce <- scater::logNormCounts(sce, log = TRUE)

        if (nrow(sce) <= varGenes) {
             topVariableGenes <- seq_len(nrow(sce))
        } else if (nrow(sce) > varGenes) {
        ## Use the top most variable genes to do rough clustering
        ## (celda_CG & Louvian graph algorithm)
            #mvTrend <- scran::trendVar(sce, use.spikes = FALSE)
            #decomposeTrend <- scran::decomposeVar(sce, mvTrend)
            #topVariableGenes <- order(decomposeTrend$bio,
            #    decreasing = TRUE)[seq(varGenes)]
            
            sce.var <- scran::modelGeneVar(sce)
            topVariableGenes <- order(sce.var$bio,
                 decreasing = TRUE)[seq(varGenes)]
        }
        countsFiltered <- as.matrix(SingleCellExperiment::counts(
            sce[topVariableGenes, ]))
        storage.mode(countsFiltered) <- "integer"

        .logMessages(
            date(),
            ".... Collapsing features into",
            L,
            "modules",
            logfile = logfile,
            append = TRUE,
            verbose = verbose
        )
        ## Celda clustering using recursive module splitting
        if (L < nrow(countsFiltered)) {
            initialModuleSplit <- recursiveSplitModule(countsFiltered,
                initialL = L, maxL = L, perplexity = FALSE, verbose = FALSE)
            initialModel <- subsetCeldaList(initialModuleSplit, list(L = L))
            fm <- factorizeMatrix(countsFiltered, initialModel, type = "counts")
            fm <- fm$counts$cell
            rm(initialModuleSplit)
            rm(initialModel)
        } else {
            fm <- countsFiltered
        }

        .logMessages(
            date(),
            ".... Reducing dimensionality with UMAP",
            logfile = logfile,
            append = TRUE,
            verbose = verbose
        )
        ## Louvan graph-based method to reduce dimension into 2 cluster
        nNeighbors <- min(15, ncol(countsFiltered))
        resUmap <- uwot::umap(t(sqrt(fm)), n_neighbors = nNeighbors,
            min_dist = 0.01, spread = 1)
        rm(fm)

        .logMessages(
            date(),
            " .... Determining cell clusters with DBSCAN (Eps=",
            dbscanEps,
            ")",
            sep = "",
            logfile = logfile,
            append = TRUE,
            verbose = verbose
        )
        # Use dbSCAN on the UMAP to identify broad cell types
        totalClusters <- 1
        while (totalClusters <= 1 & dbscanEps > 0) {
          resDbscan <- dbscan::dbscan(resUmap, dbscanEps)
          dbscanEps <- dbscanEps - (0.25 * dbscanEps)
          totalClusters <- length(unique(resDbscan$cluster))
        }

        return(list("z" = resDbscan$cluster,
                    "umap" = resUmap))
    }


## process varGenes
.processvarGenes <- function(varGenes) {
    if (is.null(varGenes)) {
        varGenes <- 5000
    } else {
        if (varGenes < 2 | !is.integer(varGenes)) {
            stop("Parameter 'varGenes' must be an integer and larger than 1.")
        }
    }
    return(varGenes)
}

## process dbscanEps for resolusion threshold using DBSCAN
.processdbscanEps <- function(dbscanEps) {
    if (is.null(dbscanEps)) {
        dbscanEps <- 1
    } else {
        if (dbscanEps < 0) {
            stop("Parameter 'dbscanEps' needs to be non-negative.")
        }
    }
    return(dbscanEps)
}

## process gene modules L
.processL <- function(L) {
    if (is.null(L)) {
        L <- 50
    } else {
        if (L < 2 | !is.integer(L)) {
            stop("Parameter 'L' must be an integer and larger than 1.")
        }
    }
    return(L)
}
