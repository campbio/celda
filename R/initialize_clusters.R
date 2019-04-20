.initializeCluster <- function(N,
    len,
    z = NULL,
    initial = NULL,
    fixed = NULL) {

    # If initial values are given, then they will not be randomly initialized
    if (!is.null(initial)) {
        initValues <- sort(unique(initial))
        if (length(unique(initial)) != N || length(initial) != len ||
            !all(initValues %in% seq(N))) {
                stop("'initial' needs to be a vector of length 'len'
                    containing N unique values.")
        }
        z <- as.integer(as.factor(initial))
    } else {
        z <- rep(NA, len)
    }

    # Set any values that need to be fixed during sampling
    if (!is.null(fixed)) {
        fixedValues <- sort(unique(fixed))
        if (length(fixed) != len || !all(fixedValues %in% seq(N))) {
            stop("'fixed' to be a vector of length 'len' where each entry is
                one of N unique values or NA.")
        }
        fixedIx <- !is.na(fixed)
        z[fixedIx] <- fixed[fixedIx]
        zNotUsed <- setdiff(seq(N), unique(fixed[fixedIx]))
    } else {
        zNotUsed <- seq(N)
        fixedIx <- rep(FALSE, len)
    }

    # Randomly sample remaining values
    zNa <- which(is.na(z))
    if (length(zNa) > 0) {
        z[zNa] <- sample(zNotUsed, length(zNa), replace = TRUE)
    }

    # Check to ensure each value is in the vector at least once
    missing <- setdiff(seq(N), z)
    for (i in missing) {
        ta <- sort(table(z[!fixedIx]), decreasing = TRUE)
        if (ta[1] == 1) {
            stop("'len' is not long enough to accomodate 'N' unique values")
        }
        ix <- which(z == as.integer(names(ta))[1] & !fixedIx)
        z[sample(ix, 1)] <- i
    }

    return(z)
}


.initializeSplitZ <- function(counts,
    K,
    KSubcluster = NULL,
    alpha = 1,
    beta = 1,
    minCell = 3) {

    s <- rep(1, ncol(counts))
    if (is.null(KSubcluster))
        KSubcluster <- ceiling(sqrt(K))

    # Initialize the model with KSubcluster clusters
    res <- .celda_C(
        counts,
        K = KSubcluster,
        maxIter = 20,
        zInitialize = "random",
        alpha = alpha,
        beta = beta,
        splitOnIter = -1,
        splitOnLast = FALSE,
        verbose = FALSE,
        reorder = FALSE
    )
    overallZ <- as.integer(as.factor(res@clusters$z))
    currentK <- max(overallZ)

    while (currentK < K) {
        # Determine which clusters are split-able
        # KRemaining <- K - currentK
        KPerCluster <- min(ceiling(K / currentK), KSubcluster)
        KToUse <- ifelse(KPerCluster < 2, 2, KPerCluster)

        zTa <- tabulate(overallZ, max(overallZ))

        zToSplit <- which(zTa > minCell & zTa > KToUse)
        if (length(zToSplit) > 1) {
            zToSplit <- sample(zToSplit)
        } else if (length(zToSplit) == 0) {
            break
        }

        # Cycle through each splitable cluster and split it up into
        # K.sublcusters
        for (i in zToSplit) {
            clustLabel <- .celda_C(counts[, overallZ == i, drop = FALSE],
                K = KToUse,
                zInitialize = "random",
                alpha = alpha,
                beta = beta,
                maxIter = 20,
                splitOnIter = -1,
                splitOnLast = FALSE,
                verbose = FALSE)
            tempZ <- as.integer(as.factor(clustLabel@clusters$z))

            # Reassign clusters with label > 1
            splitIx <- tempZ > 1
            ix <- overallZ == i
            newZ <- overallZ[ix]
            newZ[splitIx] <- currentK + tempZ[splitIx] - 1

            overallZ[ix] <- newZ
            currentK <- max(overallZ)

            # Ensure that the maximum number of clusters does not get too large'
            if (currentK > K + 10) {
                break
            }
        }
    }

    # Decompose counts for likelihood calculation
    p <- .cCDecomposeCounts(counts, s, overallZ, currentK)
    nS <- p$nS
    nG <- p$nG
    nM <- p$nM
    mCPByS <- p$mCPByS
    nGByCP <- p$nGByCP
    nCP <- p$nCP
    nByC <- p$nByC

    # Remove clusters 1-by-1 until K is reached
    while (currentK > K) {
        # Find second best assignment give current assignments for each cell
        probs <- .cCCalcEMProbZ(counts,
                s = s,
                z = overallZ,
                K = currentK,
                mCPByS = mCPByS,
                nGByCP = nGByCP,
                nByC = nByC,
                nCP = nCP,
                nG = nG,
                nM = nM,
                alpha = alpha,
                beta = beta,
                doSample = FALSE)
        zProb <- t(probs$probs)
        zProb[cbind(seq(nrow(zProb)), overallZ)] <- NA
        zSecond <- apply(zProb, 1, which.max)

        zTa <- tabulate(overallZ, currentK)
        zNonEmpty <- which(zTa > 0)

        # Find worst cluster by logLik to remove
        previousZ <- overallZ
        llShuffle <- rep(NA, currentK)
        for (i in zNonEmpty) {
            ix <- overallZ == i
            newZ <- overallZ
            newZ[ix] <- zSecond[ix]

            p <- .cCReDecomposeCounts(counts,
                s,
                newZ,
                previousZ,
                nGByCP,
                currentK)
            nGByCP <- p$nGByCP
            mCPByS <- p$mCPByS
            llShuffle[i] <- .cCCalcLL(mCPByS,
                    nGByCP,
                    s,
                    newZ,
                    currentK,
                    nS,
                    nG,
                    alpha,
                    beta)
            previousZ <- newZ
        }

        # Remove the cluster which had the the largest likelihood after removal
        zToRemove <- which.max(llShuffle)

        ix <- overallZ == zToRemove
        overallZ[ix] <- zSecond[ix]

        p <- .cCReDecomposeCounts(counts,
            s,
            overallZ,
            previousZ,
            nGByCP,
            currentK)
        nGByCP <- p$nGByCP[, -zToRemove, drop = FALSE]
        mCPByS <- p$mCPByS[-zToRemove, , drop = FALSE]
        overallZ <- as.integer(as.factor(overallZ))
        currentK <- currentK - 1
    }
    return(overallZ)
}


.initializeSplitY <- function(counts,
    L,
    LSubcluster = NULL,
    tempK = 100,
    beta = 1,
    delta = 1,
    gamma = 1,
    minFeature = 3) {

    if (is.null(LSubcluster))
        LSubcluster <- ceiling(sqrt(L))

    # Collapse cells to managable number of clusters
    if (!is.null(tempK) && ncol(counts) > tempK) {
        z <- .initializeSplitZ(counts, K = tempK)
        counts <- .colSumByGroup(counts, z, length(unique(z)))
    }

    # Initialize the model with KSubcluster clusters
    res <- .celda_G(counts,
        L = LSubcluster,
        maxIter = 10,
        yInitialize = "random",
        beta = beta,
        delta = delta,
        gamma = gamma,
        splitOnIter = -1,
        splitOnLast = FALSE,
        verbose = FALSE,
        reorder = FALSE
    )
    overallY <- as.integer(as.factor(res@clusters$y))
    currentL <- max(overallY)

    while (currentL < L) {
        # Determine which clusters are split-able
        yTa <- tabulate(overallY, max(overallY))
        yToSplit <- sample(which(yTa > minFeature & yTa > LSubcluster))

        if (length(yToSplit) == 0) {
            break
        }

        # Cycle through each splitable cluster and split it up into
        # LSublcusters
        for (i in yToSplit) {
            clustLabel <- .celda_G(
                counts[overallY == i, , drop = FALSE],
                L = LSubcluster,
                yInitialize = "random",
                beta = beta,
                delta = delta,
                gamma = gamma,
                maxIter = 20,
                splitOnIter = -1,
                splitOnLast = FALSE,
                verbose = FALSE)
            tempY <- as.integer(as.factor(clustLabel@clusters$y))

            # Reassign clusters with label > 1
            splitIx <- tempY > 1
            ix <- overallY == i
            newY <- overallY[ix]
            newY[splitIx] <- currentL + tempY[splitIx] - 1

            overallY[ix] <- newY
            currentL <- max(overallY)

            # Ensure that the maximum number of clusters does not get too large
            if (currentL > L + 10) {
                break
            }
        }
    }

    ## Decompose counts for likelihood calculation
    p <- .cGDecomposeCounts(counts = counts, y = overallY, L = currentL)
    nTSByC <- p$nTSByC
    nByG <- p$nByG
    nByTS <- p$nByTS
    nGByTS <- p$nGByTS
    nM <- p$nM
    nG <- p$nG
    rm(p)

    # Pre-compute lgamma values
    lgbeta <- lgamma((seq(0, max(.colSums(counts, nrow(counts),
        ncol(counts))))) + beta)
    lggamma <- lgamma(seq(0, nrow(counts) + L) + gamma)
    lgdelta <- c(NA, lgamma(seq(nrow(counts) + L) * delta))

    # Remove clusters 1-by-1 until L is reached
    while (currentL > L) {
        # Find second best assignment give current assignments for each cell
        probs <- .cGCalcGibbsProbY(
            counts = counts,
            y = overallY,
            L = currentL,
            nTSByC = nTSByC,
            nByTS = nByTS,
            nGByTS = nGByTS,
            nByG = nByG,
            nG = nG,
            beta = beta,
            delta = delta,
            gamma = gamma,
            lgbeta = lgbeta,
            lggamma = lggamma,
            lgdelta = lgdelta,
            doSample = FALSE)
        yProb <- t(probs$probs)
        yProb[cbind(seq(nrow(yProb)), overallY)] <- NA
        ySecond <- apply(yProb, 1, which.max)

        yTa <- tabulate(overallY, currentL)
        yNonEmpty <- which(yTa > 0)

        # Find worst cluster by logLik to remove
        previousY <- overallY
        llShuffle <- rep(NA, currentL)
        for (i in yNonEmpty) {
            ix <- overallY == i
            newY <- overallY
            newY[ix] <- ySecond[ix]

            # Move arounds counts for likelihood calculation
            p <- .cGReDecomposeCounts(counts,
                    newY,
                    previousY,
                    nTSByC,
                    nByG,
                    currentL)
            nTSByC <- p$nTSByC
            nGByTS <- p$nGByTS
            nByTS <- p$nByTS
            llShuffle[i] <- .cGCalcLL(nTSByC,
                    nByTS,
                    nByG,
                    nGByTS,
                    nM,
                    nG,
                    currentL,
                    beta,
                    delta,
                    gamma)
            previousY <- newY
        }

        # Remove the cluster which had the the largest likelihood after removal
        yToRemove <- which.max(llShuffle)

        ix <- overallY == yToRemove
        overallY[ix] <- ySecond[ix]

        # Move around counts and remove module
        p <- .cGReDecomposeCounts(counts,
            overallY,
            previousY,
            nTSByC,
            nByG,
            currentL)
        nTSByC <- p$nTSByC[-yToRemove, , drop = FALSE]
        nGByTS <- p$nGByTS[-yToRemove]
        nByTS <- p$nByTS[-yToRemove]
        overallY <- as.integer(as.factor(overallY))
        currentL <- currentL - 1
    }
    return(overallY)
}
