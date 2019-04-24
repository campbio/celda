# .cCCalcLL = function(mCPByS, nGByCP, s, z, K, nS, nG, alpha, beta)
.cCSplitZ <- function(counts,
    mCPByS,
    nGByCP,
    nCP,
    s,
    z,
    K,
    nS,
    nG,
    alpha,
    beta,
    zProb,
    maxClustersToTry = 10,
    minCell = 3) {

    ## Identify clusters to split
    zTa <- tabulate(z, K)
    zToSplit <- which(zTa >= minCell)
    zNonEmpty <- which(zTa > 0)

    if (length(zToSplit) == 0) {
        m <- paste0(date(),
            " .... Cluster sizes too small. No additional splitting was",
            " performed.")
        return(list(z = z,
            mCPByS,
            nGByCP,
            nCP = nCP,
            message = m))
    }

    ## Loop through each split-able Z and perform split
    clustSplit <- lapply(zToSplit, function(x) {
        clusters(.celda_C(counts[, z == x],
            K = 2,
            zInitialize = "random",
            maxIter = 5,
            splitOnIter = -1,
            splitOnLast = FALSE,
            verbose = FALSE))$z
    })

    clustSplit <- vector("list", K)
    for (i in zToSplit) {
        clustLabel <- .celda_C(
            counts[, z == i],
            K = 2,
            zInitialize = "random",
            maxIter = 5,
            splitOnIter = -1,
            splitOnLast = FALSE,
            verbose = FALSE
        )
        clustSplit[[i]] <- clusters(clustLabel)$z
    }

    ## Find second best assignment give current assignments for each cell
    zProb[cbind(seq(nrow(zProb)), z)] <- NA
    zSecond <- apply(zProb, 1, which.max)

    ## Set up initial variables
    zSplit <- matrix(NA,
        nrow = length(z),
        ncol = length(zToSplit) * maxClustersToTry)
    zSplitLl <- rep(NA, times = length(zToSplit) * maxClustersToTry)
    zSplitLl[1] <- .cCCalcLL(mCPByS, nGByCP, s, z, K, nS, nG, alpha, beta)
    zSplit[, 1] <- z

    ## Select worst clusters to test for reshuffling
    previousZ <- z
    llShuffle <- rep(NA, K)
    for (i in zNonEmpty) {
        ix <- z == i
        newZ <- z
        newZ[ix] <- zSecond[ix]

        p <- .cCReDecomposeCounts(counts, s, newZ, previousZ, nGByCP, K)
        nGByCP <- p$nGByCP
        mCPByS <- p$mCPByS
        llShuffle[i] <- .cCCalcLL(mCPByS, nGByCP, s, z, K, nS, nG, alpha, beta)
        previousZ <- newZ
    }
    zToShuffle <- utils::head(order(llShuffle, decreasing = TRUE, na.last = NA),
        n = maxClustersToTry)

    pairs <- c(NA, NA)
    splitIx <- 2
    for (i in zToShuffle) {
        otherClusters <- setdiff(zToSplit, i)

        for (j in otherClusters) {
            newZ <- z

            ## Assign cluster i to the next most similar cluster (excluding
            ## cluster j)
            ## as defined above by the correlation
            ixToMove <- z == i
            newZ[ixToMove] <- zSecond[ixToMove]

            ## Split cluster j according to the clustering defined above
            ixToSplit <- z == j
            newZ[ixToSplit] <- ifelse(clustSplit[[j]] == 1, j, i)

            p <- .cCReDecomposeCounts(counts, s, newZ, previousZ, nGByCP, K)
            nGByCP <- p$nGByCP
            mCPByS <- p$mCPByS

            ## Calculate likelihood of split
            zSplitLl[splitIx] <- .cCCalcLL(mCPByS,
                nGByCP,
                s,
                z,
                K,
                nS,
                nG,
                alpha,
                beta)
            zSplit[, splitIx] <- newZ
            splitIx <- splitIx + 1L
            previousZ <- newZ

            pairs <- rbind(pairs, c(i, j))
        }
    }

    select <- which.max(zSplitLl)

    if (select == 1) {
        m <- paste0(date(), " .... No additional splitting was performed.")
    } else {
        m <- paste0(date(),
            " .... Cluster ",
            pairs[select, 1],
            " was reassigned and cluster ",
            pairs[select, 2],
            " was split in two."
        )
    }

    p <- .cCReDecomposeCounts(counts, s, zSplit[, select], previousZ, nGByCP, K)
    return(list(z = zSplit[, select],
        mCPByS = p$mCPByS,
        nGByCP = p$nGByCP,
        nCP = p$nCP,
        message = m))
}


# .cCGCalcLL = function(K, L, mCPByS, nTSByCP, nByG, nByTS, nGByTS,
# nS, nG, alpha, beta, delta, gamma)
.cCGSplitZ <- function(counts,
    mCPByS,
    nTSByC,
    nTSByCP,
    nByG,
    nByTS,
    nGByTS,
    nCP,
    s,
    z,
    K,
    L,
    nS,
    nG,
    alpha,
    beta,
    delta,
    gamma,
    zProb,
    maxClustersToTry = 10,
    minCell = 3) {

    ## Identify clusters to split
    zTa <- tabulate(z, K)
    zToSplit <- which(zTa >= minCell)
    zNonEmpty <- which(zTa > 0)

    if (length(zToSplit) == 0) {
        m <- paste0(date(),
            " .... Cluster sizes too small. No additional splitting was",
            " performed.")
        return(list(z = z,
            mCPByS = mCPByS,
            nTSByCP = nTSByCP,
            nCP = nCP,
            message = m))
    }

    ## Loop through each split-able Z and perform split
    clustSplit <- vector("list", K)
    for (i in zToSplit) {
        clustLabel <- .celda_C(counts[, z == i],
            K = 2,
            zInitialize = "random",
            maxIter = 5,
            splitOnIter = -1,
            splitOnLast = FALSE,
            verbose = FALSE)
        clustSplit[[i]] <- clusters(clustLabel)$z
    }

    ## Find second best assignment give current assignments for each cell
    zProb[cbind(seq(nrow(zProb)), z)] <- NA
    zSecond <- apply(zProb, 1, which.max)

    ## Set up initial variables
    zSplit <- matrix(NA,
        nrow = length(z),
        ncol = length(zToSplit) * maxClustersToTry)
    zSplitLl <- rep(NA, ncol = length(zToSplit) * maxClustersToTry)
    zSplitLl[1] <- .cCGCalcLL(K,
        L,
        mCPByS,
        nTSByCP,
        nByG,
        nByTS,
        nGByTS,
        nS,
        nG,
        alpha,
        beta,
        delta,
        gamma)
    zSplit[, 1] <- z

    ## Select worst clusters to test for reshuffling
    previousZ <- z
    llShuffle <- rep(NA, K)
    for (i in zNonEmpty) {
        ix <- z == i
        newZ <- z
        newZ[ix] <- zSecond[ix]

        p <- .cCReDecomposeCounts(nTSByC, s, newZ, previousZ, nTSByCP, K)
        nTSByCP <- p$nGByCP
        mCPByS <- p$mCPByS
        llShuffle[i] <- .cCGCalcLL(K,
            L,
            mCPByS,
            nTSByCP,
            nByG,
            nByTS,
            nGByTS,
            nS,
            nG,
            alpha,
            beta,
            delta,
            gamma)
        previousZ <- newZ
    }
    zToShuffle <- utils::head(order(llShuffle, decreasing = TRUE, na.last = NA),
        n = maxClustersToTry)


    pairs <- c(NA, NA)
    splitIx <- 2
    for (i in zToShuffle) {
        otherClusters <- setdiff(zToSplit, i)

        for (j in otherClusters) {
            newZ <- z

            ## Assign cluster i to the next most similar cluster (excluding
            ## cluster j)
            ## as defined above by the correlation
            ixToMove <- z == i
            newZ[ixToMove] <- zSecond[ixToMove]

            ## Split cluster j according to the clustering defined above
            ixToSplit <- z == j
            newZ[ixToSplit] <- ifelse(clustSplit[[j]] == 1, j, i)

            p <- .cCReDecomposeCounts(nTSByC, s, newZ, previousZ, nTSByCP, K)
            nTSByCP <- p$nGByCP
            mCPByS <- p$mCPByS

            ## Calculate likelihood of split
            zSplitLl[splitIx] <- .cCGCalcLL(K,
                L,
                mCPByS,
                nTSByCP,
                nByG,
                nByTS,
                nGByTS,
                nS,
                nG,
                alpha,
                beta,
                delta,
                gamma)
            zSplit[, splitIx] <- newZ
            splitIx <- splitIx + 1L
            previousZ <- newZ

            pairs <- rbind(pairs, c(i, j))
        }
    }

    select <- which.max(zSplitLl)

    if (select == 1) {
        m <- paste0(date(), " .... No additional splitting was performed.")
    } else {
        m <- paste0(date(),
            " .... Cluster ",
            pairs[select, 1],
            " was reassigned and cluster ",
            pairs[select, 2],
            " was split in two.")
    }

    p <- .cCReDecomposeCounts(nTSByC,
        s,
        zSplit[, select],
        previousZ,
        nTSByCP,
        K)
    return(list(z = zSplit[, select],
        mCPByS = p$mCPByS,
        nTSByCP = p$nGByCP,
        nCP = p$nCP,
        message = m))
}


.cCGSplitY <- function(counts,
    y,
    mCPByS,
    nGByCP,
    nTSByC,
    nTSByCP,
    nByG,
    nByTS,
    nGByTS,
    nCP,
    s,
    z,
    K,
    L,
    nS,
    nG,
    alpha,
    beta,
    delta,
    gamma,
    yProb,
    maxClustersToTry = 10,
    KSubclusters = 10,
    minCell = 3) {

    #########################
    ## First, the cell dimension of the original matrix will be reduced by
    ## splitting each z cluster into 'KSubclusters'.
    #########################

    ## This will not be as big as the original matrix (which can take a lot of
    ## time to process with large number of cells), but not as small as the
    ## 'nGByCP' with current z assignments

    zTa <- tabulate(z, K)
    zNonEmpty <- which(zTa > 0)
    tempZ <- rep(0, length(z))
    currentTopZ <- 0
    for (i in zNonEmpty) {
        ix <- z == i
        if (zTa[i] <= KSubclusters) {
            tempZ[ix] <- seq(currentTopZ + 1, currentTopZ + zTa[i])
        } else {
            clustLabel <- .celda_C(counts[, z == i],
                K = KSubclusters,
                zInitialize = "random",
                maxIter = 5,
                splitOnIter = -1,
                splitOnLast = FALSE,
                verbose = FALSE)
            tempZ[ix] <- clusters(clustLabel)$z + currentTopZ
        }
        currentTopZ <- max(tempZ, na.rm = TRUE)
    }

    ## Decompose counts according to new/temp z labels
    tempNGByCP <- .colSumByGroup(counts, group = tempZ, K = currentTopZ)

    #########################
    ## Second, different y splits will be estimated and tested
    #########################

    ## Identify clusters to split
    yTa <- tabulate(y, L)
    yToSplit <- which(yTa >= minCell)
    yNonEmpty <- which(yTa > 0)

    if (length(yToSplit) == 0) {
        m <- paste0(date(),
            " .... Cluster sizes too small. No additional splitting was",
            " performed.")
        return(list(y = y,
            mCPByS = mCPByS,
            nTSByCP = nTSByCP,
            nCP = nCP,
            message = m))
    }

    ## Loop through each split-able Z and perform split
    clustSplit <- vector("list", L)
    for (i in yToSplit) {
        clustLabel <- .celda_G(tempNGByCP[y == i, ],
            L = 2,
            yInitialize = "random",
            maxIter = 5,
            splitOnIter = -1,
            splitOnLast = FALSE,
            verbose = FALSE)
        clustSplit[[i]] <- clusters(clustLabel)$y
    }

    ## Find second best assignment give current assignments for each cell
    yProb[cbind(seq(nrow(yProb)), y)] <- NA
    ySecond <- apply(yProb, 1, which.max)

    ## Set up initial variables
    ySplit <- matrix(NA,
        nrow = length(y),
        ncol = length(yToSplit) * maxClustersToTry)
    ySplitLl <- rep(NA, ncol = length(yToSplit) * maxClustersToTry)
    ySplitLl[1] <- .cCGCalcLL(K,
        L,
        mCPByS,
        nTSByCP,
        nByG,
        nByTS,
        nGByTS,
        nS,
        nG,
        alpha,
        beta,
        delta,
        gamma)
    ySplit[, 1] <- y

    ## Select worst clusters to test for reshuffling
    previousY <- y
    llShuffle <- rep(NA, L)
    for (i in yNonEmpty) {
        ix <- y == i
        newY <- y
        newY[ix] <- ySecond[ix]

        p <- .cGReDecomposeCounts(nGByCP, newY, previousY, nTSByCP, nByG, L)
        nTSByCP <- p$nTSByC
        nByTS <- p$nByTS
        nGByTS <- p$nGByTS

        llShuffle[i] <- .cCGCalcLL(K,
            L,
            mCPByS,
            nTSByCP,
            nByG,
            nByTS,
            nGByTS,
            nS,
            nG,
            alpha,
            beta,
            delta,
            gamma)
        previousY <- newY
    }
    yToShuffle <- utils::head(order(llShuffle, decreasing = TRUE, na.last = NA),
        n = maxClustersToTry)

    pairs <- c(NA, NA)
    splitIx <- 2
    for (i in yToShuffle) {
        otherClusters <- setdiff(yToSplit, i)

        for (j in otherClusters) {
            newY <- y

            ## Assign cluster i to the next most similar cluster (excluding
            ## cluster j)
            ## as defined above by the correlation
            ixToMove <- y == i
            newY[ixToMove] <- ySecond[ixToMove]

            ## Split cluster j according to the clustering defined above
            ixToSplit <- y == j
            newY[ixToSplit] <- ifelse(clustSplit[[j]] == 1, j, i)

            p <- .cGReDecomposeCounts(nGByCP, newY, previousY, nTSByCP, nByG, L)
            nTSByCP <- p$nTSByC
            nByTS <- p$nByTS
            nGByTS <- p$nGByTS

            ## Calculate likelihood of split
            ySplitLl[splitIx] <- .cCGCalcLL(
                K,
                L,
                mCPByS,
                nTSByCP,
                nByG,
                nByTS,
                nGByTS,
                nS,
                nG,
                alpha,
                beta,
                delta,
                gamma)
            ySplit[, splitIx] <- newY
            splitIx <- splitIx + 1L
            previousY <- newY

            pairs <- rbind(pairs, c(i, j))
        }
    }

    select <- which.max(ySplitLl)

    if (select == 1) {
        m <- paste0(date(), " .... No additional splitting was performed.")
    } else {
        m <- paste0(date(),
            " .... Cluster ",
            pairs[select, 1],
            " was reassigned and cluster ",
            pairs[select, 2],
            " was split in two.")
    }

    p <- .cGReDecomposeCounts(nGByCP,
        ySplit[, select],
        previousY,
        nTSByCP,
        nByG,
        L)
    return(list(y = ySplit[, select],
        nTSByCP = p$nTSByC,
        nByTS = p$nByTS,
        nGByTS = p$nGByTS,
        message = m))
}


# .cGCalcLL = function(nTSByC, nByTS, nByG, nGByTS, nM, nG, L, beta, delta,
# gamma) {
.cGSplitY <- function(counts,
    y,
    nTSByC,
    nByTS,
    nByG,
    nGByTS,
    nM,
    nG,
    L,
    beta,
    delta,
    gamma,
    yProb,
    minFeature = 3,
    maxClustersToTry = 10) {

    ## Identify clusters to split
    yTa <- table(factor(y, levels = seq(L)))
    yToSplit <- which(yTa >= minFeature)
    yNonEmpty <- which(yTa > 0)

    if (length(yToSplit) == 0) {
        m <- paste0(date(),
            " .... Cluster sizes too small. No additional splitting was",
            " performed.")
        return(list(y = y,
            nTSByC = nTSByC,
            nByTS = nByTS,
            nGByTS = nGByTS,
            message = m))
    }

    ## Loop through each split-able y and find best split
    clustSplit <- vector("list", L)
    for (i in yToSplit) {
        clustLabel <- .celda_G(counts[y == i, ],
            L = 2,
            yInitialize = "random",
            maxIter = 5,
            splitOnIter = -1,
            splitOnLast = FALSE,
            verbose = FALSE)
        clustSplit[[i]] <- clusters(clustLabel)$y
    }

    ## Find second best assignment give current assignments for each cell
    yProb[cbind(seq(nrow(yProb)), y)] <- NA
    ySecond <- apply(yProb, 1, which.max)

    ## Set up initial variables
    ySplit <- matrix(NA,
        nrow = length(y),
        ncol = length(yToSplit) * maxClustersToTry)
    ySplitLl <- rep(NA, ncol = length(yToSplit) * maxClustersToTry)
    ySplitLl[1] <- .cGCalcLL(nTSByC,
        nByTS,
        nByG,
        nGByTS,
        nM,
        nG,
        L,
        beta,
        delta,
        gamma)
    ySplit[, 1] <- y

    ## Select worst clusters to test for reshuffling
    llShuffle <- rep(NA, L)
    previousY <- y
    for (i in yNonEmpty) {
        ix <- y == i
        newY <- y
        newY[ix] <- ySecond[ix]
        p <- .cGReDecomposeCounts(counts, newY, previousY, nTSByC, nByG, L)
        llShuffle[i] <- .cGCalcLL(p$nTSByC,
            p$nByTS,
            nByG,
            p$nGByTS,
            nM,
            nG,
            L,
            beta,
            delta,
            gamma)
        previousY <- newY
    }
    yToShuffle <- utils::head(order(llShuffle, decreasing = TRUE, na.last = NA),
        n = maxClustersToTry)

    pairs <- c(NA, NA)
    splitIx <- 2
    for (i in yToShuffle) {
        otherClusters <- setdiff(yToSplit, i)

        for (j in otherClusters) {
            newY <- y

            ## Assign cluster i to the next most similar cluster (excluding
            ## cluster j)
            ## as defined above by the spearman correlation
            ixToMove <- y == i
            newY[ixToMove] <- ySecond[ixToMove]

            ## Split cluster j according to the clustering defined above
            ixToSplit <- y == j
            newY[ixToSplit] <- ifelse(clustSplit[[j]] == 1, j, i)

            ## Calculate likelihood of split
            p <- .cGReDecomposeCounts(counts, newY, previousY, nTSByC, nByG, L)
            ySplitLl[splitIx] <- .cGCalcLL(p$nTSByC,
                p$nByTS,
                nByG,
                p$nGByTS,
                nM,
                nG,
                L,
                beta,
                delta,
                gamma)
            ySplit[, splitIx] <- newY
            splitIx <- splitIx + 1L
            previousY <- newY

            pairs <- rbind(pairs, c(i, j))
        }
    }

    select <- which.max(ySplitLl)

    if (select == 1) {
        m <- paste0(date(), " .... No additional splitting was performed.")
    } else {
        m <- paste0(date(),
            " .... Cluster ",
            pairs[select, 1],
            " was reassigned and cluster ",
            pairs[select, 2],
            " was split in two.")
    }

    p <- .cGReDecomposeCounts(counts,
        ySplit[, select],
        previousY,
        nTSByC,
        nByG,
        L)
    return(list(y = ySplit[, select],
        nTSByC = p$nTSByC,
        nByTS = p$nByTS,
        nGByTS = p$nGByTS,
        message = m))
}
