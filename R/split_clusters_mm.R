# .cCCalcLL = function(mCPByS, nGByCP, s, z, K, nS, nG, alpha, beta)
.cCSplitZMM <- function(counts,
    counts2,
    mCPByS,
    nGByCP,
    nGByCP2,
    nCP,
    nCP2
    s,
    z,
    K,
    nS,
    nG,
    nG2,
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
            mCPByS = mCPByS,
            nGByCP = nGByCP,
            nGByCP2 = nGByCP2,
            nCP = nCP,
            nCP2 = nCP2,
            message = m))
    }

    ## Loop through each split-able Z and perform split
    clustSplit <- vector("list", K)
    for (i in zToSplit) {
        clustLabel <- .celda_C_mm(
            counts = counts[, z == i],
            counts2 = counts2[, z == i],
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
    zSplitLl <- rep(NA, times = length(zToSplit) * maxClustersToTry)
    zSplitLl[1] <- .cCCalcLLMM(mCPByS = mCPByS,
        nGByCP = nGByCP,
        nGByCP2 = nGByCP2,
        K = K,
        nS = nS,
        nG = nG,
        nG2 = nG2,
        alpha = alpha,
        beta = beta)
    zSplit[, 1] <- z

    ## Select worst clusters to test for reshuffling
    previousZ <- z
    llShuffle <- rep(NA, K)
    for (i in zNonEmpty) {
        ix <- z == i
        newZ <- z
        newZ[ix] <- zSecond[ix]

        p <- .cCReDecomposeCountsMM(counts,
            counts2,
            s,
            newZ,
            previousZ,
            nGByCP,
            nGByCP2,
            K)
        mCPByS <- p$mCPByS
        nGByCP <- p$nGByCP
        nGByCP2 <- p$nGByCP2

        llShuffle[i] <- .cCCalcLLMM(mCPByS = mCPByS,
            nGByCP = nGByCP,
            nGByCP2 = nGByCP2,
            K = K,
            nS = nS,
            nG = nG,
            nG2 = nG2,
            alpha = alpha,
            beta = beta)
        previousZ <- newZ
    }
    zToShuffle <- utils::head(order(llShuffle,
        decreasing = TRUE,
        na.last = NA),
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

            p <- .cCReDecomposeCountsMM(counts,
                counts2,
                s,
                newZ,
                previousZ,
                nGByCP,
                nGByCP2,
                K)
            mCPByS <- p$mCPByS
            nGByCP <- p$nGByCP
            nGByCP2 <- p$nGByCP2

            ## Calculate likelihood of split
            zSplitLl[splitIx] <- .cCCalcLLMM(mCPByS = mCPByS,
                nGByCP = nGByCP,
                nGByCP2 = nGByCP2,
                K = K,
                nS = nS,
                nG = nG,
                nG2 = nG2,
                alpha = alpha,
                beta = beta)
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

    p <- .cCReDecomposeCountsMM(counts = counts,
        counts2 = counts2,
        s = s,
        z = zSplit[, select],
        previousZ = previousZ,
        nGByCP = nGByCP,
        nGByCP2 = nGByCP2,
        K = K)
    return(list(z = zSplit[, select],
        mCPByS = p$mCPByS,
        nGByCP = p$nGByCP,
        nGByCP2 = p$nGByCP2,
        nCP = p$nCP,
        nCP2 = p$nCP2,
        message = m))
}
