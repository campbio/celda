.initializeCluster <- function(N,
    len,
    z = NULL,
    initial = NULL,
    fixed = NULL,
    seed = 12345) {

    # If initial values are given, then they will not be randomly initialized
    if (!is.null(initial)) {
        init.values <- sort(unique(initial))
        if (length(unique(initial)) != N || length(initial) != len ||
            !all(init.values %in% 1:N)) {
                stop("'initial' needs to be a vector of length 'len'
                    containing N unique values.")
        }
        z <- as.numeric(as.factor(initial))
    } else {
        z <- rep(NA, len)
    }

    # Set any values that need to be fixed during sampling
    if (!is.null(fixed)) {
        fixed.values <- sort(unique(fixed))
        if (length(fixed) != len || !all(fixed.values %in% 1:N)) {
            stop("'fixed' to be a vector of length 'len' where each entry is
                one of N unique values or NA.")
        }
        fixed.ix <- !is.na(fixed)
        z[fixed.ix] <- fixed[fixed.ix]
        z.not.used <- setdiff(1:N, unique(fixed[fixed.ix]))
    } else {
        z.not.used <- 1:N
        fixed.ix <- rep(FALSE, len)
    }

    # Randomly sample remaining values
    setSeed(seed)
    z.na <- which(is.na(z))
    if (length(z.na) > 0) {
        z[z.na] <- sample(z.not.used, length(z.na), replace = TRUE)
    }

    # Check to ensure each value is in the vector at least once
    missing <- setdiff(1:N, z)
    for (i in missing) {
        ta <- sort(table(z[!fixed.ix]), decreasing = TRUE)
        if (ta[1] == 1) {
            stop("'len' is not long enough to accomodate 'N' unique values")
        }
        ix <- which(z == as.numeric(names(ta))[1] & !fixed.ix)
        z[sample(ix, 1)] <- i
    }

    return(z)
}


.initializeSplitZ <- function(counts,
    K,
    K.subcluster = NULL,
    alpha = 1,
    beta = 1,
    min.cell = 3,
    seed = 12345) {

    s <- rep(1, ncol(counts))
    if (is.null(K.subcluster))
        K.subcluster <- ceiling(sqrt(K))

    # Initialize the model with K.subcluster clusters
    res <- .celda_C(
        counts,
        K = K.subcluster,
        max.iter = 20,
        z.initialize = "random",
        alpha = alpha,
        beta = beta,
        split.on.iter = -1,
        split.on.last = FALSE,
        verbose = FALSE,
        seed = seed,
        reorder = FALSE
    )
    overall.z <- as.integer(as.factor(res@clusters$z))
    current.K <- max(overall.z)

    while (current.K < K) {
        # Determine which clusters are split-able
        K.remaining <- K - current.K
        K.per.cluster <- min(ceiling(K / current.K), K.subcluster)
        K.to.use <- ifelse(K.per.cluster < 2, 2, K.per.cluster)

        z.ta <- tabulate(overall.z, max(overall.z))
        z.to.split <- sample(which(z.ta > min.cell & z.ta > K.to.use))

        if (length(z.to.split) == 0) {
            break()
        }

        # Cycle through each splitable cluster and split it up into K.sublcusters
        for (i in z.to.split) {
            clustLabel <- .celda_C(counts[, overall.z == i, drop = FALSE],
                K = K.to.use,
                z.initialize = "random",
                alpha = alpha,
                beta = beta,
                max.iter = 20,
                split.on.iter = -1,
                split.on.last = FALSE,
                verbose = FALSE)
            temp.z <- as.integer(as.factor(clustLabel@clusters$z))

            # Reassign clusters with label > 1
            split.ix <- temp.z > 1
            ix <- overall.z == i
            new.z <- overall.z[ix]
            new.z[split.ix] <- current.K + temp.z[split.ix] - 1

            overall.z[ix] <- new.z
            current.K <- max(overall.z)

            # Ensure that the maximum number of clusters does not get too large'
            if (current.K > K + 10) {
                break()
            }
        }
    }

    # Decompose counts for likelihood calculation
    p <- cC.decomposeCounts(counts, s, overall.z, current.K)
    nS <- p$nS
    nG <- p$nG
    nM <- p$nM
    m.CP.by.S <- p$m.CP.by.S
    n.G.by.CP <- p$n.G.by.CP
    n.CP <- p$n.CP
    n.by.C <- p$n.by.C

    # Remove clusters 1-by-1 until K is reached
    while (current.K > K) {
        # Find second best assignment give current assignments for each cell
        probs <- cC.calcEMProbZ(counts,
                s = s,
                z = overall.z,
                K = current.K,
                m.CP.by.S = m.CP.by.S,
                n.G.by.CP = n.G.by.CP,
                n.by.C = n.by.C,
                n.CP = n.CP,
                nG = nG,
                nM = nM,
                alpha = alpha,
                beta = beta,
                do.sample = FALSE
            )
        z.prob <- t(probs$probs)
        z.prob[cbind(1:nrow(z.prob), overall.z)] <- NA
        z.second <- apply(z.prob, 1, which.max)

        z.ta <- tabulate(overall.z, current.K)
        z.non.empty <- which(z.ta > 0)

        # Find worst cluster by logLik to remove
        previous.z <- overall.z
        ll.shuffle <- rep(NA, current.K)
        for (i in z.non.empty) {
            ix <- overall.z == i
            new.z <- overall.z
            new.z[ix] <- z.second[ix]

            p <- cC.reDecomposeCounts(counts,
                s,
                new.z,
                previous.z,
                n.G.by.CP,
                current.K)
            n.G.by.CP <- p$n.G.by.CP
            m.CP.by.S <- p$m.CP.by.S
            ll.shuffle[i] <- cC.calcLL(m.CP.by.S,
                    n.G.by.CP,
                    s,
                    new.z,
                    current.K,
                    nS,
                    nG,
                    alpha,
                    beta)
            previous.z <- new.z
        }

        # Remove the cluster which had the the largest likelihood after removal
        z.to.remove <- which.max(ll.shuffle)

        ix <- overall.z == z.to.remove
        overall.z[ix] <- z.second[ix]

        p <- cC.reDecomposeCounts(counts,
            s,
            overall.z,
            previous.z,
            n.G.by.CP,
            current.K)
        n.G.by.CP <- p$n.G.by.CP[, -z.to.remove, drop = FALSE]
        m.CP.by.S <- p$m.CP.by.S[-z.to.remove, , drop = FALSE]
        overall.z <- as.integer(as.factor(overall.z))
        current.K <- current.K - 1
    }
    return(overall.z)
}



.initializeSplitY <- function(counts,
    L,
    L.subcluster = NULL,
    temp.K = 100,
    beta = 1,
    delta = 1,
    gamma = 1,
    min.feature = 3,
    seed = 12345) {
    if (is.null(L.subcluster))
        L.subcluster <- ceiling(sqrt(L))

    # Collapse cells to managable number of clusters
    if (!is.null(temp.K) && ncol(counts) > temp.K) {
        z <- initialize.splitZ(counts, K = temp.K, seed = seed)
        counts <- colSumByGroup(counts, z, length(unique(z)))
    }

    # Initialize the model with K.subcluster clusters
    res <- .celda_G(counts,
        L = L.subcluster,
        max.iter = 10,
        y.initialize = "random",
        beta = beta,
        delta = delta,
        gamma = gamma,
        split.on.iter = -1,
        split.on.last = FALSE,
        verbose = FALSE,
        seed = seed,
        reorder = FALSE
    )
    overall.y <- as.integer(as.factor(res@clusters$y))
    current.L <- max(overall.y)

    while (current.L < L) {
        # Determine which clusters are split-able
        y.ta <- tabulate(overall.y, max(overall.y))
        y.to.split <-
            sample(which(y.ta > min.feature &
                    y.ta > L.subcluster))

        if (length(y.to.split) == 0)
            break()

        # Cycle through each splitable cluster and split it up into L.sublcusters
        for (i in y.to.split) {
            clustLabel <- .celda_G(
                counts[overall.y == i, , drop = FALSE],
                L = L.subcluster,
                y.initialize = "random",
                beta = beta,
                delta = delta,
                gamma = gamma,
                max.iter = 20,
                split.on.iter = -1,
                split.on.last = FALSE,
                verbose = FALSE
            )
            temp.y <- as.integer(as.factor(clustLabel@clusters$y))

            # Reassign clusters with label > 1
            split.ix <- temp.y > 1
            ix <- overall.y == i
            new.y <- overall.y[ix]
            new.y[split.ix] <- current.L + temp.y[split.ix] - 1

            overall.y[ix] <- new.y
            current.L <- max(overall.y)

            # Ensure that the maximum number of clusters does not get too large
            if (current.L > L + 10) {
                break()
            }
        }
    }

    ## Decompose counts for likelihood calculation
    p <- cG.decomposeCounts(counts = counts, y = overall.y, L = current.L)
    n.TS.by.C <- p$n.TS.by.C
    n.by.G <- p$n.by.G
    n.by.TS <- p$n.by.TS
    nG.by.TS <- p$nG.by.TS
    nM <- p$nM
    nG <- p$nG
    rm(p)

    # Pre-compute lgamma values
    lgbeta <- lgamma((0:max(.colSums(counts, nrow(counts), ncol(counts)))) +
        beta)
    lggamma <- lgamma(0:(nrow(counts) + L) + gamma)
    lgdelta <- c(NA, lgamma((1:(nrow(counts) + L) * delta)))

    # Remove clusters 1-by-1 until L is reached
    while (current.L > L) {
        # Find second best assignment give current assignments for each cell
        probs <- cG.calcGibbsProbY(
            counts = counts,
            y = overall.y,
            L = current.L,
            n.TS.by.C = n.TS.by.C,
            n.by.TS = n.by.TS,
            nG.by.TS = nG.by.TS,
            n.by.G = n.by.G,
            nG = nG,
            beta = beta,
            delta = delta,
            gamma = gamma,
            lgbeta = lgbeta,
            lggamma = lggamma,
            lgdelta = lgdelta,
            do.sample = FALSE)
        y.prob <- t(probs$probs)
        y.prob[cbind(1:nrow(y.prob), overall.y)] <- NA
        y.second <- apply(y.prob, 1, which.max)

        y.ta <- tabulate(overall.y, current.L)
        y.non.empty <- which(y.ta > 0)

        # Find worst cluster by logLik to remove
        previous.y <- overall.y
        ll.shuffle <- rep(NA, current.L)
        for (i in y.non.empty) {
            ix <- overall.y == i
            new.y <- overall.y
            new.y[ix] <- y.second[ix]

            # Move arounds counts for likelihood calculation
            p <- cG.reDecomposeCounts(counts,
                    new.y,
                    previous.y,
                    n.TS.by.C,
                    n.by.G,
                    current.L)
            n.TS.by.C <- p$n.TS.by.C
            nG.by.TS <- p$nG.by.TS
            n.by.TS <- p$n.by.TS
            ll.shuffle[i] <- cG.calcLL(n.TS.by.C,
                    n.by.TS,
                    n.by.G,
                    nG.by.TS,
                    nM,
                    nG,
                    current.L,
                    beta,
                    delta,
                    gamma)
            previous.y <- new.y
        }

        # Remove the cluster which had the the largest likelihood after removal
        y.to.remove <- which.max(ll.shuffle)

        ix <- overall.y == y.to.remove
        overall.y[ix] <- y.second[ix]

        # Move around counts and remove module
        p <- cG.reDecomposeCounts(counts,
            overall.y,
            previous.y,
            n.TS.by.C,
            n.by.G,
            current.L)
        n.TS.by.C <- p$n.TS.by.C[-y.to.remove, , drop = FALSE]
        nG.by.TS <- p$nG.by.TS[-y.to.remove]
        n.by.TS <- p$n.by.TS[-y.to.remove]
        overall.y <- as.integer(as.factor(overall.y))
        current.L <- current.L - 1
    }
    return(overall.y)
}