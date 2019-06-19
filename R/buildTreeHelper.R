# This function generates the decision tree be recursively separating classes.
.generateTreeList <- function(
    features,
    class,
    oneoffMetric,
    threshold,
    reuseFeatures,
    consectutiveOneoff = TRUE) {

    # Initialize Tree
    treeLevel <- tree <- list()

    # Initialize the first split
    treeLevel[[1]] <- list()

    # Generate the first split at the first level
    treeLevel[[1]] <- .wrapSplitHybrid(
        features,
        class,
        threshold,
        oneoffMetric
    )

    # Add set of features used at this split
    treeLevel[[1]]$fUsed <- unlist(lapply(
        treeLevel[[1]][names(treeLevel[[1]]) != "statUsed"],
        function(X) {
            X$featureName
        }))

    # Initialize split directions
    treeLevel[[1]]$dirs <- 1

    # Add split list as first level
    tree[[1]] <- treeLevel

    # Initialize tree depth
    mDepth <- 1

    # Build tree until all leafs are of a single cluster
    while (length(unlist(treeLevel)) > 0) {

        # Create list of branches on this level
        outList <- lapply(treeLevel, function(split, features, class) {

            # Check for consecutive oneoff
            tryOneoff <- TRUE
            if (!consectutiveOneoff & split$statUsed == "OO") {
                tryOneoff <- FALSE
            }

            # If length(split == 4) than this split is binary node
            if (length(split) == 4 & length(split[[1]]$group1Consensus) > 1) {

                # Create branch from this split.
                branch1 <- .wrapBranchHybrid(
                split[[1]]$group1,
                features, class,
                split$fUsed,
                threshold,
                reuseFeatures,
                oneoffMetric,
                tryOneoff)

                if (!is.null(branch1)) {

                    # Add feature to list of used features.
                    branch1$fUsed <- c(split$fUsed, unlist(lapply(
                        branch1[names(branch1) != "statUsed"],
                        function(X) {
                            X$featureName
                        })))

                    # Add the split direction (always 1 when splitting group 1)
                    branch1$dirs <- c(split$dirs, 1)
                }
            } else {
                branch1 <- NULL
            }

            # If length(split == 4) than this split is binary node
            if (length(split) == 4 & length(split[[1]]$group2Consensus) > 1) {

                # Create branch from this split
                branch2 <- .wrapBranchHybrid(
                    split[[1]]$group2,
                    features,
                    class,
                    split$fUsed,
                    threshold,
                    reuseFeatures,
                    oneoffMetric,
                    tryOneoff)

                if (!is.null(branch2)) {

                    # Add feature to list of used features.
                    branch2$fUsed <- c(split$fUsed, unlist(lapply(
                        branch2[names(branch2) != "statUsed"],
                        function(X) {
                        X$featureName
                    })))

                # Add the split direction (always 2 when splitting group 2)
                branch2$dirs <- c(split$dirs, 2)
                }

            # If length(split > 4) than this split is more than 2 edges
            # In this case group 1 will always denote leaves.
            } else if (length(split) > 4) {

                # Get samples that are never in group 1 in this split
                group1Samples <- unique(unlist(lapply(
                split[!names(split) %in% c("statUsed", "fUsed", "dirs")],
                    function(X) {
                    X$group1
                })))
                group2Samples <- unique(unlist(lapply(
                    split[!names(split) %in% c("statUsed", "fUsed", "dirs")],
                    function(X) {
                        X$group2
                })))
                group2Samples <- group2Samples[!group2Samples %in%
                    group1Samples]

                # Check that there is still more than one class
                group2Classes <- levels(droplevels(
                    class[rownames(features) %in% group2Samples]))
                if (length(group2Classes) > 1) {

                    # Create branch from this split
                    branch2 <- .wrapBranchHybrid(
                        group2Samples,
                        features,
                        class,
                        split$fUsed,
                        threshold,
                        reuseFeatures,
                        oneoffMetric,
                        tryOneoff)

                    if (!is.null(branch2)) {

                        # Add multiple features
                        branch2$fUsed <- c(split$fUsed, unlist(lapply(
                            branch2[names(branch2) != "statUsed"],
                            function(X) {
                                X$featureName
                            })))

                        # Instead of 2, this direction is 1 + the num. splits
                        branch2$dirs <- c(split$dirs,
                            sum(!names(split) %in%
                                c("statUsed", "fUsed", "dirs")) + 1)
                    }
                } else {
                    branch2 <- NULL
                }
            } else {
                branch2 <- NULL
            }

            # Combine these branches
            outBranch <- list(branch1, branch2)

            # Only keep non-null branches
            outBranch <- outBranch[!unlist(lapply(outBranch, is.null))]
            if (length(outBranch) > 0) {
                return(outBranch)
            } else {
                return(NULL)
            }
        }, features, class)

        # Unlist outList so is one list per 'treeLevel'
        treeLevel <- unlist(outList, recursive = F)

        # Increase tree depth
        mDepth <- mDepth + 1

        # Add this level to the tree
        tree[[mDepth]] <- treeLevel
    }
    return(tree)
}


# Wrapper to subset the feature and class set for each split
.wrapBranchHybrid <- function(
    groups,
    features,
    class,
    fUsed,
    threshold = 0.95,
    reuseFeatures = FALSE,
    oneoffMetric,
    tryOneoff) {

    # Subset for branch to run split
    gKeep <- rownames(features) %in% groups

    # Remove used features?
    if (reuseFeatures) {
        fSub <- features[gKeep, ]
    } else {
        fSub <- features[gKeep, !colnames(features) %in% fUsed, drop = FALSE]
    }

    # Drop levels (class that are no longer in)
    cSub <- droplevels(class[gKeep])

    # If multiple columns in fSub run split, else return null
    if (ncol(fSub) > 1) {
        return(.wrapSplitHybrid(fSub, cSub, threshold, oneoffMetric, tryOneoff))
    } else {
        return(NULL)
    }
}

# Wrapper function to perform split metrics
.wrapSplitHybrid <- function(features,
    class,
    threshold = 0.95,
    oneoffMetric,
    tryOneoff = TRUE) {

    # Get best one-2-one splits
    ## Use modified f1 or pairwise auc?
    if (tryOneoff) {
        if (oneoffMetric == "modified F1") {
            splitMetric <- .splitMetricModF1
        } else {
            splitMetric <- .splitMetricPairwiseAUC
        }
        splitStats <- .splitMetricRecursive(
            features,
            class,
            splitMetric = splitMetric)
        splitStats <- splitStats[splitStats >= threshold]
        statUsed <- "OO"
    } else {
        splitStats <- integer(0)
    }


    # If no one-2-one split meets threshold, run semi-supervised clustering
    if (length(splitStats) == 0) {
        splitMetric <- .splitMetricIGpIGd
        splitStats <- .splitMetricRecursive(features,
            class,
            splitMetric = splitMetric)[1] # Use top
    statUsed <- "IG"
    }

    # Get split for best gene
    splitList <- lapply(
        names(splitStats),
        .getSplit,
        splitStats,
        features,
        class,
        splitMetric)

    # Combine feature rules when same group1 class arises

    if (length(splitList) > 1) {
        group1Vec <- unlist(lapply(
            splitList, function(X) {
                X$group1Consensus
            }), recursive = F)

        splitList <- lapply(
            unique(group1Vec),
            function(group1, splitList, group1Vec) {

                # Get subset with same group1
                splitListSub <- splitList[group1Vec == group1]

                # Get feature, value, and stat for these
                splitFeature <- unlist(lapply(
                    splitListSub,
                    function(X) {
                        X$featureName
                    }))
                splitValue <- unlist(lapply(
                    splitListSub,
                    function(X) {
                        X$value
                    }))
                splitStat <- unlist(lapply(
                    splitListSub,
                    function(X) {
                        X$stat
                    }))

                # Create a single object and add these
                splitSingle <- splitListSub[[1]]
                splitSingle$featureName <- splitFeature
                splitSingle$value <- splitValue
                splitSingle$stat <- splitStat

                return(splitSingle)
        }, splitList, group1Vec)
    }
    names(splitList) <- unlist(lapply(
        splitList,
        function(X) {
            paste(X$featureName, collapse = ";")
        }))

    # Add statUsed
    splitList$statUsed <- statUsed

    return(splitList)
}

# Recursively run split metric on every feature
.splitMetricRecursive <- function(features, class, splitMetric) {
    splitStats <- vapply(
        colnames(features),
        function(feat, features, class, splitMetric) {
            splitMetric(feat, class, features, rPerf = TRUE)
        }, features, class, splitMetric, FUN.VALUE = double(1))
    names(splitStats) <- colnames(features)
    splitStats <- sort(splitStats, decreasing = TRUE)

    return(splitStats)
}

# Run pairwise AUC metirc on single feature
#' @importFrom pROC auc roc coords
.splitMetricPairwiseAUC <- function(feat, class, features, rPerf = FALSE) {

    # Get current feature
    currentFeature <- features[, feat]

    # Get unique classes
    classUnique <- sort(unique(class))

    # Do one-to-all to determine top cluster
    # For each class K1 determine best AUC
    auc1toAll <- vapply(classUnique, function(k1, class, currentFeature) {

        # Set value to k1
        classK1 <- as.numeric(class == k1)

        # Get AUC value
        aucK1 <- pROC::auc(pROC::roc(classK1, currentFeature, direction = "<"))

        # Return
        return(aucK1)
    }, class, currentFeature, FUN.VALUE = double(1))

    # Get class with best AUC (Class with generally highest values)
    classMax <- as.character(classUnique[which.max(auc1toAll)])

    # Get other classes
    classRest <- as.character(classUnique[classUnique != classMax])

    # for each second cluster k2
    aucFram <- as.data.frame(do.call(rbind, lapply(
        classRest,
        function(k2, k1, class, currentFeature) {

            # keep cells in k1 or k2 only
            obsKeep <- class %in% c(k1, k2)
            currentFeatureSubset <- currentFeature[obsKeep]

            # update cluster assignments
            currentClusters <- class[obsKeep]

            # label cells whether they belong to k1 (0 or 1)
            currentLabels <- as.integer(currentClusters == k1)

            # get AUC value for this feat-cluster pair
            rocK2 <- pROC::roc(currentLabels, currentFeatureSubset,
                direction = "<")
            aucK2 <- rocK2$auc
            coordK2 <- pROC::coords(rocK2, "best", ret = "threshold")[1]

            # Concatenate vectors
            statK2 <- c(threshold = coordK2, auc = aucK2)

            return(statK2)
        }, classMax, class, currentFeature)))

    # Get Min Value
    aucMin <- min(aucFram$auc)

    # Get indices where this AUC occurs
    aucMinIndices <- aucFram$auc == aucMin

    # Use maximum value if there are ties
    aucValue <- max(aucFram[, "threshold"])

    # Return performance or value?
    if (rPerf) {
        return(aucMin)
    } else {
        return(aucValue)
    }
}


# Run modified F1 metric on single feature
.splitMetricModF1 <- function(feat, class, features, rPerf = FALSE) {

    # Get number of samples
    len <- length(class)

    # Get Values
    featValues <- features[, feat]

    # Get order of values
    ord <- order(featValues, decreasing = TRUE)

    # Get sorted class and values
    featValuesSort <- featValues[ord]
    classSort <- class[ord]

    # Keep splits of the data where the class changes
    keep <- c(
        classSort[seq(1, (len - 1))] != classSort[seq(2, (len))] &
            featValuesSort[seq(1, (len - 1))] != featValuesSort[seq(2, (len))],
        FALSE)

    # Create data.matrix
    X <- model.matrix(~ 0 + classSort)

    # Get cumulative sums
    sRCounts <- apply(X, 2, cumsum)

    # Keep only values where the class changes
    sRCounts <- sRCounts[keep, , drop = FALSE]
    featValuesKeep <- featValuesSort[keep]

    # Number of each class
    Xsum <- colSums(X)

    # Remove impossible splits (No class has > 50% of there samples on one side)
    sRProbs <- sRCounts %*% diag(Xsum^-1)
    sKeepPossible <- rowSums(sRProbs >= 0.5) > 0 & rowSums(sRProbs < 0.5) > 0

    # Remove anything after a full prob (Doesn't always happen)
    maxCheck <- min(c(which(apply(sRProbs, 1, max) == 1), nrow(sRProbs)))
    sKeepCheck <- seq(1, nrow(sRProbs)) %in% seq(1, maxCheck)

    # Combine logical vectors
    sKeep <- sKeepPossible & sKeepCheck

    if (sum(sKeep) > 0) {

        # Remove these if they exist
        sRCounts <- sRCounts[sKeep, , drop = FALSE]
        featValuesKeep <- featValuesKeep[sKeep]

        # Get left counts
        sLCounts <- t(Xsum - t(sRCounts))

        # Calculate the harmonic mean of Sens, Prec, and Worst Alt Sens
        statModF1 <- vapply(
            seq(nrow(sRCounts)),
            function(i, Xsum, sRCounts, sLCounts) {

                # Right Side
                sRRowSens <- sRCounts[i, ] / Xsum # Right sensitivities
                sRRowPrec <- sRCounts[i, ] / sum(sRCounts[i, ]) # Right prec
                sRRowF1 <- 2 * (sRRowSens * sRRowPrec) / (sRRowSens + sRRowPrec)
                sRRowF1[is.nan(sRRowF1)] <- 0 # Get right F1
                bestF1Ind <- which.max(sRRowF1) # Which is the best?
                bestSens <- sRRowSens[bestF1Ind] # The corresponding sensitivity
                bestPrec <- sRRowPrec[bestF1Ind] # The corresponding precision

                # Left Side
                sLRowSens <- sLCounts[i, ] / Xsum # Get left sensitivities
                worstSens <- min(sLRowSens[-bestF1Ind]) # Get the worst

                # Get harmonic mean of best sens, best prec, and worst sens
                HMout <- (3 * bestSens * bestPrec * worstSens) /
                (bestSens * bestPrec + bestPrec * worstSens +
                    bestSens * worstSens)

                return(HMout)
            }, Xsum, sRCounts, sLCounts, FUN.VALUE = double(1))

            # Get Max Value
            ModF1Max <- max(statModF1)

            # Get indices where this value occurs (use minimum row)
            ModF1Index <- which.max(statModF1)

            # Get value at this point
            ValueCeiling <- featValuesKeep[ModF1Index]
            ValueWhich <- which(featValuesSort == ValueCeiling)
            ModF1Value <- mean(
                c(featValuesSort[ValueWhich], featValuesSort[ValueWhich + 1]))
    } else {
        ModF1Max <- 0
        ModF1Value <- NA
    }

    if (rPerf) {
        return(ModF1Max)
    } else {
        return(ModF1Value)
    }
}

# Run Information Gain (probability + density) on a single feature
.splitMetricIGpIGd <- function(feat, class, features, rPerf = FALSE) {

    # Get number of samples
    len <- length(class)

    # Get Values
    featValues <- features[, feat]

    # Get order of values
    ord <- order(featValues, decreasing = TRUE)

    # Get sorted class and values
    featValuesSort <- featValues[ord]
    classSort <- class[ord]

    # Keep splits of the data where the class changes
    keep <- c(
        classSort[seq(1, (len - 1))] != classSort[seq(2, (len))] &
            featValuesSort[seq(1, (len - 1))] != featValuesSort[seq(2, (len))],
        FALSE)

    # Create data.matrix
    X <- model.matrix(~ 0 + classSort)

    # Get cumulative sums
    sRCounts <- apply(X, 2, cumsum)

    # Keep only values where the class changes
    sRCounts <- sRCounts[keep, , drop = FALSE]
    featValuesKeep <- featValuesSort[keep]

    # Number of each class
    Xsum <- colSums(X)

    # Remove impossible splits
    sRProbs <- sRCounts %*% diag(Xsum^-1)
    sKeep <- rowSums(sRProbs >= 0.5) > 0 & rowSums(sRProbs < 0.5) > 0

    if (sum(sKeep) > 0) {

        # Remove these if they exist
        sRCounts <- sRCounts[sKeep, , drop = FALSE]
        featValuesKeep <- featValuesKeep[sKeep]

        # Get left counts
        sLCounts <- t(Xsum - t(sRCounts))

        # Multiply them to get probabilities
        sRProbs <- t(t(sRCounts) %*%
            diag(rowSums(sRCounts)^-1, nrow = nrow(sRCounts)))
        sLProbs <- t(t(sLCounts) %*%
            diag(rowSums(sLCounts)^-1, nrow = nrow(sLCounts)))

        # Multiply them by there log
        sRTrans <- sRProbs * log(sRProbs)
        sRTrans[is.na(sRTrans)] <- 0
        sLTrans <- sLProbs * log(sLProbs)
        sLTrans[is.na(sLTrans)] <- 0

        # Get entropies
        HSR <- -rowSums(sRTrans)
        HSL <- -rowSums(sLTrans)

        # Get overall probabilities and entropy
        nProbs <- colSums(X) / len
        HS <- -sum(nProbs * log(nProbs))

        # Get split proporions
        sProps <- rowSums(sRCounts) / nrow(X)

        # Get information gain (Probability)
        IGprobs <- HS - (sProps * HSR + (1 - sProps) * HSL)
        IGprobs[is.nan(IGprobs)] <- 0
        IGprobsQuantile <- IGprobs / max(IGprobs)
        IGprobsQuantile[is.nan(IGprobsQuantile)] <- 0

        # Get proportions at each split
        classProps <- sRCounts %*% diag(Xsum^-1)
        classSplit <- classProps >= 0.5

        # Initialize information gain density vector
        splitIGdensQuantile <- rep(0, nrow(classSplit))

        # Get unique splits of the data
        classSplitUnique <- unique(classSplit)
        classSplitUnique <- classSplitUnique[!rowSums(classSplitUnique) %in%
            c(0, ncol(classSplitUnique)), , drop = FALSE]

        # Get density information gain
        if (nrow(classSplitUnique) > 0) {

            # Get log(determinant of full matrix)
            DET <- .psdet(stats::cov(features))

            # Information gain of every observation
            IGdens <- apply(
                classSplitUnique,
                1,
                .infoGainDensity,
                X,
                features,
                DET)

            names(IGdens) <- apply(
                classSplitUnique * 1,
                1,
                function(X) {
                    paste(X, collapse = "")
                })

            IGdens[is.nan(IGdens) | IGdens < 0] <- 0
            IGdensQuantile <- IGdens / max(IGdens)
            IGdensQuantile[is.nan(IGdensQuantile)] <- 0

            # Get ID of each class split
            splitsIDs <- apply(
                classSplit * 1,
                1,
                function(x) {
                    paste(x, collapse = "")
                })

            # Append information gain density vector
            for (ID in names(IGdens)) {
                splitIGdensQuantile[splitsIDs == ID] <- IGdensQuantile[ID]
            }
        }

        # Add this to the other matrix
        IG <- IGprobsQuantile + splitIGdensQuantile

        # Get IG(probabilty) of maximum value
        IGreturn <- IGprobs[which.max(IG)[1]]

        # Get maximum value
        maxVal <- featValuesKeep[which.max(IG)]
        wMax <- max(which(featValuesSort == maxVal))
        IGvalue <- mean(c(featValuesSort[wMax], featValuesSort[wMax + 1]))
    } else {
        IGreturn <- 0
        IGvalue <- NA
    }

    # Report maximum ID or value at maximum IG
    if (rPerf) {
        return(IGreturn)
    } else {
        return(IGvalue)
    }
}

# Function to find pseudo-determinant
.psdet <- function(x) {
    if (sum(is.na(x)) == 0) {
        svalues <- zapsmall(svd(x)$d)
        sum(log(svalues[svalues > 0]))
    } else {
        0
    }
}

# Function to calculate density information gain
.infoGainDensity <- function(splitVector, X, features, DET) {

    # Get Subsets of the feature matrix
    sRFeat <- features[as.logical(
        rowSums(X[, splitVector, drop = F])), , drop = F]
    sLFeat <- features[as.logical(
        rowSums(X[, !splitVector, drop = F])), , drop = F]

    # Get pseudo-determinant of covariance matrices
    DETR <- .psdet(cov(sRFeat))
    DETL <- .psdet(cov(sLFeat))

    # Get relative sizes
    sJ <- nrow(features)
    sJR <- nrow(sRFeat)
    sJL <- nrow(sLFeat)

    IUout <- 0.5 * (DET - (sJR / sJ * DETR + sJL / sJ * DETL))

    return(IUout)
}

# Wrapper function for getting split statistics
.getSplit <- function(feat, splitStats, features, class, splitMetric) {
    stat <- splitStats[feat]
    splitVal <- splitMetric(feat, class, features, rPerf = FALSE)
    featValues <- features[, feat]

    # Get classes split to one node
    node1Class <- class[featValues > splitVal]

    # Get proportion of each class at each node
    group1Prop <- table(node1Class) / table(class)
    group2Prop <- 1 - group1Prop

    # Get class consensus
    group1Consensus <- names(group1Prop)[group1Prop >= 0.5]
    group2Consensus <- names(group1Prop)[group1Prop < 0.5]

    # Get group samples
    group1 <- rownames(features)[class %in% group1Consensus]
    group2 <- rownames(features)[class %in% group2Consensus]

    # Get class vector
    group1Class <- droplevels(class[class %in% group1Consensus])
    group2Class <- droplevels(class[class %in% group2Consensus])

    return(list(
        featureName = feat,
        value = splitVal,
        stat = stat,

        group1 = group1,
        group1Class = group1Class,
        group1Consensus = group1Consensus,
        group1Prop = c(group1Prop),

        group2 = group2,
        group2Class = group2Class,
        group2Consensus = group2Consensus,
        group2Prop = c(group2Prop)
    ))
}

# Function to annotate alternate split of a soley downregulated terminal nodes
.addAlternativeSplit <- function(tree, features, class) {

    # Unlist decsision decision tree
    DecTree <- unlist(tree, recursive = F)

    # Get leaves
    groupList <- lapply(DecTree, function(split) {

        # Remove directions
        split <- split[!names(split) %in% c("statUsed", "fUsed", "dirs")]

        # Get groups
        group1 <- unique(unlist(lapply(
            split,
            function(node) {
                node$group1Consensus
            })))
        group2 <- unique(unlist(lapply(
            split,
            function(node) {
                node$group2Consensus
            })))

        return(list(
            group1 = group1,
            group2 = group2
        ))
    })

    # Get vector of each group
    group1Vec <- unique(unlist(lapply(groupList, function(g) g$group1)))
    group2Vec <- unique(unlist(lapply(groupList, function(g) g$group2)))

    # Get group that is never up-regulated
    group2only <- group2Vec[!group2Vec %in% group1Vec]

    # Check whether there are solely downregulated splits
    AltSplitInd <- which(unlist(lapply(groupList, function(g, group2only) {
        group2only %in% g$group2
        }, group2only)))

    if (length(AltSplitInd) > 0) {

        AltDec <- max(which(unlist(lapply(groupList, function(g, group2only) {
            group2only %in% g$group2
        }, group2only))))

        # Get split
        downSplit <- DecTree[[AltDec]]
        downNode <- downSplit[[1]]

        # Get classes to rerun
        branchClasses <- names(downNode$group1Prop)

        # Get samples from these classes and features from this cluster
        sampKeep <- class %in% branchClasses
        featKeep <- !colnames(features) %in% downSplit$fUsed

        # Subset class and features
        cSub <- droplevels(class[sampKeep])
        fSub <- features[sampKeep, featKeep, drop = F]

        # Get best alternative split
        altStats <- do.call(rbind, lapply(
            colnames(fSub),
            function(feat, splitMetric, features, class, cInt) {
                Val <- splitMetric(feat, cSub, fSub, rPerf = F)

                # Get node1 classes
                node1Class <- class[features[, feat] > Val]

                # Get sensitivity/precision/altSens
                Sens <- sum(node1Class == cInt) / sum(class == cInt)
                Prec <- mean(node1Class == cInt)

                # Get Sensitivity of Alternate Classes
                AltClasses <- unique(class)[unique(class) != cInt]
                AltSizes <- vapply(
                    AltClasses,
                    function(cAlt, class) {
                        sum(class == cAlt)
                    }, class, FUN.VALUE = double(1))
                AltWrong <- vapply(
                    AltClasses,
                    function(cAlt, node1Class) {
                        sum(node1Class == cAlt)
                    }, node1Class, FUN.VALUE = double(1))
                AltSens <- min(1 - (AltWrong / AltSizes))

                # Get harmonic mean
                HM <- (3 * Sens * Prec * AltSens) /
                    (Sens * Prec + Prec * AltSens + Sens * AltSens)
                HM[is.nan(HM)] <- 0

                # Return
                return(data.frame(
                    feat = feat,
                    val = Val,
                    stat = HM,
                    stringsAsFactors = F))
            }, .splitMetricModF1, fSub, cSub, group2only))
        altStats <- altStats[order(altStats$stat, decreasing = TRUE), ]

        # Get alternative splits
        splitStats <- altStats$stat[1]
        names(splitStats) <- altStats$feat[1]
        altSplit <- .getSplit(
            altStats$feat[1],
            splitStats,
            fSub,
            cSub,
            .splitMetricModF1)

        # Check that this split out the group2 of interest
        if (length(altSplit$group1Consensus) == 1) {

            # Add it to split
            downSplit[[length(downSplit) + 1]] <- altSplit
            names(downSplit)[length(downSplit)] <- paste0(altStats$feat[1], "+")
            downSplit <- downSplit[c(
                which(!names(downSplit) %in% c("statUsed", "fUsed", "dirs")),
                which(names(downSplit) %in% c("statUsed", "fUsed", "dirs")))]

            # Get index of split to add it to
            branchLengths <- unlist(lapply(tree, length))
            branchCum <- cumsum(branchLengths)
            wBranch <- min(which(branchCum >= AltDec))
            wSplit <- which(seq(
                (branchCum[(wBranch - 1)] + 1),
                    branchCum[wBranch]) == AltDec)

            # Add it to decision tree
            tree[[wBranch]][[wSplit]] <- downSplit
        } else {
            cat("No non-ambiguous rule to separate", group2only, "from",
                branchClasses, ". No alternative split added.")
        }
    } else {
        print("No solely down-regulated cluster to add alternative split.")
    }

    return(tree)
}
