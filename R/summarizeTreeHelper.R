# Function to summarize and format tree list output by .generateTreeList
.summarizeTree <- function(tree, features, class) {

# Format tree into dendrogram object
    dendro <- .convertToDendrogram(tree, class)

    # Map classes to features
    class2features <- .mapClass2features(tree, features, class)

    # Get performance of the tree on training samples
    perfList <- .getPerformance(class2features$rules, features, class)

    return(list(
        rules = class2features$rules,
        dendro = dendro,
        summaryMatrix = class2features$c2fMatrix,
        prediction = perfList$prediction,
        performance = perfList$performance
    ))
}

# Function to reformat raw tree ouput to a dendrogram
.convertToDendrogram <- function(tree, class) {

    # Unlist decision tree (one element for each split)
    DecTree <- unlist(tree, recursive = F)

    # Name split by gene and threshold
    splitNames <- lapply(DecTree, function(split) {

        # Remove non-split elements
        dirs <- paste0(split$dirs, collapse = "_")
        split <- split[!names(split) %in% c("statUsed", "fUsed", "dirs")]

        # Get set of features and values for each
        featuresplits <- lapply(split, function(node) {
            nodeFeature <- node$featureName
            nodeStrings <- paste(nodeFeature, collapse = ";")
        })

        # Get split directions
        names(featuresplits) <- paste(
            dirs,
            seq(length(featuresplits)),
            sep = "_")

        return(featuresplits)
    })
    splitNames <- unlist(splitNames)
    names(splitNames) <- sub("1_", "", names(splitNames))

    # Get Stat Used
    statUsed <- unlist(lapply(
        DecTree, function(split) {
        split$statUsed
    }))
    statRep <- unlist(lapply(
        DecTree,
        function(split) {
            length(split[!names(split) %in% c("statUsed", "fUsed", "dirs")])
        }))
    statUsed <- unlist(lapply(
        seq(length(statUsed)),
        function(i) {
            rep(statUsed[i], statRep[i])
        }))
    names(statUsed) <- names(splitNames)

    # Create Matrix of results
    mat <- matrix(0, nrow = length(DecTree), ncol = length(unique(class)))
    colnames(mat) <- unique(class)
    for (i in seq(1, length(DecTree))) {

        # If only one split than ezpz
        split <- DecTree[[i]]
        split <- split[!names(split) %in% c("statUsed", "fUsed", "dirs")]
        if (length(split) == 1) {
            mat[i, split[[1]]$group1Consensus] <- 1
            mat[i, split[[1]]$group2Consensus] <- 2

            # Otherwise we need to assign > 2 splits for different higher groups
        } else {
            # Get classes in group 1
            group1classUnique <- unique(lapply(
                 split,
                 function(X) {
                    X$group1Consensus
                }))
           group1classVec <- unlist(group1classUnique)

            # Get classes always in group 2
            group2classUnique <- unique(unlist(lapply(
                split,
                function(X) {
                    X$group2Consensus
                })))
            group2classUnique <- group2classUnique[!group2classUnique %in%
                group1classVec ]

            # Assign
            for (j in seq(length(group1classUnique))) {
                mat[i, group1classUnique[[j]]] <- j
            }
            mat[i, group2classUnique ] <- j + 1
        }
    }

    ## Collapse matrix to get set of direction to include in dendrogram
    matCollapse <- sort(apply(
    mat,
        2,
        function(x) {
            paste(x[x != 0], collapse = "_")
        }))
    matUnique <- unique(matCollapse)

    # Get branchlist
    bList <- c()
    j <- 1
    for (i in seq(max(ncharX(matUnique)))) {
        sLength <- matUnique[ncharX(matUnique) >= i]
        sLength <- unique(subUnderscore(sLength, i))
        for (k in sLength) {
         bList[j] <- k
         j <- j + 1
        }
    }

    # Initialize dendrogram list
    val <- max(ncharX(matUnique)) + 1
    dendro <- list()
    attributes(dendro) <- list(
        members = length(matCollapse),
        classLabels = unique(class),
        height = val,
        midpoint = (length(matCollapse) - 1) / 2,
        label = NULL,
        name = NULL)

    for (i in bList) {

        # Add element
        iSplit <- unlist(strsplit(i, "_"))
        iPaste <- paste0("dendro",
            paste(paste0("[[", iSplit, "]]"), collapse = ""))
            eval(parse(
                text =
                paste0(iPaste, "<-list()")
            ))

        # Add attributes
        classLabels <- names(
        matCollapse[subUnderscore(matCollapse, ncharX(i)) == i])
        members <- length(classLabels)

        # Add height, set to one if leaf
        height <- val - ncharX(i)

        # Check that this isn't a terminal split
        if (members == 1) {
            height <- 1
        }

        # Add labels and stat used
        if (i %in% names(splitNames)) {
            lab <- splitNames[i]
            statUsedI <- statUsed[i]
        } else {
            lab <- NULL
            statUsedI <- NULL
        }
        att <- list(
            members = members,
            classLabels = classLabels,
            edgetext = lab,
            height = height,
            midpoint = (members - 1) / 2,
            label = lab,
            statUsed = statUsedI,
            name = i)
        eval(parse(
                text = paste0("attributes(", iPaste, ") <- att")))

        # Add leaves
        leaves <- matCollapse[matCollapse == i]
        if (length(leaves) > 0) {
            for (l in seq(1, length(leaves))) {

                # Add element
                lPaste <- paste0(iPaste, "[[", l, "]]")
                eval(parse(
                    text = paste0(lPaste, "<-list()")))

                # Add attributes
                members <- 1
                leaf <- names(leaves)[l]
                height <- 0
                att <- list(
                    members = members,
                    classLabels = leaf,
                    height = height,
                    label = leaf,
                    leaf = TRUE,
                    name = i)
                eval(parse(
                    text = paste0("attributes(", lPaste, ") <- att")))
            }
        }
    }
    class(dendro) <- "dendrogram"
    return(dendro)
}

# Function to calculate the number of non-underscore characters in a string
ncharX <- function(x) unlist(lapply(strsplit(x, "_"), length))

# Function to subset a string of characters seperated by underscores
subUnderscore <- function(x, n) unlist(lapply(
    strsplit(x, "_"),
    function(y) {
        paste(y[seq(n)], collapse = "_")
    }))

# Function to calculate performance statistics
.getPerformance <- function(rules, features, class) {

    # Get classification accuracy, balanced accurecy, and per class sensitivity
    ## Get predictions
    votes <- getDecisions(rules, t(features))
    votes[is.na(votes)] <- "MISSING"

    ## Calculate accuracy statistics and per class sensitivity
    class <- as.character(class)
    acc <- mean(votes == as.character(class))
    classCorrect <- vapply(
        unique(class),
        function(x) {
            sum(votes == x & class == x)
        }, FUN.VALUE = double(1))
    classCount <- c(table(class))[unique(class)]
    sens <- classCorrect / classCount

    ## Calculate balanced accuracy
    balacc <- mean(sens)

    ## Calculate per class and mean precision
    voteCount <- c(table(votes))[unique(class)]
    prec <- classCorrect / voteCount
    meanPrecision <- mean(prec)

    ## Add performance metrics
    performance <- list(
        accuracy = acc,
        balAcc = balacc,
        meanPrecision = meanPrecision,
        correct = classCorrect,
        sizes = classCount,
        sensitivity = sens,
        precision = prec
    )

    return(list(
        prediction = votes,
        performance = performance
    ))
}

# Create matrix of classes and features combinations
.mapClass2features <- function(tree, features, class) {

    # Create empty matrix
    c2fMatrix <- matrix(0, nrow = length(unique(class)), ncol = ncol(features))
    rownames(c2fMatrix) <- sort(unique(class))
    colnames(c2fMatrix) <- colnames(features)

    # Get class to feature indices
    class2featuresIndices <- do.call(rbind, lapply(
        seq(length(tree)),
        function(i) {
            treeLevel <- tree[[i]]
            c2fsub <- as.data.frame(do.call(rbind, lapply(
                treeLevel,
                function(split) {
                    # Keep track of stat used for rule list
                    statUsed <- split$statUsed

                    # Keep only split information
                    split <- split[!names(split) %in%
                    c("statUsed", "fUsed", "dirs")]

                    # Create data frame of split rules
                    edgeFram <- do.call(rbind, lapply(split, function(edge) {

                        # Create data.frame of groups, split-dirs, feature IDs
                        groups <- c(edge$group1Consensus, edge$group2Consensus)
                        sdir <- c(
                            rep(1, length(edge$group1Consensus)),
                            rep(-1, length(edge$group2Consensus)))
                            feat <- edge$featureName
                            val <- edge$value
                            stat <- edge$stat
                            data.frame(
                                class = rep(groups, length(feat)),
                                feature = rep(feat, each = length(groups)),
                                direction = rep(sdir, length(feat)),
                                value = rep(val, each = length(groups)),
                                stat = rep(stat, each = length(groups)),
                                stringsAsFactors = F
                            )
                        }))

                    # Add stat used
                    edgeFram$statUsed <- statUsed

                    return(edgeFram)

                })))
        c2fsub$level <- i
        return(c2fsub)
    }))
    rownames(class2featuresIndices) <- NULL

    # Add levels to matrix
    for (i in seq(nrow(class2featuresIndices))) {
        c2fMatrix[class2featuresIndices[i, "class"],
        class2featuresIndices[i, "feature"]] <-
            class2featuresIndices[i, "level"] *
                class2featuresIndices[i, "direction"]
    }

# Generate list of rules for each class
    rules <- lapply(levels(class), function(cl, class2featuresIndices) {

        class2featuresIndices[class2featuresIndices$class == cl,
        colnames(class2featuresIndices) != "class"]

    }, class2featuresIndices)
    names(rules) <- levels(class)

    return(list(
        c2fMatrix = c2fMatrix,
        rules = rules))
}
