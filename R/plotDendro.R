#' @title Plots dendrogram of `buildTreeHybrid` output
#' @description Generates a dendrogram of the rules and performance
#' (optional) of the decision tree generates by `buildTreeHybrid`.
#' @param counts List object. The output of `celda::buildTreeHybrid`.
#' @param classLabel A character value. The name of a label to which draw
#'  the path and rules. If NULL (default), the rules for every cluster is shown.
#' @param addSensPrec Logical. Print training sensitivities and precisions
#'  for each cluster below leaf label? Default is TRUE.
#' @param leafSize A numeric value. Size of text below each leaf. Default is 24.
#' @param boxSize A numeric value. Size of rule labels. Default is 7.
#' @param boxColor A character value. Color of rule labels. Default is `black`.
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' plotHeatmap(celdaCGSim$counts,
#'     z = clusters(celdaCGMod)$z, y = clusters(celdaCGMod)$y)
#' @return A ggplot2 object
#' @import ggplot2
#' @import ggdendro
#' @import dendextend
#' @export

plotDendro <- function(decisionTree,
                    classLabel = NULL,
                    addSensPrec = TRUE,
                    leafSize = 24,
                    boxSize = 7,
                    boxColor = "black"
) {

    # Get necessary elements
    dendro <- decisionTree$dendro
    
    # Get performance information (training or CV based)
    performance <- decisionTree$performance

    # Create vector of per class performance
    perfVec <- paste(performance$sizes,
        format(round(performance$sensitivity, 2), nsmall = 2),
        format(round(performance$precision, 2), nsmall = 2),
        sep = "\n"
    )
    names(perfVec) <- names(performance$sensitivity)

    # Get dendrogram segments
    dendSegs <- dendro_data(dendro, type = "rectangle")$segments

    # Get necessary coordinates to add labels to
    # These will have y > 1
    dendSegs <- unique(dendSegs[dendSegs$y > 1, c("x", "y", "yend", "xend")])
    
    # Labeled splits will be vertical (x != xend) or
    # Length 0 (x == xend & y == yend)
    dendSegsAlt <- dendSegs[
        dendSegs$x != dendSegs$xend |
            (dendSegs$x == dendSegs$xend & dendSegs$y == dendSegs$yend), 
        c("x", "xend", "y")]
    colnames(dendSegsAlt)[1] <- "xalt"

    # Label names will be at nodes, these will
    # Occur at the end of segments
    segs <- as.data.frame(get_nodes_xy(dendro))
    colnames(segs) <- c("xend", "yend")
    
    # As label and which stat was used
    # Labels will stack
    segs$label <- gsub(";", "\n", get_nodes_attr(dendro, "label"))
    segs$statUsed <- get_nodes_attr(dendro, "statUsed")
    
    # If highlighting a class label, remove non-class specific rules
    if (!is.null(classLabel)) {
        if (!classLabel %in% rownames(decisionTree$summaryMatrix)) {
            stop("classLabel not a valid class ID.")
        }
        dendro <- .highlightClassLabel(dendro, classLabel)
        keepLabel <- get_nodes_attr(dendro, "keepLabel")
        keepLabel[is.na(keepLabel)] <- FALSE
        segs$label[!keepLabel] <- NA
    }
    
    # Remove non-labelled nodes &
    # leaf nodes (yend == 0)
    segs <- segs[!is.na(segs$label) & segs$yend != 0, ]

    # Merge to full set of coordinates
    dendSegsLabelled <- merge(dendSegs, segs)
    
    # Remove duplicated labels
    dendSegsLabelled <- dendSegsLabelled[order(dendSegsLabelled$y, 
        decreasing = T),]
    dendSegsLabelled <- dendSegsLabelled[
       !duplicated(dendSegsLabelled[,c("xend", "x", "yend", "label", "statUsed")]),]
    
    # Merge with alternative x-coordinates for alternative split
    dendSegsLabelled <- merge(dendSegsLabelled, dendSegsAlt)

    # Order by height and coordinates
    dendSegsLabelled <- dendSegsLabelled[order(
        dendSegsLabelled$x),]

    # Find information gain splits
    igSplits <- dendSegsLabelled$statUsed == "IG" &
        !duplicated(dendSegsLabelled[, c("xalt", "y")])

    # Set xend for IG splits
    dendSegsLabelled$xend[igSplits] <- dendSegsLabelled$xalt[igSplits]
    
    # Set y for non-IG splits
    dendSegsLabelled$y[!igSplits] <- dendSegsLabelled$y[!igSplits] - 0.2

    # Get index of leaf labels
    leafLabels <- get_leaves_attr(dendro, "label")

    # Add sensitivity and precision measurements
    if (addSensPrec) {
        leafLabels <- paste(leafLabels, perfVec[leafLabels], sep = "\n")
    }

    # Create plot of dendrogram
    suppressMessages(dendroP <- ggdendrogram(dendro) +
        geom_label(
            data = dendSegsLabelled,
            aes(
                x = xend,
                y = y,
                label = label),
            size = boxSize,
            label.size = 1,
            fontface = "bold",
            vjust = 1,
            nudge_y = 0.1,
            color = boxColor) +
        theme_bw() +
        scale_x_reverse(breaks = seq(length(leafLabels)),
            label = leafLabels) +
        scale_y_continuous(expand = c(0, 0)) +
        theme(
            panel.grid.major.y = element_blank(),
            legend.position = "none",
            panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.border = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_text(hjust = 0.5,
                size = leafSize,
                family = "mono",
                vjust = -1),
            axis.text.y = element_blank()
        ))

    # Increase line width slightly for aesthetic purposes
    dendroP$layers[[2]]$aes_params$size <- 1.3

    return(dendroP)
}

# Function to reformat the dendrogram to draw path to a specific class
.highlightClassLabel <- function(dendro, classLabel) {

    # Reorder dendrogram
    flag <- TRUE
    bIndexString <- ""

    # Get branch
    branch <- eval(parse(text = paste0("dendro", bIndexString)))

    while (flag) {

        # Get attributes
        att <- attributes(branch)

        # Get split with the label of interest
        labList <- lapply(branch, function(split) attributes(split)$classLabels)
        wSplit <- which(unlist(lapply(
            labList,
            function(vec) {
                classLabel %in% vec
            })))

        # Keep labels for this branch
        branch <- lapply(branch, function(edge) {
            attributes(edge)$keepLabel <- TRUE
            return(edge)
        })

        # Make a dendrogram class again
        class(branch) <- "dendrogram"
        attributes(branch) <- att

        # Add branch to dendro
        eval(parse(text = paste0("dendro", bIndexString, "<- branch")))

        # Create new bIndexString
        bIndexString <- paste0(bIndexString, "[[", wSplit, "]]")

        # Get branch
        branch <- eval(parse(text = paste0("dendro", bIndexString)))

        # Add flag
        flag <- attributes(branch)$members > 1
    }

    return(dendro)
}
