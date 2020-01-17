#' @title Plots dendrogram of `findMarkers` output
#' @description Generates a dendrogram of the rules and performance
#' (optional) of the decision tree generates by `findMarkers`.
#' @param tree List object. The output of `celda::findMarkers`.
#' @param classLabel A character value. The name of a label to which draw
#'  the path and rules. If NULL (default), the rules for every cluster is shown.
#' @param addSensPrec Logical. Print training sensitivities and precisions
#'  for each cluster below leaf label? Default is FALSE.
#' @param maxFeaturePrint A numeric value. Maximum number of feature IDs to print
#'  at a given node. Default is 4.  
#' @param leafSize A numeric value. Size of text below each leaf. Default is 24.
#' @param boxSize A numeric value. Size of rule labels. Default is 7.
#' @param boxColor A character value. Color of rule labels. Default is `black`.
#' @examples
#' # Generate simulated single-cell dataset using celda 
#' sim_counts <- celda::simulateCells("celda_CG", K = 4, L = 10, G = 100)
#' # Celda clustering into 5 clusters & 10 modules
#' cm <- celda_CG(sim_counts$counts, K=4, L=10, verbose=FALSE)
#' # Get features matrix and cluster assignments
#' factorized <- factorizeMatrix(sim_counts$counts, cm)
#' features <- factorized$proportions$cell
#' class <- clusters(cm)$z
#' # Generate Decision Tree
#' DecTree <- findMarkersTree(features,class,threshold = 1)
#' # Plot dendrogram
#' plotDendro(DecTree)
#' 
#' @return A ggplot2 object
#' @import ggplot2
#' @importFrom ggdendro dendro_data ggdendrogram
#' @importFrom dendextend get_nodes_xy get_nodes_attr get_leaves_attr
#' @export
plotDendro <- function(tree,
                       classLabel = NULL,
                       addSensPrec = FALSE,
                       maxFeaturePrint = 4,
                       leafSize = 10,
                       boxSize = 2,
                       boxColor = "black") {
    
    # Get necessary elements
    dendro <- tree$dendro
    
    # Get performance information (training or CV based)
    if(addSensPrec){
        performance <- tree$performance
        
        # Create vector of per class performance
        perfVec <- paste(performance$sizes,
                         format(round(performance$sensitivity, 2), nsmall = 2),
                         format(round(performance$precision, 2), nsmall = 2),
                         sep = "\n"
        )
        names(perfVec) <- names(performance$sensitivity)
    }
    
    # Get dendrogram segments
    dendSegs <- ggdendro::dendro_data(dendro, type = "rectangle")$segments
    
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
    segs <- as.data.frame(dendextend::get_nodes_xy(dendro))
    colnames(segs) <- c("xend", "yend")
    
    # Add labels to nodes
    segs$label <- gsub(";", "\n", dendextend::get_nodes_attr(dendro, "label"))
    
    # Subset for max
    segs$label <- sapply(segs$label, function(lab, maxFeaturePrint) {
        loc <- gregexpr("\n", lab)[[1]][maxFeaturePrint]
        if(!is.na(loc)) {
            lab <- substr(lab, 1, loc-1)
        }
        return(lab)
    }, maxFeaturePrint)
    
    segs$statUsed <- dendextend::get_nodes_attr(dendro, "statUsed")
    
    # If highlighting a class label, remove non-class specific rules
    if (!is.null(classLabel)) {
        if (!classLabel %in% names(tree$rules)) {
            stop("classLabel not a valid class ID.")
        }
        dendro <- .highlightClassLabel(dendro, classLabel)
        keepLabel <- dendextend::get_nodes_attr(dendro, "keepLabel")
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
                                               decreasing = T), ]
    dendSegsLabelled <- dendSegsLabelled[
        !duplicated(dendSegsLabelled[,
                                     c("xend", "x", "yend",
                                       "label", "statUsed")]), ]
    
    # Merge with alternative x-coordinates for alternative split
    dendSegsLabelled <- merge(dendSegsLabelled, dendSegsAlt)
    
    # Order by height and coordinates
    dendSegsLabelled <- dendSegsLabelled[order(dendSegsLabelled$x), ]
    
    # Find information gain splits
    igSplits <- dendSegsLabelled$statUsed == "Split" &
        !duplicated(dendSegsLabelled[, c("xalt", "y")])
    
    # Set xend for IG splits
    dendSegsLabelled$xend[igSplits] <- dendSegsLabelled$xalt[igSplits]
    
    # Set y for non-IG splits
    dendSegsLabelled$y[!igSplits] <- dendSegsLabelled$y[!igSplits] - 0.2
    
    # Get index of leaf labels
    leafLabels <- dendextend::get_leaves_attr(dendro, "label")
    
    # Adjust leaf labels if there are metacluster labels
    if(!is.null(tree$metaclusterLabels)){
        leafLabels <- regmatches(leafLabels,
                                 regexpr(pattern = "(?<=\\().*?(?=\\)$)",
                                         leafLabels, perl = TRUE))
    }
    
    # Add sensitivity and precision measurements
    if (addSensPrec) {
        leafLabels <- paste(leafLabels, perfVec[leafLabels], sep = "\n")
        leafAngle <- 0
        leafHJust <- 0.5
        leafVJust <- -1
    } else {
        leafAngle <- 90
        leafHJust <- 1
        leafVJust <- 0.5
    }

    # Create plot of dendrogram
    suppressMessages(dendroP <- ggdendro::ggdendrogram(dendro) +
                         ggplot2::geom_label(
                             data = dendSegsLabelled,
                             ggplot2::aes(x = xend, y = y, label = label),
                             size = boxSize,
                             label.size = 1,
                             fontface = "bold",
                             vjust = 1,
                             nudge_y = 0.1,
                             color = boxColor) +
                         ggplot2::theme_bw() +
                         ggplot2::scale_x_reverse(breaks =
                                                      seq(length(leafLabels)),
                                                  label = leafLabels) +
                         ggplot2::scale_y_continuous(expand = c(0, 0)) +
                         ggplot2::theme(
                             panel.grid.major.y = ggplot2::element_blank(),
                             legend.position = "none",
                             panel.grid.minor.y = ggplot2::element_blank(),
                             panel.grid.minor.x = ggplot2::element_blank(),
                             panel.grid.major.x = ggplot2::element_blank(),
                             panel.border = ggplot2::element_blank(),
                             axis.title = ggplot2::element_blank(),
                             axis.ticks = ggplot2::element_blank(),
                             axis.text.x = ggplot2::element_text(
                                 hjust = leafHJust,
                                 angle = leafAngle,
                                 size = leafSize,
                                 family = "Palatino", 
                                 face = "bold",
                                 vjust = leafVJust),
                             axis.text.y = ggplot2::element_blank()
                         )
                     )
    
    # Check if need to add metacluster labels
    if(!is.null(tree$metaclusterLabels)){
        #store metacluster labels to add
        newLabels <- unique(tree$branchPoints$top_level$metacluster)
        
        #adjust labels for metaclusters of size one
        newLabels <- unlist(lapply(newLabels, function(curMeta){
            if(substr(curMeta, nchar(curMeta), nchar(curMeta)) == ")"){
                return(gsub(pattern = "\\(.*\\)$", replacement = "",
                            x = curMeta))
            }
            else{
                return(curMeta)
            }
        }))
        
        # Create table for metacluster labels
        metaclusterText <- dendSegsLabelled[dendSegsLabelled$y ==
                                                max(dendSegsLabelled$y), 
                                            c("xend", "y", "label")]
        metaclusterText$label <- newLabels
        
        # Add metacluster labels to top of plot
        dendroP <- dendroP +
            geom_text(data = metaclusterText, 
                      aes(x = xend, y = y,
                          label = label, fontface=2), 
                      angle = 90, 
                      nudge_y = 0.5,
                      family = "Palatino", 
                      size = leafSize/3) 
        
        #adjust coordinates of plot to show labels
        dendroP <- dendroP + coord_cartesian(ylim = 
                                                 c(0, 
                                                   max(dendSegsLabelled$y+1)))
        
    }
    
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
        labList <- lapply(branch, function(split)
            attributes(split)$classLabels)
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
