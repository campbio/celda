#' @title Generate decision tree from single-cell clustering output.
#' @description Uses decision tree procudure to generate a set of rules for each
#'  cell cluster defined by a single-cell clustering.  Splits are determined by
#'  one of two metrics at each split: a one-off metric to determine rules for
#'  identifying clusters by a single feature, and a balanced metric to determine
#'  rules for identifying sets of similar clusters.
#' @param features A L(features) by N(samples) numeric matrix.
#' @param class A vector of K label assignemnts.
#' @param cellTypes List where each element is a cell type and all the clusters within that cell type (i.e. subtypes). 
#' @param oneoffMetric A character string. What one-off metric to run, either `modified F1` or  `pairwise AUC`.
#' @param threshold A numeric value. The threshold for the oneoff metric to use between 0 and 1, 0.95 by default. Smaller values will result is more one-off splits.
#' @param reuseFeatures Logical.  Whether or not a feature can be used more than once on the same cluster. Default is TRUE.
#' @param altSplit Logical. Whether or not to force a marker for clusters that are solely defined by the absence of markers. Defsult is TRUE
#' @param consecutiveOneoff Logical. Whether or not to allow one-off splits at consecutive brances. Default is TRUE
#' @return A named list with five elements.
#' \itemize{
#'   \item rules - A named list with one `data.frame` for every label. Each
#' `data.frame` has five columns and gives the set of rules for disinguishing each label.
#'   \itemize{
#'    \item feature - Feature identifier.
#'    \item direction - Relationship to feature value, -1 if less than, 1 if greater than.
#'    \item value - The feature value which defines the decision boundary
#'    \item stat - The performance value returned by the splitting metric for this split.
#'    \item statUsed - Which performance metric was used. "IG" if information gain and "OO" if one-off.
#'    \item level - The level of the tree at which is rule was defined. 1 is the level of the first split of the tree.
#'   }
#'  \item dendro - A dendrogram object of the decision tree output
#'  \item prediction - A character vector of label of predictions of the training data using the final model. "MISSING" if label prediction was ambiguous.
#'  \item performance - A named list denoting the training performance of the model.
#'  \itemize{
#'   \item accuracy - (number correct/number of samples) for the whole set of
#'  samples.
#'   \item balAcc - mean sensitivity across all labels
#'   \item meanPrecision - mean precision across all labels
#'   \item correct - the number of correct predictions of each label
#'   \item sizes - the number of actual counts of each label
#'   \item sensitivity - the sensitivity of the prediciton of each label.
#'   \item precision - the precision of the prediciton of each label.
#'  }
#' }
#' @examples
#' library(M3DExampleData)
#' counts <- M3DExampleData::Mmus_example_list$data
#' # Subset 500 genes for fast clustering
#' counts <- as.matrix(counts[1501:2000, ])
#' # Cluster genes ans samples each into 10 modules
#' cm <- celda_CG(counts = counts, L = 10, K = 5, verbose = FALSE)
#' # Get features matrix and cluster assignments
#' factorized <- factorizeMatrix(counts, cm)
#' features <- factorized$proportions$cell
#' class <- clusters(cm)$z
#' # Generate Decision Tree
#' decTree <- findMarkers(features,
#'     class,
#'     oneoffMetric = "modified F1",
#'     threshold = 1,
#'     consecutiveOneoff = FALSE)
#'
#' # Plot dendrogram
#' plotDendro(decTree)
#' @import magrittr
#' @export
findMarkers <- function(features,
                        class,
                        cellTypes,
                        oneoffMetric = c("modified F1", "pairwise AUC"),
                        threshold = 0.95,
                        reuseFeatures = FALSE,
                        altSplit = TRUE,
                        consecutiveOneoff = TRUE) {
  
  if (ncol(features) != length(class)) 
    stop("Number of columns of features must equal length of class")
  
  if (any(is.na(class))) 
    stop("NA class values")
  
  if (any(is.na(features))) 
    stop("NA feature values")
  
  # Match the oneoffMetric argument
  oneoffMetric <- match.arg(oneoffMetric)
  
  # Transpose features
  features <- t(features)
  
  #If no detailed cell types are provided
  if(!hasArg(cellTypes)){
    
    print('Building tree...')
    
    # Set class to factor
    class <- as.factor(class)
    
    # Generate list of tree levels
    tree <- .generateTreeList(
      features,
      class,
      oneoffMetric,
      threshold,
      reuseFeatures,
      consecutiveOneoff)
    
    # Add alternative node for the solely down-regulated leaf
    if (altSplit) 
      tree <- .addAlternativeSplit(tree, features, class)
    
    print('Computing performance metrics...')
    
    # Format tree output for plotting and generate summary statistics
    DTsummary <- .summarizeTree(tree, features, class)
    
    return(DTsummary)
  }
  
  # If detailed cell types are provided
  else{

    # Check that cell types match class labels
    if(mean(unlist(cellTypes) %in% unique(class)) != 1) 
      stop("Provided cell types and class labels do not match. Please check the 'cellTypes' argument.")
    
    # Create vector with cell type class labels
    newLabels <- class
    for (i in names(cellTypes)) {
        newLabels[newLabels %in% cellTypes[[i]]] <- i
    }

    # Update cell subtype labels
    subtypeLabels <- newLabels
    subtypeLabels[subtypeLabels %in% names(cellTypes)] <- paste0(
        subtypeLabels[subtypeLabels %in% names(cellTypes)],
        "(",
        class[subtypeLabels %in% names(cellTypes)],
        ")"
    )
    
    ## Create tree for cell types ##
    print('Building tree for all cell types...')
    tree <- .generateTreeList(features, as.factor(newLabels), oneoffMetric, threshold, reuseFeatures, consecutiveOneoff)
    tree <- list(
      rules = .mapClass2features(tree, features, as.factor(newLabels))$rules,
      dendro = .convertToDendrogram(tree, as.factor(newLabels))
    )
    
    #Store tree's dendrogram in a separate variable
    dendro <- tree$dendro
    
    #Find which cell types have more than one cluster
    largeCellTypes <- names(cellTypes[lengths(cellTypes) > 1])
    
    #Create separate trees for each cell type with more than one cluster
    newTrees <- lapply(unique(newLabels), function(cellType){
      
      #Create tree for this cell type only
      if (cellType %in% largeCellTypes) {
          
        # Print current status
        print(paste('Building tree for cell type:', cellType))
        
        # Remove used features
        featUse <- colnames(features)
        if (!reuseFeatures) {
          featUse <- featUse[!featUse %in% tree$rules[[cellType]]$feature]
        }
        
        # Create new Tree
        newTree <- .generateTreeList(features[newLabels == cellType, featUse], as.factor(subtypeLabels[newLabels == cellType]), oneoffMetric, threshold, reuseFeatures, consecutiveOneoff)
        newTree <- list(
          rules = .mapClass2features(newTree, features[newLabels == cellType,], as.factor(subtypeLabels[newLabels == cellType]))$rules,
          dendro = .convertToDendrogram(newTree, as.factor(subtypeLabels[newLabels == cellType]))
        )
        
        #Adjust 'rules' table for new tree
        newTree$rules <- lapply(newTree$rules, function(rules){
          rules$level <- rules$level + max(tree$rules[[cellType]]$level)
          rules <- rbind(tree$rules[[cellType]], rules)
        })
        
        return(newTree)
      }
      else{
          
        # Adjust name of class if it contains only one cluster
        names(tree$rules)[which(names(tree$rules) == cellType)] <<- paste0(cellType,'(', unlist(cellTypes[cellType]) ,')')
        return()
      }
    })
    #Remove empty trees for single-cluster cell types
    newTrees <- newTrees[lengths(newTrees) > 0]
    names(newTrees) <- unique(newLabels)[unique(newLabels) %in% largeCellTypes]
    
    #Find indices of cell type nodes in tree
    indices <- lapply(unique(newLabels)[unique(newLabels) %in% largeCellTypes], function(cellType) {
      # Initialize sub trees, indices string, and flag
      dendSub <- dendro
      index <- ""
      flag <- TRUE
      
      while (flag) {
        # Get the edge with the class of interest
        whEdge <- which(unlist(lapply(dendSub, function(edge) cellType %in% attributes(edge)$classLabels)))
        
        # Add this as a string
        index <- paste0(index, "[[", whEdge, "]]")
        
        # Move to this branch
        dendSub <- eval(parse(text = paste0("dendro", index)))
        
        # Is this the only class in that branch
        flag <- length(attributes(dendSub)$classLabels) > 1
      }
      
      return(index)
    })
    names(indices) <- unique(newLabels)[unique(newLabels) %in% largeCellTypes]
    
    #Add each cell type tree
    for(cellType in names(newTrees)){
      
      #Get current tree
      cellTypeDendro <- newTrees[[cellType]]$dendro
      
      #Nudge nodes upward to make room for new tree
      dendro <- dendrapply(dendro, function(node) {
        if(attributes(node)$height > 1)
          attributes(node)$height <- attributes(node)$height + attributes(cellTypeDendro)$height - 1
        return(node)
      })
      
      #Adjust labels, member count, and midpoint of nodes
      dendro <- dendrapply(dendro, function(node){
        #Check if in right branch
        if(cellType %in% as.character(attributes(node)$classLabels)){
          #assign new labels replacing old one (i.e. replace cell type label with subtype labels)
          attributes(node)$classLabels <- as.character(attributes(node)$classLabels) %>%
            .[. != cellType] %>%
            c(., unique(subtypeLabels)[grep(cellType, unique(subtypeLabels))])
          
          #assign new member count for this branch
          attributes(node)$members <- length(attributes(node)$classLabels)
          
          #assign new midpoint for this branch
          attributes(node)$midpoint <- (attributes(node)$members - 1) / 2
        }
        return(node)
      })
      
      #Replace label at new tree's branch point
      branchPointLabel <- attributes(eval(parse(text = paste0("dendro", indices[[cellType]]))))$label
      branchPointStatUsed <- attributes(eval(parse(text = paste0("dendro", indices[[cellType]]))))$statUsed
      if(!is.null(branchPointLabel)) {
        attributes(cellTypeDendro)$label <- branchPointLabel
        attributes(cellTypeDendro)$statUsed <- branchPointStatUsed
      }
      
      #Add new tree to original tree
      eval(parse(text = paste0("dendro", indices[[cellType]], " <- cellTypeDendro")))
      
      #Append new tree's 'rules' tables to original tree
      tree$rules <- append(tree$rules, newTrees[[cellType]]$rules)
      
      #Remove old tree's rules
      tree$rules <- tree$rules[-which(names(tree$rules) == cellType)]
    }
    
    #Set final tree dendro
    tree$dendro <- dendro
    
    #Get performance metrics
    print('Computing performance metrics...')
    perfList <- .getPerformance(tree$rules, features, as.factor(subtypeLabels))
    tree$prediction <- perfList$prediction
    tree$performance <- perfList$performance
    
    return(tree)
  }
}
