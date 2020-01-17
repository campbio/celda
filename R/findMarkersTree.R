#' @title Generate decision tree from single-cell clustering output.
#' @description Uses decision tree procedure to generate a set of rules for each
#'  cell cluster defined by single-cell clustering. Splits are determined by
#'  one of two metrics at each split: a one-off metric to determine rules for
#'  identifying clusters by a single feature, and a balanced metric to determine
#'  rules for identifying sets of similar clusters.
#' @param features A L (features) by N (samples) numeric matrix.
#' @param class A vector of K label assignments.
#' @param oneoffMetric A character string. What one-off metric to run, either
#'  `modified F1` or  `pairwise AUC`. No default value.
#' @param metaclusters List where each element is a metacluster (e.g. known
#' cell type) and all the clusters within that metacluster (e.g. subtypes).
#' @param featureLabels Optional. A vector of feature assignments. Useful when
#'  using clusters of features (e.g. gene modules) and user wishes to expand
#'  tree results to individual features (e.g. score individual genes within
#'  marker gene modules).
#' @param counts Optional. Numeric counts matrix. Useful whenusing clusters
#'  of features (e.g. gene modules) and user wishes to expandtree results to
#'  individual features (e.g. score individual genes within marker gene
#'  modules). Row names should be individual feature names.
#' @param celda Optional. A celda_CG or celda_C object.
#' @param seurat Optional. A seurat object. Note that the seurat functions
#' RunPCA and FindClusters must have been run on the object.
#' @param threshold Numeric. The threshold for the oneoff metric to use
#'  between 0 and 1. Smaller values will result is more one-off
#'  splits. Default is 0.90.
#' @param reuseFeatures Logical.  Whether or not a feature can be used more than
#'  once on the same cluster. Default is TRUE.
#' @param altSplit Logical. Whether or not to force a marker for clusters that
#'  are solely defined by the absence of markers. Default is FALSE.
#' @param consecutiveOneoff Logical. Whether or not to allow one-off splits at
#'  consecutive brances. Default is FALSE.
#'  @param autoMetaclusters. Logical. Whether to use automatically-identified
#'  metaclusters for creating the tree. Default is TRUE.
#'  @param seed Numeric. Seed used to enable reproducible UMAP results in case
#'  parameter `meta` is set to TRUE. Default is 12345.
#' @return A named list with six elements.
#' \itemize{
#'   \item rules - A named list with one `data.frame` for every label. Each
#' `data.frame` has five columns and gives the set of rules for disinguishing
#'  each label.
#'   \itemize{
#'    \item feature - Feature identifier.
#'    \item direction - Relationship to feature value, -1 if less than, 1 if
#'  greater than.
#'    \item stat - The performance value returned by the splitting metric for
#'  this split.
#'    \item statUsed - Which performance metric was used. "Split" if information
#'  gain and "One-off" if one-off.
#'    \item level - The level of the tree at which is rule was defined. 1 is the
#'  level of the first split of the tree.
#'    \item metacluster - Optional. If metaclusters were used, the metacluster
#'     this rule is applied to.
#'   }
#'  \item dendro - A dendrogram object of the decision tree output
#'  \item classLabels - A vector of the class labels used in the model
#'  \item metaclusterLabels - A vector of the metacluster labels
#'   used in the model
#'  \item prediction - A character vector of label of predictions of the
#'  training data using the final model. "MISSING" if label prediction was
#'  ambiguous.
#'  \item performance - A named list denoting the training performance of the
#'  model.
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
#' @import magrittr
#' @importFrom methods hasArg
#' @import dbscan
#' @import uwot
#' @import pROC
#' @import withr
#' @export
findMarkersTree <- function(features,
                            class,
                            oneoffMetric = c("modified F1", "pairwise AUC"),
                            metaclusters,
                            featureLabels,
                            counts,
                            celda,
                            seurat,
                            threshold = 0.90,
                            reuseFeatures = FALSE,
                            altSplit = FALSE,
                            consecutiveOneoff = FALSE,
                            autoMetaclusters = TRUE,
                            seed = 12345) {
    
    if(methods::hasArg(celda)){
        #check that counts matrix is provided 
        if(!methods::hasArg(counts)){
            stop("Please provide counts matrix in addition to celda object.")
        }
        
        #factorize matrix (proportion of each module in each cell)
        features <- celda::factorizeMatrix(counts, celda)$proportions$cell
        
        #get class labels
        class <- celda@clusters$z
        
        #get feature labels
        featureLabels <- paste0('L',celda@clusters$y)
    }
    else if(methods::hasArg(seurat)){
        #check that counts matrix is provided 
        if(!methods::hasArg(counts)){
            stop("Please provide counts matrix in addition to seurat object.")
        }
        
        #get class labels
        class <- as.character(Idents(seurat))
        
        #get feature labels
        featureLabels <- unlist(apply(seurat@reductions$pca@feature.loadings,1,
                                      function(x)
                                          return(names(x)[which(x==max(x))])))
        
        #sum counts for each PC in each cell
        features <- matrix(unlist(lapply(unique(featureLabels), function(pc){
            colSums(counts[featureLabels==pc,])
        })),
        ncol = length(class), byrow = TRUE,
        dimnames = list(unique(featureLabels),colnames(counts)))
        
        #normalize column-wise (i.e. convert counts to proportions)
        features <- apply(features, 2, function(x) x/sum(x))
    }
    
    if (ncol(features) != length(class)) {
        stop("Number of columns of features must equal length of class")
    }
    
    if (any(is.na(class))) {
        stop("NA class values")
    }
    
    if (any(is.na(features))){
        stop("NA feature values")
    }
    
    # Match the oneoffMetric argument
    oneoffMetric <- match.arg(oneoffMetric)
    
    # Transpose features
    features <- t(features)
    
    # If no detailed cell types are provided
    if (!methods::hasArg(metaclusters) & (!autoMetaclusters)) {
        
        message("Building tree...")
        
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
        if (altSplit) {
            tree <- .addAlternativeSplit(tree, features, class)
        }
        
        message("Computing performance metrics...")
        
        # Format tree output for plotting and generate summary statistics
        DTsummary <- .summarizeTree(tree, features, class)
        
        # Remove confusing 'value' column
        DTsummary$rules <- lapply(DTsummary$rules, function(x) {
            x["value"] <- NULL; x })
        
        # Add column to each rules table which specifies its class
        DTsummary$rules <- mapply(cbind,
                                  "class"=as.character(names(DTsummary$rules)),
                                  DTsummary$rules, SIMPLIFY = FALSE)
        
        # Generate table for each branch point in the tree
        DTsummary$branchPoints <- .createBranchPoints(DTsummary$rules)
        
        # Add class labels to output
        DTsummary$classLabels <- class
        
        return(DTsummary)
    } else {
        ## If metaclusters are provided or to be identified
        
        #consecutive one-offs break the code(tricky to find 1st balanced split)
        if(consecutiveOneoff){
            stop("Cannot use metaclusters if consecutive one-offs are allowed.
                 Please set the consecutiveOneoff parameter to FALSE.")
        }
        
        # Check if need to identify metaclusters
        if(autoMetaclusters & !methods::hasArg(metaclusters)){
            message("Identifying metaclusters...")
            
            #if seurat object then use seurat's UMAP parameters
            if(methods::hasArg(seurat)){
                suppressMessages(seurat <- RunUMAP(seurat, dims = 1:ncol(seurat@reductions$pca@feature.loadings)))
                umap <- seurat@reductions$umap@cell.embeddings
            }
            else{
                if(is.null(seed)){
                    umap <- uwot::umap(t(sqrt(t(features))),
                                       n_neighbors=15, min_dist = 0.01,
                                       spread = 1, n_sgd_threads = 1)
                }
                else{
                    withr::with_seed(
                        seed,
                        umap <- uwot::umap(t(sqrt(t(features))),
                                           n_neighbors = 15, min_dist = 0.01,
                                           spread = 1, n_sgd_threads = 1))
                }
            }
            # dbscan to find metaclusters
            dbscan <- dbscan::dbscan(umap, eps = 1)
            
            # place each population in the correct metacluster
            mapping <-
                unlist(lapply(sort(as.integer(unique(class))),
                              function(population) {
                                  #get indexes of occurences of this population
                                  indexes <- which(class == population)
                                  
                                  #get corresponding metaclusters
                                  metaIndices <- dbscan$cluster[indexes]
                                  
                                  #return corresponding metacluster with majority vote
                                  return(names(sort(table(metaIndices),decreasing=TRUE)[1]))
                              }))
            
            # create list which will contain subtypes of each metacluster
            metaclusters <- vector(mode = "list")
            
            # fill in list of populations for each metacluster
            for(i in unique(mapping)){
                metaclusters[[i]] <- 
                    sort(as.integer(unique(class)))[which(mapping == i)]
            }
            names(metaclusters) <- paste0("M", unique(mapping))
            
            message(paste("Identified",length(metaclusters), "metaclusters"))
        }
        
        # Check that cell types match class labels
        if (mean(unlist(metaclusters) %in% unique(class)) != 1) {
            stop("Provided cell types do not match class labels. ",
                 "Please check the 'metaclusters' argument.")
        }
        
        # Create vector with metacluster labels
        metaclusterLabels <- class
        for (i in names(metaclusters)) {
            metaclusterLabels[metaclusterLabels %in% metaclusters[[i]]] <- i
        }
        
        #create temporary variables for top-level tree
        tmpThreshold <- threshold
        
        #create list to store split off classes at each threshold
        markerThreshold <- list()
        
        # Create top-level tree
        
        #while there is still a balanced split at the top-level
        while(TRUE){
            #create tree
            message("Building top-level tree across all metaclusters...")
            tree <-
                .generateTreeList(
                    features,
                    as.factor(metaclusterLabels),
                    oneoffMetric,
                    tmpThreshold,
                    reuseFeatures,
                    consecutiveOneoff
                )
            
            #store clusters with markers at current threshold
            topLevel <- tree[[1]][[1]]
            if(topLevel$statUsed == "One-off"){
                markerThreshold[[as.character(tmpThreshold)]] <- unlist(
                    lapply(topLevel[1:(length(topLevel)-3)], function(marker){
                        return(marker$group1Consensus)
                    })
                )
            }
            
            #if no more balanced split 
            if(length(tree) == 1){
                #if all clusters have positive markers
                if(length(tree[[1]][[1]]) == (length(metaclusters)+3)){
                    break
                }
                else{
                    #decrease threshold by 10%
                    tmpThreshold <- tmpThreshold*0.9
                    message("Decreasing classifier threshold to ",tmpThreshold)
                    next
                }
            }
            #still balanced split
            else{
                #get up-regulated clusters at first balanced split
                upClass <- tree[[2]][[1]][[1]]$group1Consensus
                
                #if only 2 clusters at the balanced split then merge them
                if((length(upClass) == 1) &&
                   (length(tree[[2]][[1]][[1]]$group2Consensus) == 1)){
                    upClass <- c(upClass, tree[[2]][[1]][[1]]$group2Consensus)
                }
                
                #update metacluster label of each cell
                tmpMeta <- metaclusterLabels
                tmpMeta[tmpMeta %in% upClass] <-
                    paste(upClass, sep="", collapse="+")
                
                
                #create top-level tree again
                tmpTree <-
                    .generateTreeList(
                        features,
                        as.factor(tmpMeta),
                        oneoffMetric,
                        tmpThreshold,
                        reuseFeatures,
                        consecutiveOneoff
                    )
                
                #if new tree still has balanced split/no markers for some
                if((length(tmpTree) > 1) ||
                   (length(tree[[1]][[1]]) != (length(metaclusters) + 3))) {
                    
                    #decrease threshold by 10%
                    tmpThreshold <- tmpThreshold*0.9
                    message("Decreasing classifier threshold to ",tmpThreshold)
                }
                else{
                    #set final metacluster labels to new set of clusters
                    metaclusterLabels <- tmpMeta
                    
                    #set final tree to current tree
                    tree <- tmpTree
                    
                    ##update 'metaclusters' (list of metaclusters)
                    #get celda clusters in these metaclusters
                    newMetacluster <- unlist(metaclusters[upClass]) 
                    #remove old metaclusters
                    metaclusters[upClass] <- NULL 
                    #add new metacluster to list of metaclusters
                    metaclusters[paste(upClass, sep="", collapse="+")] <-
                        list(unname(newMetacluster))
                    
                    break
                }
            }
        }
        
        #re-format output
        finalTree <- tree
        tree <- list(rules = .mapClass2features(finalTree, features,
                                                as.factor(metaclusterLabels),
                                                topLevelMeta = TRUE)$rules)
        
        #keep markers at first threshold they reached only
        markersToRemove <- c()
        for(thresh in names(markerThreshold)){
            thresholdClasses <- markerThreshold[[thresh]]
            for(cl in thresholdClasses){
                curRules <- tree$rules[[cl]]
                lowMarkerIndices <- which(curRules$direction==1 & 
                                              curRules$stat<as.numeric(thresh))
                if(length(lowMarkerIndices)>0){
                    markersToRemove <- c(markersToRemove,
                                         curRules[lowMarkerIndices,'feature'])
                }
            }
        }
        tree$rules <- lapply(tree$rules, function(rules){
            return(rules[!rules$feature %in% markersToRemove,])
        })
        
        #store final set of top-level markers
        topLevelMarkers <- unlist(lapply(tree$rules, function(cluster){
            markers <- cluster[cluster$direction==1, "feature"]
            return(paste(markers, collapse = ";"))
        }))
        
        #create tree dendrogram
        tree$dendro <-
            .convertToDendrogram(finalTree, as.factor(metaclusterLabels),
                                 splitNames = topLevelMarkers)
        
        #add metacluster label to rules table
        for(metacluster in names(tree$rules)){
            tree$rules[[metacluster]]$metacluster <- metacluster
        }
        
        # Store tree's dendrogram in a separate variable
        dendro <- tree$dendro
        
        # Find which metaclusters have more than one cluster
        largeMetaclusters <- names(metaclusters[lengths(metaclusters) > 1])
        
        # Update subtype labels for large metaclusters
        subtypeLabels <- metaclusterLabels
        subtypeLabels[subtypeLabels %in% largeMetaclusters] <- paste0(
            subtypeLabels[subtypeLabels %in% largeMetaclusters],
            "(",
            class[subtypeLabels %in% largeMetaclusters],
            ")"
        )
        
        # Update metaclusters list
        for(metacluster in names(metaclusters)){
            subtypes <- metaclusters[metacluster]
            subtypes <- lapply(subtypes, function(subtype){
                paste0(metacluster,"(",subtype,")")
            })
            metaclusters[metacluster] <- subtypes
        }
        
        # Create separate trees for each cell type with more than one cluster
        newTrees <- lapply(largeMetaclusters, function(metacluster){
            
            # Print current status
            message("Building tree for metacluster ", metacluster)
            
            # Remove used features
            featUse <- colnames(features)
            if (!reuseFeatures) {
                tmpRules <- tree$rules[[metacluster]]
                featUse <-
                    featUse[!featUse %in% 
                                tmpRules[tmpRules$direction == 1, "feature"]]
            }
            
            # Create new tree
            newTree <-
                .generateTreeList(
                    features[metaclusterLabels == metacluster,featUse],
                    as.factor(subtypeLabels[metaclusterLabels == metacluster]),
                    oneoffMetric,
                    threshold,
                    reuseFeatures,
                    consecutiveOneoff
                )
            newTree <- list(
                rules = .mapClass2features(newTree,
                                           features[metaclusterLabels
                                                    == metacluster,],
                                           as.factor(subtypeLabels[
                                               metaclusterLabels== metacluster]
                                           ))$rules,
                dendro = .convertToDendrogram(newTree,
                                              as.factor(
                                                  subtypeLabels[
                                                      metaclusterLabels == 
                                                          metacluster]))
            )
            
            # Adjust 'rules' table for new tree
            newTree$rules <- lapply(newTree$rules, function(rules){
                rules$level <- rules$level +
                    max(tree$rules[[metacluster]]$level)
                rules$metacluster <- metacluster
                rules <- rbind(tree$rules[[metacluster]], rules)
            })
            
            return(newTree)
        })
        names(newTrees) <- largeMetaclusters
        
        # Fix max depth in original tree
        if(length(newTrees) > 0){
            maxDepth <- max(unlist(lapply(newTrees, function(newTree) {
                lapply(newTree$rules, function(ruleDF) {
                    ruleDF$level
                })
            })))
            addDepth <- maxDepth - attributes(dendro)$height
            
            dendro <- dendrapply(dendro, function(node, addDepth) {
                if (attributes(node)$height > 1) {
                    attributes(node)$height <- attributes(node)$height +
                        addDepth + 1
                }
                return(node)
            }, addDepth)
        }
        
        # Find indices of cell type nodes in tree
        indices <- lapply(largeMetaclusters,
                          function(metacluster) {
                              # Initialize sub trees, indices string, and flag
                              dendSub <- dendro
                              index <- ""
                              flag <- TRUE
                              
                              while (flag) {
                                  # Get the edge with the class of interest
                                  whEdge <- which(unlist(
                                      lapply(dendSub,
                                             function(edge)
                                                 metacluster %in%
                                                 attributes(edge)$classLabels)
                                  ))
                                  
                                  # Add this as a string
                                  index <- paste0(index, "[[", whEdge, "]]")
                                  
                                  # Move to this branch
                                  dendSub <- eval(parse(text = paste0(
                                      "dendro", index)))
                                  
                                  # Is this the only class in that branch
                                  flag <- length(
                                      attributes(dendSub)$classLabels) > 1
                              }
                              
                              return(index)
                          }
        )
        names(indices) <- largeMetaclusters
        
        # Add each cell type tree
        for (metacluster in largeMetaclusters) {
            
            # Get current tree
            metaclusterDendro <- newTrees[[metacluster]]$dendro
            
            # Adjust labels, member count, and midpoint of nodes
            dendro <- dendrapply(dendro, function(node){
                # Check if in right branch
                if (metacluster %in%
                    as.character(attributes(node)$classLabels)) {
                    # Replace cell type label with subtype labels
                    attributes(node)$classLabels <-
                        as.character(attributes(node)$classLabels) %>%
                        .[. != metacluster] %>%
                        c(., unique(subtypeLabels)[grep(metacluster,
                                                        unique(subtypeLabels))
                                                   ])
                    
                    # Assign new member count for this branch
                    attributes(node)$members <-
                        length(attributes(node)$classLabels)
                    
                    # Assign new midpoint for this branch
                    attributes(node)$midpoint <-
                        (attributes(node)$members - 1) / 2
                }
                return(node)
            })
            
            # Replace label at new tree's branch point
            branchPointAttr <- attributes(eval(parse(
                text = paste0("dendro", indices[[metacluster]]))))
            branchPointLabel <- branchPointAttr$label
            branchPointStatUsed <- branchPointAttr$statUsed
            
            if (!is.null(branchPointLabel)) {
                attributes(metaclusterDendro)$label <- branchPointLabel
                attributes(metaclusterDendro)$statUsed <- branchPointStatUsed
            }
            
            # Fix height
            indLoc <- gregexpr("\\[\\[", indices[[metacluster]])[[1]]
            indLoc <- indLoc[length(indLoc)]
            parentIndexString <- substr(indices[[metacluster]],
                                        0,
                                        indLoc - 1)
            parentHeight <- attributes(eval(parse(
                text = paste0("dendro", parentIndexString))))$height
            metaclusterHeight <- attributes(metaclusterDendro)$height
            metaclusterDendro <- dendrapply(metaclusterDendro,
                                            function(node, parentHeight,
                                                     metaclusterHeight) {
                                                if (attributes(node)$height > 1){
                                                    attributes(node)$height <-
                                                        parentHeight - 1 - 
                                                        (metaclusterHeight -
                                                             attributes(node)$height)
                                                }
                                                return(node)
                                            }, parentHeight, metaclusterHeight)
            
            # Add new tree to original tree
            eval(parse(text = paste0(
                "dendro", indices[[metacluster]], " <- metaclusterDendro")))
            
            # Append new tree's 'rules' tables to original tree
            tree$rules <- append(tree$rules, newTrees[[metacluster]]$rules,
                                 after = which(names(tree$rules)==metacluster))
            
            # Remove old tree's rules
            tree$rules <- tree$rules[-which(names(tree$rules) == metacluster)]
        }
        
        # Set final tree dendro
        tree$dendro <- dendro
        
        # Get performance statistics
        message("Computing performance statistics...")
        perfList <- .getPerformance(tree$rules,
                                    features,
                                    as.factor(subtypeLabels))
        tree$prediction <- perfList$prediction
        tree$performance <- perfList$performance
        
        # Remove confusing 'value' column
        tree$rules <- lapply(tree$rules, function(x) { x["value"] <- NULL; x })
        
        #add column to each rules table which specifies its class
        tree$rules <- mapply(cbind, "class" = as.character(names(tree$rules)),
                             tree$rules, SIMPLIFY = FALSE)
        
        #create branch points table
        branchPoints <- .createBranchPoints(tree$rules, largeMetaclusters, metaclusters)
        
        #collapse all rules tables into one large table
        collapsed <- do.call("rbind", tree$rules)
        
        #get top-level rules
        topLevelRules <- collapsed[collapsed$level == 1,]
        
        #add 'class' column
        topLevelRules$class <- topLevelRules$metacluster
        
        #add to branch point list
        branchPoints[["top_level"]] <- topLevelRules
        
        #check if need to expand features to gene-level
        if(methods::hasArg(featureLabels) && methods::hasArg(counts)){
            
            message("Computing scores for individual genes...")
            
            #make sure feature labels match those in the tree
            if(!all(unique(collapsed$feature) %in% unique(featureLabels))){
                m<-"Provided feature labels don't match those in count matrix."
                stop(m)
            }
            
            #iterate over branch points
            branchPoints <- lapply(branchPoints, function(branch){
                
                #iterate over unique features
                featAUC <- lapply(unique(branch$feature), .getGeneAUC, branch,
                                  subtypeLabels, metaclusterLabels,
                                  featureLabels)
                
                #update branch table after merging genes data
                return(do.call("rbind", featAUC))
            })
            
            #simplify top-level in rules tables to only up-regulated markers
            tree$rules <- lapply(tree$rules, function(rule){
                return(rule[-intersect(
                    which(rule$level==1), which(rule$direction==(-1))), ])
            })
            
            ##add gene-level info to rules tables
            #collapse branch points tables into one
            collapsedBranches <- do.call("rbind", branchPoints)
            collapsedBranches$class <- as.character(collapsedBranches$class)
            
            #loop over rules tables and get relevant info
            tree$rules <- lapply(tree$rules, function(class){
                #initialize table to return
                toReturn <- data.frame(NULL)
                
                #loop over rows of this class
                for(i in 1:nrow(class)){
                    #extract relevant genes from branch points tables
                    genesAUC <- collapsedBranches[
                        collapsedBranches$feature==class$feature[i] &
                            collapsedBranches$level==class$level[i] &
                            collapsedBranches$class==class$class[i],]
                    
                    #don't forget top-level
                    if(class$level[i]==1){
                        genesAUC <- collapsedBranches[
                            collapsedBranches$feature==class$feature[i] &
                                collapsedBranches$level==class$level[i] &
                                collapsedBranches$class==class$metacluster[i],]
                    }
                    
                    #merge table
                    toReturn <- rbind(toReturn, genesAUC)
                }
                return(toReturn)
            })
            
            #remove table row names
            tree$rules <- lapply(tree$rules, function(t){
                rownames(t) <- NULL
                return(t)
            })
            
            #add feature labels to output
            tree$featureLabels <- featureLabels
            
        }
        
        #simplify top-level branch point to save memory
        branchPoints$top_level <- branchPoints$top_level[
            branchPoints$top_level$direction==1,]
        branchPoints$top_level <- branchPoints$top_level[
            !duplicated(branchPoints$top_level),]
        
        #remove branch points row names
        branchPoints <- lapply(branchPoints, function(br){
            rownames(br) <- NULL
            return(br)
        })
        
        #add to tree
        tree$branchPoints <- branchPoints
        
        #return class labels
        tree$classLabels <- subtypeLabels
        tree$metaclusterLabels <- metaclusterLabels
        
        # Final return
        return(tree)
    }
}

#helper function to create table for each branch point in the tree
.createBranchPoints <- function(rules, largeMetaclusters, metaclusters){
    # First step differs if metaclusters were used
    
    if(methods::hasArg(metaclusters) && (length(largeMetaclusters)>0)){
        #iterate over metaclusters and add the rules for each level
        branchPoints <- lapply(largeMetaclusters, function(metacluster){
            #get names of subtypes
            subtypes <- metaclusters[[metacluster]]
            
            #collapse rules tables of subtypes
            subtypeRules <- do.call("rbind", rules[subtypes])
            
            #get rules at each level
            levels <- lapply(2:max(subtypeRules$level), function(level){
                return(subtypeRules[subtypeRules$level == level,])
            })
            names(levels) <- paste0(metacluster,"_level_",
                                    1:(max(subtypeRules$level)-1))
            
            return(levels)
        })
        branchPoints <- unlist(branchPoints, recursive = FALSE)
    }
    else{
        #collapse all rules into one table
        collapsed <- do.call("rbind", rules)
        
        #subset rules at each level
        branchPoints <- lapply(1:max(collapsed$level), function(level){
            return(collapsed[collapsed$level==level,])
        })
        names(branchPoints) <- paste0("level_",1:max(collapsed$level))
    }
    
    #split each level into its branch points
    branchPoints <- lapply(branchPoints, function(level){
        #check if need to split
        firstFeat <- level$feature[1]
        firstStat <- level$stat[1]
        if(setequal(level[level$feature==firstFeat & level$stat==firstStat,
                          "class"],
                    unique(level$class))){
            return(level)
        }
        
        #initialize lists for new tables
        bSplits <- NA
        oSplits <- NA
        
        #get balanced split rows by themselves
        balS <- level[level$statUsed=="Split",]
        
        #return table for each unique value of 'stat'
        if(nrow(balS) > 0){
            #get unique splits (based on stat)
            unS <- unique(balS$stat)
            
            #return table for each unique split
            bSplits <- lapply(unS, function(s){
                balS[balS$stat==s,]
            })
        }
        
        #get one-off rows by themselves
        oneS <- level[level$statUsed=="One-off",]
        
        if(nrow(oneS) > 0){
            #check if need to split
            firstFeat <- oneS$feature[1]
            if(setequal(oneS[oneS$feature==firstFeat,"class"],
                        unique(oneS$class))){
                oSplits <- oneS
            }
            
            #get class groups for each marker
            markers <- oneS[oneS$direction==1, "feature"]
            groups <- unique(unlist(lapply(markers, function(m){
                return(paste(
                    as.character(oneS[oneS$feature==m,"class"]),
                    collapse = " "))
            })))
            
            #return table for each class group
            oSplits <- lapply(groups, function(x){
                gr <- unlist(strsplit(x, split = " "))
                oneS[as.character(oneS$class) %in% gr,]
            })
        }
        
        #rename new tables
        if(is.list(bSplits)){
            names(bSplits) <- paste0("split_",
                                     LETTERS[length(bSplits):1])
        }
        if(is.list(oSplits)){
            names(oSplits) <- paste0("one-off_",
                                     LETTERS[length(oSplits):1])
        }
        
        #return 2 sets of table
        toReturn <- list(oSplits, bSplits)
        toReturn <- toReturn[!is.na(toReturn)]
        toReturn <- unlist(toReturn, recursive = FALSE)
        return(toReturn)
    }
    )
    
    #adjust for new tables
    branchPoints <- lapply(branchPoints, function(br){
        if(inherits(br, "list")){
            return(br)
        }
        else{
            return(list(br))
        }
    })
    branchPoints <- unlist(branchPoints, recursive = FALSE)
    #replace dots in names of new branches with underscores
    names(branchPoints) <- gsub(pattern = "\\.([^\\.]*)$",
                                replacement = "_\\1",
                                names(branchPoints))
    
    return(branchPoints)
}

#helper function to get AUC for individual genes within feature
.getGeneAUC <- function(marker, table, subtypeLabels,
                        metaclusterLabels, featureLabels){
    #get up-regulated & down-regulated classes for this feature
    upClass <-
        as.character(table[table$feature == marker &
                               table$direction == 1, "class"])
    downClasses <-
        as.character(table[table$feature == marker &
                               table$direction == (-1), "class"])
    
    #subset counts matrix
    if(table$level[1] > 1){
        subCounts <- counts[,which(subtypeLabels %in% c(upClass, downClasses))]
    }
    else{
        subCounts <- counts[,which(metaclusterLabels %in%
                                       c(upClass, downClasses))]
    }
    
    #subset class labels
    if(table$level[1] > 1){
        subLabels <- subtypeLabels[which(subtypeLabels %in%
                                             c(upClass, downClasses))]
    }
    else{
        subLabels <- metaclusterLabels[which(metaclusterLabels %in%
                                                 c(upClass, downClasses))]
    }
    
    #set label to 0 if not class of interest
    subLabels <- as.numeric(subLabels %in% upClass)
    
    #get individual features within this marker
    markers <- rownames(counts)[which(featureLabels == marker)]
    
    #get one-vs-all AUC for each gene
    auc <- unlist(lapply(markers, function(markerGene){
        as.numeric(pROC::auc(pROC::roc(subLabels, subCounts[markerGene,],
                                       direction = "<", quiet = TRUE)))
    }))
    names(auc) <- markers
    
    #sort by AUC
    auc <- sort(auc, decreasing = TRUE)
    
    #create table for this marker
    featTable <- table[table$feature == marker,]
    featTable <- featTable[rep(seq_len(nrow(featTable)), each=length(auc)),]
    featTable$gene <- rep(names(auc), length(c(upClass, downClasses)))
    featTable$geneAUC <- rep(auc, length(c(upClass, downClasses)))
    
    #return table for merging with main table
    return(featTable)
}
