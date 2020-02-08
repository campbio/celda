#' @title Generate marker decision tree from single-cell clustering output
#' @description Create a decision tree that identifies gene markers for given
#'  cell populations. The algorithm uses a decision tree procedure to generate
#'  a set of rules for each cell cluster defined by single-cell clustering.
#'  Splits are determined by one of two metrics at each split: a one-off metric
#'  to determine rules for identifying clusters by a single feature, and a
#'  balanced metric to determine rules for identifying sets of similar clusters.
#' @param features features-by-samples numeric matrix, e.g. counts matrix.
#' @param class Vector of cell cluster labels.
#' @param oneoffMetric A character string. What one-off metric to run, either
#'  `modified F1` or `pairwise AUC`. Default is 'modified F1'.
#' @param metaclusters List where each element is a metacluster (e.g. known
#' cell type) and all the clusters within that metacluster (e.g. subtypes).
#' @param featureLabels  Vector of feature assignments, e.g. which cluster
#'  does each gene belong to? Useful when using clusters of features
#'  (e.g. gene modules or Seurat PCs) and user wishes to expand tree results
#'  to individual features (e.g. score individual genes within marker gene modules).
#' @param counts Numeric counts matrix. Useful when using clusters
#'  of features (e.g. gene modules) and user wishes to expand tree results to
#'  individual features (e.g. score individual genes within marker gene
#'  modules). Row names should be individual feature names.
#' @param celda A \emph{celda_CG} or \emph{celda_C} object.
#'  Counts matrix has to be provided as well.
#' @param seurat A seurat object. Note that the seurat functions
#' \emph{RunPCA} and \emph{FindClusters} must have been run on the object.
#' @param threshold Numeric between 0 and 1. The threshold for the oneoff
#'  metric. Smaller values will result in more one-off splits. Default is 0.90.
#' @param reuseFeatures Logical. Whether or not a feature can be used more than
#'  once on the same cluster. Default is TRUE.
#' @param altSplit Logical. Whether or not to force a marker for clusters that
#'  are solely defined by the absence of markers. Default is TRUE.
#' @param consecutiveOneoff Logical. Whether or not to allow one-off splits at
#'  consecutive brances. Default is FALSE.
#' @param autoMetaclusters. Logical. Whether to identify metaclusters prior to
#'  creating the tree based on the distance between clusters in a UMAP 
#'  dimensionality reduction projection. A metacluster is simply a large
#'  cluster that includes several clusters within it. Default is TRUE.
#' @param seed Numeric. Seed used to enable reproducible UMAP results
#'  for identifying metaclusters. Default is 12345.
#' @return A named list with six elements:
#' \itemize{
#'   \item rules - A named list with one data frame for every label. Each
#'  data frame has five columns and gives the set of rules for disinguishing
#'  each label.
#'   \itemize{
#'    \item feature - Marker feature, e.g. marker gene name.
#'    \item direction - Relationship to feature value. -1 if cluster is
#'    down-regulated for this feature, 1 if cluster is up-regulated.
#'    \item stat - The performance value returned by the splitting metric for
#'  this split.
#'    \item statUsed - Which performance metric was used. "Split" if information
#'  gain and "One-off" if one-off.
#'    \item level - The level of the tree at which is rule was defined. 1 is the
#'  level of the first split of the tree.
#'    \item metacluster - Optional. If metaclusters were used, the metacluster
#'     this rule is applied to.
#'   }
#'  \item dendro - A dendrogram object of the decision tree output. Plot with 
#'  plotDendro()
#'  \item classLabels - A vector of the class labels used in the model, i.e.
#'   cell cluster labels.
#'  \item metaclusterLabels - A vector of the metacluster labels
#'   used in the model
#'  \item prediction - A character vector of label of predictions of the
#'  training data using the final model. "MISSING" if label prediction was
#'  ambiguous.
#'  \item performance - A named list denoting the training performance of the
#'  model:
#'  \itemize{
#'   \item accuracy - (number correct/number of samples) for the whole set of
#'  samples.
#'   \item balAcc - mean sensitivity across all clusters
#'   \item meanPrecision - mean precision across all clusters
#'   \item correct - the number of correct predictions of each cluster
#'   \item sizes - the number of actual counts of each cluster
#'   \item sensitivity - the sensitivity of the prediciton of each cluster
#'   \item precision - the precision of the prediciton of each cluster
#'  }
#' }
#' @examples
#' # Generate simulated single-cell dataset using celda 
#' sim_counts <- celda::simulateCells("celda_CG", K = 4, L = 10, G = 100)
#' 
#' # Celda clustering into 5 clusters & 10 modules
#' cm <- celda_CG(sim_counts$counts, K=5, L=10, verbose=FALSE)
#' 
#' # Get features matrix and cluster assignments
#' factorized <- factorizeMatrix(sim_counts$counts, cm)
#' features <- factorized$proportions$cell
#' class <- clusters(cm)$z
#' 
#' # Generate Decision Tree
#' DecTree <- findMarkersTree(features, class)
#' 
#' # Plot dendrogram
#' plotDendro(DecTree)
#' 
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
                            altSplit = TRUE,
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
        #get counts matrix from seurat object
        counts <- as.matrix(seurat@assays$RNA@data)
        
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
    
    # If no detailed cell types are provided or to be identified
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
        # If metaclusters are provided or to be identified
        
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
        
        # Rename metaclusters with just one cluster
        oneCluster <- names(metaclusters[lengths(metaclusters) == 1])
        if(length(oneCluster) > 0){
            oneClusterIndices <- which(metaclusterLabels %in% oneCluster)
            metaclusterLabels[oneClusterIndices] <- 
                paste0(metaclusterLabels[oneClusterIndices],"(",
                       class[oneClusterIndices],")")
            names(metaclusters[lengths(metaclusters) == 1]) <- 
                paste0(names(metaclusters[lengths(metaclusters) == 1]), "(",
                       unlist(metaclusters[lengths(metaclusters) == 1]),")")
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
            
            # Add alternative node for the solely down-regulated leaf
            tree <- .addAlternativeSplit(tree, features,
                                         as.factor(metaclusterLabels))
            
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
                
                # Add alternative node for the solely down-regulated leaf
                tmpTree <- .addAlternativeSplit(tmpTree, features,
                                                as.factor(tmpMeta))
                
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
                if(length(lowMarkerIndices)>0 & length(which(curRules$direction==1))>1){
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
            
            # Add alternative node for the solely down-regulated leaf
            if (altSplit) {
                newTree <-
                    .addAlternativeSplit(newTree,
                                         features[metaclusterLabels == metacluster, featUse],
                                         as.factor(subtypeLabels[metaclusterLabels == metacluster]))
            }
            
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
                    labels <- attributes(node)$classLabels
                    labels <- as.character(labels)
                    labels <- labels[labels != metacluster]
                    labels <- c(labels, unique(subtypeLabels)
                                [grep(metacluster,unique(subtypeLabels))])
                    attributes(node)$classLabels <- labels
                    
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
                                  featureLabels, counts)
                
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
        
        #adjust subtype labels
        branchPoints <- lapply(branchPoints, function(br){
            br$class <- as.character(br$class)
            br$class[grepl("\\(.*\\)",br$class)] <- regmatches(
                br$class[grepl("\\(.*\\)",br$class)],
                regexpr(pattern = "(?<=\\().*?(?=\\)$)",
                        br$class[grepl("\\(.*\\)",br$class)],
                        perl = TRUE))
            
            br$metacluster <- as.character(br$metacluster)
            br$metacluster[grepl("\\(.*\\)",br$metacluster)] <-
                gsub("\\(.*\\)", "",
                     br$metacluster[grepl("\\(.*\\)",br$metacluster)])
            
            return(br)
        })
        #adjust subtype labels
        tree$rules <- suppressWarnings(lapply(tree$rules, function(r){
            r$class <- as.character(r$class)
            r$class[grepl("\\(.*\\)",r$class)] <- regmatches(
                r$class[grepl("\\(.*\\)",r$class)],
                regexpr(pattern = "(?<=\\().*?(?=\\)$)",
                        r$class[grepl("\\(.*\\)",r$class)],
                        perl = TRUE))
            
            r$metacluster[grepl("\\(.*\\)",r$metacluster)] <-
                gsub("\\(.*\\)", "",
                     r$metacluster[grepl("\\(.*\\)",r$metacluster)])
            return(r)
        }))
        
        
        #add to tree
        tree$branchPoints <- branchPoints
        
        #return class labels
        tree$classLabels <- regmatches(subtypeLabels,
                                       regexpr(pattern = "(?<=\\().*?(?=\\)$)",
                                               subtypeLabels, perl = TRUE))
        
        tree$metaclusterLabels <- metaclusterLabels
        tree$metaclusterLabels[grepl("\\(.*\\)",metaclusterLabels)] <- 
            gsub("\\(.*\\)", "",
                 metaclusterLabels[grepl("\\(.*\\)",metaclusterLabels)])
        
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
                        metaclusterLabels, featureLabels, counts){
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

# This function generates the decision tree by recursively separating classes.
.generateTreeList <- function(
    features,
    class,
    oneoffMetric,
    threshold,
    reuseFeatures,
    consecutiveOneoff = FALSE) {
    
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
            if (!consecutiveOneoff & split$statUsed == "One-off") {
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
        statUsed <- "One-off"
    } else {
        splitStats <- integer(0)
    }
    
    
    # If no one-2-one split meets threshold, run semi-supervised clustering
    if (length(splitStats) == 0) {
        splitMetric <- .splitMetricIGpIGd
        splitStats <- .splitMetricRecursive(features,
                                            class,
                                            splitMetric = splitMetric)[1] # Use top
        statUsed <- "Split"
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
        aucK1 <- pROC::auc(pROC::roc(classK1, currentFeature, direction = "<", quiet = TRUE))
        
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
            rocK2 <- pROC::roc(currentLabels, currentFeatureSubset,direction = "<", quiet=TRUE)
            aucK2 <- rocK2$auc
            coordK2 <- pROC::coords(rocK2, "best", ret = "threshold", transpose = TRUE)[1]
            
            # Concatenate vectors
            statK2 <- c(threshold = coordK2, auc = aucK2)
            
            return(statK2)
        }, classMax, class, currentFeature)))
    
    # Get Min Value
    aucMin <- min(aucFram$auc)
    
    # Get indices where this AUC occurs
    aucMinIndices <- which(aucFram$auc == aucMin)
    
    # Use maximum value if there are ties
    aucValue <- max(aucFram$threshold)
    
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
    
    # Unlist decsision tree
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
            names(downSplit)[length(downSplit)] <- altStats$feat[1]
            downSplit <- downSplit[c(
                which(!names(downSplit) %in% c("statUsed", "fUsed", "dirs")),
                which(names(downSplit) %in% c("statUsed", "fUsed", "dirs")))]
            
            # Get index of split to add it to
            branchLengths <- unlist(lapply(tree, length))
            branchCum <- cumsum(branchLengths)
            wBranch <- min(which(branchCum >= AltDec))
            if(wBranch==1){
                wSplit <- 1
            }
            else{
                wSplit <- which(seq(
                    (branchCum[(wBranch - 1)] + 1),
                    branchCum[wBranch]) == AltDec)
            }
            
            # Add it to decision tree
            tree[[wBranch]][[wSplit]] <- downSplit
        } else {
            cat("No non-ambiguous rule to separate", group2only, "from",
                branchClasses, ". No alternative split added.")
        }
    } else {
        #  print("No solely down-regulated cluster to add alternative split.")
    }
    
    return(tree)
}

#' @title Gets cluster estimates using rules generated by
#'  `celda::findMarkersTree`
#' @description Get decisions for a matrix of features. Estimate cell
#'  cluster membership using feature matrix input.
#' @param rules List object. The `rules` element from  `findMarkersTree`
#'  output. Returns NA if cluster estimation was ambiguous.
#' @param features A L(features) by N(samples) numeric matrix.
#' @return A character vector of label predicitions.

getDecisions <- function(rules, features) {
    features <- t(features)
    votes <- apply(features, 1, .predictClass, rules)
    return(votes)
}

# Function to predict class from list of rules
.predictClass <- function(samp, rules){
    
    # Initilize possible classes and level
    classes <- names(rules)
    level <- 1
    
    # Set maximum levele possible to prevent infinity run
    maxLevel <- max(unlist(lapply(rules, function(ruleSet) {
        ruleSet$level
    })))
    
    while (length(classes) > 1 & level <= maxLevel) {
        
        # Get possible classes
        clLogical <- unlist(lapply(classes, function(cl, rules, level, samp) {
            
            # Get the rules for this class
            ruleClass <- rules[[cl]]
            
            # Get the rules for this level
            ruleClass <- ruleClass[ruleClass$level == level, , drop = FALSE]
            
            # Subset class for the features at this level
            ruleClass$sample <- samp[ruleClass$feature]
            
            # For multiple direction == 1, use one with the top stat
            if (sum(ruleClass$direction == 1) > 1){
                ruleClass <- ruleClass[order(
                    ruleClass$direction
                    , decreasing = T), ]
                ruleClass <- ruleClass[c(which.max(
                    ruleClass$stat[ruleClass$direction == 1]),
                    which(ruleClass$direction == -1)), , drop = FALSE]
            }
            
            # Check for followed rules
            ruleClass$check <- ruleClass$sample >= ruleClass$value
            ruleClass$check[ruleClass$direction == -1] <- !ruleClass$check[
                ruleClass$direction == -1]
            
            # Check that all rules were followed
            ruleFollowed <- mean(
                ruleClass$check & ruleClass$direction == 1) > 0 |
                mean(ruleClass$check) == 1
            
            return(ruleFollowed)
            
        }, rules, level, samp))
        
        # Subset possible classes
        classes <- classes[clLogical]
        
        # Add level
        level <- level + 1
    }
    
    # Return if only one class selected
    if (length(classes) == 1) {
        return(classes)
    } else {
        return(NA)
    }
}

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
        prediction = perfList$prediction,
        performance = perfList$performance
    ))
}

# Function to reformat raw tree ouput to a dendrogram
.convertToDendrogram <- function(tree, class, splitNames = NULL) {
    
    # Unlist decision tree (one element for each split)
    DecTree <- unlist(tree, recursive = F)
    
    if(is.null(splitNames)){
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
    }
    else{
        names(splitNames) <- 1:(length(DecTree[[1]])-3)
    }
    
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
                                                       group1classVec]
            
            # Assign
            for (j in seq(length(group1classUnique))) {
                mat[i, group1classUnique[[j]]] <- j
            }
            mat[i, group2classUnique] <- j + 1
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

# Create rules of classes and features sequences
.mapClass2features <- function(tree, features, class, topLevelMeta = FALSE) {
    
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
    
    # Generate list of rules for each class
    if(topLevelMeta){
        orderedClass <- unique(
            class2featuresIndices[class2featuresIndices$direction==1,"class"])
    }
    else{
        orderedClass <- levels(class)
    }
    
    rules <- lapply(orderedClass, function(cl, class2featuresIndices) {
        
        class2featuresIndices[class2featuresIndices$class == cl,
                              colnames(class2featuresIndices) != "class"]
        
    }, class2featuresIndices)
    names(rules) <- orderedClass
    
    return(list(
        rules = rules))
}

#' @title Plots dendrogram of \emph{findMarkersTree} output
#' @description Generates a dendrogram of the rules and performance
#' (optional) of the decision tree generated by findMarkersTree().
#' @param tree List object. The output of findMarkersTree()
#' @param classLabel A character value. The name of a specific label to draw
#'  the path and rules. If NULL (default), the tree for all clusters is shown.
#' @param addSensPrec Logical. Print training sensitivities and precisions
#'  for each cluster below leaf label? Default is FALSE.
#' @param maxFeaturePrint Numeric value. Maximum number of markers to print
#'  at a given split. Default is 4.  
#' @param leafSize Numeric value. Size of text below each leaf. Default is 24.
#' @param boxSize Numeric value. Size of rule labels. Default is 7.
#' @param boxColor Character value. Color of rule labels. Default is black.
#' @examples
#' # Generate simulated single-cell dataset using celda 
#' sim_counts <- celda::simulateCells("celda_CG", K = 4, L = 10, G = 100)
#' 
#' # Celda clustering into 5 clusters & 10 modules
#' cm <- celda_CG(sim_counts$counts, K=5, L=10, verbose=FALSE)
#' 
#' # Get features matrix and cluster assignments
#' factorized <- factorizeMatrix(sim_counts$counts, cm)
#' features <- factorized$proportions$cell
#' class <- clusters(cm)$z
#' 
#' # Generate Decision Tree
#' DecTree <- findMarkersTree(features,class,threshold = 1)
#' 
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
        perfVec <- paste0("Sens. ",
                          format(round(performance$sensitivity, 2), nsmall = 2),
                          "\n Prec. ",
                          format(round(performance$precision, 2), nsmall = 2)
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
        leafLabels <- paste(leafLabels, perfVec, sep = "\n")
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
                             ggplot2::aes(x = dendSegsLabelled$xend,
                                          y = dendSegsLabelled$y,
                                          label = dendSegsLabelled$label),
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


#' @title Generate heatmap for a marker decision tree
#' @description Creates heatmap for a specified branch point in a marker tree.
#' @param tree A decision tree from CELDA's \emph{findMarkersTree} function.
#' @param counts Numeric matrix. Gene-by-cell counts matrix.
#' @param branchPoint Character. Name of branch point to plot heatmap for.
#' Name should match those in \emph{tree$branchPoints}.
#' @param featureLabels List of feature cluster assignments. Length should
#' be equal to number of rows in counts matrix, and formatting should match
#' that used in \emph{findMarkersTree()}. Required when using clusters
#' of features and not previously provided to \emph{findMarkersTree()}
#' @param topFeatures Integer. Number of genes to plot per marker module.
#' Genes are sorted based on their AUC for their respective cluster. 
#' Default is 10.
#' @param silent Logical. Whether to avoid plotting heatmap to screen.
#' Default is FALSE.
#' @return A heatmap visualizing the counts matrix for the cells and genes at
#' the specified branch point.
#' @examples
#' # Generate simulated single-cell dataset using celda 
#' sim_counts <- celda::simulateCells("celda_CG", K = 4, L = 10, G = 100)
#' 
#' # Celda clustering into 5 clusters & 10 modules
#' cm <- celda_CG(sim_counts$counts, K=5, L=10, verbose=FALSE)
#' 
#' # Get features matrix and cluster assignments
#' factorized <- factorizeMatrix(sim_counts$counts, cm)
#' features <- factorized$proportions$cell
#' class <- clusters(cm)$z
#' 
#' # Generate Decision Tree
#' DecTree <- findMarkersTree(features,class,threshold = 1)
#' 
#' # Plot example heatmap
#' plotMarkerHeatmap(DecTree, sim_counts$counts, branchPoint = "top_level",
#' featureLabels = paste0("L",clusters(cm)$y))
#' 
#' @export
plotMarkerHeatmap <- function(tree, counts, branchPoint, featureLabels,
                              topFeatures = 10, silent = FALSE){
    
    #get branch point to plot
    branch <- tree$branchPoints[[branchPoint]]
    
    #check that user entered valid branch point name
    if(is.null(branch)){
        stop("Invalid branch point. 
         Branch point name should match one of those in tree$branchPoints.")
    }
    
    #convert counts matrix to matrix (e.g. from dgCMatrix)
    counts <- as.matrix(counts)
    
    #get marker features
    marker <- unique(branch$feature)
    
    #add feature labels
    if('featureLabels' %in% names(tree)){
        featureLabels <- tree$featureLabels
    }
    
    #check that feature labels are provided
    if(missing(featureLabels) &
       !('featureLabels' %in% names(tree)) &
       (sum(marker %in% rownames(counts)) != length(marker))){
        stop("Please provide feature labels, i.e. gene cluster labels")
    }
    else{
        if(missing(featureLabels) &
           !('featureLabels' %in% names(tree)) &
           (sum(marker %in% rownames(counts)) == length(marker))){
            featureLabels == rownames(counts)
        }
    }
    
    #make sure feature labels match the table
    if(!all(branch$feature %in% featureLabels)){
        stop("Provided feature labels don't match those in the tree.
             Please check the feature names in the tree's rules' table.")
    }
    
    #if top-level in metaclusters tree
    if(branchPoint == "top_level"){
        #get unique metaclusters
        metaclusters <- unique(branch$metacluster)
        
        #list which will contain final set of genes for heatmap
        whichFeatures <- c()
        
        #loop over unique metaclusters
        for(meta in metaclusters){
            
            #subset table
            curMeta <- branch[branch$metacluster==meta,]
            
            #if we have gene-level info in the tree
            if("gene" %in% names(branch)){
                #sort by gene AUC score
                curMeta <- curMeta[order(curMeta$geneAUC, decreasing = TRUE),]
                
                #get genes
                genes <- unique(curMeta$gene)
                
                #keep top N features
                genes <- head(genes, topFeatures)
                
                #get gene indices
                markerGenes <- which(rownames(counts) %in% genes)
                
                #get features with non-zero variance to avoid clustering error
                markerGenes <- .removeZeroVariance(counts, 
                                                   cells = which(
                                                       tree$metaclusterLabels %in%
                                                           unique(curMeta$metacluster)),
                                                   markers = markerGenes)
                
                #add to list of features
                whichFeatures <- c(whichFeatures, markerGenes)
            }
            else{
                #current markers
                curMarker <- unique(curMeta$feature)
                
                #get marker gene indices
                markerGenes <- which(featureLabels %in% curMarker)
                
                #get features with non-zero variance to avoid error
                markerGenes <- .removeZeroVariance(counts, 
                                                   cells = which(
                                                       tree$metaclusterLabels %in%
                                                           unique(curMeta$metacluster)),
                                                   markers = markerGenes)
                
                #add to list of features
                whichFeatures <- c(whichFeatures, markerGenes)
            }
        }
        
        #order the metaclusters by size
        colOrder <- data.frame(groupName = names(
            sort(table(tree$metaclusterLabels), decreasing = T)), 
            groupIndex = seq_along(unique(tree$metaclusterLabels)))
        
        #order the markers for metaclusters
        allMarkers <- setNames(as.list(colOrder$groupName), colOrder$groupName)
        allMarkers <- lapply(allMarkers, function(x){
            unique(branch[branch$metacluster==x,"feature"])
        })
        rowOrder <- data.frame(groupName = unlist(allMarkers),
                               groupIndex = seq_along(unlist(allMarkers)))
        toRemove <- which(!rowOrder$groupName %in% featureLabels[whichFeatures])
        if(length(toRemove) > 0){
            rowOrder <- rowOrder[-toRemove,]
        }
        
        #sort cells according to metacluster size
        x <- tree$metaclusterLabels
        y <- colOrder$groupName
        sortedCells <- seq(ncol(counts))[order(match(x,y))]
        
        #create heatmap with only the markers
        return(
            plotHeatmap(counts = counts, z = tree$metaclusterLabels,
                        y = featureLabels, featureIx=whichFeatures, 
                        cellIx = sortedCells, showNamesFeature = TRUE,
                        main = "Top-level", silent = silent,
                        treeheightFeature = 0, colGroupOrder = colOrder,
                        rowGroupOrder = rowOrder, treeheightCell = 0)
        )
    }
    
    #if balanced split
    if(branch$statUsed[1] == "Split"){
        #keep entries for balanced split only (in case of alt. split)
        split <- branch$feature[1]
        branch <- branch[branch$feature == split,]
        
        #get up-regulated and down-regulated classes
        upClasses <- unique(branch[branch$direction==1, "class"])
        downClasses <- unique(branch[branch$direction==(-1), "class"])
        
        #re-order cells to keep up and down separate on the heatmap
        reorderedCells <- c((which(tree$classLabels %in% upClasses)
                             [order(tree$classLabels[
                                 tree$classLabels %in% upClasses])]),
                            (which(tree$classLabels %in% downClasses)
                             [order(tree$classLabels[
                                 tree$classLabels %in% downClasses])]))
        
        #cell annotation based on split
        cellAnno <- data.frame(split = rep("Down-regulated", ncol(counts)),
                               stringsAsFactors = FALSE)
        cellAnno$split[which(tree$classLabels %in% upClasses)] <- "Up-regulated"
        rownames(cellAnno) <- colnames(counts)
        
        #if we have gene-level info in the tree
        if(("gene" %in% names(branch))){
            #get genes
            genes <- unique(branch$gene)
            
            #keep top N features
            genes <- head(genes, topFeatures)
            
            #get gene indices
            whichFeatures <- which(rownames(counts) %in% genes)
            
            #get features with non-zero variance to avoid error
            whichFeatures <- .removeZeroVariance(counts, 
                                                 cells = which(
                                                     tree$classLabels %in%
                                                         unique(branch$class)),
                                                 markers = whichFeatures)
            
            #create heatmap with only the split feature and split classes
            return(plotHeatmap(counts = counts, z = tree$classLabels,
                               y=featureLabels, featureIx=whichFeatures,
                               cellIx = reorderedCells, clusterCell = FALSE,
                               showNamesFeature = TRUE, main = branchPoint, 
                               silent = silent, treeheightFeature = 0, 
                               treeheightCell = 0, annotationCell = cellAnno))
        }
        else{
            
            #get features with non-zero variance to avoid error
            whichFeatures <- .removeZeroVariance(counts, cells = reorderedCells, 
                                                 markers = which(
                                                     featureLabels==branch$feature[1]))
            
            #create heatmap with only the split feature and split classes
            return(plotHeatmap(counts = counts, z = tree$classLabels,
                               y=featureLabels, featureIx=whichFeatures,
                               cellIx = reorderedCells, clusterCell = FALSE,
                               showNamesFeature = TRUE, main = branchPoint,
                               silent = silent, treeheightFeature = 0,
                               treeheightCell = 0, annotationCell = cellAnno))
        }
        
    }
    
    #if one-off split
    if(branch$statUsed[1] == "One-off"){
        
        #get unique classes
        classes <- unique(branch$class)
        
        #list which will contain final set of genes for heatmap
        whichFeatures <- c()
        
        #loop over unique classes
        for(class in classes){
            
            #subset table
            curClass <- branch[branch$class==class & branch$direction==1,]
            
            #if we have gene-level info in the tree
            if(("gene" %in% names(branch))){
                #get genes
                genes <- unique(curClass$gene)
                
                #keep top N features
                genes <- head(genes, topFeatures)
                
                #get gene indices
                markerGenes <- which(rownames(counts) %in% genes)
                
                #get features with non-zero variance to avoid error
                markerGenes <- .removeZeroVariance(counts, 
                                                   cells = which(
                                                       tree$classLabels %in%
                                                           unique(curClass$class)),
                                                   markers = markerGenes)
                
                #add to list of features
                whichFeatures <- c(whichFeatures, markerGenes)
            }
            else{
                
                #get features with non-zero variance to avoid error
                markerGenes <- .removeZeroVariance(
                    counts,
                    cells = which(tree$classLabels %in%
                                      unique(curClass$class)),
                    markers = which(featureLabels %in%
                                        unique(curClass$feature))
                )
                
                #add to list of features
                whichFeatures <- c(whichFeatures, markerGenes)
            }
        }
        
        #order the clusters such that up-regulated come first
        colOrder <- data.frame(groupName = unique(
            branch[order(branch$direction, decreasing = T),"class"]),
            groupIndex = seq_along(unique(branch$class)))
        
        #order the markers for clusters
        allMarkers <- setNames(as.list(colOrder$groupName), colOrder$groupName)
        allMarkers <- lapply(allMarkers, function(x){
            unique(branch[branch$class==x & branch$direction==1,"feature"])
        })
        rowOrder <- data.frame(groupName = unlist(allMarkers),
                               groupIndex = seq_along(unlist(allMarkers)))
        toRemove <- which(!rowOrder$groupName %in% featureLabels[whichFeatures])
        if(length(toRemove) > 0){
            rowOrder <- rowOrder[-toRemove,]
        }
        
        #sort cells according to metacluster size
        x <- tree$classLabels#[tree$classLabels %in% unique(branch$class)]
        y <- colOrder$groupName
        sortedCells <- seq(ncol(counts))[order(match(x,y))]
        sortedCells <- sortedCells[seq(sum(tree$classLabels %in% classes))]
        
        #create heatmap with only the split features and split classes
        return(plotHeatmap(counts = counts, z = tree$classLabels, 
                           y = featureLabels, featureIx=whichFeatures, 
                           cellIx = sortedCells,
                           showNamesFeature = TRUE, main = branchPoint, 
                           silent = silent, treeheightFeature = 0,
                           colGroupOrder = colOrder, rowGroupOrder = rowOrder, 
                           treeheightCell = 0))
    }
    
}

#helper function to identify zero-variance genes in a counts matrix
.removeZeroVariance <- function(counts, cells, markers){
    #subset counts matrix
    counts <- counts[, cells]
    
    #scale rows
    counts <- t(scale(t(counts)))
    
    #get indices of genes which have NA
    zeroVarianceGenes <- which(!complete.cases(counts))
    
    #find overlap between zero-variance genes and marker genes
    zeroVarianceMarkers <- intersect(zeroVarianceGenes, markers)
    
    #return indices of marker genes without zero-variance
    if(length(zeroVarianceMarkers) > 0)
        return(markers[-which(markers %in% zeroVarianceMarkers)])
    else
        return(markers)
}

