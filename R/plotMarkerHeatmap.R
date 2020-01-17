#' @title Generate heatmap for a marker decision tree.
#' @description Creates heatmap for a specified branch point in a marker tree.
#' @param tree A decision tree from CELDA's findMarkersTree function.
#' @param counts Numeric matrix. Gene-by-cell counts matrix.
#' @param branchPoint Character. Name of branch point to plot heatmap for.
#' Name should match those in tree$branchPoints.
#' @param featureLabels List of feature cluster assignments. Length should
#' be equal to number of rows in counts matrix, and formatting should match
#' that used in findMarkersTree(). Required when using clusters
#' of features and not previously provided to findMarkersTree()
#' @param topFeatures Integer. Number of genes to plot per marker module.
#' Genes are sorted based on their AUC for their respective cluster. 
#' Default is 10.
#' @param silent Logical. Whether to avoid plotting heatmap to screen.
#' Default is FALSE.
#' @return A heatmap visualizing the counts matrix for the cells and genes at
#' the specified branch point.
#' #' @examples
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
#' # Plot example heatmap
#' plotMarkerHeatmap(DecTree, sim_counts$counts, branchPoint = "level_2",
#' featureLabels = paste0("L",clusters(cm)$y))
#' 
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
  
  #check if we can get individual genes from tree
  if(!("gene" %in% names(branch))){
    message("NOTE: Unable to find scores for individual features within feature
            categories. Re-run findMarkersTree() with the `counts` and
            `featureLabels` parameters set to score individual features.")
  }
  
  #if top-level in metaclusters tree
  if(branchPoint == "top_level"){
    #get unique metaclusters
    metaclusters <- unique(branch$class)
    
    #list which will contain final set of genes for heatmap
    whichFeatures <- c()
    
    #loop over unique metaclusters
    for(meta in metaclusters){
      
      #subset table
      curMeta <- branch[branch$class==meta,]
      
      #if we have gene-level info in the tree
      if(("gene" %in% names(branch))){
        #sort by gene AUC score
        curMeta <- curMeta[order(curMeta$geneAUC, decreasing = TRUE),]
        
        #get genes
        genes <- unique(curMeta$gene)
        
        #keep top N features
        genes <- head(genes, topFeatures)
        
        #get gene indices
        markerGenes <- which(rownames(counts) %in% genes)
        
        #get features with non-zero variance to avoid error
        markerGenes <- .removeZeroVariance(counts, 
                                           cells = which(
                                             tree$metaclusterLabels %in%
                                               unique(curMeta$class)),
                                           markers = markerGenes)
        
        #add to list of features
        whichFeatures <- c(whichFeatures, markerGenes)
      }
      else{
        #get marker gene indices
        markerGenes <- which(featureLabels == marker)
        
        #get features with non-zero variance to avoid error
        markerGenes <- .removeZeroVariance(counts, 
                                           cells = which(
                                             tree$metaclusterLabels %in%
                                               unique(curMarker$class)),
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
      unique(branch[branch$class==x,"feature"])
    })
    rowOrder <- data.frame(groupName = unlist(allMarkers),
                           groupIndex = seq_along(unlist(allMarkers)))
    rowOrder <- rowOrder[-which(
      !rowOrder$groupName %in% tree$featureLabels[whichFeatures]),]
    
    #create heatmap with only the markers
    return(plotHeatmap(counts = counts, z = tree$metaclusterLabels,
                              y = tree$featureLabels, featureIx=whichFeatures,
                              showNamesFeature = TRUE, main = "Top-level",
                              silent = silent, treeheightFeature = 0,
                              colGroupOrder = colOrder, rowGroupOrder = rowOrder,
                              treeheightCell = 0))
  }
  
  #if balanced split
  if(branch$statUsed[1] == "Split"){
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
                                y=tree$featureLabels, featureIx=whichFeatures,
                                cellIx = reorderedCells, clusterCell = FALSE,
                                showNamesFeature = TRUE, main = branchPoint, 
                                silent = silent, treeheightFeature = 0, 
                                treeheightCell = 0))
    }
    else{
      #get features with non-zero variance to avoid error
      whichFeatures <- .removeZeroVariance(counts, cells = reorderedCells, 
                                           markers = which(
                                             featureLabels==branch$feature[1]))
      
      #create heatmap with only the split feature and split classes
      return(plotHeatmap(counts = counts, z = tree$classLabels,
                                y=tree$featureLabels, featureIx=whichFeatures,
                                cellIx = reorderedCells, clusterCell = FALSE,
                                showNamesFeature = TRUE, main = branchPoint,
                                silent = silent, treeheightFeature = 0,
                                treeheightCell = 0))
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
    rowOrder <- rowOrder[-which(
      !rowOrder$groupName %in% tree$featureLabels[whichFeatures]),]
    
    #create heatmap with only the split features and split classes
    return(plotHeatmap(counts = counts, z = tree$classLabels, 
                       y = tree$featureLabels, featureIx=whichFeatures, 
                       cellIx = which(tree$classLabels
                                      %in% unique(branch$class)),
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
