## Barplot that shows the percentage of cells representing given marker genes
## input: 

#' Barplot that shows the percentage of cells representing given marker genes
#' @param  counts 
#' @param  z
#' @param  geneMarkers, a dataframe of marker genes which corresponding to their cell types  
#' @param  threshold, default as 
#' @param  color, default as "red3"
#' @param  textLabelSize, default as 3

celdaMarkerPlot = function(counts, z, geneMarkers, threshold = 0, color = "red3", textLabelSize=3 ){
    nC_CTbyZ = .celdabarplot(counts, z, geneMarkers, threshold)

    nC_CTbyZ.melt = reshape2::melt(nC_CTbyZ, varnames = c("cellType", "z"), value.name = "percent")

    plt = ggplot2::ggplot(nC_CTbyZ.melt, ggplot2::aes( x = z, y = percent * 100) ) + 
			geom_bar(stat = "identity", fill = color ) +
		  geom_text(aes(x = z, y= percent * 100 + 5 ,label = paste0(round(percent,4) * 100, "%") ) , size = textLabelSize ) +
			facet_grid(. ~ cellType ) + 
			theme(panel.background=element_rect(fill="white", color="grey"),
						panel.grid = element_line("grey"), legend.position="none",
						legend.key=element_rect(fill="white", color="white"),
						panel.grid.minor = element_blank(),
						panel.grid.major = element_blank(),
						text=element_text(size=10),
						axis.text.x = element_text(size=8, angle = 45, hjust = 1),
						axis.text.y=element_text(size=9),
						legend.key.size = grid::unit(8, "mm"),
						legend.text = element_text(size=10),
						strip.text.x = element_text(size = 10))
		return(plt)
}

.celdabarplot = function(counts, z, geneMarkers, threshold = 4){

  rNames = rownames(counts)
	gNames = geneMarkers[, "geneMarkers"]

  # Only use intersected genes 
	intersectedGenes = intersect(rNames, gNames)

	if (length(intersectedGenes) == 0) {
		stop("Cannot find marker genes in 'counts'. Make sure the name of the marker genes in 'geneMarker' consistent with row-names in 'count'.")
	}

	sub_counts = counts[rNames %in% intersectedGenes, ]
	geneMarkers = geneMarkers[gNames %in% intersectedGenes, ]

	geneMarkers = .geneMarkerProcess(geneMarkers) 

	if (is.null(nrow(sub_counts))) {
		# There is only one gene after subsetting 
		# which resulted into a vector, instead of a matrix
		nC_CTbyZ = (sub_counts > threshold) * 1

	} else {
       # collapse gene by cell matrix into cell-type-marker by cell-type matrix,
       # with each element being the number of cells in that cell-type presented any marker in that cell type
       nC_CTbyZ = collapseRowByGeneMarker( sub_counts = sub_counts, genePresented = geneMarkers, z = z, threshold = threshold)
	}

	return(nC_CTbyZ)
}
  


collapseRowByGeneMarker = function(sub_counts, genePresented, z, threshold = 1) {
	K = max(z)

  nTC = length(levels(genePresented[, "cellType"]))
	nZ = length(unique(z))
	nCTbyZ = matrix(0, nrow=nTC,  ncol=nZ)

	# convert matrix into dgCMatrix if it is not
	mtx = as( sub_counts > threshold, "dgCMatrix")
	ij_pair = Ringo::nonzero(mtx)
	i_celltype = plyr::mapvalues(ij_pair[,"row"], from=genePresented[,"geneMarkers"], to=genePresented[,"cellType"])

	if (nrow(genePresented) == length(unique(genePresented[, "geneMarkers"]))){
	# When each gene is only specified in ONE cell type

	} else {
	# When at least one gene is specified as marker for multiple cell types

		duplicateTF = duplicated(genePresented[, "geneMarkers"])
    duplicatedMarker = genePresented[duplicateTF, ]

    for ( r in 1:nrow(duplicatedMarker) ) {
			  celltype = duplicatedMarker[ r, "cellType"]
		    gene = duplicatedMarker[ r, "geneMarkers"]

			  duplicated_pair = ij_pair[, "row"] == gene
        ij_pair = rbind(ij_pair, ij_pair[duplicated_pair, ])
				i_celltype = c(i_celltype, celltype)
		}
	}

    CTnames = levels(factor(genePresented[, "cellName"]))
		Cnames = colnames(sub_counts)
		ng_CTbyC = sparseMatrix( i = i_celltype, j = ij_pair[, "col"], x = 1, giveCsparse=TRUE, dimnames = list(CTnames, Cnames), dims=c(nTC, ncol(sub_counts)) )
		binary_CTbyC = as((ng_CTbyC > 0) * 1, "matrix") 
		storage.mode(binary_CTbyC) = "integer"
		nC_CTbyZ = celda:::.colSumByGroup(binary_CTbyC, z, K)
		rownames(nC_CTbyZ) = rownames(ng_CTbyC)
		colnames(nC_CTbyZ) = 1:K

		pct_CTbyZ = sweep(nC_CTbyZ, MARGIN=2, STATS=table(z), FUN="/")
	#return (nC_CTbyZ)
	return(pct_CTbyZ)
}



.geneMarkerProcess = function( geneMarkers ){
    # geneMarkers should be a dataframe,  column names should be `cellType` and `geneMarkers`
	  # `cellType` should be factor + integer
	  # `geneMarkers` should be factor + integer
        geneMarkers[, "cellName"] = geneMarkers[, "cellType"]
				geneMarkers[, "cellType"] = factor(geneMarkers[, "cellType"])
				levels(geneMarkers[, "cellType"]) = 1:length(levels(geneMarkers[, "cellType"]))

				geneMarkers[, "geneName"] = geneMarkers[, "geneMarkers"]
				geneMarkers[, "geneMarkers"] = factor(geneMarkers[, "geneMarkers"])
				levels(geneMarkers[, "geneMarkers"]) = 1:length(levels(geneMarkers[, "geneMarkers"]))

    return(geneMarkers)
}

# 
#geneMarkers should be a dataframe
#counts = matrix(1:70, nrow=7, dimnames=list(1:7, NULL))
#geneMarkers = data.frame( cellType = c(rep("Tcells", 3), rep("Bcells", 3), "DC"), geneMarkers = 1:7) # string as factor
#z = c(rep(1, 4), rep(2,4), rep(3,2))
#a = .celdabarplot( counts = counts, z = z, geneMarkers = geneMarkers)
#plt = celdaMarkerPlot(counts = counts, z = z, geneMarkers = geneMarkers)
#plt
