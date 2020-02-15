## Barplot that shows the percentage of cells representing given marker genes
## input: 

#' Barplot that shows the percentage of cells representing given marker genes
#' @param  counts 
#' @param  z
#' @param  cellType
#' @param  geneMarkers, a list of marker genes which corresponding to their cell types  

.celdabarplot = function(counts, z, cellType, geneMarkers){

	rowNames = rownames(counts)

	# Subset count matrix by gene names 
	gene_names = unique(geneMarkers[, "geneMarkers"])
	sub_counts = counts[ rowNames %in% gene_names,   ]
	if (nrow(sub_counts) == 0) {
		stop("Cannot find marker genes in 'counts'. Make sure the name of the marker genes in 'geneMarker' consistent with row-names in 'count'.")
	}

	if (is.null(nrow(sub_counts))) {
		# There is only one gene after subsetting 
		# which resulted into a vector, instead of a matrix

		# collapse gene by cell matrix into cell-type-marker by cell-type matrix
		nGenesByC = sub_counts

	} else {
       # collapse gene by cell matrix into cell-type-marker by cell-type matrix 
       # with each element being the number of cells in that cell-type presented any marker in that cell type
       nGenesByC = collapseRowByGeneMarker( counts = sub_counts, geneMarkers = geneMarkers, genesProvided = genesProvided )
	}

  





#geneMarkers should be a dataframe
geneMarkers = data.frame( cellType = c(rep("Tcells", 3), rep("Bcells", 3), "DC"), geneMarkers = 1:7) # string as factor
geneMarkers[, "cellTypeFactr"] = geneMarkers["cellType"]
levels(geneMarkers[, "cellTypeFactr"]) = 1:3

collapseRowByGeneMarker = function(sub_counts, geneMarkers, z, threshold = 1) {
	genePresented = geneMarkers[geneMarkers[, "geneMarkers"] %in% rownames(sub_counts),  ]  # only use marker genes that presented in the count matrix

	names(genePresented) 

  nTC = length(unique(genePresented[, "cellType"]))
	nZ = length(unique(z))
	nCTbyZ = matrix(0, nrow=nTC ncol=nZ)

	if (length(genePresented) == length(unique(genePresented))){
	# When each gene is only specified in ONE cell type
		rowLabel = plyr::mapvalues(rownames(sub_counts), from = genePresented, to = names(genePresented) )
		# convert matrix into dgCMatrix if it is not
		mtx = as( sub_counts > threshold, "dgCMatrix")

		ij_pair = Ringo::nonzero(mtx)
		i_celltype = plyr::mapvalues(ij_pair[,"row"], from=geneMarkers[,"geneMarkers"], to=geneMarkers[,"cellType"])
		#mtx_binary= sparseMatrix( i = ij_pair[, "row"], j = ij_pair[, "col"], x = 1, giveCsparse=TRUE, dimnames = list( rownames(sub_counts), colnames(sub_counts)))
		ng_CTbyC = sparseMatrix( i = i_celltype, j = ij_pair[, "col"], x = 1, giveCsparse=TRUE, dimnames = list( rownames(sub_counts), colnames(sub_counts)), dims=c(nTC, ncol(sub_counts)) )
		ng_CTbyC = as(ng_CTbyC, "matrix") 
		storage.mode(ng_CTbyC) = "integer"
		ng_CTbyZ = .colSumByGroup(ng_CTbyC, z)


	} else {
	# When at least one gene is specified as marker for multiple cell types

	}

	return (ng_CTbyZ)
}


cdf = function( name, a ){ data.frame("celltype"= name, "geneMarkers" = a[name]) }


.geneMarkerProcess = function( geneMarkers ){
	cellTypes = names(geneMarkers)
	if (length(unique(cellTypes)) != length(cellTypes) ) {
		stop("Some cell type(s) are listed twice in provided 'geneMarkers' list")
	}
	return(unlist(geneMarkers))
}
