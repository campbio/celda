#' @title Plots percentage of cells with their cell type markers
#'
#' @description Barplot that shows the percentage of cells within subpopulations
#' with detectable levels of given marker genes
#' @param  counts Matrix of counts
#' @param  z a vector. Cell cluster label, same length as the column numbers of the count matrix
#' @param  geneMarkers a dataframe of marker genes which corresponding to their cell types
#' @param  threshold default as 1.
#' @param  color default as "red3"
#' @param  textLabelSize default as 3
#' @param  precision default as 2. Precision of percentage of cells showing the marker gene shown on the barplot
#' @export

decontXMarkerPlot <- function(counts, z, geneMarkers, threshold = 1, color = "red3", textLabelSize = 3, precision = 2) {
  z_names <- levels(factor(z))
  z <- .processZ(z)

  pct_CTbyZ <- .calCellPct(counts, z, geneMarkers, threshold)
  colnames(pct_CTbyZ) <- z_names

  pct_CTbyZ.melt <- reshape2::melt(pct_CTbyZ, varnames = c("cellType", "z"), value.name = "percent")

  plt <- ggplot2::ggplot(pct_CTbyZ.melt, ggplot2::aes(x = z, y = pct_CTbyZ.melt$percent * 100)) +
    ggplot2::geom_bar(stat = "identity", fill = color) +
    ggplot2::geom_text(aes(x = z, y = pct_CTbyZ.melt$percent * 100 + 5, label = paste0(round(pct_CTbyZ.melt$percent, precision) * 100, "%")), size = textLabelSize) +
    ggplot2::xlab("Cluster") +
    ggplot2::ylab("Percentage of cells expressing cell-type\nspecific markers") +
    ggplot2::facet_grid(. ~ cellType) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = "grey"),
      panel.grid = ggplot2::element_line("grey"), legend.position = "none",
      legend.key = ggplot2::element_rect(fill = "white", color = "white"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      text = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 9),
      legend.key.size = grid::unit(8, "mm"),
      legend.text = ggplot2::element_text(size = 10),
      strip.text.x = ggplot2::element_text(size = 10)
    )
  return(plt)
}


.calCellPct <- function(counts, z, geneMarkers, threshold = 1) {
  rNames <- rownames(counts)
  gNames <- geneMarkers[, "geneMarkers"]

  # Only use intersected marker genes
  intersectedGenes <- intersect(rNames, gNames)

  if (length(intersectedGenes) == 0) {
    stop("Cannot find marker genes in 'counts'. Make sure the name of the marker genes in 'geneMarker' consistent with row-names in 'count'.")
  }

  sub_counts <- counts[rNames %in% intersectedGenes, ]
  geneMarkers <- geneMarkers[gNames %in% intersectedGenes, ]

  geneMarkers <- .geneMarkerProcess(geneMarkers, orders = rownames(sub_counts))

  if (is.null(nrow(sub_counts))) {
    # When there is only one gene --> a vector
    pct_CTbyZ <- .vctrRowProjectColCollapse(sub_vector = sub_counts, genePresented = geneMarkers, z = z, threshold = threshold)
  } else {
    # When multiple genes --> matrix
    pct_CTbyZ <- .mtxRowProjectColCollapse(sub_counts = sub_counts, genePresented = geneMarkers, z = z, threshold = threshold)
  }

  return(pct_CTbyZ)
}


# Collapse (1-)Gene x Cell count vector into cell-type-marker by cluster matrix
# with each element being the cells% in that cluster shows at least one marker in that cell type
.vctrRowProjectColCollapse <- function(sub_vector, genePresented, z, threshold = 1) {
  K <- max(z)
  binary_2byZ <- table(sub_vector >= threshold, z)
  CTnames <- genePresented[, "cellType"]
  pct_CTbyZ <- sweep(binary_2byZ, MARGIN = 2, STATS = table(z), FUN = "/")[rownames(binary_2byZ) == TRUE, ]

  if (length(CTnames) > 1) {
    # Multiple cell types have this marker
    pct_CTbyZ <- matrix(rep(pct_CTbyZ, length(CTnames)), nrow = length(CTnames), dimnames = list(CTnames, 1:K))
  } else {
    # Only one cell type has this marker
    pct_CTbyZ <- matrix(pct_CTbyZ, nrow = 1, dimnames = list(CTnames, 1:K))
  }
  return(pct_CTbyZ)
}

# Collapse Gene x Cell count matrix into cell-type-marker by cluster matrix
# with each element being the cells% in that cluster shows at least one marker in that cell type
.mtxRowProjectColCollapse <- function(sub_counts, genePresented, z, threshold = 1) {
  K <- max(z)

  nTC <- length(levels(genePresented[, "cellType"]))
  nZ <- length(unique(z))
  nCTbyZ <- matrix(0, nrow = nTC, ncol = nZ)

  # convert matrix into dgCMatrix if it is not
  ij_pair <- nonzero(sub_counts)
	index_filter = ij_pair[["val"]] >= threshold
	iR = ij_pair[["row"]] = ij_pair[["row"]][index_filter]
	iC = ij_pair[["col"]] = ij_pair[["col"]][index_filter]
  i_celltype <- plyr::mapvalues(iR, from = genePresented[, "geneMarkers"], to = genePresented[, "cellType"])

  if (nrow(genePresented) == length(unique(genePresented[, "geneMarkers"]))) {
    # When each gene qs only specified in ONE cell type
  } else {
    # When at least one gene is specified as marker for multiple cell types

    duplicateTF <- duplicated(genePresented[, "geneMarkers"]) 
    # assume if a gene has shown more than once,
    # it is because this gene is also a marker for another cell type
		# TODO: make sure there is no same gene same celltype pair has shown more than once

    duplicatedMarker <- genePresented[duplicateTF, ]

    for (r in 1:nrow(duplicatedMarker)) {
      celltype <- duplicatedMarker[r, "cellType"]
      gene <- duplicatedMarker[r, "geneMarkers"]

      duplicated_pair <- ij_pair[["row"]] == gene
			iR <- c(iR, ij_pair[["row"]][duplicated_pair])
			iC <- c(iC, ij_pair[["col"]][duplicated_pair]) 
      i_celltype <- c(i_celltype, celltype)
    }
  }

  CTnames <- levels(factor(genePresented[, "cellName"]))
  Cnames <- colnames(sub_counts)
  ng_CTbyC <- Matrix::sparseMatrix(i = i_celltype, j = iC, x = 1, giveCsparse = TRUE, dimnames = list(CTnames, Cnames), dims = c(nTC, ncol(sub_counts)))
  binary_CTbyC <- methods::as((ng_CTbyC > 0) * 1, "matrix")
  storage.mode(binary_CTbyC) <- "integer"
  nC_CTbyZ <- .colSumByGroup(binary_CTbyC, z, K)
  rownames(nC_CTbyZ) <- rownames(ng_CTbyC)
  colnames(nC_CTbyZ) <- 1:K

  pct_CTbyZ <- sweep(nC_CTbyZ, MARGIN = 2, STATS = table(z), FUN = "/")

  return(pct_CTbyZ)
}


# geneMarkers should be a dataframe,  w/ 2 column names being `cellType` and `geneMarkers`
# convert both `cellType` and `geneMarkers` are factors with levels being integer
.geneMarkerProcess <- function(geneMarkers, orders = NULL) {
  geneMarkers[, "cellName"] <- geneMarkers[, "cellType"]
  geneMarkers[, "cellType"] <- factor(geneMarkers[, "cellType"])
  levels(geneMarkers[, "cellType"]) <- 1:length(levels(geneMarkers[, "cellType"]))

  geneMarkers[, "geneName"] <- geneMarkers[, "geneMarkers"]
  if (is.null(orders)) {
    geneMarkers[, "geneMarkers"] <- factor(geneMarkers[, "geneMarkers"])
  } else {
    geneMarkers[, "geneMarkers"] <- factor(geneMarkers[, "geneMarkers"], levels = orders)
  }
  levels(geneMarkers[, "geneMarkers"]) <- 1:length(levels(geneMarkers[, "geneMarkers"]))

  return(geneMarkers)
}

# Convert z to be factor with levels being integer
.processZ <- function(z) {
  z <- factor(z)
  levels(z) <- seq_len(length(levels(z)))
  z <- as.integer(z)
  return(z)
}

