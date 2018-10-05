#' @title Gene expression markers for cell clusters using MAST
#' @description Finds markers (differentially expressed genes) for cell clusters
#'    using MAST: a flexible statistical framework for assessing transcriptional
#'    changes and characterizing heterogeneity in single-cell RNA sequencing data
#'    (Finak et al, Genome Biology, 2015)
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda model. Options are "celda_C" or "celda_CG". 
#' @param c1 Integer vector. Cell populations to include in group 1 for the differential expression analysis. 
#' @param c2 Integer vector. Cell populations to include in group 2 for the differential expression analysis. If NULL, everything in c1 is compared to all other clusters. Default NULL. 
#' @param only.pos Logical. Whether to only return markers with positive log2fc. Default FALSE.

#' @param log2fc.threshold Numeric. A number greater than 0 that specifies the absolute log2 fold change threshold. Only features with absolute value above this threshold will be returned.
#' @param fdr.threshold Numeric. A number between 0 and 1 that specifies the false discovery rate (FDR) threshold. Only features below this threshold will be returned.
#' @return Data frame containing a ranked list (based on the absolute value of log2fc) of putative markers,
#'    and associated statistics (p-value, log2fc and FDR).
#' @export
#' @import data.table
differentialExpression <- function(counts, celda.mod, c1, c2 = NULL, only.pos = FALSE, log2fc.threshold = NULL, fdr.threshold = 1) {
  if (is.null(counts)) {
    stop("'counts' should be a numeric count matrix")
  }
  if (is.null(celda.mod) || is.null(celda.mod$z)){
    stop("'celda.mod' should be an object of class celda_C or celda_CG")
  }
  if (is.null(c1)) {
    stop("'c1' should be a numeric vector of cell cluster(s)")
  }
  compareCountMatrix(counts, celda.mod)
  
  if (is.null(c2)){
    c2 <- sort(setdiff(unique(celda.mod$z),c1))
  }
  if (length(c1) > 1){
    cells1 <-
      celda.mod$names$column[which(celda.mod$z %in% c1)]
  }
  else{
    cells1 <-
      celda.mod$names$column[which(celda.mod$z == c1)]
  }
  if (length(c2) > 1){
    cells2 <-
      celda.mod$names$column[which(celda.mod$z %in% c2)]
  }
  else{
    cells2 <-
      celda.mod$names$column[which(celda.mod$z == c2)]
  }
  mat <- counts[,c(cells1,cells2)]
  log_normalized_mat <- normalizeCounts(mat, normalize="cpm", transformation.fun=log1p)
  cdat <-
    data.frame(
      wellKey = c(cells1, cells2),
      condition = c(rep("c1", length(cells1)), rep("c2", length(cells2))),
      ngeneson = rep("", (length(cells1) + length(cells2))),
      stringsAsFactors = FALSE
    )
  sca <- MAST::FromMatrix(log_normalized_mat, cdat)
  cdr2 <- colSums(SummarizedExperiment::assay(sca) > 0)
  SummarizedExperiment::colData(sca)$cngeneson <- scale(cdr2)
  cond <- factor(SummarizedExperiment::colData(sca)$condition)
  cond <- relevel(cond, "c2")
  SummarizedExperiment::colData(sca)$condition <- cond
  zlmCond <- MAST::zlm( ~ condition + cngeneson, sca)
  summaryCond <- MAST::summary(zlmCond, doLRT = 'conditionc1')
  summaryDt <- summaryCond$datatable
  fcHurdle <-
    merge(summaryDt[contrast == 'conditionc1' &
                      component == 'H', .(primerid, `Pr(>Chisq)`)],
          summaryDt[contrast == 'conditionc1' &
                      component == 'logFC', .(primerid, coef, ci.hi, ci.lo)], by = 'primerid')
  
  fcHurdle[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
  ###Some genes aren't outputted because log2FC gets NaN if one or both clusters have 0 counts for a gene
  ###and then they're discarded because NaN !> 0
  if(is.null(log2fc.threshold)){
   fcHurdleSig <-fcHurdle  
  }else{
  fcHurdleSig <-
    merge(fcHurdle[fdr < fdr.threshold &
                     abs(coef) > log2fc.threshold], data.table::as.data.table(GenomicRanges::mcols(sca)), by = 'primerid')
    if(only.pos){
      fcHurdleSig <- fcHurdleSig[which(fcHurdleSig$log2fc > 0),]
    }
  }
  fcHurdleSig <- fcHurdleSig[,-c(4,5)]
  names(fcHurdleSig)[c(1,3)] <- c("Gene", "log2fc")
  fcHurdleSig <- fcHurdleSig[order(abs(fcHurdleSig$log2fc), decreasing = TRUE),]
  return(fcHurdleSig)
}

