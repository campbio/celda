#' @title Differential expression for cell subpopulations using MAST
#' @description Uses MAST to find differentially expressed features for specified cell subpopulations.
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Celda object of class `celda_C` or `celda_CG`. 
#' @param c1 Integer vector. Cell populations to include in group 1 for the differential expression analysis. 
#' @param c2 Integer vector. Cell populations to include in group 2 for the differential expression analysis. If NULL, the clusters in the c1 group are compared to all other clusters. Default NULL. 
#' @param only.pos Logical. Whether to only return markers with positive log2 fold change. Default FALSE.
#' @param log2fc.threshold Numeric. A number greater than 0 that specifies the absolute log2 fold change threshold. Only features with absolute value above this threshold will be returned. If NULL, this filter will not be applied. Default NULL.
#' @param fdr.threshold Numeric. A number between 0 and 1 that specifies the false discovery rate (FDR) threshold. Only features below this threshold will be returned. Default 1. 
#' @return Data frame containing MAST results including statistics such as p-value, log2 fold change, and FDR.
#' @examples
#' cluster.diffexp.res = differentialExpression(celda.CG.sim$counts, celda.CG.mod, c1=c(1,2))
#' @export
#' @import data.table
#' @import plyr
differentialExpression <- function(counts, celda.mod, c1, c2 = NULL, only.pos = FALSE, log2fc.threshold = NULL, fdr.threshold = 1) {
  if (!is.matrix(counts)) {
    stop("'counts' should be a numeric count matrix")
  }
  if (!(methods::is(celda.mod, "celda_C") || methods::is(celda.mod, "celda_CG"))){
    stop("'celda.mod' should be an object of class celda_C or celda_CG")
  }
  if (is.null(c1)) {
    stop("'c1' should be a numeric vector of cell cluster(s)")
  }
  compareCountMatrix(counts, celda.mod)
  
  if (is.null(c2)){
    c2 <- sort(setdiff(unique(celda.mod@clusters$z),c1))
  }
  if (length(c1) > 1){
    cells1 <-
      celda.mod@names$column[which(celda.mod@clusters$z %in% c1)]
  }else{
    cells1 <-
      celda.mod@names$column[which(celda.mod@clusters$z == c1)]
  }
  if (length(c2) > 1){
    cells2 <-
      celda.mod@names$column[which(celda.mod@clusters$z %in% c2)]
  }else{
    cells2 <-
      celda.mod@names$column[which(celda.mod@clusters$z == c2)]
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
  sca <- suppressMessages(MAST::FromMatrix(log_normalized_mat, cdat))
  cdr2 <- colSums(SummarizedExperiment::assay(sca) > 0)
  SummarizedExperiment::colData(sca)$cngeneson <- scale(cdr2)
  cond <- factor(SummarizedExperiment::colData(sca)$condition)
  cond <- stats::relevel(cond, "c2")
  SummarizedExperiment::colData(sca)$condition <- cond
  zlmCond <- MAST::zlm( ~ condition + cngeneson, sca)
  summaryCond <- MAST::summary(zlmCond, doLRT = 'conditionc1')
  summaryDt <- summaryCond$datatable
  contrast <- component <- primerid <- coef <- ci.hi <- ci.lo <- `Pr(>Chisq)` <- fdr <- NULL # Avoid NSE notes in check
  fcHurdle <-
    merge(summaryDt[contrast == 'conditionc1' &
                      component == 'H', .(primerid, `Pr(>Chisq)`)],
          summaryDt[contrast == 'conditionc1' &
                      component == 'logFC', .(primerid, coef, ci.hi, ci.lo)], by = 'primerid')
  
  fcHurdle[, fdr := stats::p.adjust(`Pr(>Chisq)`, 'fdr')]
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
  names(fcHurdleSig)[c(1,2,3,4)] <- c("Gene","Pvalue","Log2_FC","FDR")
  fcHurdleSig <- fcHurdleSig[order(fcHurdleSig$Pvalue, decreasing = FALSE),]
  return(fcHurdleSig)
}

