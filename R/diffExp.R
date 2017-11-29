#' @title Gene expression markers for cell clusters using MAST
#' @description Finds markers (differentially expressed genes) for cell clusters
#'    using MAST: a flexible statistical framework for assessing transcriptional
#'    changes and characterizing heterogeneity in single-cell RNA sequencing data
#'    (Finak et al, Genome Biology, 2015)
#' 
#' @param counts A numeric count matrix.
#' @param celda.mod An object of class celda_C or celda_CG.
#' @param c1 Numeric. Cell cluster(s) to define markers for.
#' @param c2 Numeric. Second cell cluster(s) for comparison. If NULL (default) - use all
#'    other cells for comparison.
#' @param only.pos Logical. Only return markers with positive log2fc (FALSE by default).
#' @param log2fc.threshold Numeric. Only return markers whose the absolute values of log2fc
#'    are greater than this threshold (0 by default).
#' @param fdr.threshold Numeric. Only return markers whose false discovery rates (FDRs) are less
#'    than this threshold (1 by default).
#' @return Data frame containing a ranked list (based on the absolute value of log2fc) of putative markers,
#'    and associated statistics (p-value, log2fc and FDR).
#' @export
diffExp_MAST <- function(counts, celda.mod, c1, c2 = NULL, only.pos = FALSE, log2fc.threshold = 0, fdr.threshold = 1) {
  if (is.null(counts)) {
    stop("'counts' should be a numeric count matrix")
  }
  if (is.null(celda.mod) || is.null(celda.mod$z)){
    stop("'celda.mod' should be an object of class celda_C or celda_CG")
  }
  if (is.null(c1)) {
    stop("'c1' should be a numeric vector of cell cluster(s)")
  }
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
  log_normalized_mat <- log2(normalizeCounts(mat) + 1)
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
  fcHurdleSig <-
    merge(fcHurdle[fdr < fdr.threshold &
                     abs(coef) > log2fc.threshold], data.table::as.data.table(GenomicRanges::mcols(sca)), by = 'primerid')
  fcHurdleSig <- fcHurdleSig[,-c(4,5)]
  names(fcHurdleSig)[c(1,3)] <- c("Gene", "log2fc")
  if(only.pos){
    fcHurdleSig <- fcHurdleSig[which(fcHurdleSig$log2fc > 0),]
  }
  fcHurdleSig <- fcHurdleSig[order(abs(fcHurdleSig$log2fc), decreasing = TRUE),]
  return(fcHurdleSig)
}

#' @title Gene expression markers for cell clusters using ZingeR and edgeR
#' @description Finds markers (differentially expressed genes) for cell clusters
#'    using zingeR: unlocking RNA-seq tools for zero-inflation and single cell 
#'    applications (Berge et al, bioRxiv, 2017)
#'
#' @param counts A numeric count matrix.
#' @param celda.mod An object of class celda_C or celda_CG.
#' @param c1 Numeric. Cell cluster(s) to define markers for.
#' @param c2 Numeric. Second cell cluster(s) for comparison. If NULL (default) - use all
#'    other cells for comparison.
#' @param correct_gene_detection Logical. Should the number of detected genes for
#'    each cell be a covariate in the zero-inflated negative binomial model? (FALSE by default)
#' @param only.pos Logical. Only return markers with positive log2fc (FALSE by default).
#' @param fdr.threshold Numeric. Only return markers whose false discovery rates (FDRs) are less
#'    than this threshold (1 by default).
#' @param max.iter Numeric. The number of iterations for the EM-algorithm (200 by default).
#' @return Data frame containing a ranked list (based on the absolute value of log2FC) of putative markers,
#'    and associated statistics (log2FC, logCPM, LR, PValue and FDR).
#' @export 
diffExp_zingeR_edgeR <- function(counts, celda.mod, c1, c2 = NULL, correct_gene_detection = FALSE, only.pos = FALSE, fdr.threshold = 1, max.iter = 200) {
  if (is.null(c1)) {
    stop("'c1' should be a numeric vector of cell cluster(s)")
  }
  if (is.null(c2)){
    c2 <- sort(setdiff(unique(celda.mod$z),c1))
  }
  filtered_counts <- counts[,celda.mod$z %in% c(c1,c2)]
  d=edgeR::DGEList(filtered_counts)
  d=suppressWarnings(edgeR::calcNormFactors(d))
  grp <- celda.mod$z[celda.mod$z %in% c(c1,c2)]
  if(length(c2) > 1){
    grp[which(grp %in% c2)] <- 0
  }else{
    grp[which(grp == c2)] <- 0
  }
  if(length(c1) > 1){
    grp[which(grp %in% c1)] <- 1
  }else{
    grp[which(grp == c1)] <- 1
  }
  if(correct_gene_detection){
    gene_counts_each_cell <- (sapply(1:length(celda.mod$z), function(i){sum(counts[,i] > 0)}))[celda.mod$z %in% c(c1,c2)]  
    design=model.matrix(~grp+gene_counts_each_cell)
  }else{
    design=model.matrix(~grp)
  }
  weights <- zingeR::zeroWeightsLS(counts=d$counts, design=design, maxit=max.iter, normalization="TMM")
  d$weights <- weights
  d=edgeR::estimateDisp(d, design)
  fit=edgeR::glmFit(d,design)
  lrt=zingeR::glmWeightedF(fit,coef=2, independentFiltering = TRUE)
  signif_genes <- lrt$table[lrt$table$padjFilter < fdr.threshold,]
  colnames(signif_genes)[5] <- "FDR"
  if(only.pos){
    signif_genes <- signif_genes[which(signif_genes$logFC > 0),]
  }
  signif_genes <- signif_genes[order(abs(signif_genes$logFC), decreasing = TRUE),]
  return(na.omit(signif_genes)) # remove rows with NAs
}