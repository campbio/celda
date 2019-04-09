#' @title Differential expression for cell subpopulations using MAST
#' @description Uses MAST to find differentially expressed features for
#'  specified cell subpopulations.
#' @param counts Integer matrix. Rows represent features and columns represent
#'  cells. This matrix should be the same as the one used to generate
#'  `celdaMod`.
#' @param celdaMod Celda object of class `celda_C` or `celda_CG`.
#' @param c1 Integer vector. Cell populations to include in group 1 for the
#'  differential expression analysis.
#' @param c2 Integer vector. Cell populations to include in group 2 for the
#'  differential expression analysis. If NULL, the clusters in the c1 group are
#'  compared to all other clusters. Default NULL.
#' @param onlyPos Logical. Whether to only return markers with positive log2
#'  fold change. Default FALSE.
#' @param log2fcThreshold Numeric. A number greater than 0 that specifies the
#'  absolute log2 fold change threshold. Only features with absolute value
#'  above this threshold will be returned. If NULL, this filter will not be
#'  applied. Default NULL.
#' @param fdrThreshold Numeric. A number between 0 and 1 that specifies the
#'  false discovery rate (FDR) threshold. Only features below this threshold
#'  will be returned. Default 1.
#' @return Data frame containing MAST results including statistics such as
#'  p-value, log2 fold change, and FDR.
#' @examples
#' library(SummarizedExperiment)
#' clusterDiffexpRes = differentialExpression(celdaCGSim$counts,
#'     celdaCGMod, c1 = c(1, 2))
#' @export
#' @import data.table plyr
#' @rawNamespace import(MAST, except = c(combine))
#' @rawNamespace import(SummarizedExperiment, except = c(shift, rowRanges))
differentialExpression <- function(counts,
    celdaMod,
    c1,
    c2 = NULL,
    onlyPos = FALSE,
    log2fcThreshold = NULL,
    fdrThreshold = 1) {

    if (!is.matrix(counts)) {
        stop("'counts' should be a numeric count matrix")
    }
    if (!(methods::is(celdaMod, "celda_C") ||
            methods::is(celdaMod, "celda_CG"))) {
        stop("'celdaMod' should be an object of class celda_C or celda_CG")
    }
    if (is.null(c1)) {
        stop("'c1' should be a numeric vector of cell cluster(s)")
    }
    compareCountMatrix(counts, celdaMod)

    if (is.null(c2)) {
        c2 <- sort(setdiff(unique(celdaMod@clusters$z), c1))
    }
    if (length(c1) > 1) {
        cells1 <- celdaMod@names$column[which(celdaMod@clusters$z %in% c1)]
    } else {
        cells1 <- celdaMod@names$column[which(celdaMod@clusters$z == c1)]
    }

    if (length(c2) > 1) {
        cells2 <- celdaMod@names$column[which(celdaMod@clusters$z %in% c2)]
    } else {
        cells2 <- celdaMod@names$column[which(celdaMod@clusters$z == c2)]
    }

    mat <- counts[, c(cells1, cells2)]
    log_normalized_mat <- normalizeCounts(mat,
        normalize = "cpm",
        transformationFun = log1p)
    cdat <- data.frame(wellKey = c(cells1, cells2),
        condition = c(rep("c1", length(cells1)), rep("c2", length(cells2))),
        ngeneson = rep("", (length(cells1) + length(cells2))),
        stringsAsFactors = FALSE)

    sca <- suppressMessages(MAST::FromMatrix(log_normalized_mat, cdat))
    cdr2 <- colSums(SummarizedExperiment::assay(sca) > 0)
    SummarizedExperiment::colData(sca)$cngeneson <- scale(cdr2)
    cond <- factor(SummarizedExperiment::colData(sca)$condition)
    cond <- stats::relevel(cond, "c2")
    SummarizedExperiment::colData(sca)$condition <- cond
    zlmCond <- MAST::zlm(~ condition + cngeneson, sca)
    summaryCond <- MAST::summary(zlmCond, doLRT = "conditionc1")
    summaryDt <- summaryCond$datatable
    contrast <-
        component <-
        primerid <-
        coef <-
        ci.hi <-
        ci.lo <-
        `Pr(>Chisq)` <-
        fdr <- NULL # Avoid NSE notes in check

    fcHurdle <- merge(summaryDt[contrast == "conditionc1" &
            component == "H", .(primerid, `Pr(>Chisq)`)],
        summaryDt[contrast == "conditionc1" & component == "logFC",
            .(primerid, coef, ci.hi, ci.lo)], by = "primerid")

    fcHurdle[, fdr := stats::p.adjust(`Pr(>Chisq)`, "fdr")]
    ### Some genes aren't outputted because log2FC gets NaN if one or both
    ### clusters have 0 counts for a gene
    ### and then they're discarded because NaN !> 0
    if (is.null(log2fcThreshold)) {
        fcHurdleSig <- fcHurdle
    } else {
        fcHurdleSig <- merge(fcHurdle[fdr < fdrThreshold &
                abs(coef) > log2fcThreshold],
            data.table::as.data.table(GenomicRanges::mcols(sca)),
            by = 'primerid')
        if (onlyPos) {
            fcHurdleSig <- fcHurdleSig[which(fcHurdleSig$log2fc > 0), ]
        }
    }
    fcHurdleSig <- fcHurdleSig[, -c(4, 5)]
    names(fcHurdleSig)[c(1, 2, 3, 4)] <- c("Gene", "Pvalue", "Log2_FC", "FDR")
    fcHurdleSig <- fcHurdleSig[order(fcHurdleSig$Pvalue, decreasing = FALSE), ]

    return(fcHurdleSig)
}
