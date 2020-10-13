#' @title Differential expression for cell subpopulations using MAST
#' @description Uses MAST to find differentially expressed features for
#'  specified cell subpopulations.
#' @param x A numeric \link{matrix} of counts or a
#'  \linkS4class{SingleCellExperiment}
#'  with the matrix located in the assay slot under \code{useAssay}.
#'  Rows represent features and columns represent cells. Must contain cluster
#'  labels in \code{celdaClusters(x, altExpName = altExpName)} if \code{x} is a
#'  \linkS4class{SingleCellExperiment} object.
#' @param useAssay A string specifying which \link[SummarizedExperiment]{assay}
#'  slot to use if \code{x} is a
#'  \link[SingleCellExperiment]{SingleCellExperiment} object. Default "counts".
#' @param altExpName The name for the \link[SingleCellExperiment]{altExp} slot
#'  to use. Default "featureSubset".
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
#' @param ... Ignored. Placeholder to prevent check warning.
#' @return Data frame containing MAST results including statistics such as
#'  p-value, log2 fold change, and FDR.
#' @export
#' @rawNamespace import(data.table, except = c(melt, shift))
#' @importFrom MAST FromMatrix
#' @importFrom MAST zlm
#' @importFrom MAST summary
#' @importFrom S4Vectors mcols
#' @importFrom plyr .
#' @import SummarizedExperiment
setGeneric("differentialExpression", function(x, ...) {
    standardGeneric("differentialExpression")})


#' @rdname differentialExpression
#' @examples
#' data(sceCeldaCG)
#' clusterDiffexpRes <- differentialExpression(sceCeldaCG, c1 = c(1, 2))
#' @export
setMethod("differentialExpression",
    signature(x = "SingleCellExperiment"),
    function(x,
        useAssay = "counts",
        altExpName = "featureSubset",
        c1,
        c2 = NULL,
        onlyPos = FALSE,
        log2fcThreshold = NULL,
        fdrThreshold = 1) {

        if (is.null(c1)) {
            stop("'c1' should be a numeric vector of cell cluster(s)")
        }

        cdiff <- setdiff(c1, celdaClusters(x, altExpName = altExpName))

        if (length(cdiff) > 0) {
            warning("cluster ", cdiff, "in 'c1' does not exist in",
                " 'celdaClusters(x, altExpName = altExpName)'!")
            if (length(cdiff) == length(c1)) {
                stop("All clusters in 'c1' does not exist in",
                    " 'celdaClusters(x, altExpName = altExpName)'!")
            }
        }

        altExp <- SingleCellExperiment::altExp(x, altExpName)
        counts <- SummarizedExperiment::assay(x, i = useAssay)

        if (is.null(c2)) {
            c2 <- sort(setdiff(unique(celdaClusters(x,
                altExpName = altExpName)), c1))
        }
        if (length(c1) > 1) {
            cells1 <- SummarizedExperiment::colData(altExp)$colnames[
                which(celdaClusters(x, altExpName = altExpName) %in% c1)]
        } else {
            cells1 <- SummarizedExperiment::colData(altExp)$colnames[
                which(celdaClusters(x, altExpName = altExpName) == c1)]
        }

        if (length(c2) > 1) {
            cells2 <- SummarizedExperiment::colData(altExp)$colnames[
                which(celdaClusters(x, altExpName = altExpName) %in% c2)]
        } else {
            cells2 <- SummarizedExperiment::colData(altExp)$colnames[
                which(celdaClusters(x, altExpName = altExpName) == c2)]
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
                data.table::as.data.table(S4Vectors::mcols(sca)),
                by = "primerid")
            if (onlyPos) {
                fcHurdleSig <- fcHurdleSig[which(fcHurdleSig$log2fc > 0), ]
            }
        }
        fcHurdleSig <- fcHurdleSig[, -c(4, 5)]
        names(fcHurdleSig)[c(1, 2, 3, 4)] <- c("Gene", "Pvalue", "Log2_FC",
            "FDR")
        fcHurdleSig <- fcHurdleSig[order(fcHurdleSig$Pvalue,
            decreasing = FALSE), ]

        return(fcHurdleSig)
    }
)


#' @rdname differentialExpression
#' @examples
#' data(celdaCGSim, celdaCGMod)
#' clusterDiffexpRes <- differentialExpression(celdaCGSim$counts,
#'   celdaCGMod,
#'   c1 = c(1, 2))
#' @export
setMethod("differentialExpression",
    signature(x = "matrix"),
    function(x,
        celdaMod,
        c1,
        c2 = NULL,
        onlyPos = FALSE,
        log2fcThreshold = NULL,
        fdrThreshold = 1) {

        if (!(methods::is(celdaMod, "celda_C") ||
                methods::is(celdaMod, "celda_CG"))) {
            stop("'celdaMod' should be an object of class celda_C or celda_CG")
        }
        if (is.null(c1)) {
            stop("'c1' should be a numeric vector of cell cluster(s)")
        }

        cdiff <- setdiff(c1, celdaClusters(celdaMod)$z)

        if (length(cdiff) > 0) {
            warning("cluster ", cdiff, "in 'c1' does not exist in",
                " 'celdaClusters(celdaMod)$z'!")
            if (length(cdiff) == length(c1)) {
                stop("All clusters in 'c1' does not exist in",
                    " 'celdaClusters(celdaMod)$z'!")
            }
        }

        compareCountMatrix(x, celdaMod)

        if (is.null(c2)) {
            c2 <- sort(setdiff(unique(celdaClusters(celdaMod)$z), c1))
        }
        if (length(c1) > 1) {
            cells1 <- matrixNames(celdaMod)$column[which(
                celdaClusters(celdaMod)$z %in% c1
            )]
        } else {
            cells1 <- matrixNames(celdaMod)$column[which(
                celdaClusters(celdaMod)$z == c1
            )]
        }

        if (length(c2) > 1) {
            cells2 <- matrixNames(celdaMod)$column[which(
                celdaClusters(celdaMod)$z %in% c2
            )]
        } else {
            cells2 <- matrixNames(celdaMod)$column[which(
                celdaClusters(celdaMod)$z == c2
            )]
        }

        mat <- x[, c(cells1, cells2)]
        log_normalized_mat <- normalizeCounts(mat,
            normalize = "cpm",
            transformationFun = log1p
        )
        cdat <- data.frame(
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
            summaryDt[
                contrast == "conditionc1" & component == "logFC",
                .(primerid, coef, ci.hi, ci.lo)
                ],
            by = "primerid"
        )

        fcHurdle[, fdr := stats::p.adjust(`Pr(>Chisq)`, "fdr")]
        ### Some genes aren't outputted because log2FC gets NaN if one or both
        ### clusters have 0 counts for a gene
        ### and then they're discarded because NaN !> 0
        if (is.null(log2fcThreshold)) {
            fcHurdleSig <- fcHurdle
        } else {
            fcHurdleSig <- merge(fcHurdle[fdr < fdrThreshold &
                    abs(coef) > log2fcThreshold],
                data.table::as.data.table(S4Vectors::mcols(sca)),
                by = "primerid"
            )
            if (onlyPos) {
                fcHurdleSig <- fcHurdleSig[which(fcHurdleSig$log2fc > 0), ]
            }
        }
        fcHurdleSig <- fcHurdleSig[, -c(4, 5)]
        names(fcHurdleSig)[c(1, 2, 3, 4)] <- c("Gene", "Pvalue", "Log2_FC",
            "FDR")
        fcHurdleSig <- fcHurdleSig[order(fcHurdleSig$Pvalue,
            decreasing = FALSE), ]

        return(fcHurdleSig)
    }
)
