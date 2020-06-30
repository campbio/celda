#' @title t-Distributed Stochastic Neighbor Embedding (t-SNE) dimension
#'  reduction for celda \code{sce} object
#' @description Embeds cells in two dimensions using \link[Rtsne]{Rtsne} based
#'  on a celda model. For celda_C \code{sce} objects, PCA on the normalized
#'  counts is used to reduce the number of features before applying t-SNE. For
#'  celda_CG and celda_G \code{sce} objects, tSNE is run on module
#'  probabilities to reduce the number of features instead of using PCA.
#'  Module probabilities are square-root transformed before applying tSNE.
#' @param sce A \linkS4class{SingleCellExperiment} object
#'  returned by \link{celda_C}, \link{celda_G}, or \link{celda_CG}.
#' @param useAssay A string specifying which \link[SummarizedExperiment]{assay}
#'  slot to use. Default "counts".
#' @param maxCells Integer. Maximum number of cells to plot. Cells will be
#'  randomly subsampled if \code{ncol(counts) > maxCells}. Larger numbers of
#'  cells requires more memory. If \code{NULL}, no subsampling will be
#'  performed. Default \code{NULL}.
#' @param minClusterSize Integer. Do not subsample cell clusters below this
#'  threshold. Default 100.
#' @param initialDims Integer. PCA will be used to reduce the dimensionality
#'  of the dataset. The top 'initialDims' principal components will be used
#'  for tSNE. Default 20.
#' @param modules Integer vector. Determines which feature modules to use for
#'  tSNE. If \code{NULL}, all modules will be used. Default \code{NULL}.
#' @param perplexity Numeric. Perplexity parameter for tSNE. Default 20.
#' @param maxIter Integer. Maximum number of iterations in tSNE generation.
#'  Default 2500.
#' @param normalize Character. Passed to \link{normalizeCounts} in
#'  normalization step. Divides counts by the library sizes for each
#'  cell. One of 'proportion', 'cpm', 'median', or 'mean'. 'proportion' uses
#'  the total counts for each cell as the library size. 'cpm' divides the
#'  library size of each cell by one million to produce counts per million.
#'  'median' divides the library size of each cell by the median library size
#'  across all cells. 'mean' divides the library size of each cell by the mean
#'  library size across all cells.
#' @param scaleFactor Numeric. Sets the scale factor for cell-level
#'  normalization. This scale factor is multiplied to each cell after the
#'  library size of each cell had been adjusted in \code{normalize}. Default
#'  \code{NULL} which means no scale factor is applied.
#' @param transformationFun Function. Applys a transformation such as 'sqrt',
#'  'log', 'log2', 'log10', or 'log1p'. If \code{NULL}, no transformation will
#'  be applied. Occurs after applying normalization and scale factor. Default
#'  \code{NULL}.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @return \code{sce} with t-SNE coordinates
#'  (columns "celda_tSNE1" & "celda_tSNE2") added to
#'  \code{\link[SingleCellExperiment]{reducedDim}(sce, "celda_tSNE")}.
#' @export
setGeneric("celdaTsne",
    function(sce, ...) {
        standardGeneric("celdaTsne")
    })


#' @rdname celdaTsne
#' @examples
#' data(sceCeldaCG)
#' tsneRes <- celdaTsne(sceCeldaCG)
#' @export
setMethod("celdaTsne", signature(sce = "SingleCellExperiment"),
    function(sce,
        useAssay = "counts",
        maxCells = NULL,
        minClusterSize = 100,
        initialDims = 20,
        modules = NULL,
        perplexity = 20,
        maxIter = 2500,
        normalize = "proportion",
        scaleFactor = NULL,
        transformationFun = sqrt,
        seed = 12345) {

        if (is.null(seed)) {
            res <- .celdaTsne(sce = sce,
                useAssay = useAssay,
                maxCells = maxCells,
                minClusterSize = minClusterSize,
                initialDims = initialDims,
                modules = modules,
                perplexity = perplexity,
                maxIter = maxIter,
                normalize = normalize,
                scaleFactor = scaleFactor,
                transformationFun = transformationFun)
        } else {
            with_seed(seed,
                res <- .celdaTsne(sce = sce,
                    useAssay = useAssay,
                    maxCells = maxCells,
                    minClusterSize = minClusterSize,
                    initialDims = initialDims,
                    modules = modules,
                    perplexity = perplexity,
                    maxIter = maxIter,
                    normalize = normalize,
                    scaleFactor = scaleFactor,
                    transformationFun = transformationFun))
        }

        SingleCellExperiment::reducedDim(sce, "celda_tSNE") <- res
        return(sce)
    })


.celdaTsne <- function(sce,
    useAssay,
    maxCells,
    minClusterSize,
    initialDims,
    modules,
    perplexity,
    maxIter,
    normalize,
    scaleFactor,
    transformationFun) {

    celdaMod <- celdaModel(sce)

    if (celdaMod == "celda_C") {
        res <- .celdaTsneC(sce = sce,
            useAssay = useAssay,
            maxCells = maxCells,
            minClusterSize = minClusterSize,
            initialDims = initialDims,
            perplexity = perplexity,
            maxIter = maxIter,
            normalize = normalize,
            scaleFactor = scaleFactor,
            transformationFun = transformationFun)
    } else if (celdaMod == "celda_CG") {
        res <- .celdaTsneCG(sce = sce,
            useAssay = useAssay,
            maxCells = maxCells,
            minClusterSize = minClusterSize,
            initialDims = initialDims,
            modules = modules,
            perplexity = perplexity,
            maxIter = maxIter,
            normalize = normalize,
            scaleFactor = scaleFactor,
            transformationFun = transformationFun)
    } else if (celdaMod == "celda_G") {
        res <- .celdaTsneG(sce = sce,
            useAssay = useAssay,
            maxCells = maxCells,
            minClusterSize = minClusterSize,
            initialDims = initialDims,
            modules = modules,
            perplexity = perplexity,
            maxIter = maxIter,
            normalize = normalize,
            scaleFactor = scaleFactor,
            transformationFun = transformationFun)
    } else {
        stop("S4Vectors::metadata(sce)$celda_parameters$model must be",
            " one of 'celda_C', 'celda_G', or 'celda_CG'")
    }
    return(res)
}


.celdaTsneC <- function(sce,
    useAssay,
    maxCells,
    minClusterSize,
    initialDims,
    perplexity,
    maxIter,
    normalize,
    scaleFactor,
    transformationFun) {

    preparedCountInfo <- .prepareCountsForDimReductionCeldaC(sce = sce,
        useAssay = useAssay,
        maxCells = maxCells,
        minClusterSize = minClusterSize,
        normalize = normalize,
        scaleFactor = scaleFactor,
        transformationFun = transformationFun)

    res <- .calculateTsne(preparedCountInfo$norm,
        perplexity = perplexity,
        maxIter = maxIter,
        doPca = TRUE,
        initialDims = initialDims)

    final <- matrix(NA, nrow = ncol(sce), ncol = 2)
    final[preparedCountInfo$cellIx, ] <- res
    rownames(final) <- colnames(sce)
    colnames(final) <- c("celda_tSNE1", "celda_tSNE2")
    return(final)
}


.celdaTsneCG <- function(sce,
    useAssay,
    maxCells,
    minClusterSize,
    initialDims,
    modules,
    perplexity,
    maxIter,
    normalize,
    scaleFactor,
    transformationFun) {

    preparedCountInfo <- .prepareCountsForDimReductionCeldaCG(sce = sce,
        useAssay = useAssay,
        maxCells = maxCells,
        minClusterSize = minClusterSize,
        modules = modules,
        normalize = normalize,
        scaleFactor = scaleFactor,
        transformationFun = transformationFun)
    norm <- preparedCountInfo$norm
    res <- .calculateTsne(norm,
        doPca = FALSE,
        perplexity = perplexity,
        maxIter = maxIter,
        initialDims = initialDims)
    final <- matrix(NA, nrow = ncol(sce), ncol = 2)
    final[preparedCountInfo$cellIx, ] <- res
    rownames(final) <- colnames(sce)
    colnames(final) <- c("celda_tSNE1", "celda_tSNE2")
    return(final)
}


.celdaTsneG <- function(sce,
    useAssay,
    maxCells,
    minClusterSize,
    initialDims,
    modules,
    perplexity,
    maxIter,
    normalize,
    scaleFactor,
    transformationFun) {

    preparedCountInfo <- .prepareCountsForDimReductionCeldaG(sce = sce,
        useAssay = useAssay,
        maxCells = maxCells,
        minClusterSize = minClusterSize,
        modules = modules,
        normalize = normalize,
        scaleFactor = scaleFactor,
        transformationFun = transformationFun)
    res <- .calculateTsne(preparedCountInfo$norm,
        perplexity = perplexity,
        maxIter = maxIter,
        doPca = FALSE,
        initialDims = initialDims)
    final <- matrix(NA, nrow = ncol(sce), ncol = 2)
    final[preparedCountInfo$cellIx, ] <- res
    rownames(final) <- colnames(sce)
    colnames(final) <- c("celda_tSNE1", "celda_tSNE2")
    return(final)
}


# Run the t-SNE algorithm for dimensionality reduction
# @param norm Normalized count matrix.
# @param perplexity Numeric vector. Determines perplexity for tsne. Default 20.
# @param maxIter Numeric vector. Determines iterations for tsne. Default 1000.
# @param doPca Logical. Whether to perform
# dimensionality reduction with PCA before tSNE.
# @param initialDims Integer. Number of dimensions from PCA to use as
# input in tSNE. Default 50.
#' @importFrom Rtsne Rtsne
.calculateTsne <- function(norm,
    perplexity,
    maxIter,
    doPca,
    initialDims) {

    res <- Rtsne::Rtsne(
        norm,
        pca = doPca,
        max_iter = maxIter,
        perplexity = perplexity,
        check_duplicates = FALSE,
        is_distance = FALSE,
        initial_dims = initialDims)$Y

    return(res)
}
