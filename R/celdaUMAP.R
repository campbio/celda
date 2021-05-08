#' @title Uniform Manifold Approximation and Projection (UMAP) dimension
#'  reduction for celda \code{sce} object
#' @description Embeds cells in two dimensions using \link[uwot]{umap} based on
#'  a celda model. For celda_C \code{sce} objects, PCA on the normalized counts
#'  is used to reduce the number of features before applying UMAP. For celda_CG
#'  \code{sce} object, UMAP is run on module probabilities to reduce the number
#'  of features instead of using PCA. Module probabilities are square-root
#'  transformed before applying UMAP.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object
#'  returned by \link{celda_C}, \link{celda_G}, or \link{celda_CG}.
#' @param useAssay A string specifying which \link{assay}
#'  slot to use. Default "counts".
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default "featureSubset".
#' @param maxCells Integer. Maximum number of cells to plot. Cells will be
#'  randomly subsampled if \code{ncol(sce) > maxCells}. Larger numbers of cells
#'  requires more memory. If NULL, no subsampling will be performed.
#'  Default NULL.
#' @param minClusterSize Integer. Do not subsample cell clusters below this
#'  threshold. Default 100.
#' @param modules Integer vector. Determines which features modules to use for
#'  UMAP. If NULL, all modules will be used. Default NULL.
#' @param seed Integer. Passed to \link[withr]{with_seed}. For reproducibility,
#'  a default value of 12345 is used. If NULL, no calls to
#'  \link[withr]{with_seed} are made.
#' @param nNeighbors The size of local neighborhood used for
#'  manifold approximation. Larger values result in more global
#'  views of the manifold, while smaller values result in more
#'  local data being preserved. Default 30.
#'  See \link[uwot]{umap} for more information.
#' @param minDist The effective minimum distance between embedded points.
#'  Smaller values will result in a more clustered/clumped
#'  embedding where nearby points on the manifold are drawn
#'  closer together, while larger values will result on a more
#'  even dispersal of points. Default 0.75.
#'  See \link[uwot]{umap} for more information.
#' @param spread The effective scale of embedded points. In combination with
#'  \code{min_dist}, this determines how clustered/clumped the
#'   embedded points are. Default 1. See \link[uwot]{umap} for more information.
#' @param pca Logical. Whether to perform
#'  dimensionality reduction with PCA before UMAP. Only works for celda_C
#'  \code{sce} objects.
#' @param initialDims Integer. Number of dimensions from PCA to use as
#'  input in UMAP. Default 50. Only works for celda_C \code{sce} objects.
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
#' @param cores Number of threads to use. Default 1.
#' @param ... Additional parameters to pass to \link[uwot]{umap}.
#' @return \code{sce} with UMAP coordinates
#'  (columns "celda_UMAP1" & "celda_UMAP2") added to
#'  \code{\link{reducedDim}(sce, "celda_UMAP")}.
#' @export
setGeneric("celdaUmap",
    function(sce,
        useAssay = "counts",
        altExpName = "featureSubset",
        maxCells = NULL,
        minClusterSize = 100,
        modules = NULL,
        seed = 12345,
        nNeighbors = 30,
        minDist = 0.75,
        spread = 1,
        pca = TRUE,
        initialDims = 50,
        normalize = "proportion",
        scaleFactor = NULL,
        transformationFun = sqrt,
        cores = 1,
        ...) {

        standardGeneric("celdaUmap")
    })


#' @rdname celdaUmap
#' @examples
#' data(sceCeldaCG)
#' umapRes <- celdaUmap(sceCeldaCG)
#' @export
setMethod("celdaUmap", signature(sce = "SingleCellExperiment"),
    function(sce,
        useAssay = "counts",
        altExpName = "featureSubset",
        maxCells = NULL,
        minClusterSize = 100,
        modules = NULL,
        seed = 12345,
        nNeighbors = 30,
        minDist = 0.75,
        spread = 1,
        pca = TRUE,
        initialDims = 50,
        normalize = "proportion",
        scaleFactor = NULL,
        transformationFun = sqrt,
        cores = 1,
        ...) {

        if (is.null(seed)) {
            sce <- .celdaUmap(sce = sce,
                useAssay = useAssay,
                altExpName = altExpName,
                maxCells = maxCells,
                minClusterSize = minClusterSize,
                modules = modules,
                seed = seed,
                nNeighbors = nNeighbors,
                minDist = minDist,
                spread = spread,
                pca = pca,
                initialDims = initialDims,
                normalize = normalize,
                scaleFactor = scaleFactor,
                transformationFun = transformationFun,
                cores = cores,
                ...)
        } else {
            with_seed(seed,
                sce <- .celdaUmap(sce = sce,
                    useAssay = useAssay,
                    altExpName = altExpName,
                    maxCells = maxCells,
                    minClusterSize = minClusterSize,
                    modules = modules,
                    seed = seed,
                    nNeighbors = nNeighbors,
                    minDist = minDist,
                    spread = spread,
                    pca = pca,
                    initialDims = initialDims,
                    normalize = normalize,
                    scaleFactor = scaleFactor,
                    transformationFun = transformationFun,
                    cores = cores,
                    ...))
        }
        return(sce)
    })


.celdaUmap <- function(sce,
    useAssay,
    altExpName,
    maxCells,
    minClusterSize,
    modules,
    seed,
    nNeighbors,
    minDist,
    spread,
    pca,
    initialDims,
    cores,
    normalize,
    scaleFactor,
    transformationFun,
    ...) {

    celdaMod <- celdaModel(sce, altExpName = altExpName)
    altExp <- SingleCellExperiment::altExp(sce, altExpName)

    if (celdaMod == "celda_C") {
        res <- .celdaUmapC(sce = altExp,
            useAssay = useAssay,
            maxCells = maxCells,
            minClusterSize = minClusterSize,
            nNeighbors = nNeighbors,
            minDist = minDist,
            spread = spread,
            pca = pca,
            initialDims = initialDims,
            normalize = normalize,
            scaleFactor = scaleFactor,
            transformationFun = transformationFun,
            cores = cores,
            ...)
    } else if (celdaMod == "celda_CG") {
        res <- .celdaUmapCG(sce = altExp,
            useAssay = useAssay,
            maxCells = maxCells,
            minClusterSize = minClusterSize,
            modules = modules,
            seed = seed,
            nNeighbors = nNeighbors,
            minDist = minDist,
            spread = spread,
            normalize = normalize,
            scaleFactor = scaleFactor,
            transformationFun = transformationFun,
            cores = cores,
            ...)
    } else if (celdaMod == "celda_G") {
        res <- .celdaUmapG(sce = altExp,
            useAssay = useAssay,
            maxCells = maxCells,
            minClusterSize = minClusterSize,
            modules = modules,
            seed = seed,
            nNeighbors = nNeighbors,
            minDist = minDist,
            spread = spread,
            normalize = normalize,
            scaleFactor = scaleFactor,
            transformationFun = transformationFun,
            cores = cores,
            ...)
    } else {
        stop("S4Vectors::metadata(altExp(sce, altExpName))$",
            "celda_parameters$model must be",
            " one of 'celda_C', 'celda_G', or 'celda_CG'")
    }
    SingleCellExperiment::reducedDim(altExp, "celda_UMAP") <- res
    SingleCellExperiment::altExp(sce, altExpName) <- altExp
    return(sce)
}


.celdaUmapC <- function(sce,
    useAssay,
    maxCells,
    minClusterSize,
    nNeighbors,
    minDist,
    spread,
    pca,
    initialDims,
    normalize,
    scaleFactor,
    transformationFun,
    cores,
    ...) {

    preparedCountInfo <- .prepareCountsForDimReductionCeldaC(sce = sce,
        useAssay = useAssay,
        maxCells = maxCells,
        minClusterSize = minClusterSize,
        normalize = normalize,
        scaleFactor = scaleFactor,
        transformationFun = transformationFun)
    umapRes <- .calculateUmap(preparedCountInfo$norm,
        nNeighbors = nNeighbors,
        minDist = minDist,
        spread = spread,
        pca = pca,
        initialDims = initialDims,
        cores = cores,
        ...
    )

    final <- matrix(NA, nrow = ncol(sce), ncol = 2)
    final[preparedCountInfo$cellIx, ] <- umapRes
    rownames(final) <- colnames(sce)
    colnames(final) <- c("celda_UMAP1", "celda_UMAP2")
    return(final)
}


.celdaUmapCG <- function(sce,
    useAssay,
    maxCells,
    minClusterSize,
    modules,
    seed,
    nNeighbors,
    minDist,
    spread,
    normalize,
    scaleFactor,
    transformationFun,
    cores,
    ...) {

    preparedCountInfo <- .prepareCountsForDimReductionCeldaCG(sce = sce,
        useAssay = useAssay,
        maxCells = maxCells,
        minClusterSize = minClusterSize,
        modules = modules,
        normalize = normalize,
        scaleFactor = scaleFactor,
        transformationFun = transformationFun)
    umapRes <- .calculateUmap(preparedCountInfo$norm,
        nNeighbors = nNeighbors,
        minDist = minDist,
        spread = spread,
        cores = cores,
        ...)

    final <- matrix(NA, nrow = ncol(sce), ncol = 2)
    final[preparedCountInfo$cellIx, ] <- umapRes
    rownames(final) <- colnames(sce)
    colnames(final) <- c("celda_UMAP1", "celda_UMAP2")
    return(final)
}


.celdaUmapG <- function(sce,
    useAssay,
    maxCells,
    minClusterSize,
    modules,
    seed,
    nNeighbors,
    minDist,
    spread,
    normalize,
    scaleFactor,
    transformationFun,
    cores,
    ...) {

    preparedCountInfo <- .prepareCountsForDimReductionCeldaG(sce = sce,
        useAssay = useAssay,
        maxCells = maxCells,
        minClusterSize = minClusterSize,
        modules = modules,
        normalize = normalize,
        scaleFactor = scaleFactor,
        transformationFun = transformationFun)
    umapRes <- .calculateUmap(preparedCountInfo$norm,
        nNeighbors = nNeighbors,
        minDist = minDist,
        spread = spread,
        cores = cores,
        ...)

    final <- matrix(NA, nrow = ncol(sce), ncol = 2)
    final[preparedCountInfo$cellIx, ] <- umapRes
    rownames(final) <- colnames(sce)
    colnames(final) <- c("celda_UMAP1", "celda_UMAP2")
    return(final)
}


# Run the UMAP algorithm for dimensionality reduction
# @param norm Normalized count matrix.
# @param nNeighbors The size of local neighborhood used for
#   manifold approximation. Larger values result in more global
#   views of the manifold, while smaller values result in more
#   local data being preserved. Default 30.
#    See `?uwot::umap` for more information.
# @param minDist The effective minimum distance between embedded points.
#    Smaller values will result in a more clustered/clumped
#    embedding where nearby points on the manifold are drawn
#    closer together, while larger values will result on a more
#    even dispersal of points. Default 0.2.
#    See `?uwot::umap` for more information.
# @param spread The effective scale of embedded points. In combination with
#    'min_dist', this determines how clustered/clumped the
#    embedded points are. Default 1.
#    See `?uwot::umap` for more information.
# @param pca Logical. Whether to perform
# dimensionality reduction with PCA before UMAP.
# @param initialDims Integer. Number of dimensions from PCA to use as
# input in UMAP. Default 50.
# @param cores Number of threads to use. Default 1.
# @param ... Other parameters to pass to `uwot::umap`.
#' @import uwot
.calculateUmap <- function(norm,
    nNeighbors = 30,
    minDist = 0.75,
    spread = 1,
    pca = FALSE,
    initialDims = 50,
    cores = 1,
    ...) {
    if (isTRUE(pca)) {
        doPCA <- initialDims
    } else {
        doPCA <- NULL
    }

    res <- uwot::umap(norm,
        n_neighbors = nNeighbors,
        min_dist = minDist, spread = spread,
        n_threads = cores, n_sgd_threads = 1, pca = doPCA, ...
    )
    return(res)
}
