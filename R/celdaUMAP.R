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
#' @param useAssay A string specifying which \link[SummarizedExperiment]{assay}
#'  slot to use. Default "counts".
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
#'   manifold approximation. Larger values result in more global
#'   views of the manifold, while smaller values result in more
#'   local data being preserved. Default 30.
#'   See \link[uwot]{umap} for more information.
#' @param minDist The effective minimum distance between embedded points.
#'   Smaller values will result in a more clustered/clumped
#'   embedding where nearby points on the manifold are drawn
#'   closer together, while larger values will result on a more
#'   even dispersal of points. Default 0.75.
#'   See \link[uwot]{umap} for more information.
#' @param spread The effective scale of embedded points. In combination with
#'  \code{min_dist}, this determines how clustered/clumped the
#'   embedded points are. Default 1. See \link[uwot]{umap} for more information.
#' @param pca Logical. Whether to perform
#' dimensionality reduction with PCA before UMAP. Only works for celda_C
#'  \code{sce} objects.
#' @param initialDims Integer. Number of dimensions from PCA to use as
#' input in UMAP. Default 50. Only works for celda_C \code{sce} objects.
#' @param cores Number of threads to use. Default 1.
#' @param ... Additional parameters to pass to \link[uwot]{umap}.
#' @return \code{sce} with UMAP coordinates
#'  (columns "celda_UMAP1" & "celda_UMAP2") added to
#'  \code{\link[SingleCellExperiment]{reducedDim}(sce, "celda_UMAP")}.
#' @export
setGeneric("celdaUmap",
    function(sce, ...) {
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
        maxCells = NULL,
        minClusterSize = 100,
        modules = NULL,
        seed = 12345,
        nNeighbors = 30,
        minDist = 0.75,
        spread = 1,
        pca = TRUE,
        initialDims = 50,
        cores = 1,
        ...) {

        if (is.null(seed)) {
            res <- .celdaUmap(sce = sce,
                useAssay = useAssay,
                maxCells = maxCells,
                minClusterSize = minClusterSize,
                modules = modules,
                seed = seed,
                nNeighbors = nNeighbors,
                minDist = minDist,
                spread = spread,
                pca = pca,
                initialDims = initialDims,
                cores = cores,
                ...)
        } else {
            with_seed(seed,
                res <- .celdaUmap(sce = sce,
                    useAssay = useAssay,
                    maxCells = maxCells,
                    minClusterSize = minClusterSize,
                    modules = modules,
                    seed = seed,
                    nNeighbors = nNeighbors,
                    minDist = minDist,
                    spread = spread,
                    pca = pca,
                    initialDims = initialDims,
                    cores = cores,
                    ...))
        }

        SingleCellExperiment::reducedDim(sce, "celda_UMAP") <- res
        return(sce)
    })


.celdaUmap <- function(sce,
    useAssay,
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
    ...) {

    celdaMod <- celdaModel(sce)

    if (celdaMod == "celda_C") {
        res <- .celdaUmapC(sce = sce,
            useAssay = useAssay,
            maxCells = maxCells,
            minClusterSize = minClusterSize,
            nNeighbors = nNeighbors,
            minDist = minDist,
            spread = spread,
            pca = pca,
            initialDims = initialDims,
            cores = cores,
            ...)
    } else if (celdaMod == "celda_CG") {
        res <- .celdaUmapCG(sce = sce,
            useAssay = useAssay,
            maxCells = maxCells,
            minClusterSize = minClusterSize,
            modules = modules,
            seed = seed,
            nNeighbors = nNeighbors,
            minDist = minDist,
            spread = spread,
            cores = cores,
            ...)
    } else if (celdaMod == "celda_G") {
        res <- .celdaUmapG(sce = sce,
            useAssay = useAssay,
            maxCells = maxCells,
            minClusterSize = minClusterSize,
            modules = modules,
            seed = seed,
            nNeighbors = nNeighbors,
            minDist = minDist,
            spread = spread,
            cores = cores,
            ...)
    } else {
        stop("S4Vectors::metadata(sce)$celda_parameters$model must be",
            " one of 'celda_C', 'celda_G', or 'celda_CG'")
    }
    return(res)

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
    cores,
    ...) {

    preparedCountInfo <- .prepareCountsForDimReductionCeldaC(sce,
        useAssay,
        maxCells,
        minClusterSize)
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
    cores,
    ...) {

    preparedCountInfo <- .prepareCountsForDimReductionCeldaCG(sce,
        useAssay,
        maxCells,
        minClusterSize,
        modules)
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


.celdaUmapG <- function(sce = sce,
    useAssay = useAssay,
    maxCells = maxCells,
    minClusterSize = minClusterSize,
    modules = modules,
    seed = seed,
    nNeighbors = nNeighbors,
    minDist = minDist,
    spread = spread,
    cores = cores,
    ...) {

    preparedCountInfo <- .prepareCountsForDimReductionCeldaG(sce = sce,
        useAssay = useAssay,
        maxCells = maxCells,
        minClusterSize = minClusterSize,
        modules = modules)
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
