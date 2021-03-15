
#' @title Convert old celda model object to \code{SCE} object
#' @description Convert a old celda model object (\code{celda_C},
#'  \code{celda_G}, or \code{celda_CG} object) to a
#'  \linkS4class{SingleCellExperiment} object containing celda model
#'  information in \code{metadata} slot. Counts matrix is stored in the
#'  \code{"counts"} assay slot in \code{assays}.
#' @param celdaModel A \code{celdaModel} or \code{celdaList} object generated
#'  using older versions of \code{celda}.
#' @param counts A numeric \link{matrix} of counts used to generate
#'  \code{celdaModel}. Dimensions and MD5 checksum will be checked by
#'  \link{compareCountMatrix}.
#' @param useAssay A string specifying the name of the
#'  \link{assay} slot to use. Default "counts".
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default "featureSubset".
#' @param ... Ignored. Placeholder to prevent check warning.
#' @return A \linkS4class{SingleCellExperiment} object. Function
#'  parameter settings are stored in the \link{metadata}
#'  \code{"celda_parameters"} slot.
#'  Columns \code{celda_sample_label} and \code{celda_cell_cluster} in
#'  \link{colData} contain sample labels and celda cell
#'  population clusters. Column \code{celda_feature_module} in
#'  \link{rowData} contain feature modules.
#' @export
setGeneric("celdatosce", function(celdaModel, counts, ...) {
    standardGeneric("celdatosce")})


#' @rdname celdatosce
#' @examples
#' data(celdaCMod, celdaCSim)
#' sce <- celdatosce(celdaCMod, celdaCSim$counts)
#' @export
setMethod("celdatosce",
    signature(celdaModel = "celda_C"),
    function(celdaModel,
        counts,
        useAssay = "counts",
        altExpName = "featureSubset") {

        compareCountMatrix(counts, celdaModel, errorOnMismatch = FALSE)

        ls <- list()
        ls[[useAssay]] <- counts
        sce <- SingleCellExperiment::SingleCellExperiment(assays = ls)
        SingleCellExperiment::altExp(sce, altExpName) <- sce
        xClass <- "matrix"

        altExp <- .createSCEceldaC(celdaCMod = celdaModel,
            sce = SingleCellExperiment::altExp(sce, altExpName),
            xClass = xClass,
            useAssay = useAssay,
            algorithm = NULL,
            stopIter = NULL,
            maxIter = NULL,
            splitOnIter = NULL,
            splitOnLast = NULL,
            nchains = NULL,
            zInitialize = NULL,
            zInit = NULL,
            logfile = NULL,
            verbose = NULL)
        SingleCellExperiment::altExp(sce, altExpName) <- altExp
        return(sce)
    }
)


#' @rdname celdatosce
#' @examples
#' data(celdaGMod, celdaGSim)
#' sce <- celdatosce(celdaGMod, celdaGSim$counts)
#' @export
setMethod("celdatosce",
    signature(celdaModel = "celda_G"),
    function(celdaModel,
        counts,
        useAssay = "counts",
        altExpName = "featureSubset") {

        compareCountMatrix(counts, celdaModel, errorOnMismatch = FALSE)

        ls <- list()
        ls[[useAssay]] <- counts
        sce <- SingleCellExperiment::SingleCellExperiment(assays = ls)
        SingleCellExperiment::altExp(sce, altExpName) <- sce
        xClass <- "matrix"

        altExp <- .createSCEceldaG(celdaGMod = celdaModel,
            sce = SingleCellExperiment::altExp(sce, altExpName),
            xClass = xClass,
            useAssay = useAssay,
            stopIter = NULL,
            maxIter = NULL,
            splitOnIter = NULL,
            splitOnLast = NULL,
            nchains = NULL,
            yInitialize = NULL,
            yInit = NULL,
            logfile = NULL,
            verbose = NULL)
        SingleCellExperiment::altExp(sce, altExpName) <- altExp
        return(sce)
    }
)


#' @rdname celdatosce
#' @examples
#' data(celdaCGMod, celdaCGSim)
#' sce <- celdatosce(celdaCGMod, celdaCGSim$counts)
#' @export
setMethod("celdatosce",
    signature(celdaModel = "celda_CG"),
    function(celdaModel,
        counts,
        useAssay = "counts",
        altExpName = "featureSubset") {

        compareCountMatrix(counts, celdaModel, errorOnMismatch = FALSE)

        ls <- list()
        ls[[useAssay]] <- counts
        sce <- SingleCellExperiment::SingleCellExperiment(assays = ls)
        SingleCellExperiment::altExp(sce, altExpName) <- sce
        xClass <- "matrix"

        altExp <- .createSCEceldaCG(celdaCGMod = celdaModel,
            sce = SingleCellExperiment::altExp(sce, altExpName),
            xClass = xClass,
            useAssay = useAssay,
            algorithm = NULL,
            stopIter = NULL,
            maxIter = NULL,
            splitOnIter = NULL,
            splitOnLast = NULL,
            nchains = NULL,
            zInitialize = NULL,
            yInitialize = NULL,
            zInit = NULL,
            yInit = NULL,
            logfile = NULL,
            verbose = NULL)
        SingleCellExperiment::altExp(sce, altExpName) <- altExp
        return(sce)
    }
)


#' @rdname celdatosce
#' @examples
#' data(celdaCGGridSearchRes, celdaCGSim)
#' sce <- celdatosce(celdaCGGridSearchRes, celdaCGSim$counts)
#' @export
setMethod("celdatosce",
    signature(celdaModel = "celdaList"),
    function(celdaModel,
        counts,
        useAssay = "counts",
        altExpName = "featureSubset") {

        compareCountMatrix(counts, celdaModel, errorOnMismatch = FALSE)

        ls <- list()
        ls[[useAssay]] <- counts
        sce <- SingleCellExperiment::SingleCellExperiment(assays = ls)
        SingleCellExperiment::altExp(sce, altExpName) <- sce
        xClass <- "matrix"
        model <- celdaModel@celdaGridSearchParameters$model
        paramsTest <- celdaModel@celdaGridSearchParameters$paramsTest
        paramsFixed <-
            celdaModel@celdaGridSearchParameters$paramsFixed
        maxIter <- celdaModel@celdaGridSearchParameters$maxIter
        nchains <- celdaModel@celdaGridSearchParameters$nchains
        cores <- celdaModel@celdaGridSearchParameters$cores
        bestOnly <- celdaModel@celdaGridSearchParameters$bestOnly
        perplexity <- celdaModel@celdaGridSearchParameters$perplexity
        verbose <- celdaModel@celdaGridSearchParameters$verbose
        logfilePrefix <-
            celdaModel@celdaGridSearchParameters$logfilePrefix

        altExp <- .createSCEceldaGridSearch(celdaList = celdaModel,
            sce = SingleCellExperiment::altExp(sce, altExpName),
            xClass = xClass,
            useAssay = useAssay,
            model = model,
            paramsTest = paramsTest,
            paramsFixed = paramsFixed,
            maxIter = maxIter,
            seed = NULL,
            nchains = nchains,
            cores = cores,
            bestOnly = bestOnly,
            perplexity = perplexity,
            verbose = verbose,
            logfilePrefix = logfilePrefix)
        SingleCellExperiment::altExp(sce, altExpName) <- altExp
        return(sce)
    }
)
