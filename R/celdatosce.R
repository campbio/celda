
#' @title Convert old celda model object to \code{SCE} object
#' @description Convert a old celda model object (\code{celda_C},
#'  \code{celda_G}, or \code{celda_CG} object) to a
#'  \linkS4class{SingleCellExperiment} object containing celda model
#'  information in \code{metadata} slot. Counts matrix is stored in the
#'  \code{"counts"} assay slot in \code{assays}.
#' @param celdaModel A celda model object generated using older versions of
#'  \code{celda}.
#' @param counts A numeric \link{matrix} of counts used to generate
#'  \code{celdaModel}. Dimensions and MD5 checksum will be checked by
#'  \link{compareCountMatrix}.
#' @return A \linkS4class{SingleCellExperiment} object. Function
#'  parameter settings are stored in the \link[S4Vectors]{metadata}
#'  \code{"celda_parameters"} slot.
#'  Columns \code{celda_sample_label} and \code{celda_cell_cluster} in
#'  \link[SummarizedExperiment]{colData} contain sample labels and celda cell
#'  population clusters. Column \code{celda_feature_module} in
#'  \link[SummarizedExperiment]{rowData} contain feature modules.

#' @export
setGeneric("celdatosce", function(celdaModel, counts) {
    standardGeneric("celdatosce")})


#' @rdname celdatosce
#' @examples
#' data(celdaCMod, celdaCSim)
#' sce <- celdatosce(celdaCMod, celdaCSim$counts)
#' @export
setMethod("celdatosce",
    signature(celdaModel = "celda_C"),
    function(celdaModel, counts) {
        compareCountMatrix(counts, celdaModel, errorOnMismatch = FALSE)

        xClass <- "matrix"
        useAssay <- NULL
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = counts))

        sce <- .createSCEceldaC(celdaCMod = celdaModel,
            sce = sce,
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
    function(celdaModel, counts) {
        compareCountMatrix(counts, celdaModel, errorOnMismatch = FALSE)

        xClass <- "matrix"
        useAssay <- NULL
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = counts))

        sce <- .createSCEceldaG(celdaGMod = celdaModel,
            sce = sce,
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
    function(celdaModel, counts) {
        compareCountMatrix(counts, celdaModel, errorOnMismatch = FALSE)

        xClass <- "matrix"
        useAssay <- NULL
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = counts))

        sce <- .createSCEceldaCG(celdaCGMod = celdaModel,
            sce = sce,
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
        return(sce)
    }
)
