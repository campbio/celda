#' @title Generate an HTML report for celda_CG
#' @name reportceldaCG
#' @description \code{reportCeldaCGRun} will run \link{recursiveSplitModule} and
#'   \link{recursiveSplitCell} to find the number of modules (\code{L}) and the
#'   number of cell populations (\code{K}). A final \link{celda_CG} model will
#'   be selected from \link{recursiveSplitCell}. After a \link{celda_CG} model
#'   has been fit, \code{reportCeldaCGPlotResults} can be used to create an HTML
#'   report for visualization and exploration of the \link{celda_CG} model
#'   results. Some of the plotting and feature selection functions require the
#'   installation of the Bioconductor package \code{singleCellTK}.
#' @param sce A \linkS4class{SingleCellExperiment} with the matrix located in
#'   the assay slot under \code{useAssay}. Rows represent features and columns
#'   represent cells.
#' @param L Integer. Final number of feature modules. See \code{celda_CG} for
#'   more information.
#' @param K Integer. Final number of cell populations. See \code{celda_CG} for
#'   more information.
#' @param sampleLabel Vector or factor. Denotes the sample label for each cell
#'   (column) in the count matrix.
#' @param altExpName The name for the \link{altExp} slot to use. Default
#'   \code{"featureSubset"}.
#' @param useAssay A string specifying which \link{assay} slot to use. Default
#'   \code{"counts"}.
#' @param initialL Integer. Minimum number of modules to try. See
#'   \link{recursiveSplitModule} for more information. Defailt \code{10}.
#' @param maxL Integer. Maximum number of modules to try. See
#'   \link{recursiveSplitModule} for more information. Default \code{150}.
#' @param initialK Integer. Initial number of cell populations to try.
#' @param maxK Integer. Maximum number of cell populations to try.
#' @param minCell Integer. Minimum number of cells required for feature
#'   selection. See \link{selectFeatures} for more information. Default
#'   \code{3}.
#' @param minCount Integer. Minimum number of counts required for feature
#'   selection. See \link{selectFeatures} for more information. Default
#'   \code{3}.
#' @param maxFeatures Integer. Maximum number of features to include. If the
#'   number of features after filtering for \code{minCell} and \code{minCount}
#'   are greater than \code{maxFeature}, then Seurat's VST function is used to
#'   select the top variable features. Default \code{5000}.
#' @param reducedDimName Character. Name of the reduced dimensional object to be
#'   used in 2-D scatter plots throughout the report. Default \code{celda_UMAP}.
#' @param features Character vector.  Expression of these features will be
#'   displayed on a reduced dimensional plot defined by \code{reducedDimName}.
#'   If \code{NULL}, then no plotting of features on a reduced dimensinoal plot
#'   will be performed. Default \code{NULL}.
#' @param displayName Character. The name to use for display in scatter plots
#'   and heatmaps. If \code{NULL}, then the rownames of the \code{sce} object
#'   will be used. This can also be set to the name of a column in the row data
#'   of the \code{sce}. Default \code{NULL}. data. Default \code{"rownames"}.
#' @param cellAnnot Character vector. The cell-level annotations to display on
#'   the reduced dimensional plot. These variables should be present in the
#'   column data of the \code{sce} object. Default \code{NULL}.
#' @param cellAnnotLabel Character vector. Additional cell-level annotations
#'   to display on the reduced dimensional plot. Variables will be treated
#'   as categorial and labels for each group will be placed on the plot.
#'   These variables should be present in the column data of the \code{sce}
#'   object. Default \code{NULL}.
#' @param exactMatch Boolean. Whether to only identify exact matches or to
#'   identify partial matches using \code{\link{grep}}. Default \code{FALSE}.
#' @param moduleFilePrefix Character. The features in each module will be
#' written to a a csv file starting with this name. If \code{NULL}, then no
#' file will be written. Default \code{"module_features"}.
#' @param output_file Character. Prefix of the html file. Default
#'   \code{"CeldaCG_ResultReport"}.
#' @param output_sce_prefix Character. The \code{sce} object with
#'   \code{celda_CG} results will be saved to an \code{.rds} file starting with
#'   this prefix. Default \code{celda_cg}.
#' @param output_dir Character. Path to save the html file. Default \code{.}.
#' @param pdf Boolean. Whether to create PDF versions of each plot in addition
#'   to PNGs. Default \code{FALSE}.
#' @param showSetup Boolean. Whether to show the setup code at the beginning.
#'   Default \code{TRUE}.
#' @param showSession Boolean. Whether to show the session information at the
#'   end. Default \code{TRUE}.
#' @return .html file
#' @rdname reportceldaCG
#' @examples
#' data(sceCeldaCG)
#' \donttest{
#' library(SingleCellExperiment)
#' sceCeldaCG$sum <- colSums(counts(sceCeldaCG))
#' sceCeldaCG <- reportCeldaCGRun(sceCeldaCG,
#'               initialL=5, maxL=20, initialK=5, maxK=20, L=10, K=5)
#' reportCeldaCGPlotResults(sce = sceCeldaCG, reducedDimName = "celda_UMAP",
#'                          features = c("Gene_1", "Gene_100"),
#'                          cellAnnot="sum")
#' }
NULL

#' @rdname reportceldaCG
#' @export
reportCeldaCGRun <-
  function(sce,
           L,
           K,
           sampleLabel = NULL,
           altExpName = "featureSubset",
           useAssay = "counts",
           initialL = 10,
           maxL = 150,
           initialK = 5,
           maxK = 50,
           minCell = 3,
           minCount = 3,
           maxFeatures = 5000,
           output_file = "CeldaCG_RunReport",
           output_sce_prefix = "celda_cg",
           output_dir = ".",
           pdf = FALSE,
           showSession = TRUE) {
    sceFile <-
      file.path(normalizePath(output_dir),
                paste0(output_sce_prefix, ".rds"))

    rmarkdown::render(
      system.file("rmarkdown/CeldaCG_Run.Rmd", package = "celda"),
      params = list(
        sce = sce,
        L = L,
        K = K,
        sampleLabel = sampleLabel,
        altExpName = altExpName,
        useAssay = useAssay,
        initialL = initialL,
        maxL = maxL,
        initialK = initialK,
        maxK = maxK,
        minCell = minCell,
        minCount = minCount,
        maxFeatures = maxFeatures,
        sceFile = sceFile,
        pdf = isTRUE(pdf),
        showSession = isTRUE(showSession)
      ),
      output_file = output_file,
      output_dir = output_dir,
      intermediates_dir = output_dir,
      knit_root_dir = output_dir
    )

    if (!is.null(output_sce_prefix)) {
      if (file.exists(sceFile)) {
        sce <- readRDS(sceFile)
        invisible(sce)
      } else {
        warning(
          "The file '",
          sceFile,
        "' could not be found. The SCE with celda_CG results was not reloaded."
        )
      }
    }
  }




#' @rdname reportceldaCG
#' @export
reportCeldaCGPlotResults <-
  function(sce,
           reducedDimName,
           features = NULL,
           displayName = NULL,
           altExpName = "featureSubset",
           useAssay = "counts",
           cellAnnot = NULL,
           cellAnnotLabel = NULL,
           exactMatch = TRUE,
           moduleFilePrefix = "module_features",
           output_file = "CeldaCG_ResultReport",
           output_dir = ".",
           pdf = FALSE,
           showSetup = TRUE,
           showSession = TRUE) {

    moduleFileName <- NULL
    if (!is.null(moduleFilePrefix)) {
      moduleFileName <-
        file.path(normalizePath(output_dir),
                  paste0(moduleFilePrefix, ".csv"))
    }

    rmarkdown::render(
      system.file("rmarkdown/CeldaCG_PlotResults.Rmd", package = "celda"),
      params = list(
        sce = sce,
        altExpName = altExpName,
        useAssay = useAssay,
        reducedDimName = reducedDimName,
        features = features,
        displayName = displayName,
        cellAnnot = cellAnnot,
        cellAnnotLabel = cellAnnotLabel,
        exactMatch = isTRUE(exactMatch),
        moduleFileName = moduleFileName,
        pdf = isTRUE(pdf),
        showSetup = isTRUE(showSetup),
        showSession = isTRUE(showSession)
      ),
      output_file = output_file,
      output_dir = output_dir,
      intermediates_dir = output_dir,
      knit_root_dir = output_dir
    )
  }
