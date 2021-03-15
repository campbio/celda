#' @title Generate a report for celda_CG results
#' @description After a celda_CG model has been fitted, this function will
#' create an html report which can be used to visualize and explore the results. 
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object
#'   returned by \link{celda_CG}.
#' @param reducedDimName Character. Name of the reduced dimensional object to be used in 2-D scatter plots throughout the report.
#' @param features Character vector.  Expression of these features will be
#'   displayed on a reduced dimensional plot defined by \code{reducedDimName}. If \code{NULL}, then no plotting of features
#'   on a reduced dimensinoal plot will be performed. Default \code{NULL}.
#' @param displayName Character. The name to use for display in scatter plots
#'   and heatmaps. This can be \code{"rownames"} or the name of a column in the row
#'   data. Default \code{"rownames"}.
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default "featureSubset".
#' @param useAssay A string specifying which \link{assay}
#'  slot to use. Default \code{"counts"}.
#' @param cellAnnot Character vector. The cell-level annotations to display on the reduced dimensional plot. These variables should be present in the column data of the \code{sce} object.
#' @param exactMatch Boolean. Whether to only identify exact matches
#' or to identify partial matches using \code{\link{grep}}. Default \code{FALSE}.
#' @param output_file Character. Prefix of the html file. Default \code{"CeldaCG_ResultReport"}.
#' @param output_dir Character. Path to save the html file. Default \code{.}.
#' @param pdf Boolean. Whether to create PDF versions of each plot in addition to PNGs. Default \code{FALSE}.
#' @param showSetup Boolean. Whether to show the setup code at the beginning. Default \code{TRUE}.
#' @param showSession Boolean. Whether to show the session information at the end. Default \code{TRUE}.
#' @return .html file
#' @examples
#' data(sceCeldaCG)
#' \dontrun{
#' reportCeldaCG_PlotResults(sce = sceCeldaCG, features = c("Gene_1", "Gene_100"))
#' }
#' @export
reportCeldaCG_PlotResults <-
  function(sce,
           reducedDimName,
           features = NULL,
           displayName = "rownames",
           altExpName = "featureSubset",
           useAssay = "counts",
           cellAnnot = NULL,
           exactMatch = TRUE,
           output_file = "CeldaCG_ResultReport",
           output_dir = ".",
           pdf = FALSE,
           showSetup = TRUE,
           showSession = TRUE) {
    rmarkdown::render(
      system.file("rmarkdown/CeldaCG_PlotResults.rmd", package = "celda"),
      params = list(
        sce = sce,
        altExpName = altExpName,
        useAssay = useAssay,
        reducedDimName = reducedDimName,
        features = features,
        displayName = displayName,
        cellAnnot = cellAnnot,
        exactMatch = isTRUE(exactMatch),
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
