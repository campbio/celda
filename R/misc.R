#' @title Celda models
#' @description List of available Celda models with correpsonding descriptions.
#' @export
#' @examples
#' celda()
#' @return None
celda <- function() {
    message(
        "celda_C: Clusters the columns of a count matrix containing",
        " single-cell data into K subpopulations."
    )
    message(
        "celda_G: Clusters the rows of a count matrix containing",
        " single-cell data into L modules."
    )
    message(
        "celda_CG: Clusters the rows and columns of a count matrix",
        " containing single-cell data into L modules and K subpopulations,",
        " respectively."
    )
    message(
        "celdaGridSearch: Run Celda with different combinations of",
        " parameters and multiple chains in parallel."
    )
}


#' @title Get celda model from a celda
#'  \link[SingleCellExperiment]{SingleCellExperiment} object
#' @description Return the celda model for \code{sce} returned by
#'  \link{celda_C}, \link{celda_G} or \link{celda_CG}.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object
#'  returned by \link{celda_C}, \link{celda_G}, or \link{celda_CG}.
#' @return Character. The celda model. Can be one of "celda_C", "celda_G", or
#'  "celda_CG".
#' @examples
#' data(sceCeldaCG)
#' celdaModel(sceCeldaCG)
#' @export
setGeneric("celdaModel",
    function(sce) {
        standardGeneric("celdaModel")
    })
#' @rdname celdaModel
#' @export
setMethod("celdaModel",
    signature(sce = "SingleCellExperiment"),
    function(sce) {
        tryCatch(
            if (S4Vectors::metadata(sce)$celda_parameters$model %in%
                    c("celda_C", "celda_G", "celda_CG")) {
                return(S4Vectors::metadata(sce)$celda_parameters$model)
            } else {
                stop("S4Vectors::metadata(sce)$celda_parameters$model must be",
                    " one of 'celda_C', 'celda_G', or 'celda_CG'")
            },
            error = function(e) {
                message("S4Vectors::metadata(sce)$celda_parameters$model must",
                    " exist! Try running celda model (celda_C, celda_CG, or",
                    " celda_G) first.")
                stop(e)
            })
    })


#' @title Get perplexity for every model in a celdaList
#' @description Returns perplexity for each model in a celdaList as calculated
#'  by `perplexity().`
#' @param celdaList An object of class celdaList.
#' @return List. Contains one celdaModel object for each of the parameters
#'  specified in the `runParams()` of the provided celda list.
#' @examples
#' data(celdaCGGridSearchRes)
#' celdaCGGridModelPerplexities <- celdaPerplexity(celdaCGGridSearchRes)
#' @export
setGeneric(
  "celdaPerplexity",
  function(celdaList) {
    standardGeneric("celdaPerplexity")
  }
)
#' @title Get perplexity for every model in a celdaList
#' @description Returns perplexity for each model in a celdaList as calculated
#'  by `perplexity().`
#' @param celdaList An object of class celdaList.
#' @return List. Contains one celdaModel object for each of the parameters
#'  specified in the `runParams()` of the provided celda list.
#' @examples
#' data(celdaCGGridSearchRes)
#' celdaCGGridModelPerplexities <- celdaPerplexity(celdaCGGridSearchRes)
#' @export
setMethod("celdaPerplexity",
  signature = c(celdaList = "celdaList"),
  function(celdaList) {
    celdaList@perplexity
  }
)


#' @title Append two celdaList objects
#' @description Returns a single celdaList representing the combination of two
#'  provided celdaList objects.
#' @return A celdaList object. This object contains all resList entries and
#'  runParam records from both lists.
#' @param list1 A celda_list object
#' @param list2 A celda_list object to be joined with list_1
#' @examples
#' data(celdaCGGridSearchRes)
#' appendedList <- appendCeldaList(
#'   celdaCGGridSearchRes,
#'   celdaCGGridSearchRes
#' )
#' @importFrom methods new
#' @export
appendCeldaList <- function(list1, list2) {
  if (!is.element("celdaList", class(list1)) |
    !is.element("celdaList", class(list2))) {
    stop("Both parameters to appendCeldaList must be of class celdaList.")
  }
  if (!(countChecksum(list1) == countChecksum(list2))) {
    warning(
      "Provided lists have different countChecksums and may have",
      " been generated from different count matrices. Using checksum",
      " from first list..."
    )
  }
  newList <- methods::new(
    "celdaList",
    runParams = rbind(runParams(list1), runParams(list2)),
    resList = c(resList(list1), resList(list2)),
    countChecksum = countChecksum(list1),
    perplexity = matrix(nrow = 0, ncol = 0)
  )
  return(newList)
}


#' @title Get the MD5 hash of the count matrix from the celdaList
#' @description Returns the MD5 hash of the count matrix used to generate the
#'  celdaList.
#' @param celdaList An object of class celdaList.
#' @return A character string of length 32 containing the MD5 digest of
#'  the count matrix.
#' @examples
#' data(celdaCGGridSearchRes)
#' countChecksum <- countChecksum(celdaCGGridSearchRes)
#' @export
setGeneric(
  "countChecksum",
  function(celdaList) {
    standardGeneric("countChecksum")
  }
)
#' @title Get the MD5 hash of the count matrix from the celdaList
#' @description Returns the MD5 hash of the count matrix used to generate the
#'  celdaList.
#' @param celdaList An object of class celdaList.
#' @return A character string of length 32 containing the MD5 digest of
#'  the count matrix.
#' @examples
#' data(celdaCGGridSearchRes)
#' countChecksum <- countChecksum(celdaCGGridSearchRes)
#' @export
setMethod("countChecksum",
  signature = c(celdaList = "celdaList"),
  function(celdaList) {
    celdaList@countChecksum
  }
)
