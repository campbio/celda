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
