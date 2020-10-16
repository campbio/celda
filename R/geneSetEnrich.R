#' @title Gene set enrichment
#' @description Identify and return significantly-enriched terms for each gene
#'  module in a Celda object or a \linkS4class{SingleCellExperiment} object.
#'  Performs gene set enrichment analysis for Celda
#'  identified modules using the \link[enrichR]{enrichr}.
#' @author Ahmed Youssef, Zhe Wang
#' @param x A numeric \link{matrix} of counts or a
#'  \linkS4class{SingleCellExperiment}
#'  with the matrix located in the assay slot under \code{useAssay}.
#'  Rows represent features and columns represent cells. Rownames of the
#'  matrix or \linkS4class{SingleCellExperiment} object should be gene names.
#' @param useAssay A string specifying which \link{assay}
#'  slot to use if \code{x} is a
#'  \linkS4class{SingleCellExperiment} object. Default "counts".
#' @param altExpName The name for the \link{altExp} slot
#'  to use. Default "featureSubset".
#' @param celdaModel Celda object of class \code{celda_G} or \code{celda_CG}.
#' @param databases Character vector. Name of reference database. Available
#'  databases can be viewed by \link[enrichR]{listEnrichrDbs}.
#' @param fdr False discovery rate (FDR). Numeric. Cutoff value for adjusted
#'  p-value, terms with FDR below this value are considered significantly
#'  enriched.
#' @param ... Ignored. Placeholder to prevent check warning.
#' @return List of length 'L' where each member contains the significantly
#'  enriched terms for the corresponding module.
#' @importFrom enrichR enrichr
#' @importFrom enrichR listEnrichrDbs
#' @export
setGeneric("geneSetEnrich", function(x, ...) {
    standardGeneric("geneSetEnrich")})


#' @rdname geneSetEnrich
#' @examples
#' library(M3DExampleData)
#' counts <- M3DExampleData::Mmus_example_list$data
#' # subset 500 genes for fast clustering
#' counts <- counts[seq(1501, 2000), ]
#' # cluster genes into 10 modules for quick demo
#' sce <- celda_G(x = as.matrix(counts), L = 10, verbose = FALSE)
#' gse <- geneSetEnrich(sce,
#'   databases = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018"))
#' @export
setMethod("geneSetEnrich",
    signature(x = "SingleCellExperiment"),
    function(x,
        useAssay = "counts",
        altExpName = "featureSubset",
        databases,
        fdr = 0.05) {

        altExp <- SingleCellExperiment::altExp(x, e = altExpName)

        # initialize list with one entry for each gene module
        modules <- vector("list",
            length = S4Vectors::metadata(altExp)$celda_parameters$L)

        # create dataframe with gene-module associations
        genes <- data.frame(gene = rownames(altExp),
            module = celdaModules(x, altExpName = altExpName))

        # iterate over each module, get genes in that module, add to list
        for (i in seq_len(S4Vectors::metadata(altExp)$celda_parameters$L)) {
            modules[[i]] <- as.character(genes[genes$module == i, "gene"])
        }

        # enrichment analysis
        enrichment <- lapply(modules, function(module) {
            invisible(utils::capture.output(table <- enrichR::enrichr(
                genes = module,
                databases = databases)))
            table <- Reduce(f = rbind, x = table)
            table[table$Adjusted.P.value < fdr, "Term"]
        })

        # return results as a list
        return(enrichment)
    }
)


#' @rdname geneSetEnrich
#' @export
setMethod("geneSetEnrich",
    signature(x = "matrix"),
    function(x, celdaModel, databases, fdr = 0.05) {
        # check for correct celda object
        if (!(class(celdaModel) %in% c("celda_G", "celda_CG"))) {
            stop(
                "No gene modules in celda object. ",
                "Please provide object of class celda_G or celda_CG."
            )
        }

        # initialize list with one entry for each gene module
        modules <- vector("list", length = params(celdaModel)$L)

        # create dataframe with gene-module associations
        genes <- data.frame(gene = rownames(x),
            module = celdaClusters(celdaModel)$y)

        # iterate over each module, get genes in that module, add to list
        for (i in seq_len(params(celdaModel)$L)) {
            modules[[i]] <- as.character(genes[genes$module == i, "gene"])
        }

        # enrichment analysis
        enrichment <- lapply(modules, function(module) {
            invisible(utils::capture.output(table <- enrichR::enrichr(
                genes = module,
                databases = databases
            )))
            table <- Reduce(f = rbind, x = table)
            table[table$Adjusted.P.value < fdr, "Term"]
        })

        # return results as a list
        return(enrichment)
    }
)
