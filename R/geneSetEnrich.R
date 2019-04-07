#' @title Gene set enrichment
#' @description Identify and return significantly-enriched terms for each gene
#'  module in a Celda object. Performs gene set enrichment analysis for Celda
#'  identified modules using the enrichR package.
#' @author Ahmed Youssef
#' @param counts Integer count matrix. Rows represent genes and columns
#'  represent cells. Row names of the matrix should be gene names.
#' @param celdaModel Celda object of class `celda_G` or `celda_CG`.
#' @param databases Character vector. Name of reference database. Available
#'  databases can be viewed by running \code{enrichR::listEnrichrDbs()}.
#' @param fdr False discovery rate (FDR). Numeric. Cutoff value for adjusted
#'  p-value, terms with FDR below this value are considered significantly
#'  enriched.
#' @return List of length 'L' where each member contains the significantly
#'  enriched terms for the corresponding module.
#' @examples
#' library(M3DExampleData)
#' counts <- M3DExampleData::Mmus_example_list$data
#' #subset 100 genes for fast clustering
#' counts <- counts[seq(1200, 2000), ]
#' #cluster genes into 10 modules for quick demo
#' cm <- celda_G(counts = as.matrix(counts), L = 10, verbose = FALSE)
#' gse <- geneSetEnrich(counts,
#'     cm,
#'     databases = c('GO_Biological_Process_2018','GO_Molecular_Function_2018'))
#' @export
geneSetEnrich <- function(counts, celdaModel, databases, fdr = 0.05) {
    #check for correct celda object
    if (!(class(celdaModel) %in% c("celda_G", "celda_CG"))) {
        stop("No gene modules in celda object. ",
            "Please provide object of class celda_G or celda_CG.")
    }

    #initialize list with one entry for each gene module
    modules <- vector("list", length = celdaModel@params$L)

    #create dataframe with gene-module associations
    genes <- data.frame(gene = rownames(counts), module = celdaModel@clusters$y)

    #iterate over each module, get genes in that module, add to list
    for (i in seq_len(celdaModel@params$L)) {
        modules[[i]] <- as.character(genes[genes$module == i, "gene"])
    }

    #enrichment analysis
    enrichment <- lapply(modules, function(module) {
        invisible(utils::capture.output(table <- enrichR::enrichr(
            genes = module,
            databases = databases)))
        table <- Reduce(f = rbind, x = table)
        table[table$Adjusted.P.value < fdr, "Term"]
    })

    #return results as a list
    return(enrichment)
}
