## This script defines a function that performs gene enrichment analysis for celda-identified modules using the enrichR tool ##

#' @title Gene set enrichment
#' @description Identify and return significantly-enriched terms for each gene module in a Celda object
#' @author Ahmed Youssef
#' 
#' @param counts Integer matrix. Rows represent genes and columns represent cells. 
#' @param celda Celda object of class `celda_G` or `celda_CG`. 
#' @param databases Character. Name of reference database.
#' @param fdr Numeric. Cutoff value for adjusted p-value, terms below this value are considered significantly-enriched.
#' @return List of length 'L' where each member contains the significantly-enriched terms for corresponding module. 
#' @example celdaEnrichment(counts, celda, databases = c('GO_Biological_Process_2018','GO_Molecular_Function_2018'))

celdaEnrichment <- function(counts, celda, databases, fdr = 0.05){
  #check for correct celda object
  if(!(class(celda) %in% c('celda_G', 'celda_CG')))
    stop('No gene modules in celda object. Please provide object of class celda_G or celda_CG.')
  
  #initialize list with one entry for each gene module
  modules <- vector("list", length = celda@params$L)
  
  #create dataframe with gene-module associations
  genes <- data.frame(gene = rownames(counts), module = celda@clusters$y)
  
  #iterate over each module, get genes in that module, add to list
  for(i in seq_len(celda@params$L))
    modules[[i]] <- as.character(genes[genes$module==i,'gene'])
  
  #enrichment analysis
  enrichment <- lapply(modules, function(module){
    invisible(capture.output(table <- enrichR::enrichr(genes = module, databases = databases)))
    table <- Reduce(f = rbind, x = table)
    table[table$Adjusted.P.value < fdr, 'Term']
  })
  
  #return results as a list
  return(enrichment)
}