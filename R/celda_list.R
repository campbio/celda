################################################################################
# S3 Methods for celda_list objects                                            #
################################################################################


#' Select models from a list of results generated in a celda run.
#' 
#' Convenience function for picking out specific models from a celda_list object.
#' Models can be selected by various parameters, most importantly the K/L parameters (number of cell
#'  clusters / number of feature clusters). 
#' 
#' @param celda.list Object of class "celda_list". An object containing celda models returned from `celdaGridSearch()`.
#' @param K The K parameter for the desired model in the results list. Matches all K by default. Accepts ranges.
#' @param L The L parameter for the desired model in the results list. Matches all L by default. Accepts ranges.
#' @param chain The desired chain(s) for the specified model, for the specified K/L. Matches all chains by default. Accepts ranges.
#' @param index The index of the desired model in the run.parameters in the provided celda_list. Overrides all other parameters if provided. Defaults to NULL.
#' @return A celda model object matching the provided parameters, or a list of celda model objects if multiple models were matched (of class "celda_C", "celda_G", "celda_CG" accordingly), or NA if one is not found.
#' @examples
#' celda.mods = celdaGridSearch(celda::pbmc_select, model="celda_CG", K.to.test=c(10,15), 
#'                       L.to.test=c(50, 100), max.iter=2, nchains=1)
#' desired.mod = filterCeldaList(celda.mods, K=10, L=100)
#' @export
filterCeldaList = function(celda.list, params=NULL, index=NULL) {
  if (!isTRUE(class(celda.list)[1] == "celda_list")) stop("celda.list parameter was not of class celda_list.")

  # If user provides index / range of indices to select, override all else
  if (!is.null(index)) {
    if(!is.numeric(index) | length(index) != 1) stop("'index' needs to be an integer of length 1.") 
    ix = which(celda.list$run.params$index == index)
    if(length(ix) == 0) {
      stop(paste0("The index '", index, "' was not found in the 'run.params' within the given 'celda.list'."))
    }
    
    return(celda.list$res.list[[ix]])
    
  } else if(!is.null(params)) {
  
    ## Check for bad parameter names
    if(!all(names(params) %in% colnames(celda.list$run.params))) {
      bad.params = setdiff(names(params), colnames(celda.list$run.params))
      stop(paste0("The following elements in 'params' are not columns in celda.list$run.params: ", paste(bad.params, collapse=","), "."))
    }
    
    ## Subset 'run.params' based on items in 'params'
    new.run.params = celda.list$run.params
    for(i in names(params)) {
      new.run.params = subset(new.run.params, new.run.params[,i] %in% params[[i]])
      
      if(nrow(new.run.params) == 0) {
        stop("No runs matched the criteria given in 'params'. Check 'celda.list$run.params' for complete list of parameters used to generate 'celda.list'.")
      }
    }
    
    ## Get index of selected models, subset celda.list, and return
    ix = match(new.run.params$index, celda.list$run.params$index)
    if(length(ix) == 1) {
      return(celda.list$res.list[[ix]])
    } else {
      celda.list$run.params = as.data.frame(new.run.params)
      celda.list$res.list = celda.list$res.list[ix]
      return(celda.list)
    }
  } else {
    stop("One of 'params' or 'index' must be given.")
  }
   
  return(celda.list$res.list[filtered.run.params$index])
}


#' Select the best model from a celda_list by final log-likelihood.
#' 
#' This function returns the celda model from a celda_list object containing
#' the maxiumim final log-likelihood. K, L, or combination of K and L must be 
#' provided for celda_list objects containing celda_C, celda_G, and celda_CG 
#' models, respectively.
#' 
#' @param celda.list Object of class "celda_list". An object containing celda models returned from `celdaGridSearch()`.
#' @return The celda model object with the highest finalLogLik attribute, meeting any K/L criteria provided
#' @examples
#' celda.mods = celdaGridSearch(celda::pbmc_select, model="celda_CG", K.to.test=c(5,10), 
#'                       L=c(10,20,30), max.iter=2, nchains=1)
#' best.mod.k5.l20 = selectBestModel(celda.mods, K=5, L=20)
#' @export
selectBestModel = function(celda.list) {
  if (!isTRUE(class(celda.list)[1] == "celda_list")) stop("celda.list parameter was not of class celda_list.")
 
  group = setdiff(colnames(celda.list$run.params), c("index", "chain", "log_likelihood"))
  new.run.params = celda.list$run.params %>% group_by(.dots=group) %>% slice(which.max(log_likelihood))
  
  ix = match(new.run.params$index, celda.list$run.params$index)
  if(nrow(new.run.params) == 1) {
    return(celda.list$res.list[[ix]])
  } else {
    celda.list$run.params = as.data.frame(new.run.params)
    celda.list$res.list = celda.list$res.list[ix]
    return(celda.list)
  }
}


validateGetModelParams = function(celda.list, K, L, chain) {
  if (class(celda.list)[1] != "celda_list") stop("First argument to getModel() should be an object of class 'celda_list'")
  
  if ((is.null(K) | is.null(L)) & celda.list$content.type == "celda_CG") {
    stop("Both K and L parameters needed for subsetting celda_CG result lists")
  }

  if (is.null(K) & celda.list$content.type == "celda_C") {
    stop("K parameter needed when subsetting celda_C result lists")
  }

  if (is.null(L) & celda.list$content.type == "celda_G") {
    stop("L parameter needed when subsetting celda_G result lists")
  }
  
  if (!is.null(K)){
    if(!(K %in% celda.list$run.params$K)) {
      stop("Provided K was not profiled in the provided celda_list object")
    } 
  }
  
  if (!is.null(L)){ 
    if(!(L %in% celda.list$run.params$L)) {
      stop("Provided L was not profiled in the provided celda_list object")
    }
  }
}
