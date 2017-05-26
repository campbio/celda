#' Identify features with the highest influence on clustering
#' 
#' @param fm factorized matrix.
#' @param n Maximum number of items returned for each entry. 
#' @param margin 1 for rows, 2 for columns.
#' @param threshold only include entries in the matrix that is greader than the threshold.
#' @param decreasing Logical; specifying if rank should be decreasing. Default to be TRUE. 
#' @export
topRank = function(fm, n=25, margin=2, threshold=0, decreasing=TRUE) {
  if(is.null(threshold) | is.na(threshold)) {
    threshold = min(m) - 1 
  }
  
  ## Function to sort values in a vector and return 'n' top results
  ## If there are not 'n' top results above 'thresh', then the
  ## number of entries in 'v' that are above 'thresh' will be returned
  topFunction = function(v, n, thresh) {
    v.above.thresh = sum(v > thresh)
    n.to.select = min(v.above.thresh, n)
    
    h = NA
    if(n.to.select > 0) {
      h = head(order(v, decreasing=decreasing), n.to.select)
    }  
    return(h)
  }

  ## Parse top ranked indices from matrix
  top.ix = apply(fm, margin, topFunction, thresh=threshold, n=n)
  
  ## Convert to list if apply converted to a matrix because all
  ## elements had the same length
  if(is.matrix(top.ix)) {
    top.ix = lapply(1:ncol(top.ix), function(i) top.ix[,i]) 
    names(top.ix) = dimnames(fm)[[margin]]
  }  
  
  ## Parse names from returned margin
  opposite.margin = ifelse(margin - 1 > 0, margin-1, length(dim(fm)))
  top.names = NULL
  names.to.parse = dimnames(fm)[[opposite.margin]]
  if(!is.null(names.to.parse) & all(!is.na(top.ix))) {
    top.names = lapply(1:length(top.ix),
    			function(i) {
    			  ifelse(is.na(top.ix[[i]]), NA, names.to.parse[top.ix[[i]]])
    			})
    names(top.names) = names(top.ix)
  }
  
  return(list(index=top.ix, names=top.names))
}


