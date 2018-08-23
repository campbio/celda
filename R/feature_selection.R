#' Identify features with the highest influence on clustering
#' 
#' @param matrix Numeric matrix. 
#' @param n Integer. Maximum number of items above `threshold` returned for each ranked row or column.  
#' @param margin Integer. Dimension of `matrix` to rank, with 1 for rows, 2 for columns. Default 2. 
#' @param threshold Numeric. Only return ranked rows or columns in the matrix that are above this threshold. Default 0. 
#' @param decreasing Logical. Specifies if the rank should be decreasing. Default TRUE.  
#' @export
topRank = function(matrix, n=25, margin=2, threshold=0, decreasing=TRUE) {
  if(is.null(threshold) | is.na(threshold)) {
    threshold = min(matrix) - 1 
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
  top.ix = base::apply(matrix, margin, topFunction, thresh=threshold, n=n)
  
  ## Convert to list if apply converted to a matrix because all
  ## elements had the same length
  if(is.matrix(top.ix)) {
    top.ix = lapply(1:ncol(top.ix), function(i) top.ix[,i]) 
    names(top.ix) = dimnames(matrix)[[margin]]
  }  
  
  ## Parse names from returned margin
  opposite.margin = ifelse(margin - 1 > 0, margin-1, length(dim(matrix)))
  top.names = NULL
  names.to.parse = dimnames(matrix)[[opposite.margin]]
  if(!is.null(names.to.parse) & all(!is.na(top.ix))) {
    top.names = lapply(1:length(top.ix),
    			function(i) {
    			  ifelse(is.na(top.ix[[i]]), NA, names.to.parse[top.ix[[i]]])
    			})
    names(top.names) = names(top.ix)
  }
  
  return(list(index=top.ix, names=top.names))
}
