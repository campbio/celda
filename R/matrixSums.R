#' @useDynLib celda _rowSumByGroup
rowSumByGroup <- function(x, group, L) {
  if(any(group > L)) {
    stop("An entry in 'group' is greater than L.")
  }

  group <- factor(group, levels=1:L)
  res <- .Call("_rowSumByGroup", x, group)
  return(res)
}

#' @useDynLib celda _rowSumByGroupChange
rowSumByGroupChange <- function(x, px, group, pgroup, L) {
  if(any(group > L)) {
    stop("An entry in 'group' is greater than L.")
  }
  if(any(pgroup > L)) {
    stop("An entry in 'pgroup' is greater than L.")
  }
  
  group <- factor(group, levels=1:L)
  pgroup <- factor(pgroup, levels=1:L)  
  res <- .Call("_rowSumByGroupChange", x, px, group, pgroup)
  return(res)
}

#' @useDynLib celda _colSumByGroup
colSumByGroup <- function(x, group, K) {
  if(any(group > K)) {
    stop("An entry in 'group' is greater than K.")
  }
  group <- factor(group, levels=1:K)
  res <- .Call("_colSumByGroup", x, group)
  return(res)
}

#' @useDynLib celda _colSumByGroupChange
colSumByGroupChange <- function(x, px, group, pgroup, K) {
  if(any(group > K)) {
    stop("An entry in 'group' is greater than K.")
  }
  if(any(pgroup > K)) {
    stop("An entry in 'pgroup' is greater than K.")
  }
  
  group <- factor(group, levels=1:K)
  pgroup <- factor(pgroup, levels=1:K)    
  res <- .Call("_colSumByGroupChange", x, px, group, pgroup)
  return(res)
}


#' @useDynLib celda _rowSumByGroup_numeric 
rowSumByGroup.numeric <- function(x, group, L) {
  group <- factor(group, levels=1:L)
  res <- .Call("_rowSumByGroup_numeric", x, group)
  return(res)
}

#' @useDynLib celda _colSumByGroup_numeric
colSumByGroup.numeric <- function(x, group, K) {
  group <- factor(group, levels=1:K)
  res <- .Call("_colSumByGroup_numeric", x, group)
  return(res)
}
