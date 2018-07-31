#' @useDynLib celda _rowSumByGroup
rowSumByGroup <- function(x, group, L) {
  group <- factor(group, levels=1:L)
  res <- .Call("_rowSumByGroup", x, group)
  return(res)
}

#' @useDynLib celda _rowSumByGroupChange
rowSumByGroupChange <- function(x, px, group, pgroup, L) {
  group <- factor(group, levels=1:L)
  pgroup <- factor(pgroup, levels=1:L)  
  res <- .Call("_rowSumByGroupChange", x, px, group, pgroup)
  return(res)
}

#' @useDynLib celda _colSumByGroup
colSumByGroup <- function(x, group, K) {
  group <- factor(group, levels=1:K)
  res <- .Call("_colSumByGroup", x, group)
  return(res)
}

#' @useDynLib celda _colSumByGroupChange
colSumByGroupChange <- function(x, px, group, pgroup, K) {
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

#' @useDynLib celda _mvMult
mvMult <- function(x, v, by.row) {
  res <- .Call("_mvMult", x, v, by.row) 
  return(res)
}

#' @useDynLib celda _mvAdd 
mvAdd <- function(x, v, by.row) {
  res <- .Call("_mvAdd", x, v, by.row)
  return(res)
}