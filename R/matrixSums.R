#' @useDynLib celda rowSumByGroup_
rowSumByGroup <- function(x, group, L) {
  group <- factor(group, levels=1:L)
  res <- .Call("rowSumByGroup", x, group)
  return(res)
}

#' @useDynLib celda rowSumByGroupChange_
rowSumByGroupChange <- function(x, px, group, pgroup, L) {
  group <- factor(group, levels=1:L)
  pgroup <- factor(pgroup, levels=1:L)  
  res <- .Call("rowSumByGroupChange", x, px, group, pgroup)
  return(res)
}

#' @useDynLib celda colSumByGroup_
colSumByGroup <- function(x, group, K) {
  group <- factor(group, levels=1:K)
  res <- .Call("colSumByGroup", x, group)
  return(res)
}

#' @useDynLib celda colSumByGroupChange_
colSumByGroupChange <- function(x, px, group, pgroup, K) {
  group <- factor(group, levels=1:K)
  pgroup <- factor(pgroup, levels=1:K)    
  res <- .Call("colSumByGroupChange", x, px, group, pgroup)
  return(res)
}

