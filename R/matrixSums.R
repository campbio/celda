.rowSumByGroup <- function(counts, group, L) {
  if (inherits(counts, "matrix") & is.integer(counts)) {
    res <- .rowSumByGroupInteger(counts, group, L)
  } else if (inherits(counts, "matrix") & is.numeric(counts)) {
    res <- .rowSumByGroupNumeric(counts, group, L)
  } else if (inherits(counts, "dgCMatrix")) {
    res <- rowSumByGroupSparse(counts, group, L)
  } else {
    stop("'counts' must be an integer, numeric, or dgCMatrix matrix.")
  }
  return(res)
}

.rowSumByGroupChange <- function(counts, pcounts, group, pgroup, L) {
  if (inherits(counts, "matrix") & is.integer(counts)) {
    res <- .rowSumByGroupChangeInteger(counts, pcounts, group, pgroup, L)
  } else if (inherits(counts, "matrix") & is.numeric(counts)) {
    res <- .rowSumByGroupChangeNumeric(counts, pcounts, group, pgroup, L)
  } else if (inherits(counts, "dgCMatrix")) {
    res <- rowSumByGroupChangeSparse(counts, pcounts, group, pgroup, L)
  } else {
    stop("'counts' must be an integer, numeric, or dgCMatrix matrix.")
  }
  return(res)
}

.colSumByGroup <- function(counts, group, K) {
  if (inherits(counts, "matrix") & is.integer(counts)) {
    res <- .colSumByGroupInteger(counts, group, K)
  } else if (inherits(counts, "matrix") & is.numeric(counts)) {
    res <- .colSumByGroupNumeric(counts, group, K)
  } else if (inherits(counts, "dgCMatrix")) {
    res <- colSumByGroupSparse(counts, group, K)
  } else {
    stop("'counts' must be an integer, numeric, or dgCMatrix matrix.")
  }
  return(res)
}

.colSumByGroupChange <- function(counts, pcounts, group, pgroup, K) {
  if (inherits(counts, "matrix") & is.integer(counts)) {
    res <- .colSumByGroupChangeInteger(counts, pcounts, group, pgroup, K)
  } else if (inherits(counts, "matrix") & is.numeric(counts)) {
    res <- .colSumByGroupChangeNumeric(counts, pcounts, group, pgroup, K)
  } else if (inherits(counts, "dgCMatrix")) {
    res <- colSumByGroupChangeSparse(counts, pcounts, group, pgroup, K)
  } else {
    stop("'counts' must be an integer, numeric, or dgCMatrix matrix.")
  }
  return(res)
}



#' @useDynLib celda _rowSumByGroup
.rowSumByGroupInteger <- function(x, group, L) {
  group <- factor(group, levels = seq(L))
  res <- .Call("_rowSumByGroup", x, group)
  return(res)
}

#' @useDynLib celda _rowSumByGroupChange
.rowSumByGroupChangeInteger <- function(x, px, group, pgroup, L) {
  group <- factor(group, levels = seq(L))
  pgroup <- factor(pgroup, levels = seq(L))
  res <- .Call("_rowSumByGroupChange", x, px, group, pgroup)
  return(res)
}

#' @useDynLib celda _colSumByGroup
.colSumByGroupInteger <- function(x, group, K) {
  group <- factor(group, levels = seq(K))
  res <- .Call("_colSumByGroup", x, group)
  return(res)
}

#' @useDynLib celda _colSumByGroupChange
.colSumByGroupChangeInteger <- function(x, px, group, pgroup, K) {
  group <- factor(group, levels = seq(K))
  pgroup <- factor(pgroup, levels = seq(K))
  res <- .Call("_colSumByGroupChange", x, px, group, pgroup)
  return(res)
}


#' @useDynLib celda _rowSumByGroup_numeric
.rowSumByGroupNumeric <- function(x, group, L) {
  group <- factor(group, levels = seq(L))
  res <- .Call("_rowSumByGroup_numeric", x, group)
  return(res)
}

#' @useDynLib celda _colSumByGroup_numeric
.colSumByGroupNumeric <- function(x, group, K) {
  group <- factor(group, levels = seq(K))
  res <- .Call("_colSumByGroup_numeric", x, group)
  return(res)
}

#' @useDynLib celda _rowSumByGroupChange_numeric
.rowSumByGroupChangeNumeric <- function(x, px, group, pgroup, L) {
  group <- factor(group, levels = seq(L))
  pgroup <- factor(pgroup, levels = seq(L))
  res <- .Call("_rowSumByGroupChange_numeric", x, px, group, pgroup)
  return(res)
}

#' @useDynLib celda _colSumByGroupChange_numeric
.colSumByGroupChangeNumeric <- function(x, px, group, pgroup, K) {
  group <- factor(group, levels = seq(K))
  pgroup <- factor(pgroup, levels = seq(K))
  res <- .Call("_colSumByGroupChange_numeric", x, px, group, pgroup)
  return(res)
}

#' @useDynLib celda _perplexityG
.perplexityGLogPx <- function(x, phi, psi, group, L) {
  group <- factor(group, levels = seq(L))
  res <- .Call("_perplexityG", x, phi, psi, group)
  return(res)
}
