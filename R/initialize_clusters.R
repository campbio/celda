
initialize.cluster = function(N, len, z = NULL, initial = NULL, fixed = NULL, seed=12345) {

  ## If initial values are given, then they will not be randomly initialized
  if(!is.null(initial)) {
    init.values = sort(unique(initial))
    if(length(unique(initial)) != N || length(initial) != len || !all(init.values %in% 1:N)) {
      stop("'initial' needs to be a vector of length 'len' containing N unique values.")
    }
    z = as.numeric(as.factor(initial))
  } else {
    z = rep(NA, len)
  } 
  
  ## Set any values that need to be fixed during sampling
  if(!is.null(fixed)) {
    fixed.values = sort(unique(fixed))
    if(length(fixed) != len || !all(fixed.values %in% 1:N)) {
      stop("'fixed' to be a vector of length 'len' where each entry is one of N unique values or NA.")
    }
    fixed.ix = !is.na(fixed)
    z[fixed.ix] = fixed[fixed.ix]
    z.not.used = setdiff(1:N, unique(fixed[fixed.ix]))
  } else {
    z.not.used = 1:N
    fixed.ix = rep(FALSE, len)
  }  

  ## Randomly sample remaining values
  if (!is.null(seed)) {
    set.seed(seed)
  }
  z.na = which(is.na(z))
  if(length(z.na) > 0) {
    z[z.na] = sample(z.not.used, length(z.na), replace=TRUE)
  }  
  
  ## Check to ensure each value is in the vector at least once
  missing = setdiff(1:N, z)
  for(i in missing) {
    ta = sort(table(z[!fixed.ix]), decreasing = TRUE)
    if(ta[1] == 1) {
      stop("'len' is not long enough to accomodate 'N' unique values")
    }
    ix = which(z == as.numeric(names(ta))[1] & !fixed.ix)
    z[sample(ix, 1)] = i
  }
  
  return(z)
}

