sample.ll = function(ll.probs) {
  probs.sub <- exp(ll.probs - max(ll.probs))
  probs.norm <- probs.sub / sum(probs.sub)
  probs.select <- sample.int(length(probs.norm), size=1L, replace=TRUE, prob=probs.norm)
  return(probs.select)
}


cosineDist = function(x){
  x <- t(x)
  y <- (1 - cosine(x)) / 2
  return(stats::as.dist(y))
}

cosine = function(x) {
  y = x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))
  return(y)
}
  
spearmanDist = function(x){
  y = (1 - stats::cor(x, method="spearman"))/2
  return(stats::as.dist(y))
}

hellingerDist = function(x) {
  y = dist(t(sqrt(x)), method = "euclidean") * 1/sqrt(2)
  return(y)
}  


normalizeLogProbs = function(ll.probs) {
  ll.probs <- exp(sweep(ll.probs, 1, base::apply(ll.probs, 1, max), "-"))
  probs <- sweep(ll.probs, 1, rowSums(ll.probs), "/")
  return(probs)
}


#' Performs normalization, transformation, and/or scaling on a counts matrix
#' 
#' @param counts A count matrix 
#' @param normalize Character. Divides counts by the library sizes for each cell. One of "proportion", "cpm", "median", or "mean". "proportion" uses the total counts for each cell as the library size. "cpm" divides the library size of each cell by one million to produce counts per million. "median" divides the library size of each cell by the median library size across all cells.  "mean" divides the library size of each cell by the mean library size across all cells.
#' @param transformation.fun Function. Applys a transformation such as `sqrt`, `log`, `log2`, `log10`, or `log1p`. If NULL, no transformation will be applied. Occurs after normalization. Default NULL.
#' @param scale.fun Function. Scales the rows of the normalized and transformed count matrix. Default NULL.
#' @param pseudocount.normalize Numeric. Add a pseudocount to counts before normalization. Default  0. 
#' @param pseudocount.transform Numeric. Add a pseudocount to normalized counts before applying the transformation function. Adding a pseudocount can be useful before applying a log transformation. Default  0. 
#' @export
normalizeCounts = function(counts, normalize=c("proportion", "cpm", "median", "mean"),
							transformation.fun=NULL, scale.fun=NULL,
							pseudocount.normalize=0, pseudocount.transform=0) {

  normalize = match.arg(normalize)
  if(!is.null(transformation.fun) && !is.function(transformation.fun)) {
    stop("'transformation.fun' needs to be of class 'function'")
  }
  if(!is.null(scale.fun) && !is.function(scale.fun)) {
    stop("'scale.fun' needs to be of class 'function'")
  }

  # Perform normalization  
  counts = counts + pseudocount.normalize
  cs = .colSums(counts, nrow(counts), ncol(counts))
  norm = switch(normalize,
    "proportion" = sweep(counts, 2, cs, "/"),
    "cpm" = sweep(counts, 2, cs / 1e6, "/"),
    "median" = sweep(counts, 2, cs / median(cs), "/"),
    "mean" = sweep(counts, 2, cs / mean(cs), "/")
  )  
  
  if(!is.null(transformation.fun)){
    norm <- do.call(transformation.fun, list(norm + pseudocount.transform))
  }
  if(!is.null(scale.fun)) {
    norm <- t(base::apply(norm, 1, scale.fun))
  }  

  return(norm)
}
  
  
reorder.label.by.size = function(z, K) {
  z.ta = as.numeric(names(sort(table(factor(z, levels=1:K)), decreasing=TRUE)))
  
  new.z = z
  for(i in 1:length(z.ta)) {
    new.z[z == z.ta[i]] = i
  }
  return(list(new.labels=new.z, map=z.ta))
}  

reorder.labels.by.size.then.counts = function(counts, z, y, K, L) {
  z.ta = as.numeric(names(sort(table(factor(z, levels=1:K)), decreasing=TRUE)))
  
  new.z = rep(NA, length(z))
  for(i in 1:length(z.ta)) {
    new.z[z == z.ta[i]] = i
  }

  n.TS.by.C = rowsum(counts, group=y, reorder=TRUE)
  n.CP.by.TS = rowsum(t(n.TS.by.C), group=new.z, reorder=TRUE)
  TS.order = lapply(1:K, function(i) order(n.CP.by.TS[i,], decreasing=TRUE))  
  
  ## Determine the number of transcriptional states to order by in each cell population
  if(K > L) {
    num.L.per.K = rep(1:0, c(K, K-L))
  } else { 
    #num.L.per.K = rep(c(ceiling(L/K), floor(L/K)), c(ceiling(K/2),floor(K/2)))
    temp = rep(floor(L/K), K)
    num.L.per.K = temp + rep(1:0, c(L-sum(temp), K-L+sum(temp)))
  }

  ## Cycle through each K and assign the state(s) with the highest counts
  y.to.choose.from = 1:L  
  cp.with.ts = which(num.L.per.K > 0)  
  y.order = c()
  for(i in cp.with.ts) {
    ix = setdiff(TS.order[[i]], y.order)[1:num.L.per.K[i]]
    y.order = c(y.order, ix)
    y.to.choose.from = setdiff(y.to.choose.from, ix)
  }
  
  new.y = rep(NA, length(y))
  for(i in 1:length(y.order)) {
    new.y[y == y.order[i]] = i
  }

  return(list(z=new.z, y=new.y, z.map=z.ta, y.map=y.order))  
}  


#' Obtain the gene module of a gene of interest
#' 
#' This function will output the gene module of a specific gene(s) from a celda model
#'  
#' @param counts Integer matrix. Rows represent features and columns represent cells. This matrix should be the same as the one used to generate `celda.mod`.
#' @param celda.mod Model of class "celda_G" or "celda_CG".
#' @param gene Character vector. Identify feature modules for the specified feature names. 
#' @export
featureModuleLookup <- function(counts, celda.mod, gene){
  list <- list()
  for(x in 1:length(gene)){
    if(gene[x] %in% rownames(counts)){
      list[x] <- model$y[which(rownames(counts) == gene[x])]
    }else{
      list[x] <- c("The gene you selected does not exist within your data")
    }
  } 
  names(list) <- gene
  return(list)
}


#' Re-code cell cluster labels by provided mapping scheme
#' 
#' This function will re-code _cell_ cluster labels based off of a mapping provided by the user,
#' for all fields on a celda object involving cluster labels.
#' e.g. if Z (cell cluster) values range from 1-4 and the user would like all 3's switched to 1's and
#' vice versa, this function can be useful. NOTE: it is recommended that this function's results
#' aren't used to overwrite the original celda model object provided, in the event of a mis-mapping.
#' 
#' @param celda.mod Celda object of class "celda_C" or "celda_CG". 
#' @param from Numeric vector. Unique values in the range of 1:K that correspond to the original cluster labels in `celda.mod`.  
#' @param to Numeric vector. Unique values in the range of 1:K that correspond to the new cluster labels. 
#' @export
recodeClusterZ = function(celda.mod, from, to) {
  if (length(setdiff(from, to)) != 0) {
    stop("All values in 'from' must have a mapping in 'to'")
  }
  if (is.null(celda.mod$z)) {
    stop("Provided celda.mod argument does not have a z attribute")
  }
  celda.mod$z = plyr::mapvalues(celda.mod$z, from, to)
  return(celda.mod)
}
  

#' Re-code gene cluster labels by provided mapping scheme
#' 
#' This function will re-code _gene_ cluster labels based off of a mapping provided by the user,
#' for all fields on a celda object involving cluster labels.
#' e.g. if Y (gene cluster) values range from 1-4 and the user would like all 3's switched to 1's and
#' vice versa, this function can be useful. NOTE: it is recommended that this function's results
#' aren't used to overwrite the original celda model object provided, in the event of a mis-mapping.
#' 
#' @param celda.mod Celda object of class "celda_G" or "celda_CG". 
#' @param from Numeric vector. Unique values in the range of 1:L that correspond to the original cluster labels in `celda.mod`. 
#' @param to Numeric vector. Unique values in the range of 1:L that correspond to the new cluster labels. 
#' @export
recodeClusterY = function(celda.mod, from, to) {
  if (length(setdiff(from, to)) != 0) {
    stop("All values in 'from' must have a mapping in 'to'")
  }
  if (is.null(celda.mod$y)) {
    stop("Provided celda.mod argument does not have a y attribute")
  }
  celda.mod$y = plyr::mapvalues(celda.mod$y, from, to)
  return(celda.mod)
}


#' Check whether a count matrix was the one used in a given celda run
#' 
#' Ensures that the provided celda object was generated from a counts matrix
#' with similar dimensions to the one provided. 
#' 
#' Then, compare the MD5 checksum of a provided count.matrix to the count matrix
#' checksum on a celda_list object, to see if they're the same.
#' @param counts Integer matrix. Rows represent features and columns represent cells. 
#' @param celda.mod Celda model. Options available in `celda::available.models`.
#' @param error.on.mismatch Logical. Whether to stop execution in the event of a count matrix mismatch. Default TRUE.
#' @return TRUE if provided count matrix matches the one used in the celda run, FALSE otherwise. Error on FALSE if error.on.mismatch is TRUE.
#' @export
compareCountMatrix = function(counts, celda.mod, error.on.mismatch=TRUE) {
  if (length(celda.mod$y != 0) & nrow(counts) != length(celda.mod$y)) {
    stop("The provided celda object was generated from a counts matrix with a different number of genes than the one provided.")
  }
  
  if (length(celda.mod$z != 0 ) & ncol(counts) != length(celda.mod$z)) {
    stop("The provided celda object was generated from a counts matrix with a different number of cells than the one provided.")
  }
  
  celda.checksum = celda.mod$count.checksum
  count.md5 = digest::digest(counts, algo="md5")
  res = isTRUE(count.md5 == celda.checksum)
  if (res) return(TRUE)
  if (!res && error.on.mismatch) stop("There was a mismatch between the provided count matrix and the count matrix used to generate the provided celda result.")
  else if (!res && !error.on.mismatch) return(FALSE)
}


# Check whether the given celda_list's run.params attribute ordering matches it's res.list models
# 
# @param celda.list A celda_list model as generated by celda().
# @return TRUE if the celda.list's run.params K/L attributes are ordered similarly to the res.list models' K/L values, FALSE otherwise.
validateRunParams = function(celda.list) {
  if (class(celda.list)[1] != "celda_list") {
    stop("Provided object is not of class celda_list")
  }
  
  run.params = celda.list$run.params
  new.run.params = newRunParamsFromResList(celda.list)
  if (!is.null(run.params$K)) {
    if (!isTRUE(all.equal(new.run.params$K, run.params$K))) return(FALSE)
  }
  if (!is.null(run.params$L)) {
    if (!isTRUE(all.equal(new.run.params$L, run.params$L))) return(FALSE)
  }
  
  return(TRUE)
} 


newRunParamsFromResList = function(celda.list) {
  new.res.list = switch(celda.list$content.type,
                        "celda_CG" = data.frame(K=sapply(celda.list$res.list, getK), 
                                                L=sapply(celda.list$res.list, getL)),
                        "celda_C" = data.frame(K=sapply(celda.list$res.list, getK)), 
                        "celda_G" = data.frame(L=sapply(celda.list$res.list, getL))) 
  new.res.list$index = 1:nrow(new.res.list)
  return(new.res.list)
}


logMessages = function(..., sep = " ", logfile = NULL, append = FALSE) {
  if(!is.null(logfile)) {
    if(!is.character(logfile) || length(logfile) > 1) {
      stop("The log file parameter needs to be a single character string.")
    }    
    cat(paste(..., "\n", sep=sep), file=logfile, append=append)
    
  } else {
    message(paste(..., sep=sep))
  }
}





#' Generate a distinct palette for coloring different clusters
#' 
#' @param n Integer. Number of colors to generate. 
#' @param hues Character vector. Colors available from `colors()`. These will be used as the base colors for the clustering scheme in HSV. Different saturations and values will be generated for each hue. Default c("red", "cyan", "orange", "blue", "yellow", "purple", "green", "magenta").
#' @param saturation.range Numeric vector. A vector of length 2 denoting the saturation for HSV. Values must be in [0,1]. Default: c(0.25, 1).
#' @param value.range Numeric vector. A vector of length 2 denoting the range of values for HSV. Values must be in [0,1]. Default: `c(0.5, 1)`.
#' @return A vector of distinct colors that have been converted to  HEX from HSV.
#' @export
distinct_colors = function(n,
						   hues = c("red", "cyan", "orange", "blue", "yellow", "purple", "green", "magenta"),
                           saturation.range = c(0.7, 1),
                           value.range = c(0.7, 1)) {
                           
  if(!(all(hues %in% grDevices::colors()))) {
    stop("Only color names listed in the 'color' function can be used in 'hues'")
  }
  
  ## Convert R colors to RGB and then to HSV color format
  hues.hsv = grDevices::rgb2hsv(grDevices::col2rgb(hues))
  
  ## Calculate all combination of saturation/value pairs
  ## Note that low saturation with low value (i.e. high darkness) is too dark for all hues
  ## Likewise, high saturation with high value (i.e. low darkness) is hard to distinguish
  ## Therefore, saturation and value are set to be anticorrelated
  num.vs = ceiling(n / length(hues))
  s = seq(from=saturation.range[1], to=saturation.range[2], length=num.vs)
  v = seq(from=value.range[2], to=value.range[1], length=num.vs)

  ## Create all combination of hues with saturation/value pairs
  new.hsv = c()
  for(i in 1:num.vs) {
    temp = rbind(hues.hsv[1,], s[i], v[i])
    new.hsv = cbind(new.hsv, temp)  
  }

  ## Convert to hex
  col = grDevices::hsv(new.hsv[1,], new.hsv[2,], new.hsv[3,])
  
  return(col[1:n])
}


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
  set.seed(seed)
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


processCounts = function(counts) {
  if (typeof(counts) != "integer") {
    counts = round(counts)
    storage.mode(counts) = "integer"
  }
  return(counts)  
}


## Generate n random deviates from the Dirichlet function with shape parameters alpha
## Adapted from gtools v3.5
rdirichlet <- function(n, alpha) {
    l <- length(alpha);
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow=TRUE);
    
    ## Check for case where all sampled entries are zero due to round off
    ## One entry will be randomly chosen to be one
    is_zero <- rowSums(x) == 0
    assignment <- sample(1:l, size=sum(is_zero), replace=TRUE)
    x[cbind(which(is_zero), assignment)] <- 1
    
    ## Normalize
    sm <- x %*% rep(1, l);
    y <- x / as.vector(sm);
    return(y)
}


## Make sure provided sample labels are the right type, or generate some if none were provided
processSampleLabels = function(sample.label, num.cells) {
  if(is.null(sample.label)) {
    s = rep(1, num.cells)
    sample.label = s 
  } else if(is.factor(sample.label)) {
    s = as.numeric(sample.label)
  } else {
    sample.label = as.factor(sample.label)
    s = as.integer(sample.label)
  }
  
  return(s)
}
