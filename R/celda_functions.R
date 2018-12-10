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
  y = stats::dist(t(sqrt(x)), method = "euclidean") * 1/sqrt(2)
  return(y)
}  


normalizeLogProbs = function(ll.probs) {
  ll.probs <- exp(sweep(ll.probs, 1, base::apply(ll.probs, 1, max), "-"))
  probs <- sweep(ll.probs, 1, rowSums(ll.probs), "/")
  return(probs)
}


#' @title Normalization of count data
#' @description Performs normalization, transformation, and/or scaling of a counts matrix
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. 
#' @param normalize Character. Divides counts by the library sizes for each cell. One of 'proportion', 'cpm', 'median', or 'mean'. 'proportion' uses the total counts for each cell as the library size. 'cpm' divides the library size of each cell by one million to produce counts per million. 'median' divides the library size of each cell by the median library size across all cells.  'mean' divides the library size of each cell by the mean library size across all cells.
#' @param transformation.fun Function. Applys a transformation such as `sqrt`, `log`, `log2`, `log10`, or `log1p`. If NULL, no transformation will be applied. Occurs after normalization. Default NULL.
#' @param scale.fun Function. Scales the rows of the normalized and transformed count matrix. For example, 'scale' can be used to z-score normalize the rows. Default NULL.
#' @param pseudocount.normalize Numeric. Add a pseudocount to counts before normalization. Default  0. 
#' @param pseudocount.transform Numeric. Add a pseudocount to normalized counts before applying the transformation function. Adding a pseudocount can be useful before applying a log transformation. Default  0. 
#' @return Numeric Matrix. A normalized matrix.
#' @examples
#' normalized.counts = normalizeCounts(celda.CG.sim$counts, "proportion", 
#'                                     pseudocount.normalize=1)
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
  if(normalize == "proportion") {
    norm = fastNormProp(counts, pseudocount.normalize)
  } else {
    counts = counts + pseudocount.normalize
    cs = .colSums(counts, nrow(counts), ncol(counts))
    norm = switch(normalize,
      "proportion" = sweep(counts, 2, cs, "/"),
      "cpm" = sweep(counts, 2, cs / 1e6, "/"),
      "median" = sweep(counts, 2, cs / stats::median(cs), "/"),
      "mean" = sweep(counts, 2, cs / mean(cs), "/"))
  }  
 
  
  if(!is.null(transformation.fun)){
    norm <- do.call(transformation.fun, list(norm + pseudocount.transform))
  }
  if(!is.null(scale.fun)) {
    norm <- t(base::apply(norm, 1, scale.fun))
  }  

  colnames(norm) = colnames(counts)
  rownames(norm) = rownames(counts)
  return(norm)
}
  

#' @title Recode cell cluster labels 
#' @description Recode cell subpopulaton clusters using a mapping in the `from` and `to` arguments.
#' 
#' @param celda.mod Celda object of class `celda_C` or `celda_CG`. 
#' @param from Numeric vector. Unique values in the range of 1:K that correspond to the original cluster labels in `celda.mod`.  
#' @param to Numeric vector. Unique values in the range of 1:K that correspond to the new cluster labels. 
#' @return Celda object with cell subpopulation clusters, with class corresponding to that of `celda.mod`.
#' @examples
#' celda.mod.reordered.z = recodeClusterZ(celda.CG.mod, c(1, 3), c(3, 1))
#' @export
#' @export
recodeClusterZ = function(celda.mod, from, to) {
  if (length(setdiff(from, to)) != 0) {
    stop("All values in 'from' must have a mapping in 'to'")
  }
  if (is.null(celda.mod@clusters$z)) {
    stop("Provided celda.mod argument does not have a z attribute")
  }
  celda.mod@clusters$z = plyr::mapvalues(celda.mod@clusters$z, from, to)
  return(celda.mod)
}
  

#' @title Recode feature module clusters
#' @description Recode feature module clusters using a mapping in the `from` and `to` arguments.
#' 
#' @param celda.mod Celda object of class `celda_G` or `celda_CG`. 
#' @param from Numeric vector. Unique values in the range of 1:L that correspond to the original cluster labels in `celda.mod`. 
#' @param to Numeric vector. Unique values in the range of 1:L that correspond to the new cluster labels. 
#' @return Celda object with recoded feature module clusters, with class corresponding to that of `celda.mod`.
#' @examples
#' celda.mod.reordered.y = recodeClusterY(celda.CG.mod, c(1, 3), c(3, 1))
#' @export
recodeClusterY = function(celda.mod, from, to) {
  if (length(setdiff(from, to)) != 0) {
    stop("All values in 'from' must have a mapping in 'to'")
  }
  if (is.null(celda.mod@clusters$y)) {
    stop("Provided celda.mod argument does not have a y attribute")
  }
  celda.mod@clusters$y = plyr::mapvalues(celda.mod@clusters$y, from, to)
  return(celda.mod)
}


#' @title Check count matrix consistency
#' @description Checks if the counts matrix is the same one used to generate the celda model object by comparing dimensions and MD5 checksum.
#' 
#' @param counts Integer matrix. Rows represent features and columns represent cells. 
#' @param celda.mod Celda model object. 
#' @param error.on.mismatch Logical. Whether to throw an error in the event of a mismatch. Default TRUE.
#' @return Returns TRUE if provided count matrix matches the one used in the celda object and/or `error.on.mismatch=FALSE`, FALSE otherwise.
#' @examples
#'  compareCountMatrix(celda.CG.sim$counts, celda.CG.mod, 
#'                     error.on.mismatch=FALSE)
#' @export
compareCountMatrix = function(counts, celda.mod, error.on.mismatch=TRUE) {
  if (methods::.hasSlot(celda.mod, "y")) {
    if (nrow(counts) != length(celda.mod@clusters$y)) {
      stop(paste0("The provided celda object was generated from a counts matrix with a different number of features than the one provided."))
    }  
  }
  
  if (methods::.hasSlot(celda.mod, "z")) {  
    if (ncol(counts) != length(celda.mod@clusters$z)) {
      stop(paste0("The provided celda object was generated from a counts matrix with a different number of cells than the one provided."))
    }
  }
  
  celda.checksum = celda.mod@params$count.checksum
  counts = processCounts(counts)  # Checksums are generated in celdaGridSearch and model functions after processing
  count.md5 = digest::digest(counts, algo="md5")
  res = isTRUE(count.md5 == celda.checksum)
  if (res) return(TRUE)
  if (!res && error.on.mismatch) stop("There was a mismatch between the provided count matrix and the count matrix used to generate the provided celda result.")
  else if (!res && !error.on.mismatch) return(FALSE)
}



logMessages = function(..., sep = " ", logfile = NULL, append = FALSE, verbose = TRUE) {
  if(isTRUE(verbose)) {
	if(!is.null(logfile)) {
	  if(!is.character(logfile) || length(logfile) > 1) {
		stop("The log file parameter needs to be a single character string.")
	  }    
	  cat(paste(..., "\n", sep=sep), file=logfile, append=append)
	} else {
	  message(paste(..., sep=sep))
	}
  }	
}

#' @title Create a color palette
#' @description Generate a palette of `n` distinct colors. 
#' 
#' @param n Integer. Number of colors to generate. 
#' @param hues Character vector. Colors available from `colors()`. These will be used as the base colors for the clustering scheme in HSV. Different saturations and values will be generated for each hue. Default c("red", "cyan", "orange", "blue", "yellow", "purple", "green", "magenta").
#' @param saturation.range Numeric vector. A vector of length 2 denoting the saturation for HSV. Values must be in [0,1]. Default: c(0.25, 1).
#' @param value.range Numeric vector. A vector of length 2 denoting the range of values for HSV. Values must be in [0,1]. Default: `c(0.5, 1)`.
#' @return A vector of distinct colors that have been converted to  HEX from HSV.
#' @examples
#' color.pal = distinct_colors(6)  # can be used in plotting functions
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
  list <- lapply(1:num.vs, function(x){
    rbind(hues.hsv[1,], s[x], v[x])
  })
  new.hsv <- do.call(cbind,list)
  
  ## Convert to hex
  col = grDevices::hsv(new.hsv[1,], new.hsv[2,], new.hsv[3,])
  
  return(col[1:n])
}


processCounts = function(counts) {
  if (typeof(counts) != "integer") {
    counts = round(counts)
    storage.mode(counts) = "integer"
  }
  return(counts)  
}


# Perform some simple checks on the counts matrix, to ensure celda modeling
# expectations are met
validateCounts = function(counts) {
  # And each row/column of the count matrix must have at least one count
  count.row.sum = rowSums(counts)
  count.col.sum = colSums(counts)
  
  if (sum(count.row.sum == 0) > 1 | sum(count.col.sum == 0) > 1) {
    stop("Each row and column of the count matrix must have at least one count")
  }
}


## Generate n random deviates from the Dirichlet function with shape parameters alpha
## Adapted from gtools v3.5
rdirichlet <- function(n, alpha) {
    l <- length(alpha);
    x <- matrix(stats::rgamma(l * n, alpha), ncol = l, byrow=TRUE);
    
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
    sample.label = as.factor(rep("Sample_1", num.cells))
  } else {
    if(length(sample.label) != num.cells) {
      stop("'sample.label' must be the same length as the number of columns in the 'counts' matrix.")
    }
  }
  if(!is.factor(sample.label)) {
    sample.label = as.factor(sample.label)
  }
  return(sample.label)
}
