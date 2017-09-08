sample.ll = function(ll.probs) {
  probs.sub <- exp(ll.probs - max(ll.probs))
  probs.norm <- probs.sub / sum(probs.sub)
  probs.select <- sample(1:length(ll.probs), size=1, prob=probs.norm)
  return(probs.select)
}


cosineDist = function(x){
  x <- t(x)
  y <- (1 - cosine(x)) / 2
  return(as.dist(y))
}

cosine = function(x) {
  y = x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))
  return(y)
}
  
spearmanDist = function(x){
  y = (1 - cor(x, method="spearman"))/2
  return(as.dist(y))
}


stability = function(probs) {
  nStates <- ncol(probs)
  nData <- nrow(probs)
  stability <- sum(1-base::apply(probs, 1, entropy::entropy) / log(nStates)) / nData
  return(stability)
}

second.best = function(v) sort(v, decreasing=TRUE)[2]

normalizeLogProbs = function(ll.probs) {
  ll.probs <- exp(sweep(ll.probs, 1, base::apply(ll.probs, 1, max), "-"))
  probs <- sweep(ll.probs, 1, rowSums(ll.probs), "/")
  return(probs)
}


#' Normalize a counts matrix by a scalar factor
#' 
#' @param counts A count matrix 
#' @param scale.factor the scalar for the normalization 
#' @export
normalizeCounts = function(counts, scale.factor=1e6) {
  counts.norm = sweep(counts, 2, colSums(counts) / scale.factor, "/")
  return(counts.norm)
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


#' Re-code cluster label by provided mapping scheme
#' 
#' This function will re-code cluster labels based off of a mapping provided by the user,
#' for all fields on a celda object involving cluster labels.
#' e.g. if Z values range from 1-4 and the user would like all 3's switched to 1's and
#' vice versa, this function can be useful. NOTE: it is recommended that this function's results
#' aren't used to overwrite the original celda model object provided, in the event of a mis-mapping.
#' 
#' @param celda.mod An object of class celda_C, celda_G, or celda_CG
#' @param from A numeric vector of the "keys", corresponding to the "values" in the to parameter
#' @param to A numeric vector of the "values"; what each corresponding number in "from" should be mapped to
#' @export
recode_cluster_labels = function(celda.mod, from, to) {
  # Make sure that every value gets a replacement
  if (length(setdiff(from, to)) != 0) {
    stop("All values in 'from' must have a mapping in 'to'")
  }
  
  if (class(celda.mod) %in% c("celda_C", "celda_CG")) {
    celda.mod$z = plyr::mapvalues(celda.mod$z, from, to)
    celda.mod$z.probability = celda.mod$z.probability[, to]
  }
  if (class(celda.mod) %in% c("celda_G", "celda_CG")) {
    celda.mod$y = plyr::mapvalues(celda.mod$y, from, to)
    celda.mod$y.probability = celda.mod$y.probability[, to]
  }
  return(celda.mod)
}


#' Check whether a count matrix was the one used in a given celda run
#' 
#' Compare the MD5 checksum of a provided count.matrix to the count matrix
#' checksum on a celda_list object, to see if they're the same.
#' @param count.matrix A numeric matrix of counts
#' @param celda.checksum An MD5 checksum from a celda_list or celda model object (celda_C, celda_G, celda_CG)
#' @return TRUE if provided count matrix matches the one used in the celda run, FALSE otherwise
#' @export
compare_count_matrix = function(count.matrix, celda.checksum) {
  count.md5 = digest::digest(count.matrix, algo="md5")
  return(count.md5 == celda.checksum)
}


log_messages = function(..., sep = " ", logfile = NULL, append = FALSE) {
  if(!is.null(logfile)) {
    if(!is.character(logfile) || length(logfile) > 1) {
      stop("The log file parameter needs to be a single character string.")
    }    
    cat(paste(..., "\n", sep=sep), file=logfile, append=append)
    
  } else {
    message(paste(..., sep=sep))
  }
}





#' Generate a distinct palette for coloring clusters
#' 
#' @param n Integer; Number of colors to generate
#' @param saturation.range Numeric vector of length 2 with values between 0 and 1. Default: c(0.25, 1)
#' @param value.range Numeric vector of length 2 with values between 0 and 1. Default: c(0.5, 1)
#' @return A vector of distinct colors in HEX
#' @export
distinct_colors = function(n,
						   hues = c("red", "cyan", "orange", "blue", "yellow", "purple", "green", "magenta"),
                           saturation.range = c(0.4, 1),
                           value.range = c(0.4, 1)) {
                           
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



