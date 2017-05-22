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
  stability <- sum(1-apply(probs, 1, entropy::entropy) / log(nStates)) / nData
  return(stability)
}

second.best = function(v) sort(v, decreasing=TRUE)[2]

normalizeLogProbs = function(ll.probs) {
  ll.probs <- exp(sweep(ll.probs, 1, apply(ll.probs, 1, max), "-"))
  probs <- sweep(ll.probs, 1, rowSums(ll.probs), "/")
  return(probs)
}


#' normalizeCounts function 
#' @param counts A count matrix 
#' @param scale.factor the scalor for the normalization 
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
  return(new.z)
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

  return(list(z=new.z, y=new.y))  
}  
