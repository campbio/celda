sample.ll = function(ll.probs) {
  probs.sub <- exp(ll.probs - max(ll.probs))
  probs.norm <- probs.sub / sum(probs.sub)
  probs.select <- sample(1:length(ll.probs), size=1, prob=probs.norm)
  return(probs.select)
}


cosineDist = function(x){
  x <- t(x)
  y <- (1 - x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))/2)
  return(as.dist(y))
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

normalizeCounts = function(counts, scale.factor=1e6) {
  counts.norm = sweep(counts, 2, colSums(counts) * scale.factor, "/")
  return(counts.norm)
}
  
