sample.ll = function(ll.probs) {
  probs.sub = exp(ll.probs - max(ll.probs))
  probs.norm = probs.sub / sum(probs.sub)
  probs.select = sample(1:length(ll.probs), size=1, prob=probs.norm)
  return(probs.select)
}


cosineDist <- function(x){
  x = t(x)
  y = as.dist(1 - x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))) 
  return(y)
}

second.best = function(v) sort(v, decreasing=TRUE)[2]

stability = function(ll.probs) {
  require(entropy)
  nStates = ncol(ll.probs)
  nData = nrow(ll.probs)
  
  ## Normalize probabilties 
  ll.probs = exp(sweep(probs, 1, apply(ll.probs, 1, max), "-"))
  probs = sweep(ll.probs, 1, rowSums(ll.probs), "/")
  
  stability = sum(1-apply(probs, 1, entropy) / log(nStates)) / nData
  return(stability)
}
  