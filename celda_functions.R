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
