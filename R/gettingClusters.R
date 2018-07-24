#' Assists in selection of K/L parameter for downstream analysis.
#' 
#' @param celda.list celda.list object. Ouput of celda function.
#' @param matrix matrix. Counts matrix that was inputted to celda function. Do not normalize prior to using function.
#' @param iterations Numeric. Number of iterations to run the function. A higher number will generally output a smoother plot
#' @export
gettingClusters <- function(celda.list, matrix, iterations){
  #matrix <- log(normalizeCounts(matrix) + 1)
  all.max.value <- matrix(ncol = length(unique(celda.list$run.params$L)) * length(unique(celda.list$run.params$K)))
  for(times in 1:iterations){
    
    #max.value will contain the "worst" p-value upon pairwise MAST.
    max.value <- double(length = length(unique(celda.list$run.params$L)) * length(unique(celda.list$run.params$K)))
    list.index = 1
    
    #For all transcriptional states
    for(states in unique(celda.list$run.params$L)){
      #For all cell states
      
      for(x in unique(celda.list$run.params$K)){
        #Best model via log-likelihood has been shown to get closest to "true" labels
        model <- getBestModel(celda.list,K = x,L = states)
        
        #Allocate memory for all.min, which will contain all p-values for pairwise MAST comparisons
        #The maximum value of all.min after running all pairwise comparisons will be stored into max.value
        all.min <- double(length = (length(model$K) * (length(model$K) - 1))/2)
        all.min.index = 1
        
        #Make factorized matrix off of model, use population states
        factorize.matrix <- factorizeMatrix(counts = matrix, celda.mod = model)
        fm <- factorize.matrix$proportions$population.states
        #print(paste("K:",x))
        #print(paste("L:",states))
        
        
        p.values = c()
        pairs <- c()
        distance.matrix <- dist(t(fm))
        for(first in 1:(model$K-1)){
          for(second in (first+1):(model$K)){
            pairs <- c(pairs,(paste0(first,"-",second)))}}
        top.distances <- head(order(distance.matrix, decreasing = FALSE),10)
        
        #MAST will do pairwise comparison on cell subpopulation combinations that are shown to 
        #have be the closest neighbors.
        #first.k contains the first subpopulation, second.k contains the second subpopulation.
        first.k <- as.numeric(gsub(pattern = "-.*$",replacement = "",pairs[top.distances]))
        second.k <- as.numeric(gsub(pattern = ".*-",replacement = "",pairs[top.distances]))
        
        #Take all p-values from MAST
        for(z in 1:length(first.k)){
          k = first.k[z]
          y = second.k[z]
          
          cluster1 <- rmultinom(n = 250,size = 250,prob = fm[,k])
          cluster2 <- rmultinom(n = 250,size = 250,prob = fm[,y])
          
          colnames(cluster1) <- paste0("A",c(1:250))
          colnames(cluster2) <- paste0("B",c(1:250))
          mat <- cbind(cluster1,cluster2)
          
          cells1 <- colnames(cluster1)
          cells2 <- colnames(cluster2)
          log_normalized_mat <- log2(normalizeCounts(mat) + 1)
          cdat <-
            data.frame(
              wellKey = c(cells1, cells2),
              condition = c(rep("c1", length(cells1)), rep("c2", length(cells2))),
              ngeneson = rep("", (length(cells1) + length(cells2))),
              stringsAsFactors = FALSE
            )
          
          sca <- MAST::FromMatrix(log_normalized_mat, cdat)
          ##If cdr2 is the exact same for both cell types zlm will fail
          cdr2 <- colSums(SummarizedExperiment::assay(sca) > 0)
          SummarizedExperiment::colData(sca)$cngeneson <- scale(cdr2)
          cond <- factor(SummarizedExperiment::colData(sca)$condition)
          cond <- relevel(cond, "c2")
          SummarizedExperiment::colData(sca)$condition <- cond
          if(all(!is.nan(scale(cdr2)))){
            zlmCond <- suppressMessages(MAST::zlm( ~ condition + cngeneson, sca))
          }else{
            zlmCond <- suppressMessages(MAST::zlm( ~ condition, sca))
          }
          summaryCond <- suppressMessages(MAST::summary(zlmCond, doLRT = 'conditionc1')) 
          summaryDt <- summaryCond$datatable
          
          fcHurdle <-
            merge(summaryDt[contrast == 'conditionc1' &
                              component == 'H', .(primerid, `Pr(>Chisq)`)],
                  summaryDt[contrast == 'conditionc1' &
                              component == 'logFC', .(primerid, coef, ci.hi, ci.lo)], by = 'primerid')
          
          fcHurdle[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
          fcHurdle <- fcHurdle[,-c(4,5)]
          names(fcHurdle)[c(1,3)] <- c("Gene", "log2fc")
          fcHurdle <- fcHurdle[order(abs(fcHurdle$log2fc), decreasing = TRUE),]
          
          min.value <- (min(fcHurdle[,2]))
          
          all.min[all.min.index] <- min.value
          all.min.index <- all.min.index + 1
          
        }
        max.value[list.index] <- max(all.min)
        list.index = list.index + 1
      }
    }
    all.max.value <- rbind(all.max.value, max.value)
  }
  all.max.value <- all.max.value[-1,]
  
  #Create for ggplot dataframe
  Ks = rep(unique(celda.list$run.params$K), times = length(unique(celda.list$run.params$L)))
  L = rep(unique(celda.list$run.params$L), each = length(unique(celda.list$run.params$K)))
  
  if(iterations > 1){
    logmaxp <- -log(colMeans(as.matrix(all.max.value)))
  }else{
    logmaxp <- -log(all.max.value)
  }
  
  #data frame for ggplot
  df <- data.frame(Ks, logmaxp, L)
  
  colnames(df)[2] <- "negative_log"
  
  plot <- ggplot2::ggplot(data = df, ggplot2::aes(x = Ks, y = negative_log, group = L)) + 
    ggplot2::geom_line(ggplot2::aes(color = as.factor(L))) + 
    ggplot2::geom_point(ggplot2::aes(color = as.factor(L)))
  return(plot)
  
}