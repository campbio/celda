setClass("celdaModel", 
         representation(completeLogLik = "numeric", 
                        finalLogLik = "numeric",
                        seed = "numeric", 
                        count.checksum = "numeric"))

setGeneric("completeLogLik",
           function(celda.mod){ standardGeneric("completeLogLik") })
setGeneric("completeLogLik<-",
           function(celd.amod, val){ standardGeneric("completeLogLik") })
setGeneric("finalLogLik",
           function(celda.mod, val){ standardGeneric("finalLogLik") })
setGeneric("seed",
           function(celda.mod, val){ standardGeneric("seed") })
setGeneric("countChecksum",
           function(celda.mod, val){ standardGeneric("countChecksum") })

setMethod("completeLogLik",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod, val){ 
             if (!missing(val)) celda.mod@completeLogLik = val
             return(celda.mod@completeLogLik) 
           })
setMethod("finalLogLik",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod){ return(celda.mod@finalLogLik) })
setMethod("seed",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod){ return(celda.mod@seed) })
setMethod("count.checksum",
           signature=c(celda.mod="celdaModel"),
           function(celda.mod){ return(celda.mod@count.checksum) })





setClass("celdaList",
         representation(runParams = "data.frame",
                        res.list = "list",
                        count.checksum = "numeric"))
