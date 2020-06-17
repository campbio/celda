setClass("celdaModel",
    slots = c(params = "list",
        # K, L, model priors, checksum
        names = "list",
        completeLogLik = "numeric",
        finalLogLik = "numeric",
        clusters = "list")
) # z and or y

setClass("celda_C",
    representation(sampleLabel = "factor"),
    contains = "celdaModel")

setClass("celda_G", contains = "celdaModel")

setClass("celda_CG", contains = c("celda_C", "celda_G"))

setClass("celdaList",
    slots = c(runParams = "data.frame",
        resList = "list",
        countChecksum = "character",
        perplexity = "matrix",
        celdaGridSearchParameters = "list")
)
