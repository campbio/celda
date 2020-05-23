#' @title Calculate the Log-likelihood of a celda model
#' @description Calculate the log-likelihood for cell population
#'  and feature module cluster assignments on the count matrix, per celda model.
#' @param x A \linkS4class{SingleCellExperiment} object returned by
#'  \link{celda_C}, \link{celda_G}, or \link{celda_CG}, with the matrix
#'  located in the \code{useAssay} assay slot.
#'  Rows represent features and columns represent cells.
#' @param useAssay A string specifying which \link[SummarizedExperiment]{assay}
#'  slot to use. Default "counts".
#' @param celdaMod celda model object. Ignored if \code{x} is a
#'  \linkS4class{SingleCellExperiment} object.
#' @return The log-likelihood of the cluster assignment for the
#'  provided \linkS4class{SingleCellExperiment}.
#' @seealso `celda_C()` for clustering cells
#' @export
setGeneric("logLikelihood", function(x, celdaMod, ...) {
    standardGeneric("logLikelihood")
})


#' @rdname logLikelihood
#' @examples
#' data(sceCeldaC, sceCeldaCG)
#' loglikC <- logLikelihood(sceCeldaC)
#' loglikCG <- logLikelihood(sceCeldaCG)
#' @export
setMethod("logLikelihood", signature(x = "SingleCellExperiment"),
    function(x, useAssay = "counts") {

        counts <- SummarizedExperiment::assay(x, i = useAssay)
        sampleLabel <- sampleLabel(x)
        z <- celdaClusters(x)
        y <- celdaModules(x)
        K <- S4Vectors::metadata(x)$celda_parameters$K
        L <- S4Vectors::metadata(x)$celda_parameters$L
        alpha <- S4Vectors::metadata(x)$celda_parameters$alpha
        beta <- S4Vectors::metadata(x)$celda_parameters$beta
        delta = S4Vectors::metadata(x)$celda_parameters$delta
        gamma = S4Vectors::metadata(x)$celda_parameters$gamma

        if (celdaModel(x) == "celda_C") {
            ll <- .logLikelihoodcelda_C(counts = counts,
                sampleLabel = sampleLabel,
                z = z,
                K = K,
                alpha = alpha,
                beta = beta)
        } else if (celdaModel(x) == "celda_CG") {
            ll <- .logLikelihoodcelda_CG(counts = counts,
                sampleLabel = sampleLabel,
                z = z,
                y = y,
                K = K,
                L = L,
                alpha = alpha,
                beta = beta,
                delta = delta,
                gamma = gamma)
        } else if (celdaModel(x) == "celda_G") {
            ll <- .logLikelihoodcelda_G(counts = counts,
                y = y,
                L = L,
                beta = beta,
                delta = delta,
                gamma = gamma)
        } else {
            stop("S4Vectors::metadata(x)$celda_parameters$model must be",
                " one of 'celda_C', 'celda_G', or 'celda_CG'!")
        }
        return(ll)
    }
)


#' @rdname logLikelihood
#' @export
setMethod("logLikelihood", signature(x = "matrix", celdaMod = "celda_C"),
    function(x, celdaMod) {
        sampleLabel <- sampleLabel(celdaMod)
        z = celdaClusters(celdaMod)$z
        K = params(celdaMod)$K
        alpha = params(celdaMod)$alpha
        beta = params(celdaMod)$beta

        ll <- .logLikelihoodcelda_C(counts = x,
            sampleLabel = sampleLabel,
            z = z,
            K = K,
            alpha = alpha,
            beta = beta)
        return(ll)
    }
)


#' @rdname logLikelihood
#' @export
setMethod("logLikelihood", signature(x = "matrix", celdaMod = "celda_G"),
    function(x, celdaMod) {
        y <- celdaClusters(celdaMod)$y
        L <- params(celdaMod)$L
        beta = params(celdaMod)$beta
        delta = params(celdaMod)$delta
        gamma = params(celdaMod)$gamma

        ll <- .logLikelihoodcelda_G(counts = x,
            y = y,
            L = L,
            beta = beta,
            delta = delta,
            gamma = gamma)
        return(ll)
    }
)


#' @rdname logLikelihood
#' @export
setMethod("logLikelihood", signature(x = "matrix", celdaMod = "celda_CG"),
    function(x, celdaMod) {
        sampleLabel <- sampleLabel(celdaMod)
        z <- celdaClusters(celdaMod)$z
        y <- celdaClusters(celdaMod)$y
        K <- params(celdaMod)$K
        L <- params(celdaMod)$L
        alpha = params(celdaMod)$alpha
        beta = params(celdaMod)$beta
        delta = params(celdaMod)$delta
        gamma = params(celdaMod)$gamma

        ll <- .logLikelihoodcelda_CG(counts = x,
            sampleLabel = sampleLabel,
            z = z,
            y = y,
            K = K,
            L = L,
            alpha = alpha,
            beta = beta,
            delta = delta,
            gamma = gamma)
        return(ll)
    }
)


#' @title Get log-likelihood history
#' @description Retrieves the complete log-likelihood from all iterations of
#'  Gibbs sampling used to generate a celda model.
#' @param x A \linkS4class{SingleCellExperiment} object
#'  returned by \link{celda_C}, \link{celda_G}, or \link{celda_CG}, or a celda
#'  model object.
#' @return Numeric. The log-likelihood at each step of Gibbs sampling used to
#'  generate the model.
#' @export
setGeneric(
    "logLikelihoodHistory",
    function(x) {
        standardGeneric("logLikelihoodHistory")
    }
)


#' @rdname logLikelihoodHistory
#' @examples
#' data(sceCeldaCG)
#' logLikelihoodHistory(sceCeldaCG)
#' @export
setMethod("logLikelihoodHistory",
    signature(x = "SingleCellExperiment"),
    function(x) {
        cll <- S4Vectors::metadata(x)$celda_parameters$completeLogLik
        return(cll)
    }
)


#' @rdname logLikelihoodHistory
#' @examples
#' data(celdaCGMod)
#' logLikelihoodHistory(celdaCGMod)
#' @export
setMethod("logLikelihoodHistory",
    signature(x = "celdaModel"),
    function(x) {
        cll <- x@completeLogLik
        return(cll)
    }
)


#' @title Get the log-likelihood
#' @description Retrieves the final log-likelihood from all iterations of Gibbs
#'  sampling used to generate a celdaModel.
#' @return Numeric. The log-likelihood at the final step of Gibbs sampling used
#'  to generate the model.
#' @param x A \linkS4class{SingleCellExperiment} object
#'  returned by \link{celda_C}, \link{celda_G}, or \link{celda_CG}, or a celda
#'  model object.
#' @examples
#' data(celdaCGMod)
#' bestLogLikelihood(celdaCGMod)
#' @export
setGeneric(
    "bestLogLikelihood",
    function(x) {
        standardGeneric("bestLogLikelihood")
    }
)


#' @rdname bestLogLikelihood
#' @examples
#' data(sceCeldaCG)
#' bestLogLikelihood(sceCeldaCG)
#' @export
setMethod("bestLogLikelihood",
    signature(x = "celdaModel"),
    function(x) {
        fll <- S4Vectors::metadata(x)$celda_parameters$finalLogLik
        return(fll)
    }
)


#' @rdname bestLogLikelihood
#' @examples
#' data(celdaCGMod)
#' bestLogLikelihood(celdaCGMod)
#' @export
setMethod("bestLogLikelihood",
    signature(x = "celdaModel"),
    function(x) {
        fll <- x@finalLogLik
        return(fll)
    }
)
