# https://stackoverflow.com/questions/35194048/using-r-how-to-calculate-the-distance-from-one-point-to-a-line
# http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
# Kimberling, C. "Triangle Centers and Central Triangles." Congr.
# Numer. 129, 1-295, 1998.
.dist2d <- function(a, b, c) {
    v1 <- b - c
    v2 <- a - b
    m <- cbind(v1, v2)
    d <- abs(det(m)) / sqrt(sum(v1 * v1))
}

.secondDerivativeEstimate <- function(v) {
    nv <- length(v)
    res <- rep(NA, nv)
    for (i in 2:(nv - 1)) {
        res[i] <- v[i + 1] + v[i - 1] - (2 * v[i])
    }
    return(res)
}

.curveElbow <- function(var, perplexity, pval.cutoff = 0.05) {
    len <- length(perplexity)

    a <- c(var[1], perplexity[1])
    b <- c(var[len], perplexity[len])
    res <- rep(NA, len)
    for (i in seq_along(var)) {
        res[i] <- dist2d(c(var[i], perplexity[i]), a, b)
    }

    elbow <- which.max(res)
    ix <- var > var[elbow]
    perplexity.sde <- secondDerivativeEstimate(perplexity)
    perplexity.sde.sd <- stats::sd(perplexity.sde[ix], na.rm = TRUE)
    perplexity.sde.mean <-
        stats::mean(perplexity.sde[ix], na.rm = TRUE)
    perplexity.sde.pval <-
        stats::pnorm(
            perplexity.sde,
            mean = perplexity.sde.mean,
            sd = perplexity.sde.sd,
            lower.tail = FALSE
        )

    other <- which(ix & perplexity.sde.pval < pval.cutoff)
    return(list(elbow = var[elbow])) # , secondary=l[other]))
}
