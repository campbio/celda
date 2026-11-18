# plotDimReduceFeature
library(SingleCellExperiment)

set.seed(12345)
counts <- matrix(stats::rpois(30, lambda = 5),
                 nrow = 3,
                 dimnames = list(c("A", "B", "C"), paste0("cell", 1:10)))
dim1 <- stats::rnorm(10)
dim2 <- stats::rnorm(10)

# Labels that the plot shows on its facet panels, in panel order
.facetLabels <- function(g) {
  panels <- data.frame(facet = levels(g$data$facet))
  unname(unlist(g$facet$params$labeller(panels)))
}

test_that("plotDimReduceFeature keeps headers on the right panels", {
  expect_warning(
    g <- plotDimReduceFeature(counts,
                              features = c("A", "missing", "C"),
                              dim1 = dim1,
                              dim2 = dim2,
                              headers = c("hA", "hMissing", "hC")),
    regexp = "not present in 'x': missing"
  )
  expect_equal(levels(g$data$facet), c("A", "C"))
  expect_equal(.facetLabels(g), c("hA", "hC"))
})

test_that("plotDimReduceFeature keeps displayName labels on the right panels", {
  sce <- SingleCellExperiment(assays = list(counts = counts))
  rowData(sce)$symbol <- c("symA", "symB", "symC")
  altExp(sce, "featureSubset") <- sce
  reducedDim(altExp(sce, "featureSubset"), "dr") <- cbind(dim1, dim2)

  expect_warning(
    g <- plotDimReduceFeature(sce,
                              features = c("symA", "missing", "symC"),
                              reducedDimName = "dr",
                              displayName = "symbol"),
    regexp = "not present in 'x': missing"
  )
  expect_equal(levels(g$data$facet), c("A", "C"))
  expect_equal(.facetLabels(g), c("symA", "symC"))
})
