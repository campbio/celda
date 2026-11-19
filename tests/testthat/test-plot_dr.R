# plotDimReduceFeature
library(SingleCellExperiment)

# Fixed values, so no random seed is needed
countsMat <- matrix(c(1, 5, 2, 8, 3, 6, 4, 9, 7, 2,
                      3, 1, 6, 2, 9, 4, 8, 5, 2, 7,
                      6, 3, 1, 7, 2, 8, 5, 4, 9, 3),
                    nrow = 3,
                    byrow = TRUE,
                    dimnames = list(c("A", "B", "C"), paste0("cell", 1:10)))
dim1 <- seq_len(10)
dim2 <- rev(seq_len(10))

# Labels that the plot shows on its facet panels, in panel order
.facetLabels <- function(g) {
  panels <- data.frame(facet = levels(g$data$facet))
  unname(unlist(g$facet$params$labeller(panels)))
}

test_that("plotDimReduceFeature keeps headers on the right panels", {
  expect_warning(
    g <- plotDimReduceFeature(countsMat,
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
  sce <- SingleCellExperiment(assays = list(counts = countsMat))
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

test_that("plotDimReduceFeature uses headers over displayName", {
  sce <- SingleCellExperiment(assays = list(counts = countsMat))
  rowData(sce)$symbol <- c("symA", "symB", "symC")
  altExp(sce, "featureSubset") <- sce
  reducedDim(altExp(sce, "featureSubset"), "dr") <- cbind(dim1, dim2)

  expect_warning(
    g <- plotDimReduceFeature(sce,
                              features = c("symA", "missing", "symC"),
                              reducedDimName = "dr",
                              displayName = "symbol",
                              headers = c("hA", "hMissing", "hC")),
    regexp = "not present in 'x': missing"
  )
  expect_equal(levels(g$data$facet), c("A", "C"))
  expect_equal(.facetLabels(g), c("hA", "hC"))
})

test_that("plotDimReduceFeature checks headers length with displayName", {
  sce <- SingleCellExperiment(assays = list(counts = countsMat))
  rowData(sce)$symbol <- c("symA", "symB", "symC")
  altExp(sce, "featureSubset") <- sce
  reducedDim(altExp(sce, "featureSubset"), "dr") <- cbind(dim1, dim2)

  expect_error(
    plotDimReduceFeature(sce,
                         features = c("symA", "symC"),
                         reducedDimName = "dr",
                         displayName = "symbol",
                         headers = c("hA", "hB", "hC")),
    regexp = "must be the same length as .features. \\(3 and 2\\)"
  )
})

test_that("plotDimReduceFeature falls back to displayName if exactMatch FALSE", {
  sce <- SingleCellExperiment(assays = list(counts = countsMat))
  rowData(sce)$symbol <- c("symA", "symB", "symC")
  altExp(sce, "featureSubset") <- sce
  reducedDim(altExp(sce, "featureSubset"), "dr") <- cbind(dim1, dim2)

  expect_warning(
    g <- plotDimReduceFeature(sce,
                              features = c("symA", "symC"),
                              reducedDimName = "dr",
                              displayName = "symbol",
                              headers = c("hA", "hC"),
                              exactMatch = FALSE),
    regexp = "not used because .exactMatch. is FALSE"
  )
  expect_equal(.facetLabels(g), c("symA", "symC"))
})
