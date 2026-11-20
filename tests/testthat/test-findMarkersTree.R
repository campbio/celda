# findMarkersTree
library(celda)

# Build a small features-by-cells matrix with `nClass` classes of 20 cells.
# f1 marks class 1 and f2 marks class 2. f3 is high in classes 1-3 and low
# in the others, so it separates class 3 from class 4 (and above) but is not
# a one-off marker. Classes 4 and above differ only by noise.
# With `f3Class3Only = TRUE`, f3 is high only in class 3, apart from three
# class-3 cells, so its one-off statistic is about 0.94.
.markerTreeFixture <- function(nClass = 4, f3Class3Only = FALSE) {
  withr::with_seed(410, {
    class <- rep(seq_len(nClass), each = 20)
    n <- length(class)
    f1 <- ifelse(class == 1, runif(n, 5, 6), runif(n, 0, 1))
    f2 <- ifelse(class == 2, runif(n, 5, 6), runif(n, 0, 1))
    if (f3Class3Only) {
      f3 <- ifelse(class == 3, runif(n, 5, 6), runif(n, 0, 1))
      f3[which(class == 3)[seq_len(3)]] <- runif(3, 0, 1)
    } else {
      f3 <- ifelse(class %in% c(1, 2, 3), runif(n, 5, 6), runif(n, 0, 1))
    }
  })
  features <- rbind(f1 = f1, f2 = f2, f3 = f3)
  colnames(features) <- paste0("Cell_", seq_len(n))
  list(features = features, class = class)
}

# These tests supply `metaclusters` so that the top-level threshold loop is
# reached without the UMAP step of `autoMetaclusters = TRUE`, which is the
# path from issue 410. Both paths share the loop.

test_that("findMarkersTree merges metaclusters left without markers", {
  fixture <- .markerTreeFixture()
  metaclusters <- list(M1 = 1, M2 = 2, M3 = 3, M4 = 4)
  expect_warning(
    tree <- findMarkersTree(
      fixture$features,
      fixture$class,
      metaclusters = metaclusters
    ),
    regexp = "M3\\(3\\), M4\\(4\\).*no positive marker"
  )
  merged <- "M3+M4"
  expect_setequal(unique(tree$metaclusterLabels), c("M1", "M2", merged))
  expect_named(
    tree$rules,
    c("M1(1)", "M2(2)", paste0(merged, c("(3)", "(4)"))),
    ignore.order = TRUE
  )
  # The merged clusters are separated in their own subtree by f3
  subtreeRules <- tree$rules[[paste0(merged, "(3)")]]
  f3Up <- subtreeRules$feature == "f3" & subtreeRules$direction == 1
  expect_true(any(f3Up & subtreeRules$metacluster == merged))
  expect_s3_class(plotDendro(tree), "ggplot")
})

test_that("findMarkersTree merges a multi-cluster metacluster", {
  fixture <- .markerTreeFixture(nClass = 5)
  metaclusters <- list(M1 = 1, M2 = 2, M3 = 3, M4 = c(4, 5))
  expect_warning(
    tree <- findMarkersTree(
      fixture$features,
      fixture$class,
      metaclusters = metaclusters
    ),
    regexp = "M3\\(3\\), M4 have no positive marker"
  )
  merged <- "M3+M4"
  expect_setequal(unique(tree$metaclusterLabels), c("M1", "M2", merged))
  expect_named(
    tree$rules,
    c("M1(1)", "M2(2)", paste0(merged, c("(3)", "(4)", "(5)"))),
    ignore.order = TRUE
  )
})

test_that("findMarkersTree keeps a single metacluster without markers", {
  fixture <- .markerTreeFixture(nClass = 3)
  metaclusters <- list(M1 = 1, M2 = 2, M3 = 3)
  expect_warning(
    tree <- findMarkersTree(
      fixture$features,
      fixture$class,
      metaclusters = metaclusters
    ),
    regexp = "Metacluster M3\\(3\\) has no positive marker"
  )
  expect_named(tree$rules, c("M1(1)", "M2(2)", "M3(3)"), ignore.order = TRUE)
  expect_false(any(tree$rules[["M3(3)"]]$direction == 1))
  expect_s3_class(plotDendro(tree), "ggplot")
})

test_that("findMarkersTree skips ambiguous alternative splits", {
  fixture <- .markerTreeFixture()
  expect_message(
    tree <- findMarkersTree(
      fixture$features,
      fixture$class,
      autoMetaclusters = FALSE
    ),
    regexp = "No unambiguous rule to separate 3"
  )
  expect_named(tree$rules, as.character(seq_len(4)), ignore.order = TRUE)
})

test_that("findMarkersTree adds an alternative split to its own class", {
  fixture <- .markerTreeFixture(f3Class3Only = TRUE)
  expect_message(
    tree <- findMarkersTree(
      fixture$features,
      fixture$class,
      threshold = 0.95,
      autoMetaclusters = FALSE
    ),
    regexp = "No unused features to separate 4"
  )
  rules3 <- tree$rules[["3"]]
  rules4 <- tree$rules[["4"]]
  expect_true(any(rules3$feature == "f3" & rules3$direction == 1))
  expect_false(any(rules4$direction == 1))
})
