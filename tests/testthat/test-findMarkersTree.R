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
    regexp = "Metaclusters M3, M4 have no positive marker"
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
  # The root of the dendrogram holds every cluster, merged ones included
  root <- attributes(tree$dendro)
  expect_setequal(
    root$classLabels,
    c("M1(1)", "M2(2)", paste0(merged, c("(3)", "(4)")))
  )
  expect_equal(root$members, 4)
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
    regexp = "Metaclusters M3, M4 have no positive marker"
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
    regexp = "Metacluster M3 has no positive marker"
  )
  expect_named(tree$rules, c("M1(1)", "M2(2)", "M3(3)"), ignore.order = TRUE)
  expect_false(any(tree$rules[["M3(3)"]]$direction == 1))
  expect_s3_class(plotDendro(tree), "ggplot")
})

# Classes 1, 2 and 3 have 20 cells each. f1 is high in classes 1 and 2, so
# it marks them together but neither alone, and f3 marks class 3. f2 is only
# slightly higher in class 2, too weak to be a top-level marker, but enough
# to separate classes 1 and 2 from each other.
.balancedSplitFixture <- function() {
  withr::with_seed(441, {
    class <- rep(seq_len(3), each = 20)
    n <- length(class)
    f1 <- ifelse(class %in% c(1, 2), runif(n, 5, 6), runif(n, 0, 1))
    f2 <- ifelse(class == 2, runif(n, 0.6, 1.6), runif(n, 0, 1))
    f3 <- ifelse(class == 3, runif(n, 5, 6), runif(n, 0, 1))
  })
  features <- rbind(f1 = f1, f2 = f2, f3 = f3)
  colnames(features) <- paste0("Cell_", seq_len(n))
  list(features = features, class = class)
}

test_that("findMarkersTree merges metaclusters in a top-level balanced split", {
  fixture <- .balancedSplitFixture()
  metaclusters <- list(M1 = 1, M2 = 2, M3 = 3)
  expect_no_warning(
    tree <- suppressMessages(findMarkersTree(
      fixture$features,
      fixture$class,
      metaclusters = metaclusters
    ))
  )
  merged <- "M1+M2"
  expect_setequal(unique(tree$metaclusterLabels), c(merged, "M3"))
  expect_setequal(unique(tree$classLabels), c("1", "2", "3"))
  expect_named(
    tree$rules,
    c(paste0(merged, c("(1)", "(2)")), "M3(3)"),
    ignore.order = TRUE
  )
  expect_true(paste0(merged, "_level_1") %in% names(tree$branchPoints))
  # The merged clusters are separated in their own subtree by f2
  subtreeRules <- tree$rules[[paste0(merged, "(2)")]]
  f2Up <- subtreeRules$feature == "f2" & subtreeRules$direction == 1
  expect_true(any(f2Up & subtreeRules$level == 2))
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

# Classes 1 and 2 have 40 cells, class 3 has 4. f4 marks class 1 and f2
# marks class 2. f1 is high in class 1 and in two class-3 cells, so its
# best one-off split (statistic about 0.74) marks both class 1 and class 3.
.multiClassFixture <- function() {
  withr::with_seed(440, {
    class <- rep(c(1, 2, 3), c(40, 40, 4))
    n <- length(class)
    f4 <- ifelse(class == 1, runif(n, 5, 6), runif(n, 0, 1))
    f2 <- ifelse(class == 2, runif(n, 5, 6), runif(n, 0, 1))
    f1 <- ifelse(class == 1, runif(n, 5, 6), runif(n, 0, 1))
    f1[which(class == 3)[1:2]] <- c(5.5, 5.6)
  })
  features <- rbind(f1 = f1, f2 = f2, f4 = f4)
  colnames(features) <- paste0("Cell_", seq_len(n))
  list(features = features, class = class)
}

test_that("one-off splits each mark exactly one class", {
  fixture <- .multiClassFixture()
  split <- celda:::.wrapSplitHybrid(
    t(fixture$features),
    as.factor(fixture$class),
    threshold = 0.7,
    oneoffMetric = "modified F1"
  )
  # The fixture's premise: f1 passes the threshold, and its split marks two
  # classes
  f1Alone <- celda:::.getSplit(
    "f1",
    c(f1 = 1),
    t(fixture$features),
    as.factor(fixture$class),
    celda:::.splitMetricModF1
  )
  expect_identical(f1Alone$group1Consensus, c("1", "3"))
  f1Stat <- celda:::.splitMetricModF1(
    "f1",
    as.factor(fixture$class),
    t(fixture$features),
    rPerf = TRUE
  )
  expect_gte(f1Stat, 0.7)

  expect_identical(split$statUsed, "One-off")
  nodes <- split[names(split) != "statUsed"]
  nodeFeatures <- unlist(lapply(nodes, function(node) node$featureName))
  expect_false("f1" %in% nodeFeatures)
  expect_true(all(nzchar(names(nodes))))
  group1 <- lapply(nodes, function(node) node$group1Consensus)
  expect_true(all(lengths(group1) == 1))
  expect_false(anyDuplicated(unlist(group1)) > 0)
  # Every feature of a node, split on its own, marks that node's class
  for (node in nodes) {
    for (feat in node$featureName) {
      alone <- celda:::.getSplit(
        feat,
        node$stat[node$featureName == feat],
        t(fixture$features),
        as.factor(fixture$class),
        celda:::.splitMetricModF1
      )
      expect_identical(alone$group1Consensus, node$group1Consensus)
    }
  }
})

test_that("findMarkersTree works when a one-off split marks two classes", {
  fixture <- .multiClassFixture()
  tree <- findMarkersTree(
    fixture$features,
    fixture$class,
    threshold = 0.7,
    autoMetaclusters = FALSE
  )
  rules <- do.call(rbind, tree$rules)
  expect_false(anyNA(rules$feature))
  expect_true(all(nzchar(rules$feature)))
  f4Up <- rules$feature == "f4" & rules$direction == 1
  expect_true(any(rules$class == 1 & f4Up))
  expect_s3_class(suppressWarnings(plotDendro(tree)), "ggplot")
})

test_that("a balanced split is used when no one-off split marks one class", {
  fixture <- .multiClassFixture()
  split <- celda:::.wrapSplitHybrid(
    t(fixture$features["f1", , drop = FALSE]),
    as.factor(fixture$class),
    threshold = 0.7,
    oneoffMetric = "modified F1"
  )
  expect_identical(split$statUsed, "Split")
  expect_length(split[names(split) != "statUsed"], 1)
})

# Silence only pROC's deprecation of coords(transpose = TRUE), which
# .splitMetricPairwiseAUC() still uses
.withoutTransposeWarning <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("transpose=TRUE", conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}

test_that("findMarkersTree with pairwise AUC keeps one class per node", {
  fixture <- .multiClassFixture()
  split <- .withoutTransposeWarning(celda:::.wrapSplitHybrid(
    t(fixture$features),
    as.factor(fixture$class),
    threshold = 0.7,
    oneoffMetric = "pairwise AUC"
  ))
  nodes <- split[names(split) != "statUsed"]
  group1 <- lapply(nodes, function(node) node$group1Consensus)
  expect_true(all(lengths(group1) == 1))
  expect_true(all(nzchar(names(nodes))))

  tree <- .withoutTransposeWarning(findMarkersTree(
    fixture$features,
    fixture$class,
    oneoffMetric = "pairwise AUC",
    threshold = 0.7,
    autoMetaclusters = FALSE
  ))
  rules <- do.call(rbind, tree$rules)
  expect_true(all(nzchar(rules$feature)))
  expect_named(tree$rules, c("1", "2", "3"), ignore.order = TRUE)
})
