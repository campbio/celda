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

# These tests supply `metaclusters` so that the top-level tree is built
# without the UMAP step of `autoMetaclusters = TRUE`, which is the path from
# issue 410. Both paths share the top-level tree.

test_that("findMarkersTree merges metaclusters left without markers", {
  fixture <- .markerTreeFixture()
  metaclusters <- list(M1 = 1, M2 = 2, M3 = 3, M4 = 4)
  expect_message(
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
  expect_message(
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
  expect_message(
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

# Classes 1, 2 and 3 have 20 cells each. f1 marks class 1 and f2 marks
# class 2. f3 is spread over 2 to 3 in class 3 and in 11 of the 20 cells of
# class 1, and below 0.5 elsewhere. The pairwise AUC cut is below 0.5, so it
# takes in classes 1 and 3 and f3 cannot be a one-off marker. The modified
# F1 cut of an alternative split falls between 2 and 3, where it keeps over
# half of class 3 and under half of class 1, so f3 marks class 3 there.
.weakMarkerFixture <- function() {
  withr::with_seed(443, {
    class <- rep(seq_len(3), each = 20)
    n <- length(class)
    f1 <- ifelse(class == 1, runif(n, 5, 6), runif(n, 0, 1))
    f2 <- ifelse(class == 2, runif(n, 5, 6), runif(n, 0, 1))
  })
  f3 <- c(
    seq(2.02, 2.98, length.out = 11), seq(0, 0.48, length.out = 9),
    seq(0.01, 0.49, length.out = 20),
    seq(2, 3, length.out = 20)
  )
  features <- rbind(f1 = f1, f2 = f2, f3 = f3)
  colnames(features) <- paste0("Cell_", seq_len(n))
  list(features = features, class = class)
}

test_that("the top level never gets a forced marker", {
  # Issues 443 and 450: f3 could be forced on M3 as a marker, but the top
  # level keeps only real markers, whatever altSplit is
  fixture <- .weakMarkerFixture()
  metaclusters <- list(M1 = 1, M2 = 2, M3 = 3)
  for (altSplit in c(TRUE, FALSE)) {
    expect_message(
      tree <- findMarkersTree(
        fixture$features,
        fixture$class,
        oneoffMetric = "pairwise AUC",
        altSplit = altSplit,
        metaclusters = metaclusters
      ),
      regexp = "Metacluster M3 has no positive marker"
    )
    expect_false(any(tree$rules[["M3(3)"]]$direction == 1))
  }
})

test_that("the top level does not try an alternative split", {
  # Issue 443: no feature is left for M3, so an attempt would say so
  fixture <- .markerTreeFixture(nClass = 3)
  features <- fixture$features[c("f1", "f2"), ]
  for (altSplit in c(TRUE, FALSE)) {
    messages <- capture_messages(findMarkersTree(
      features,
      fixture$class,
      altSplit = altSplit,
      metaclusters = list(M1 = 1, M2 = 2, M3 = 3)
    ))
    expect_false(any(grepl("No alternative split added", messages)))
  }
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
  expect_s3_class(plotDendro(tree), "ggplot")
})

test_that("plotDendro does not use ggplot2's deprecated label.size", {
  # Show soft deprecations every time, not once per session (#445)
  withr::local_options(lifecycle_verbosity = "warning")
  fixture <- .multiClassFixture()
  tree <- findMarkersTree(
    fixture$features,
    fixture$class,
    threshold = 0.7,
    autoMetaclusters = FALSE
  )
  # expect_no_warning() skips deprecations unless their class is named
  plot <- expect_no_warning(
    plotDendro(tree),
    class = "lifecycle_warning_deprecated"
  )
  expect_s3_class(plot, "ggplot")
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

test_that("pairwise AUC uses the largest pairwise cut-off", {
  # B overlaps A the most (lowest AUC) but its best cut-off, 6.5, would let
  # all of C through; C's cut-off, 9.75, leaves all of C and most of B
  # below it (#457)
  class <- factor(rep(c("A", "B", "C"), c(4, 5, 4)))
  values <- c(10, 11, 12, 13, 0, 1, 2, 3, 10.5, 8, 8.5, 9, 9.5)
  features <- matrix(
    values,
    ncol = 1,
    dimnames = list(paste0("cell", seq_along(values)), "f1")
  )
  value <- celda:::.splitMetricPairwiseAUC("f1", class, features)
  expect_equal(value, 9.75)
  stat <- celda:::.splitMetricPairwiseAUC("f1", class, features, rPerf = TRUE)
  expect_equal(stat, 0.95)
  split <- celda:::.getSplit(
    "f1",
    c(f1 = stat),
    features,
    class,
    celda:::.splitMetricPairwiseAUC
  )
  expect_identical(split$group1Consensus, "A")
})

test_that("findMarkersTree with pairwise AUC keeps one class per node", {
  # pROC's deprecated coords(transpose = TRUE) is no longer used (#447).
  # Matching the message also catches it if pROC makes it a lifecycle
  # deprecation, which expect_no_warning() otherwise skips
  withr::local_options(lifecycle_verbosity = "warning")
  fixture <- .multiClassFixture()
  split <- expect_no_warning(
    celda:::.wrapSplitHybrid(
      t(fixture$features),
      as.factor(fixture$class),
      threshold = 0.7,
      oneoffMetric = "pairwise AUC"
    ),
    message = "transpose"
  )
  nodes <- split[names(split) != "statUsed"]
  group1 <- lapply(nodes, function(node) node$group1Consensus)
  expect_true(all(lengths(group1) == 1))
  expect_true(all(nzchar(names(nodes))))

  tree <- expect_no_warning(
    findMarkersTree(
      fixture$features,
      fixture$class,
      oneoffMetric = "pairwise AUC",
      threshold = 0.7,
      autoMetaclusters = FALSE
    ),
    message = "transpose"
  )
  rules <- do.call(rbind, tree$rules)
  expect_true(all(nzchar(rules$feature)))
  expect_named(tree$rules, c("1", "2", "3"), ignore.order = TRUE)
})

# Stop `expr` with an error once it has run for `seconds`, so that a
# threshold loop that never ends fails the test instead of hanging it
.withTimeLimit <- function(expr, seconds = 60) {
  setTimeLimit(elapsed = seconds, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf))
  expr
}

test_that("findMarkersTree merges unmarked metaclusters with reuseFeatures", {
  # Issue 438: f3 separates class 3 from class 4 in a balanced split, but
  # can never mark M3, M4 or their merge alone, which once kept the
  # threshold loop going
  fixture <- .markerTreeFixture()
  metaclusters <- list(M1 = 1, M2 = 2, M3 = 3, M4 = 4)
  expect_message(
    tree <- .withTimeLimit(findMarkersTree(
      fixture$features,
      fixture$class,
      metaclusters = metaclusters,
      reuseFeatures = TRUE
    )),
    regexp = "Metaclusters M3, M4 have no positive marker"
  )
  merged <- "M3+M4"
  expect_setequal(unique(tree$metaclusterLabels), c("M1", "M2", merged))
  expect_named(
    tree$rules,
    c("M1(1)", "M2(2)", paste0(merged, c("(3)", "(4)"))),
    ignore.order = TRUE
  )
})

# Classes 1-4 have 20 cells each. g marks classes 3 and 4 (metacluster M2)
# and h is high in classes 1 and 3, so it is the only feature left to split
# each metacluster
.twoMetaclusterFixture <- function() {
  withr::with_seed(455, {
    class <- rep(seq_len(4), each = 20)
    n <- length(class)
    g <- ifelse(class %in% c(3, 4), runif(n, 5, 6), runif(n, 0, 1))
    h <- ifelse(class %in% c(1, 3), runif(n, 5, 6), runif(n, 0, 1))
  })
  features <- rbind(g = g, h = h)
  colnames(features) <- paste0("Cell_", seq_len(n))
  list(features = features, class = class)
}

test_that("findMarkersTree splits a metacluster with one feature left", {
  fixture <- .twoMetaclusterFixture()
  expect_message(
    tree <- findMarkersTree(
      fixture$features,
      fixture$class,
      metaclusters = list(M1 = c(1, 2), M2 = c(3, 4))
    ),
    regexp = "Metacluster M1 has no positive marker"
  )
  expect_named(tree$rules, c("M1(1)", "M1(2)", "M2(3)", "M2(4)"),
               ignore.order = TRUE)
})

# Classes 1-4 have 20 cells each. f1 marks class 1, f2 marks class 2 and b
# marks classes 3 and 4 together. f3 is high in class 3 only, so it splits
# classes 3 and 4.
.namedMetaclusterFixture <- function() {
  withr::with_seed(455, {
    class <- rep(seq_len(4), each = 20)
    n <- length(class)
    f1 <- ifelse(class == 1, runif(n, 5, 6), runif(n, 0, 1))
    f2 <- ifelse(class == 2, runif(n, 5, 6), runif(n, 0, 1))
    b <- ifelse(class %in% c(3, 4), runif(n, 5, 6), runif(n, 0, 1))
    f3 <- ifelse(class == 3, runif(n, 5, 6), runif(n, 0, 1))
  })
  features <- rbind(f1 = f1, f2 = f2, b = b, f3 = f3)
  colnames(features) <- paste0("Cell_", seq_len(n))
  list(features = features, class = class)
}

test_that("findMarkersTree keeps metacluster names with parentheses", {
  # Issue 455: "(cluster)" was removed from names with a greedy regex
  fixture <- .namedMetaclusterFixture()
  metaclusters <- list("T (CD4)" = 1, "T (CD8)" = 2, "B (all)" = c(3, 4))
  expect_no_warning(
    tree <- suppressMessages(findMarkersTree(
      fixture$features,
      fixture$class,
      metaclusters = metaclusters
    ))
  )
  expect_setequal(unique(tree$metaclusterLabels), names(metaclusters))
  expect_identical(tree$classLabels, as.character(fixture$class))
  expect_named(
    tree$rules,
    c("T (CD4)(1)", "T (CD8)(2)", "B (all)(3)", "B (all)(4)"),
    ignore.order = TRUE
  )
  for (label in names(tree$rules)) {
    rules <- tree$rules[[label]]
    expect_identical(unique(rules$class),
                     sub("^.*\\((.*)\\)$", "\\1", label))
    expect_identical(unique(rules$metacluster),
                     sub("\\([^()]*\\)$", "", label))
  }
  expect_setequal(unique(tree$branchPoints$top_level$metacluster),
                  names(metaclusters))
  expect_true("B (all)_level_1" %in% names(tree$branchPoints))
  plot <- plotDendro(tree)
  plotLabels <- unlist(lapply(plot$layers, function(layer) {
    as.character(layer$data$label)
  }))
  expect_true(all(names(metaclusters) %in% plotLabels))
  # Leaves show the clusters
  leafLabels <- plot$scales$get_scales("x")$labels
  expect_setequal(leafLabels, as.character(seq_len(4)))
})

test_that("merged metaclusters keep names with parentheses", {
  fixture <- .markerTreeFixture()
  metaclusters <- list("A (x)" = 1, "B (y)" = 2, "C (z)" = 3, "D (w)" = 4)
  expect_message(
    tree <- findMarkersTree(
      fixture$features,
      fixture$class,
      metaclusters = metaclusters
    ),
    regexp = "Metaclusters C (z), D (w) have no positive marker",
    fixed = TRUE
  )
  expect_setequal(unique(tree$metaclusterLabels),
                  c("A (x)", "B (y)", "C (z)+D (w)"))
})

test_that("branch points group class labels that contain spaces", {
  # Level 1 splits classes "T (a)" and "T (b)" from "T (c)" and "T (d)".
  # Level 2 has a one-off split in each branch.
  rulesRow <- function(cl, feature, direction, statUsed, level) {
    data.frame(feature = feature, direction = direction, stat = 0.9,
               statUsed = statUsed, level = level, class = cl,
               stringsAsFactors = FALSE)
  }
  classes <- c("T (a)", "T (b)", "T (c)", "T (d)")
  rules <- lapply(seq_along(classes), function(i) {
    cl <- classes[i]
    left <- i <= 2
    rbind(
      rulesRow(cl, "f0", if (left) 1 else -1, "Split", 1),
      rulesRow(cl, if (left) "fx" else "fy", if (i %% 2 == 1) 1 else -1,
               "One-off", 2)
    )
  })
  names(rules) <- classes
  branchPoints <- celda:::.createBranchPoints(rules, character(0))
  levelTwo <- branchPoints[startsWith(names(branchPoints), "level_2")]
  groups <- lapply(levelTwo, function(br) sort(unique(br$class)))
  expect_setequal(groups, list(classes[1:2], classes[3:4]))
})

# Classes 1-12 have 10 cells each. Metacluster M1 holds classes 1 and 2,
# which m1 marks, and f1 marks class 1. Each other class k has its own
# metacluster, M2 to M11, marked by f<k>.
.manyMetaclusterFixture <- function() {
  withr::with_seed(444, {
    class <- rep(seq_len(12), each = 10)
    n <- length(class)
    features <- t(vapply(seq_len(12), function(k) {
      ifelse(class == k, runif(n, 5, 6), runif(n, 0, 1))
    }, FUN.VALUE = double(n)))
    m1 <- ifelse(class %in% c(1, 2), runif(n, 5, 6), runif(n, 0, 1))
  })
  features <- rbind(features[-2, ], m1)
  rownames(features) <- c(paste0("f", c(1, 3:12)), "m1")
  colnames(features) <- paste0("Cell_", seq_len(n))
  metaclusters <- c(list(M1 = c(1, 2)),
                    stats::setNames(as.list(3:12), paste0("M", 2:11)))
  list(features = features, class = class, metaclusters = metaclusters)
}

test_that("findMarkersTree relabels a metacluster's nodes exactly", {
  # Issue 444: M1 must not pick up the clusters of M10 and M11
  fixture <- .manyMetaclusterFixture()
  tree <- suppressMessages(findMarkersTree(
    fixture$features,
    fixture$class,
    metaclusters = fixture$metaclusters
  ))
  root <- attributes(tree$dendro)
  expect_false(anyDuplicated(root$classLabels) > 0)
  expect_equal(root$members, 12)
  expect_setequal(root$classLabels, names(tree$rules))
})

# A top-level balanced split of metaclusters A and B against C, with the
# fields .generateTreeList() adds
.balancedTopLevel <- function() {
  list(f = list(group1Consensus = c("A", "B"), group2Consensus = "C"),
       statUsed = "Split", fUsed = "f", dirs = 1)
}

test_that("metaclusters sharing a balanced split's marker are unmarked", {
  labels <- rep(c("A", "B", "C"), each = 2)
  expect_identical(
    celda:::.unmarkedMetaclusters(.balancedTopLevel(), labels),
    c("A", "B", "C")
  )
  oneOff <- list(f = list(group1Consensus = "A",
                          group2Consensus = c("B", "C")),
                 statUsed = "One-off", fUsed = "f", dirs = 1)
  expect_identical(celda:::.unmarkedMetaclusters(oneOff, labels),
                   c("B", "C"))
})

test_that("findMarkersTree explains when no feature is left for a subtree", {
  # g and h both mark M2, so neither is left to split classes 3 and 4
  withr::with_seed(455, {
    class <- rep(seq_len(4), each = 20)
    n <- length(class)
    g <- ifelse(class %in% c(3, 4), runif(n, 5, 6), runif(n, 0, 1))
    h <- ifelse(class %in% c(3, 4), runif(n, 5, 6), runif(n, 0, 1))
  })
  features <- rbind(g = g, h = h)
  colnames(features) <- paste0("Cell_", seq_len(n))
  expect_error(
    suppressWarnings(suppressMessages(findMarkersTree(
      features,
      class,
      metaclusters = list(M1 = c(1, 2), M2 = c(3, 4))
    ))),
    regexp = "No features are left to split the clusters of metacluster M2"
  )
})

# Classes 1, 2 and 3 have 20 cells each. f1 marks class 1 and f2 marks
# class 2, but no feature marks class 3. w is high in 15 class-3 cells and
# in 8 class-2 cells, so a marker forced on it for class 3 would also take
# in class-2 cells and miss the other class-3 cells.
.overlapFixture <- function() {
  withr::with_seed(450, {
    class <- rep(seq_len(3), each = 20)
    n <- length(class)
    f1 <- ifelse(class == 1, runif(n, 5, 6), runif(n, 0, 1))
    f2 <- ifelse(class == 2, runif(n, 5, 6), runif(n, 0, 1))
    wHigh <- c(rep(FALSE, 20), rep(c(TRUE, FALSE), c(8, 12)),
               rep(c(TRUE, FALSE), c(15, 5)))
    w <- ifelse(wHigh, runif(n, 2, 3), runif(n, 0, 1))
  })
  features <- rbind(f1 = f1, f2 = f2, w = w)
  colnames(features) <- paste0("Cell_", seq_len(n))
  list(features = features, class = class)
}

test_that("the top level keeps only real markers (issue 450)", {
  fixture <- .overlapFixture()
  for (altSplit in c(TRUE, FALSE)) {
    messages <- capture_messages(
      tree <- findMarkersTree(
        fixture$features,
        fixture$class,
        altSplit = altSplit,
        metaclusters = list(M1 = 1, M2 = 2, M3 = 3)
      )
    )
    expect_false(any(grepl("Decreasing classifier threshold", messages)))
    expect_true(any(grepl("Metacluster M3 has no positive marker",
                          messages)))
    topM3 <- tree$rules[["M3(3)"]]
    expect_false(any(topM3$level == 1 & topM3$direction == 1))
    # M3 is identified by the absence of f1 and f2, so no cell is missed
    expect_false(any(tree$prediction == "MISSING"))
    expect_equal(tree$performance$accuracy, 1)
  }
})

test_that("forced markers are not used for prediction", {
  # Without metaclusters, class 3 has no marker and gets one forced on f3
  fixture <- .weakMarkerFixture()
  trees <- lapply(c(TRUE, FALSE), function(altSplit) {
    suppressMessages(findMarkersTree(
      fixture$features,
      fixture$class,
      oneoffMetric = "pairwise AUC",
      altSplit = altSplit,
      autoMetaclusters = FALSE
    ))
  })
  rules3 <- trees[[1]]$rules[["3"]]
  forced <- rules3$feature == "f3" & rules3$direction == 1
  expect_true(any(forced))
  # Only the forced marker's row is marked as alternative
  expect_identical(rules3$alternative, forced)
  expect_identical(trees[[1]]$prediction, trees[[2]]$prediction)
})

test_that("getDecisions skips rules marked as alternative", {
  # Class A is low on g and has a forced marker on f; class B is high on g
  rules <- list(
    A = data.frame(feature = c("f", "g"), direction = c(1, -1),
                   value = c(1, 1), stat = c(0.6, 0.9),
                   statUsed = "One-off", level = 1,
                   alternative = c(TRUE, FALSE)),
    B = data.frame(feature = c("f", "g"), direction = c(-1, 1),
                   value = c(1, 1), stat = c(0.6, 0.9),
                   statUsed = "One-off", level = 1,
                   alternative = c(TRUE, FALSE))
  )
  features <- rbind(f = c(0, 2, 0), g = c(0, 0, 2))
  expect_identical(unname(getDecisions(rules, features)),
                   c("A", "A", "B"))
  # Rules without the column, as in older trees, use every rule
  oldRules <- lapply(rules, function(r) r[names(r) != "alternative"])
  expect_identical(unname(getDecisions(oldRules, features)),
                   c(NA, "A", "B"))
})

test_that("findMarkersTree drops metaclusters that have no marker", {
  # f1 and f2 each mark one cluster of M1, and f3 is high in M1 and in
  # half of M2, so no feature marks M1 or M2 alone
  fixture <- .markerTreeFixture()
  expect_message(
    tree <- findMarkersTree(
      fixture$features,
      fixture$class,
      metaclusters = list(M1 = c(1, 2), M2 = c(3, 4))
    ),
    regexp = "built without metaclusters"
  )
  flat <- suppressMessages(findMarkersTree(
    fixture$features,
    fixture$class,
    autoMetaclusters = FALSE
  ))
  expect_identical(tree$rules, flat$rules)
  expect_identical(tree$prediction, flat$prediction)
})

test_that("gene-level rules keep metaclusters without a marker", {
  # Each feature of .overlapFixture() becomes two genes
  fixture <- .overlapFixture()
  features <- fixture$features
  counts <- round(rbind(features * 10, features * 10 + 1))
  rownames(counts) <- c(paste0(rownames(features), "a"),
                        paste0(rownames(features), "b"))
  featureLabels <- rep(rownames(features), 2)
  # M3 alone has no marker, then M1 is the only metacluster with one
  for (metaclusters in list(list(M1 = 1, M2 = 2, M3 = 3),
                            list(M1 = 1, M23 = c(2, 3)))) {
    tree <- suppressMessages(findMarkersTree(
      features,
      fixture$class,
      metaclusters = metaclusters,
      featureLabels = featureLabels,
      counts = counts
    ))
    rulesM1 <- tree$rules[["M1(1)"]]
    expect_true(all(rulesM1$direction == 1))
    expect_setequal(rulesM1$gene, c("f1a", "f1b"))
  }
  # M3 is described by the absence of the top-level markers
  rulesM3 <- tree$rules[["M23(3)"]]
  expect_true(all(rulesM3$direction[rulesM3$level == 1] == -1))
  expect_true(any(rulesM3$level == 1))
})
