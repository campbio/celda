# recursiveSplitCell, recursiveSplitModule, and plotRPC
library(celda)
context("Testing recursiveSplit and plotRPC")

rsSim <- simulateCells("celda_CG",
  S = 1, K = 4, L = 6, G = 150, CRange = c(80, 100), seed = 7
)
rsSim <- selectFeatures(rsSim, minCount = 0, minCell = 0)
rsY <- as.integer(factor(
  rowData(altExp(rsSim))$celda_feature_module
))
rsZ <- as.integer(rsSim$celda_cell_cluster)

# celda_C models, K varies
rsCellC <- recursiveSplitCell(rsSim,
  initialK = 2, maxK = 6, tempL = 6, verbose = FALSE
)
# celda_CG models, K varies and L is fixed
rsCellCG <- recursiveSplitCell(rsSim,
  initialK = 2, maxK = 6, yInit = rsY, verbose = FALSE
)
# celda_CG models, L varies and K is fixed
rsModuleCG <- recursiveSplitModule(rsSim,
  initialL = 3, maxL = 7, zInit = rsZ, verbose = FALSE
)

kneeLayerData <- function(g) {
  layer <- g$layers[[length(g$layers)]]
  expect_true(inherits(layer$geom, "GeomText"))
  layer$data
}

expectedKnee <- function(sce, by) {
  rp <- runParams(sce)
  rp <- rp[order(rp[[by]]), ]
  .findKnee(rp[[by]], rp$mean_perplexity, nPoints = 1)
}

samePartition <- function(a, b) {
  length(unique(paste(a, b))) == length(unique(b))
}

test_that(desc = "Testing recursiveSplitModule with character zInit", {
  zChar <- paste0("pop_", rsZ)
  for (useMeta in c(FALSE, TRUE)) {
    res <- recursiveSplitModule(rsSim,
      initialL = 3, maxL = 4, zInit = zChar, useMetaCells = useMeta,
      verbose = FALSE
    )
    expect_equal(runParams(res)$L, 3:4)
  }
  # With useMetaCells = TRUE the cell partition stays as given
  for (mod in resList(res)) {
    expect_true(samePartition(celdaClusters(mod)$z, zChar))
  }
})

test_that(desc = "Testing recursiveSplitModule useMetaCells needs zInit", {
  expect_error(
    recursiveSplitModule(rsSim, useMetaCells = TRUE, verbose = FALSE),
    "'zInit' needs to be specified"
  )
})

test_that(desc = "Testing plotRPC knees for celda_C models", {
  g <- plotRPC(rsCellC, plotKnees = TRUE, knees = 1)
  knee <- kneeLayerData(g)
  expect_equal(knee$pos, expectedKnee(rsCellC, "K"))
  # The line is drawn at the knee's position on the discrete x axis
  xLevels <- ggplot2::layer_scales(g)$x$range$range
  expect_equal(xLevels[knee$nx], as.character(knee$pos))
})

test_that(desc = "Testing plotRPC knees for celdaList input", {
  celdaList <- S4Vectors::metadata(altExp(rsCellC))$celda_grid_search
  g <- plotRPC(celdaList, plotKnees = TRUE, knees = 1)
  expect_equal(kneeLayerData(g)$pos, expectedKnee(rsCellC, "K"))
})

test_that(desc = "Testing plotRPC knees use the mean of repeated values", {
  rp <- data.frame(
    K = rep(2:7, each = 2),
    mean_perplexity = rep(c(100, 60, 40, 35, 33, 32), each = 2) +
      c(-30, 30)
  )
  g <- ggplot2::ggplot()
  knee <- kneeLayerData(.addRPCKnees(g, rp, "celda_C", knees = 1))
  expect_equal(knee$pos, .findKnee(2:7, c(100, 60, 40, 35, 33, 32)))
})

test_that(desc = "Testing plotRPC knees for celda_CG models varying K", {
  expect_no_error(g <- plotRPC(rsCellCG, plotKnees = TRUE, knees = 1))
  expect_equal(kneeLayerData(g)$pos, expectedKnee(rsCellCG, "K"))
})

test_that(desc = "Testing plotRPC knees for celda_CG models varying L", {
  g <- plotRPC(rsModuleCG, plotKnees = TRUE, knees = 1)
  expect_equal(kneeLayerData(g)$pos, expectedKnee(rsModuleCG, "L"))
})

test_that(desc = "Testing plotRPC knees reject grids varying K and L", {
  data(sceCeldaCGGridSearch)
  sce <- resamplePerplexity(sceCeldaCGGridSearch)
  expect_error(
    plotRPC(sce, plotKnees = TRUE),
    "only K or only L"
  )
})

test_that(desc = "Testing plotRPC knees argument checking", {
  expect_error(
    plotRPC(rsCellC, plotKnees = TRUE, knees = 0),
    "positive integers"
  )
  expect_error(
    plotRPC(rsCellC, plotKnees = TRUE, knees = 6),
    "cannot exceed 5\\. Requested: 6"
  )
})

test_that(desc = "Testing plotRPC builds without ggplot2 warnings", {
  withr::local_options(lifecycle_verbosity = "warning")
  for (sce in list(rsCellC, rsCellCG, rsModuleCG)) {
    expect_no_warning(ggplot2::ggplot_build(plotRPC(sce, plotKnees = TRUE)))
  }
})

test_that(desc = "Testing recursiveSplitModule useMetaCells rejects NA zInit", {
  zNA <- rsZ
  zNA[1] <- NA
  expect_error(
    recursiveSplitModule(rsSim,
      initialL = 3, maxL = 4, zInit = zNA, useMetaCells = TRUE,
      verbose = FALSE
    ),
    "'zInit' contains NA"
  )
})

test_that(desc = "Testing recursiveSplitCell clears an old tempY", {
  expect_false(is.null(S4Vectors::metadata(altExp(rsCellC))$tempY))
  rerun <- recursiveSplitCell(rsCellC,
    initialK = 2, maxK = 3, yInit = rsY, verbose = FALSE
  )
  expect_null(S4Vectors::metadata(altExp(rerun))$tempY)
})
