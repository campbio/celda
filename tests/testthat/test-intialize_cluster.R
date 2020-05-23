# test_that(desc = "Test 'random' initialization for all models", {
#     simRes <- simulateCells(model = "celda_CG")
#     modelCG <- celda_CG(simRes$counts,
#         simRes$sampleLabel,
#         K = simRes$K,
#         L = simRes$L,
#         zInitialize = "random",
#         yInitialize = "random",
#         maxIter = 5,
#         splitOnLast = FALSE,
#         splitOnIter = -1)
#     expect_true(is(modelCG, "celda_CG"))
#     modelG <- celda_G(simRes$counts,
#         L = simRes$L,
#         yInitialize = "random",
#         maxIter = 5,
#         splitOnLast = FALSE,
#         splitOnIter = -1)
#     expect_true(is(modelG, "celda_G"))
#     modelC <- celda_C(simRes$counts,
#         simRes$sampleLabel,
#         K = simRes$K,
#         zInitialize = "random",
#         maxIter = 5,
#         splitOnLast = FALSE,
#         splitOnIter = -1)
#     expect_true(is(modelC, "celda_C"))
# })
#
# test_that(desc = "Testing .initializeCluster for random initialization", {
#     ## Completely random
#     z <- .initializeCluster(10, 100)
#     expect_true(length(z) == 100 & length(unique(z) == 10))
#     expect_error(z <- .initializeCluster(100, 10))
#
#     ## With all values initialized
#     initZ <- rep(seq(10), each = 10)
#     z <- .initializeCluster(10, 100, initial = initZ)
#     expect_true(all(initZ == z))
#     expect_error(z <- .initializeCluster(10, 100, initial = initZ[seq(99)]))
#     expect_error(z <- .initializeCluster(11, 100, initial = initZ))
#     expect_error(z <- .initializeCluster(10, 99, initial = initZ))
#
#     ## With only a few values initialized
#     fixedZ <- rep(NA, 100)
#     fixedZ[seq(10)] <- 1
#     z <- .initializeCluster(10, 100, fixed = fixedZ)
#     expect_true(all(z[seq(10)] == 1) &
#             length(z) == 100 & length(unique(z)) == 10)
#     expect_error(z <- .initializeCluster(10, 100, fixed = fixedZ[seq(99)]))
#     fixedZ[seq(10)] <- 11
#     expect_error(z <- .initializeCluster(10, 100, fixed = fixedZ))
# })
