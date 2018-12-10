library(celda)
context("Testing error checking in C-level matrix sum functions")

## Test internal error checking
mat = matrix(1:5, ncol=10, nrow=10)
label1 = rep(1:2, each=5)
label2 = as.factor(1:100) 
label3 = as.factor(label1)
label4 = label3
label4[1:2] = 2
label5 = as.factor(rep(1:5, each=2))

test_that(desc = "Testing rowSumByGroup",{
  expect_error(.Call("_rowSumByGroup", mat, label1))
  expect_error(.Call("_rowSumByGroup", mat, label2))
  res <- .Call("_rowSumByGroup", mat, label3)
  expect_true(all(res == rowsum(mat, label3)))
  res <- rowSumByGroup(mat, label3, 2)
  expect_true(all(res == rowsum(mat, label3))) 
})

test_that(desc = "Testing rowSumByGroupChange",{
  res <- rowsum(mat, label3)
  expect_error(.Call("_rowSumByGroupChange", mat, res, label4, label1))
  expect_error(.Call("_rowSumByGroupChange", mat, res, label4, label2))
  expect_error(.Call("_rowSumByGroupChange", mat[-1,], res, label4, label3))
  expect_error(.Call("_rowSumByGroupChange", mat[,-1], res, label4, label3))
  expect_error(.Call("_rowSumByGroupChange", mat, res, label4, label5))
  res2 <- .Call("_rowSumByGroupChange", mat, res, label4, label3)
  expect_true(all(res2 == rowsum(mat, label4)))
  res <- rowsum(mat, label3)
  res2 <- rowSumByGroupChange(mat, res, label4, label3, 2)
  expect_true(all(res2 == rowsum(mat, label4)))
})

test_that(desc = "Testing colSumByGroup",{
  expect_error(.Call("_colSumByGroup", mat, label1))
  expect_error(.Call("_colSumByGroup", mat, label2))
  res <- .Call("_colSumByGroup", mat, label3)
  expect_true(all(res == t(rowsum(t(mat), label3))))
  res <- colSumByGroup(mat, label3, 2)
  expect_true(all(res == t(rowsum(t(mat), label3))))
})

test_that(desc = "Testing colSumByGroupChange",{
  res <- t(rowsum(t(mat), label3))
  expect_error(.Call("_colSumByGroupChange", mat, res, label4, label1))
  expect_error(.Call("_colSumByGroupChange", mat, res, label4, label2))
  expect_error(.Call("_colSumByGroupChange", mat[,-1], res, label4, label3))
  expect_error(.Call("_colSumByGroupChange", mat, res, label4, label5))
  expect_error(.Call("_colSumByGroupChange", mat[-1,], res, label4, label3))
  res2 <- .Call("_colSumByGroupChange", mat, res, label4, label3)
  expect_true(all(res2 == t(rowsum(t(mat), label4))))
  res <- t(rowsum(t(mat), label3))
  res2 <- colSumByGroupChange(mat, res, label4, label3, 2)
  expect_true(all(res2 == t(rowsum(t(mat), label4))))
})

storage.mode(mat) = "numeric"
test_that(desc = "Testing rowSumByGroup.numeric",{
  expect_error(.Call("_rowSumByGroup_numeric", mat, label1))
  expect_error(.Call("_rowSumByGroup_numeric", mat, label2))
  res <- .Call("_rowSumByGroup_numeric", mat, label3)
  expect_true(all(res == rowsum(mat, label3)))
  res <- rowSumByGroup.numeric(mat, label3, 2)
  expect_true(all(res == rowsum(mat, label3)))
})

test_that(desc = "Testing colSumByGroup.numeric",{
  expect_error(.Call("_colSumByGroup_numeric", mat, label1))
  expect_error(.Call("_colSumByGroup_numeric", mat, label2))
  res <- .Call("_colSumByGroup_numeric", mat, label3)
  expect_true(all(res == t(rowsum(t(mat), label3))))
  res <- colSumByGroup.numeric(mat, label3, 2)
  expect_true(all(res == t(rowsum(t(mat), label3))))
})

# test_that(desc = "Testing fastNormProp Rcpp funtion",{
#   res <- fastNormProp(mat, 0)
#   res2 <- prop.table(mat, 2)
#   expect_true(all(res == res2))
# })
# 
# test_that(desc = "Testing fastNormPropLog Rcpp funtion",{
#   res <- fastNormPropLog(mat, 0)
#   res2 <- log(prop.table(mat, 2))
#   expect_true(all(res == res2))
# })
# 
# test_that(desc = "Testing fastNormPropSqrt Rcpp funtion",{
#   res <- fastNormPropSqrt(mat, 0)
#   res2 <- sqrt(prop.table(mat, 2))
#   expect_true(all(res == res2))
# })


