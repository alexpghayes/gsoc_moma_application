context("test-sfpca.R")

# tests:
#   - gets the same answers as MATLAB implementation
#   - input checking
#   - matches PCA in base case
#   - matches other factorizations in non-base cases?

X <- as.matrix(iris[, 1:4])

sfpca_r(X, 10, 4)
sfpca_r(X, 40, 40)

sfpca(X, 10, 4)
sfpca(X, 40, 40)

sfpca(as.matrix(iris[, 1:4]), 10, 4)

test_that("R implementation", {
  expect_equal(2 * 2, 4)
})

