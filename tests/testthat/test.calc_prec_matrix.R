context("testing function calc_prec_matrix()")

test_that("calc_prec_matrix", {

n <- 100
gamma <- 0.5
neighborhood_matrix <- Matrix::sparseMatrix(1:n, 1:n, x=rep(1,n))

precision <- calc_prec_matrix(neighborhood_matrix, gamma)

expect_equal(dim(precision)[1], n)
expect_equal(dim(precision)[1], n)
})
