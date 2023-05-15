test_that("Get errors for input x", {
  x1 = matrix(rnorm(300), ncol = 3)
  x2 = matrix(rnorm(300), ncol = 3)
  expect_error(ssMRCD(rbind(x1, x2),
                      lambda = 0.4,
                      weights = matrix(c(0,1,1,0), 2)))
})

test_that("Get errors for input lambda", {
  x1 = matrix(rnorm(300), ncol = 3)
  x2 = matrix(rnorm(300), ncol = 3)
  expect_error(ssMRCD(list(x1, x2),
                      lambda = -0.4,
                      weights = matrix(c(0,1,1,0), 2)))
})

test_that("Get errors for input weights", {
  x1 = matrix(rnorm(300), ncol = 3)
  x2 = matrix(rnorm(300), ncol = 3)
  expect_error(ssMRCD(list(x1, x2),
                      lambda = 0.4,
                      weights = matrix(c(1,1,1,1), 2)))
})

test_that("Get errors for input weights", {
  x1 = matrix(rnorm(300), ncol = 3)
  x2 = matrix(rnorm(300), ncol = 3)
  expect_error(ssMRCD(list(x1, x2),
                      lambda = 0.4,
                      weights = data.frame(matrix(c(0,1,1,0), 2))))
})

test_that("Get errors for input alpha", {
  x1 = matrix(rnorm(300), ncol = 3)
  x2 = matrix(rnorm(300), ncol = 3)
  expect_error(ssMRCD(list(x1, x2),
                      lambda = 0.4,
                      weights = matrix(c(0,1,1,0), 2),
                      alpha = 0.1))
})

test_that("Some correct calculations", {
  x1 = matrix(rnorm(300), ncol = 3)
  x2 = matrix(rnorm(300), ncol = 3)
  x3 = matrix(rnorm(600), ncol = 3)
  N = 3
  out = ssMRCD(list(x1, x2, x3),
         lambda = 0.4,
         weights = matrix(c(0,1/2,1/2,1,0,0, 1, 0, 0), 3, byrow = T),
         alpha = 0.75)
  expect_equal(out$MRCDcov[[N]] %*% out$MRCDicov[[N]], diag(1, 3))
})

