test_that("dd_poly2d_powers enumerates a total-degree basis in the documented order", {
  p <- dd_poly2d_powers(2)
  expect_equal(nrow(p), 6L)  # (2+1)(2+2)/2
  expect_equal(rowSums(p), c(0, 1, 1, 2, 2, 2))
  expect_equal(unname(p), matrix(c(
    0, 0,   # 1
    1, 0,   # x1
    0, 1,   # x2
    2, 0,   # x1^2
    1, 1,   # x1*x2
    0, 2    # x2^2
  ), ncol = 2, byrow = TRUE))
})

test_that("dd_poly_dictionary2d columns match direct power computation", {
  x1 <- c(0.5, -1, 2, 0)
  x2 <- c(1, 0.5, -2, 3)
  dict <- dd_poly_dictionary2d(x1, x2, degree = 2)
  expect_equal(dim(dict), c(4L, 6L))
  expect_equal(dict[, 1], rep(1, 4))
  expect_equal(dict[, 2], x1)
  expect_equal(dict[, 3], x2)
  expect_equal(dict[, 4], x1^2)
  expect_equal(dict[, 5], x1 * x2)
  expect_equal(dict[, 6], x2^2)
})

test_that("dd_fit_polynomial2d recovers a known polynomial exactly (threshold = 0, no noise)", {
  set.seed(3)
  x1 <- runif(300, -2, 2)
  x2 <- runif(300, -2, 2)
  # f(x1, x2) = 1.5 - 2*x1 + 0.5*x2 + 0.3*x1^2 - 0.7*x1*x2
  y <- 1.5 - 2 * x1 + 0.5 * x2 + 0.3 * x1^2 - 0.7 * x1 * x2

  fit <- dd_fit_polynomial2d(x1, x2, y, degree = 2, threshold = 0)
  powers <- dd_poly2d_powers(2)
  expected <- numeric(nrow(powers))
  expected[powers[, 1] == 0 & powers[, 2] == 0] <- 1.5
  expected[powers[, 1] == 1 & powers[, 2] == 0] <- -2
  expected[powers[, 1] == 0 & powers[, 2] == 1] <- 0.5
  expected[powers[, 1] == 2 & powers[, 2] == 0] <- 0.3
  expected[powers[, 1] == 1 & powers[, 2] == 1] <- -0.7
  expect_equal(fit$coeffs, expected, tolerance = 1e-8)

  # predict.poly2d should reproduce y exactly (noiseless data, correct degree)
  pred <- stats::predict(fit, x1, x2)
  expect_equal(pred, y, tolerance = 1e-8)
})

test_that("dd_fit_polynomial2d sparsifies small terms when thresholded", {
  set.seed(4)
  x1 <- runif(500, -1, 1)
  x2 <- runif(500, -1, 1)
  y <- 3 * x1 + 0.001 * x2 + rnorm(500, sd = 0.01)  # x2 term negligible
  fit <- dd_fit_polynomial2d(x1, x2, y, degree = 1, threshold = 0.05)
  powers <- dd_poly2d_powers(1)
  x2_idx <- which(powers[, 1] == 0 & powers[, 2] == 1)
  expect_equal(fit$coeffs[x2_idx], 0)
  x1_idx <- which(powers[, 1] == 1 & powers[, 2] == 0)
  expect_equal(fit$coeffs[x1_idx], 3, tolerance = 0.2)
})

test_that("format.poly2d produces a readable, non-crashing string", {
  fit <- new_poly2d(c(1, -2, 0, 0.5, 0, 0), degree = 2)
  s <- format(fit)
  expect_type(s, "character")
  expect_true(grepl("x1", s))
  expect_false(grepl("x2", s))  # its coefficient (index 3) is 0, so x2 term dropped

  zero_fit <- new_poly2d(rep(0, 6), degree = 2)
  expect_equal(format(zero_fit), "0")
})

test_that("dd_tune_threshold2d returns a valid poly2d with a cv_curve attribute", {
  set.seed(5)
  x1 <- runif(400, -1, 1)
  x2 <- runif(400, -1, 1)
  y <- 2 * x1 - 0.5 * x2 + rnorm(400, sd = 0.05)
  fit <- dd_tune_threshold2d(x1, x2, y, degree = 1, folds = 4)
  expect_s3_class(fit, "poly2d")
  expect_true(!is.null(attr(fit, "cv_curve")))
  expect_true(!is.null(attr(fit, "threshold")))
})
