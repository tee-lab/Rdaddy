test_that("dd_ridge_fit recovers exact coefficients for noiseless data (alpha = 0)", {
  set.seed(1)
  x <- seq(-2, 2, length.out = 200)
  true_coeffs <- c(0.5, -1.2, 0.3, 0.8)  # 1, x, x^2, x^3
  dict <- dd_poly_dictionary(x, degree = 3)
  y <- as.vector(dict %*% true_coeffs)

  fit <- dd_ridge_fit(dict, y, alpha = 0)
  expect_equal(fit, true_coeffs, tolerance = 1e-8)
})

test_that("dd_ridge_fit shrinks coefficients as alpha grows", {
  set.seed(2)
  x <- seq(-2, 2, length.out = 200)
  dict <- dd_poly_dictionary(x, degree = 2)
  y <- 1 + 2 * x + rnorm(200, sd = 0.01)

  fit_small_alpha <- dd_ridge_fit(dict, y, alpha = 0.01)
  fit_large_alpha <- dd_ridge_fit(dict, y, alpha = 1000)

  expect_lt(sum(fit_large_alpha^2), sum(fit_small_alpha^2))
})

test_that("dd_stlsq recovers a sparse polynomial and zeroes small terms", {
  set.seed(3)
  x <- seq(-1, 1, length.out = 500)
  # y = 0.9 x^3 + tiny noise; degree-3 dictionary, other coefficients should be thresholded to 0
  dict <- dd_poly_dictionary(x, degree = 3)
  true_coeffs <- c(0, 0, 0, 0.9)
  y <- as.vector(dict %*% true_coeffs) + rnorm(500, sd = 1e-4)

  fit <- dd_stlsq(dict, y, threshold = 0.05)
  expect_equal(fit$coeffs[1:3], c(0, 0, 0), tolerance = 1e-6)
  expect_equal(fit$coeffs[4], 0.9, tolerance = 1e-2)
})

test_that("dd_fit_polynomial errors clearly with too few points for the degree", {
  expect_error(dd_fit_polynomial(x = c(1, 2), y = c(1, 2), degree = 5),
               "Not enough finite data points")
})

test_that("poly1d predict/format/print work", {
  p <- new_poly1d(c(1, 2, 3))  # 1 + 2x + 3x^2
  expect_equal(predict(p, c(0, 1, 2)), c(1, 6, 17))
  expect_match(format(p), "x\\^2")
  expect_output(print(p), "x\\^2")
})

test_that("dd_tune_threshold picks a threshold and returns a valid fit", {
  set.seed(4)
  x <- seq(-1, 1, length.out = 300)
  y <- 0.5 * x + rnorm(300, sd = 0.05)
  fit <- dd_tune_threshold(x, y, degree = 3, steps = 5, folds = 3)
  expect_s3_class(fit, "poly1d")
  expect_true(is.numeric(attr(fit, "threshold")))
})
