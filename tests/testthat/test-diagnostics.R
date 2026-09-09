test_that("dd_gaussianity_test accepts genuinely Gaussian noise", {
  set.seed(1)
  x <- rnorm(5000)
  res <- dd_gaussianity_test(x)
  expect_true(res$is_gaussian)
  expect_gt(res$p_value, 0.05)
})

test_that("dd_gaussianity_test rejects strongly non-Gaussian noise", {
  set.seed(2)
  x <- rexp(5000, rate = 1)  # strongly skewed, not remotely Gaussian
  res <- dd_gaussianity_test(x)
  expect_false(res$is_gaussian)
  expect_lt(res$p_value, 0.01)
})

test_that("dd_gaussianity_test errors clearly on too little data", {
  expect_error(dd_gaussianity_test(c(1, 2, 3)), "at least 8")
})

test_that("dd_autocorrelation returns a tidy ACF and a plausible autocorr_time", {
  set.seed(3)
  n <- 5000
  x <- numeric(n)
  for (i in 2:n) x[i] <- 0.9 * x[i - 1] + rnorm(1)  # strongly autocorrelated AR(1)

  res <- dd_autocorrelation(x, lags = 200)
  expect_true(is.data.frame(res$acf))
  expect_equal(names(res$acf), c("lag", "acf"))
  expect_true(res$autocorr_time > 1)  # AR(1) with phi=0.9 decays slowly
})
