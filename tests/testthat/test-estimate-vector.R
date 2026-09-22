test_that("dd_estimate_vector gives full 2D bin coverage with no double-counting", {
  # Regression test for a real PyDaddy vector-path bug: its 2D bin membership
  # test is closed on both ends (bin_x <= x & x <= bin_x + inc_x), so a point
  # sitting exactly on a shared edge is counted in two bins at once. daddyR
  # reuses dd_bin_index()'s half-open, exactly-tiling bins for both axes, so
  # every point should land in exactly one (x1, x2) bin, with no gaps.
  set.seed(1)
  n <- 3000
  x1 <- runif(n, -1, 1)
  x2 <- runif(n, -1, 1)
  est <- dd_estimate_vector(x1, x2, t_int = 1, Dt = 1L, dt = 1L, bins = 10L)

  expect_length(est$op_x, 10)
  expect_length(est$op_y, 10)
  expect_equal(dim(est$avg_drift1), c(10L, 10L))

  # Every point (via dd_bin_index on the full series) must fall in exactly
  # one bin per axis -- no NAs, all in range.
  bx <- dd_bin_index(x1, est$edges_x)
  by <- dd_bin_index(x2, est$edges_y)
  expect_false(anyNA(bx)); expect_false(anyNA(by))
  expect_true(all(bx >= 1 & bx <= 10))
  expect_true(all(by >= 1 & by <= 10))
})

test_that("dd_estimate_vector recovers known drift/diffusion for an independent 2D OU process", {
  # Two independent OU processes: dx1 = -theta1*x1*dt + sigma1*dW1,
  # dx2 = -theta2*x2*dt + sigma2*dW2 (G12 should be ~0 since the noise
  # sources are independent).
  set.seed(7)
  theta1 <- 0.6; sigma1 <- 0.25
  theta2 <- 1.0; sigma2 <- 0.4
  dt_sim <- 0.01; n <- 150000
  x1 <- numeric(n); x2 <- numeric(n)
  for (i in 2:n) {
    x1[i] <- x1[i - 1] - theta1 * x1[i - 1] * dt_sim + sigma1 * sqrt(dt_sim) * rnorm(1)
    x2[i] <- x2[i - 1] - theta2 * x2[i - 1] * dt_sim + sigma2 * sqrt(dt_sim) * rnorm(1)
  }

  est <- dd_estimate_vector(x1, x2, t_int = dt_sim, Dt = 1L, dt = 1L, bins = 12L,
                             op_range_x = c(-1, 1), op_range_y = c(-1, 1))

  # F1 should depend (linearly, with slope -theta1) on x1 alone, i.e. be
  # roughly constant across x2 for a fixed x1 bin; check the marginal slope
  # across x1, averaging over x2 and well-sampled bins.
  well_sampled <- est$drift1_n > 30
  x1_bin_avg <- rowMeans(est$avg_drift1, na.rm = TRUE)
  x1_bin_n <- rowSums(est$drift1_n)
  central <- abs(est$op_x) < 0.6 & x1_bin_n > 100
  expect_gt(sum(central), 4)
  slope1 <- stats::coef(stats::lm(x1_bin_avg[central] ~ est$op_x[central]))[2]
  expect_equal(unname(slope1), -theta1, tolerance = 0.25)

  x2_bin_avg <- colMeans(est$avg_drift2, na.rm = TRUE)
  x2_bin_n <- colSums(est$drift2_n)
  central2 <- abs(est$op_y) < 0.6 & x2_bin_n > 100
  expect_gt(sum(central2), 4)
  slope2 <- stats::coef(stats::lm(x2_bin_avg[central2] ~ est$op_y[central2]))[2]
  expect_equal(unname(slope2), -theta2, tolerance = 0.25)

  # Diagonal diffusion should be close to sigma_i^2; cross term should be
  # small relative to the diagonal terms (independent noise sources).
  expect_equal(mean(est$avg_diff11[central, ], na.rm = TRUE), sigma1^2, tolerance = 0.3)
  expect_equal(mean(est$avg_diff22[, central2], na.rm = TRUE), sigma2^2, tolerance = 0.3)
  expect_lt(abs(mean(est$avg_diff12, na.rm = TRUE)), 0.3 * max(sigma1^2, sigma2^2))
})

test_that("dd_analyse_vector validates its inputs", {
  expect_error(dd_analyse_vector(1:5, 1:5, t = 1), "too short")
  expect_error(dd_analyse_vector(rnorm(20), rnorm(19), t = 1), "same length")
  expect_error(dd_analyse_vector(rnorm(20), letters[1:20], t = 1), "finite numeric")
})
