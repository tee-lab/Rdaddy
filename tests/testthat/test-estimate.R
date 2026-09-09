test_that("dd_bin_edges + dd_bin_index give full, non-overlapping bin coverage", {
  # This is a regression test for a real bug in PyDaddy: it computes the bin
  # width for grouping data as (max-min)/bins, but reports bin locations as
  # np.linspace(min, max, bins), whose spacing is (max-min)/(bins-1). These
  # differ whenever bins > 1, leaving small gaps (or overlaps) between
  # consecutive bins. daddyR's edges are constructed to exactly tile the
  # range, so every finite point must land in exactly one bin.
  set.seed(1)
  x <- runif(5000, -1, 1)
  edges <- dd_bin_edges(c(-1, 1), bins = 20L)
  expect_length(edges, 21)
  expect_equal(edges[1], -1)
  expect_equal(edges[21], 1)

  idx <- dd_bin_index(x, edges)
  expect_true(all(idx >= 1 & idx <= 20))
  expect_false(anyNA(idx))

  # Bin widths must all be equal (evenly tiled), unlike PyDaddy's mismatched
  # inc vs. linspace-derived spacing.
  widths <- diff(edges)
  expect_equal(widths, rep(widths[1], 20), tolerance = 1e-12)
})

test_that("dd_estimate recovers known drift/diffusion for a synthetic OU process", {
  # Simulate dx = -theta*x*dt + sigma*dW via Euler-Maruyama, then check that
  # dd_estimate's binned drift/diffusion approximate the known theta, sigma^2.
  set.seed(42)
  theta <- 0.7; sigma <- 0.3; dt_sim <- 0.01; n <- 200000
  x <- numeric(n)
  for (i in 2:n) {
    x[i] <- x[i - 1] - theta * x[i - 1] * dt_sim + sigma * sqrt(dt_sim) * rnorm(1)
  }

  est <- dd_estimate(x, t_int = dt_sim, Dt = 1L, dt = 1L, bins = 15L,
                      op_range = c(-1, 1))

  # Drift should be close to -theta * op across the central, well-sampled bins.
  central <- abs(est$op) < 0.6 & est$drift_n > 50
  expect_gt(sum(central), 5)
  slope <- stats::coef(stats::lm(est$avg_drift[central] ~ est$op[central]))[2]
  expect_equal(unname(slope), -theta, tolerance = 0.15)

  # Diffusion should be close to sigma^2 across the same bins.
  expect_equal(mean(est$avg_diff[central], na.rm = TRUE), sigma^2, tolerance = 0.15)
})

test_that("dd_analyse rejects series that are too short to analyse", {
  expect_error(dd_analyse(1:5, t = 1), "too short")
})
