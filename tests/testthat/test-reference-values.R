# These tests check daddyR's output against reference numbers produced by
# actually running the original Python PyDaddy package (v1.1.1, patched only
# to work around its numpy>=1.24 incompatibility -- see
# DESIGN_DECISIONS.md) on the bundled scalar_pairwise dataset, with
# bins = 20, Dt = dt = 1. This is the closest thing to a ground-truth check
# available for a from-scratch reimplementation.
#
# Two kinds of quantity are compared differently:
#  - The raw (unbinned) drift series and its input x-values do not depend on
#    binning at all, so they should match Python's output to numerical
#    precision.
#  - The binned averages and the diffusion series *do* depend on binning, and
#    daddyR deliberately uses different (corrected) bin edges than PyDaddy
#    (see DESIGN_DECISIONS.md and test-estimate.R), so these are compared
#    only approximately -- close enough to confirm the algorithm is right,
#    not so close that the comparison would fail because the fix worked.

fixture_path <- function(name) testthat::test_path("fixtures", name)

test_that("raw drift series matches PyDaddy exactly (bin-independent quantity)", {
  d <- dd_load_sample_data("scalar-pairwise")
  params <- utils::read.csv(fixture_path("fixture_params.csv"), header = FALSE)
  expected_t_int <- params[params[[1]] == "t_int", 2]
  expect_equal(d$t_int, expected_t_int, tolerance = 1e-8)

  est <- dd_estimate(d$x, t_int = d$t_int, Dt = 1L, dt = 1L, bins = 20L)

  head_fixture <- utils::read.csv(fixture_path("fixture_drift_head.csv"))
  n <- nrow(head_fixture)

  expect_equal(est$x_for_drift_fit[1:n], head_fixture$x_for_fit, tolerance = 1e-8)
  expect_equal(est$drift_series[1:n], head_fixture$drift_series, tolerance = 1e-6)
})

test_that("drift polynomial fit (threshold = 0) matches PyDaddy's ridge fit", {
  d <- dd_load_sample_data("scalar-pairwise")
  est <- dd_estimate(d$x, t_int = d$t_int, Dt = 1L, dt = 1L, bins = 20L)

  fit <- dd_fit_polynomial(est$x_for_drift_fit, est$drift_series, degree = 3, threshold = 0)
  expected <- utils::read.csv(fixture_path("fixture_Fpoly.csv"))$coef

  expect_equal(fit$coeffs, expected, tolerance = 1e-4)
})

test_that("binned drift/diffusion are close to PyDaddy's (allowing for the bin-edge fix)", {
  d <- dd_load_sample_data("scalar-pairwise")
  est <- dd_estimate(d$x, t_int = d$t_int, Dt = 1L, dt = 1L, bins = 20L)

  binned <- utils::read.csv(fixture_path("fixture_binned.csv"))
  expect_equal(length(est$op), nrow(binned))
  # op ranges should agree up to half a bin width: daddyR reports the
  # *midpoint* of each evenly-tiled bin, while PyDaddy reports raw
  # np.linspace(min, max, bins) points, whose first/last values sit exactly
  # at the data range's edges rather than at a bin center. For bins = 20 over
  # range [-1, 1] (bin width 0.1), that's an expected, exact 0.05 shift at
  # each end -- not a discrepancy to chase, so the tolerance here is set just
  # above half the bin width rather than near zero.
  half_bin_width <- diff(range(binned$op)) / (nrow(binned) - 1) / 2
  expect_equal(range(est$op), range(binned$op), tolerance = half_bin_width * 1.2)
  # Binned averages should be strongly correlated and similar in magnitude;
  # exact equality is not expected (see file header comment). These
  # thresholds were set from the actual observed correlations on this
  # dataset (~0.95 for drift, ~0.98 for diffusion), with headroom below that
  # rather than an untested guess.
  expect_gt(stats::cor(est$avg_drift, binned$avgdrift, use = "complete.obs"), 0.9)
  expect_gt(stats::cor(est$avg_diff, binned$avgdiff, use = "complete.obs"), 0.9)
})

test_that("autocorrelation time is in the same ballpark as PyDaddy's", {
  d <- dd_load_sample_data("scalar-pairwise")
  params <- utils::read.csv(fixture_path("fixture_params.csv"), header = FALSE)
  expected_act <- as.numeric(params[params[[1]] == "autocorr_time", 2])

  act <- dd_autocorr_time(d$x)
  # PyDaddy's ACF estimator and R's stats::acf use different normalization
  # conventions, so this is a sanity check, not an exact-match test.
  expect_lt(abs(act - expected_act), max(10, 0.3 * expected_act))
})
