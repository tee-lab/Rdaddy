test_that("dd_analyse_vector builds a well-formed vector daddy object", {
  d <- dd_load_sample_data("vector-pairwise")
  dd <- dd_analyse_vector(d$x1, d$x2, t = d$t_int, Dt = 1L, dt = 1L, bins = 15L)

  expect_s3_class(dd, "daddy")
  expect_true(dd$vector)
  expect_equal(dd$t_int, d$t_int)
  expect_length(dd$estimate$op_x, 15)
  expect_length(dd$estimate$op_y, 15)
  expect_null(dd$fits$F1)

  # print()/summary() should not error, with or without fits set
  expect_output(print(dd), "vector stochastic differential equation")
  expect_output(summary(dd), "Summary of daddy object \\(vector\\)")
})

test_that("dd_fit dispatches to the 2D path for a vector daddy object", {
  d <- dd_load_sample_data("vector-pairwise")
  dd <- dd_analyse_vector(d$x1, d$x2, t = d$t_int, Dt = 1L, dt = 1L, bins = 15L)

  dd <- dd_fit(dd, "F1", degree = 2, threshold = 0)
  dd <- dd_fit(dd, "F2", degree = 2, threshold = 0)
  dd <- dd_fit(dd, "G11", degree = 1, threshold = 0)
  dd <- dd_fit(dd, "G22", degree = 1, threshold = 0)
  dd <- dd_fit(dd, "G12", degree = 1, threshold = 0)

  expect_s3_class(dd$fits$F1, "poly2d")
  expect_s3_class(dd$fits$G12, "poly2d")
  expect_error(dd_fit(dd, "drift", degree = 2), "should be one of")
})

test_that("dd_drift/dd_diffusion/dd_export_data dispatch to vector-shaped output", {
  d <- dd_load_sample_data("vector-pairwise")
  dd <- dd_analyse_vector(d$x1, d$x2, t = d$t_int, Dt = 1L, dt = 1L, bins = 10L)

  drift_df <- dd_drift(dd)
  expect_equal(nrow(drift_df), 100L)  # 10 x 10 grid
  expect_named(drift_df, c("x1", "x2", "avg_drift1", "avg_drift2", "n"))

  diff_df <- dd_diffusion(dd)
  expect_equal(nrow(diff_df), 100L)
  expect_named(diff_df, c("x1", "x2", "avg_diff11", "avg_diff22", "avg_diff12", "n"))

  out <- dd_export_data(dd)
  expect_equal(nrow(out), 100L)
  expect_true(all(c("avg_drift1", "avg_drift2", "avg_diff11", "avg_diff22", "avg_diff12") %in% names(out)))
})

test_that("dd_gaussianity_test_vector runs on both noise components", {
  d <- dd_load_sample_data("vector-pairwise")
  dd <- dd_analyse_vector(d$x1, d$x2, t = d$t_int, Dt = 1L, dt = 1L, bins = 15L)
  res <- dd_gaussianity_test_vector(dd)
  expect_named(res, c("x1", "x2"))
  expect_true(is.numeric(res$x1$p_value))
  expect_true(is.numeric(res$x2$p_value))
})
