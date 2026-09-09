#' Autocorrelation function of a time series
#'
#' A thin, documented wrapper around [stats::acf()] returning a tidy
#' data frame, together with the estimated autocorrelation time (the lag at
#' which the ACF first drops below `1/e`).
#'
#' @param x Numeric vector.
#' @param lags Maximum lag to compute out to.
#' @return A list with `acf` (a data frame of `lag` and `acf`) and
#'   `autocorr_time` (integer).
#' @export
dd_autocorrelation <- function(x, lags = min(1000L, length(x) - 1L)) {
  x <- x[is.finite(x)]
  lags <- max(1L, min(lags, length(x) - 1L))
  acf_obj <- stats::acf(x, lag.max = lags, plot = FALSE)
  list(
    acf = data.frame(lag = 0:lags, acf = as.vector(acf_obj$acf)),
    autocorr_time = dd_autocorr_time(x, max_lag = lags)
  )
}

#' Test whether extracted noise is Gaussian
#'
#' Extracts the residual noise (the diffusion-series input, standardized) at
#' one bin of the order parameter and tests it for normality.
#'
#' @details
#' # Deviation from PyDaddy
#' PyDaddy's `GaussianTest` builds a null distribution by drawing 1,000 pairs
#' of standard-normal samples, computing "KL divergence" between each pair as
#' `sum(p * log(abs((p + eps) / (q + eps))))` applied elementwise to the raw
#' samples themselves (rather than to density estimates of `p` and `q`), and
#' checking whether the same quantity computed between the real residual and
#' a fresh normal sample falls within the 2.5-97.5 percentile band of that
#' null. This is not actually an estimate of KL divergence between
#' distributions — real KL divergence requires a density or probability mass
#' function, not raw samples plugged into the divergence formula — so the
#' resulting test statistic doesn't have a clean interpretation, and (as a
#' secondary issue) the reference implementation also calls
#' `numpy.histogram(..., normed=True)`, an argument numpy removed years ago,
#' so the function raises an error on any numpy released since 2020.
#'
#' daddyR replaces this with a standard, well-calibrated normality test: a
#' one-sample Kolmogorov-Smirnov test of the standardized residuals against
#' N(0, 1). It answers the same scientific question (is the extracted noise
#' Gaussian?) using an established procedure instead of a bespoke, non-
#' standard statistic.
#'
#' @param residual Numeric vector of residuals (e.g. `dd_estimate()`'s
#'   `diff_series`, or a subset of it at one order-parameter bin).
#' @param alpha Significance level for the test (default 0.05).
#' @return A list with `statistic`, `p_value`, and `is_gaussian` (logical,
#'   `p_value > alpha`).
#' @export
dd_gaussianity_test <- function(residual, alpha = 0.05) {
  residual <- residual[is.finite(residual)]
  if (length(residual) < 8) {
    stop("Need at least 8 finite residual values to test for Gaussianity; got ",
         length(residual), ".", call. = FALSE)
  }
  z <- scale(residual)[, 1]
  test <- suppressWarnings(stats::ks.test(z, "pnorm"))
  list(
    statistic = unname(test$statistic),
    p_value = test$p.value,
    is_gaussian = test$p.value > alpha,
    n = length(residual)
  )
}
