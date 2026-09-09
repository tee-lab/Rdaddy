#' Validate a scalar time series input
#' @keywords internal
dd_validate_series <- function(x, t) {
  if (!is.numeric(x)) stop("`data` must be a numeric vector.", call. = FALSE)
  if (any(is.infinite(x))) stop("`data` must not contain Inf/-Inf values.", call. = FALSE)
  if (length(x) < 10) stop("`data` is too short (", length(x), " points) to analyse.", call. = FALSE)

  if (length(t) > 1) {
    if (length(t) != length(x)) {
      stop("If `t` is a vector of timestamps, it must be the same length as `data`.", call. = FALSE)
    }
    t_int <- (t[length(t)] - t[1]) / (length(t) - 1)
  } else {
    if (!is.numeric(t) || t <= 0) stop("`t` must be a positive number, or a vector of timestamps.", call. = FALSE)
    t_int <- t
  }
  t_int
}

#' Estimate the autocorrelation time of a series
#'
#' Defined as the lag at which the autocorrelation function first drops to
#' `1/e`, matching the convention used throughout PyDaddy. Falls back to the
#' full lag range if the ACF never drops that low.
#'
#' @param x Numeric vector.
#' @param max_lag Maximum lag to search out to.
#' @return Integer autocorrelation time (in samples).
#' @export
dd_autocorr_time <- function(x, max_lag = min(1000L, length(x) - 1L)) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(1L)
  max_lag <- max(1L, min(max_lag, length(x) - 1L))
  acf_vals <- stats::acf(x, lag.max = max_lag, plot = FALSE, na.action = stats::na.pass)$acf[, 1, 1]
  below <- which(acf_vals < exp(-1))
  if (length(below) == 0) return(max_lag)
  max(1L, below[1] - 1L)
}

#' Auto-select a drift/diffusion timescale (Dt) from the autocorrelation time
#'
#' @details
#' # Deviation from PyDaddy
#' PyDaddy's internal `Main` class computes exactly this (`Dt =
#' ceiling(autocorrelation_time / 10)`) when `Dt` is left unset — but the
#' public, documented entry point (`pydaddy.Characterize()`) hardcodes
#' `Dt = 1` as its default, so in practice this auto-selection is unreachable
#' from the documented API. Since the auto-selection logic clearly exists and
#' is exercised elsewhere in the code (the interactive timescale-slider
#' feature), daddyR treats it as the intended default: `dd_analyse()` uses
#' `Dt = NULL` (meaning "auto") by default. Pass an explicit `Dt` to opt out.
#'
#' @param x Numeric vector (the time series, or its squared magnitude for
#'   vector data).
#' @keywords internal
dd_auto_dt <- function(x) {
  max(1L, ceiling(dd_autocorr_time(x) / 10))
}

#' Auto-select a bin count via the Freedman-Diaconis rule
#'
#' PyDaddy defines a Freedman-Diaconis auto-bin-count method
#' (`Main._autobins()`) but never calls it anywhere — the actual default is
#' an unconditional, data-independent 20 bins. daddyR keeps 20 as the
#' explicit default (so results match what PyDaddy actually produces by
#' default), but exposes this as an opt-in via `bins = "auto"`.
#'
#' @param x Numeric vector.
#' @keywords internal
dd_auto_bins <- function(x) {
  grDevices::nclass.FD(x[is.finite(x)])
}
