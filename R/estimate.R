#' Kramers-Moyal drift and diffusion estimation
#'
#' @description
#' Core numerical machinery for estimating drift and diffusion coefficients
#' from a scalar time series, following the Kramers-Moyal / Friedrich-Peinke
#' approach used throughout PyDaddy. Given a time series `x` sampled at
#' interval `t_int`, the drift at lag `Dt` is the finite-difference estimate
#' `(x[t+Dt] - x[t]) / (t_int * Dt)`, binned by the value of `x[t]`; the
#' diffusion is estimated from the squared residual of `x` after removing the
#' locally-averaged drift.
#'
#' @details
#' # Deviation from PyDaddy: consistent bin edges
#'
#' PyDaddy computes the bin width for grouping data (`inc`) as
#' `(max - min) / bins`, but separately reports the bin locations (`op`) as
#' `numpy.linspace(min, max, bins)`, whose spacing is `(max - min) / (bins -
#' 1)`. These two are different whenever `bins > 1`, which means the bins
#' used to group data internally do not exactly tile the reported `op` axis:
#' depending on rounding, this leaves narrow gaps between consecutive bins
#' (points silently excluded from every bin) or slight overlaps. daddyR fixes
#' this by constructing `bins + 1` edges that exactly span the data range,
#' using bin midpoints as the reported order-parameter values. Coverage is
#' therefore exact and every point falls in exactly one bin. This means
#' `dd_estimate()`'s binned averages will not match PyDaddy's bit-for-bit
#' (see `DESIGN_DECISIONS.md`), but the underlying unbinned drift/diffusion
#' series used for polynomial fitting are unaffected and match exactly.
#'
#' @param x Numeric vector, the (scalar) time series.
#' @param t_int Numeric scalar, time interval between observations.
#' @param Dt Integer, lag (in samples) used for the drift estimate.
#' @param dt Integer, lag (in samples) used for the diffusion estimate.
#' @param bins Integer, number of bins to use across the range of `x`
#'   (ignored if `inc` is supplied).
#' @param inc Numeric, explicit bin width to use instead of `bins`.
#' @param op_range Optional length-2 numeric vector giving the range over
#'   which to bin; defaults to `range(x, na.rm = TRUE)`.
#'
#' @return A list with components:
#'   \item{op}{Bin-center order-parameter values (length `bins`).}
#'   \item{avg_drift, avg_diff}{Binned average drift/diffusion at each `op`.}
#'   \item{drift_ebar, diff_ebar}{Standard errors of the binned averages.}
#'   \item{drift_n, diff_n}{Number of points contributing to each bin.}
#'   \item{drift_series, diff_series}{The full (unbinned) drift/diffusion
#'     series, aligned to `x_for_fit`.}
#'   \item{x_for_fit}{The `x` values (`x[1:(n - Dt)]`) that
#'     `drift_series` and `diff_series` are indexed against; use this,
#'     together with `drift_series`/`diff_series`, as the regression data for
#'     [dd_fit_polynomial()].}
#'
#' @keywords internal
dd_estimate <- function(x, t_int, Dt = 1L, dt = 1L, bins = 20L, inc = NULL,
                         op_range = NULL) {
  n <- length(x)
  stopifnot(Dt >= 1L, dt >= 1L, n > max(Dt, dt) + 1L)

  if (is.null(op_range)) {
    op_range <- range(x, na.rm = TRUE)
  }

  edges <- dd_bin_edges(op_range, bins = bins, inc = inc)
  op <- (utils::head(edges, -1) + utils::tail(edges, -1)) / 2
  n_bins <- length(op)

  # --- Drift: (x[t+Dt] - x[t]) / (t_int * Dt), indexed by x[t] ------------
  drift_series <- (x[(Dt + 1):n] - x[1:(n - Dt)]) / (t_int * Dt)
  x_drift <- x[1:(n - Dt)]

  drift_bin <- dd_bin_index(x_drift, edges)
  drift_stats <- dd_bin_stats(drift_series, drift_bin, n_bins)

  # --- Diffusion: squared residual after removing the *binned* drift -----
  # Every point x[i] has its own bin's average drift subtracted (matching
  # PyDaddy's "fast_mode" estimator, which is what its documented public
  # entry point always uses).
  x_bin_all <- dd_bin_index(x, edges)
  residual <- x - drift_stats$avg[x_bin_all] * t_int

  diff_series <- (residual[(dt + 1):n] - residual[1:(n - dt)])^2 / (t_int * dt)

  m <- n - max(Dt, dt)
  x_diff <- x[1:m]
  diff_bin <- dd_bin_index(x_diff, edges)
  diff_stats <- dd_bin_stats(diff_series[1:m], diff_bin, n_bins)

  list(
    op = op,
    edges = edges,
    avg_drift = drift_stats$avg,
    avg_diff = diff_stats$avg,
    drift_ebar = drift_stats$se,
    diff_ebar = diff_stats$se,
    drift_n = drift_stats$n,
    diff_n = diff_stats$n,
    drift_series = drift_series,
    diff_series = diff_series,
    x_for_drift_fit = x_drift,
    x_for_diff_fit = x_diff
  )
}

#' Construct bin edges spanning a range
#' @keywords internal
dd_bin_edges <- function(op_range, bins = 20L, inc = NULL) {
  lo <- op_range[1]; hi <- op_range[2]
  if (!is.null(inc)) {
    if (inc <= 0) stop("`inc` must be > 0.", call. = FALSE)
    n_bins <- max(1L, ceiling((hi - lo) / inc))
    seq(lo, by = inc, length.out = n_bins + 1)
  } else {
    if (bins < 1L) stop("`bins` must be >= 1.", call. = FALSE)
    seq(lo, hi, length.out = bins + 1)
  }
}

#' Assign values to half-open bins `[edges[i], edges[i+1])`, clamped to range
#' @keywords internal
dd_bin_index <- function(x, edges) {
  idx <- findInterval(x, edges, rightmost.closed = TRUE, all.inside = TRUE)
  idx
}

#' Per-bin mean, standard error, and count, ignoring NA/NaN/Inf
#' @keywords internal
dd_bin_stats <- function(values, bin_index, n_bins) {
  avg <- rep(NA_real_, n_bins)
  se <- rep(NA_real_, n_bins)
  cnt <- integer(n_bins)
  ok <- is.finite(values)
  split_vals <- split(values[ok], factor(bin_index[ok], levels = seq_len(n_bins)))
  for (k in seq_len(n_bins)) {
    v <- split_vals[[k]]
    cnt[k] <- length(v)
    if (cnt[k] > 0) {
      avg[k] <- mean(v)
      se[k] <- if (cnt[k] > 1) stats::sd(v) / sqrt(cnt[k]) else NA_real_
    }
  }
  list(avg = avg, se = se, n = cnt)
}
