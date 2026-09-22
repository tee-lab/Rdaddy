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

  noise_series <- (residual[(dt + 1):n] - residual[1:(n - dt)]) / sqrt(t_int * dt)
  diff_series <- noise_series^2

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
    noise_series = noise_series,
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

#' Per-2D-bin mean, standard error, and count, ignoring NA/NaN/Inf
#'
#' Like [dd_bin_stats()], but for a joint `(bin_x, bin_y)` grid. Returned
#' matrices are indexed `[bin_x, bin_y]` throughout (row = x-bin, column =
#' y-bin) -- see the "Deviation from PyDaddy" note in [dd_estimate_vector()].
#'
#' @keywords internal
dd_bin_stats_2d <- function(values, bin_x, bin_y, nx, ny) {
  avg <- matrix(NA_real_, nx, ny)
  se <- matrix(NA_real_, nx, ny)
  cnt <- matrix(0L, nx, ny)
  ok <- is.finite(values) & is.finite(bin_x) & is.finite(bin_y)
  key <- (bin_y[ok] - 1L) * nx + bin_x[ok]
  split_vals <- split(values[ok], factor(key, levels = seq_len(nx * ny)))
  for (k in seq_len(nx * ny)) {
    v <- split_vals[[k]]
    n_k <- length(v)
    if (n_k > 0) {
      avg[k] <- mean(v)
      se[k] <- if (n_k > 1) stats::sd(v) / sqrt(n_k) else NA_real_
      cnt[k] <- n_k
    }
  }
  list(avg = avg, se = se, n = cnt)
}

#' Kramers-Moyal drift and diffusion estimation for a bivariate (2D) series
#'
#' @description
#' The vector-data analogue of [dd_estimate()]: given two co-observed series
#' `x1`, `x2`, estimates the drift vector `(F1, F2)` and the (symmetric)
#' diffusion matrix `(G11, G22, G12)` on a joint 2D grid of bins over
#' `(x1, x2)`.
#'
#' @details
#' # Deviations from PyDaddy's vector estimator
#'
#' PyDaddy's own vector implementation (`SDE._vector_drift_diff`) has several
#' idiosyncrasies beyond the bin-edge issue already fixed in [dd_estimate()]:
#' fixed here too, for both axes independently.
#'
#' * **Closed-interval double counting.** PyDaddy's 2D bin-membership test is
#'   closed on both ends (`bin_x <= x & x <= bin_x + inc_x`), so a point
#'   sitting exactly on a shared edge between two bins is counted in *both*
#'   (its scalar binning, by contrast, correctly uses a half-open interval).
#'   daddyR uses the same half-open, exactly-tiling bins for both axes as
#'   [dd_estimate()] (via [dd_bin_index()]), so every point falls in exactly
#'   one 2D bin.
#' * **Inconsistent matrix orientation.** PyDaddy declares its binned-average
#'   matrices as shape `(nx, ny)` but then fills and reads them as
#'   `[y_bin, x_bin]` -- the transpose of what the declared shape suggests.
#'   A related diagnostic function elsewhere in PyDaddy indexes the *same*
#'   matrices the other way around (`[x_bin, y_bin]`), silently reading the
#'   wrong bin's average whenever a point's x- and y-bin indices differ.
#'   daddyR uses one explicit convention everywhere: every matrix returned
#'   here is indexed `[bin_x, bin_y]` (row = x-bin, column = y-bin).
#' * **Conflated timescales in the diffusion estimate.** PyDaddy's diffusion
#'   matrix formula lags the raw increment by the *drift* timescale `Dt`,
#'   but then divides by the *diffusion* timescale `dt` -- structurally
#'   different from (and only equivalent to) its own scalar formula when
#'   `Dt == dt`. daddyR instead extends the same residual-based estimator
#'   used by [dd_estimate()] to the vector case: subtract each point's own
#'   bin's average drift (scaled by `t_int` only, not `Dt`) to form a
#'   residual series, then lag *that* by `dt` for both the diagonal
#'   (`G11`, `G22`) and cross (`G12`) terms. This keeps `Dt` and `dt`
#'   independently meaningful, exactly as in the scalar case.
#' * **No separate "G21".** Because the diffusion matrix is symmetric by
#'   construction (`G12` and `G21` are computed from the same
#'   `residual1 * residual2` product), daddyR does not expose a redundant
#'   `G21` quantity to fit separately -- there is only `G12`, used for both
#'   off-diagonal entries wherever the full matrix is needed (e.g.
#'   [dd_simulate()]).
#'
#' @param x1,x2 Numeric vectors of equal length: the two co-observed series.
#' @param t_int Numeric scalar, time interval between observations.
#' @param Dt Integer, lag (in samples) used for the drift estimate.
#' @param dt Integer, lag (in samples) used for the diffusion estimate.
#' @param bins Integer, number of bins per axis (ignored on an axis where
#'   `inc_x`/`inc_y` is supplied).
#' @param inc_x,inc_y Optional explicit bin widths for `x1`/`x2`.
#' @param op_range_x,op_range_y Optional length-2 ranges to bin over;
#'   default to the observed range of `x1`/`x2`.
#'
#' @return A list with (among others) `op_x`, `op_y` (bin-center vectors),
#'   `avg_drift1`, `avg_drift2`, `avg_diff11`, `avg_diff22`, `avg_diff12`
#'   (all `nx` by `ny` matrices, indexed `[bin_x, bin_y]`), the unbinned
#'   `drift1_series`/`drift2_series`/`diff11_series`/`diff22_series`/
#'   `diff12_series` (for polynomial fitting), and `noise1`/`noise2` (the
#'   pre-squared residual noise series, suitable for [dd_gaussianity_test()]).
#' @keywords internal
dd_estimate_vector <- function(x1, x2, t_int, Dt = 1L, dt = 1L, bins = 20L,
                                inc_x = NULL, inc_y = NULL,
                                op_range_x = NULL, op_range_y = NULL) {
  n <- length(x1)
  stopifnot(length(x2) == n, Dt >= 1L, dt >= 1L, n > max(Dt, dt) + 1L)

  if (is.null(op_range_x)) op_range_x <- range(x1, na.rm = TRUE)
  if (is.null(op_range_y)) op_range_y <- range(x2, na.rm = TRUE)

  edges_x <- dd_bin_edges(op_range_x, bins = bins, inc = inc_x)
  edges_y <- dd_bin_edges(op_range_y, bins = bins, inc = inc_y)
  op_x <- (utils::head(edges_x, -1) + utils::tail(edges_x, -1)) / 2
  op_y <- (utils::head(edges_y, -1) + utils::tail(edges_y, -1)) / 2
  nx <- length(op_x); ny <- length(op_y)

  # --- Drift: per-component Kramers-Moyal first moment, joint-binned -------
  drift1_series <- (x1[(Dt + 1):n] - x1[1:(n - Dt)]) / (t_int * Dt)
  drift2_series <- (x2[(Dt + 1):n] - x2[1:(n - Dt)]) / (t_int * Dt)
  x1_drift <- x1[1:(n - Dt)]
  x2_drift <- x2[1:(n - Dt)]

  bx_drift <- dd_bin_index(x1_drift, edges_x)
  by_drift <- dd_bin_index(x2_drift, edges_y)
  drift1_stats <- dd_bin_stats_2d(drift1_series, bx_drift, by_drift, nx, ny)
  drift2_stats <- dd_bin_stats_2d(drift2_series, bx_drift, by_drift, nx, ny)

  # --- Diffusion matrix: residuals after removing each point's own bin's --
  # --- average drift, then dt-lagged (see "Deviations" above) -------------
  bx_all <- dd_bin_index(x1, edges_x)
  by_all <- dd_bin_index(x2, edges_y)
  res1 <- x1 - drift1_stats$avg[cbind(bx_all, by_all)] * t_int
  res2 <- x2 - drift2_stats$avg[cbind(bx_all, by_all)] * t_int

  d1 <- res1[(dt + 1):n] - res1[1:(n - dt)]
  d2 <- res2[(dt + 1):n] - res2[1:(n - dt)]
  noise1 <- d1 / sqrt(t_int * dt)
  noise2 <- d2 / sqrt(t_int * dt)
  diff11_series <- noise1^2
  diff22_series <- noise2^2
  diff12_series <- noise1 * noise2

  m <- n - max(Dt, dt)
  x1_diff <- x1[1:m]
  x2_diff <- x2[1:m]
  bx_diff <- dd_bin_index(x1_diff, edges_x)
  by_diff <- dd_bin_index(x2_diff, edges_y)
  diff11_stats <- dd_bin_stats_2d(diff11_series[1:m], bx_diff, by_diff, nx, ny)
  diff22_stats <- dd_bin_stats_2d(diff22_series[1:m], bx_diff, by_diff, nx, ny)
  diff12_stats <- dd_bin_stats_2d(diff12_series[1:m], bx_diff, by_diff, nx, ny)

  list(
    op_x = op_x, op_y = op_y, edges_x = edges_x, edges_y = edges_y,
    avg_drift1 = drift1_stats$avg, avg_drift2 = drift2_stats$avg,
    avg_diff11 = diff11_stats$avg, avg_diff22 = diff22_stats$avg,
    avg_diff12 = diff12_stats$avg,
    drift1_ebar = drift1_stats$se, drift2_ebar = drift2_stats$se,
    diff11_ebar = diff11_stats$se, diff22_ebar = diff22_stats$se,
    diff12_ebar = diff12_stats$se,
    drift1_n = drift1_stats$n, drift2_n = drift2_stats$n,
    diff11_n = diff11_stats$n, diff22_n = diff22_stats$n, diff12_n = diff12_stats$n,
    drift1_series = drift1_series, drift2_series = drift2_series,
    diff11_series = diff11_series, diff22_series = diff22_series,
    diff12_series = diff12_series,
    noise1 = noise1, noise2 = noise2,
    x1_for_drift_fit = x1_drift, x2_for_drift_fit = x2_drift,
    x1_for_diff_fit = x1_diff, x2_for_diff_fit = x2_diff
  )
}
