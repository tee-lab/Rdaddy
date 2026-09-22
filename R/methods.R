# ggplot2/patchwork are Suggests, not Imports (see DESCRIPTION), so plot.daddy()
# below only calls them via `pkg::fn()` after an explicit requireNamespace()
# check, and this file uses bare column names in aes() (standard tidy
# evaluation) rather than `.data$col`, so as not to require an importFrom
# binding to a package that might not be installed. The globalVariables()
# call below just tells R CMD check that these bare names are intentional
# (data-frame columns referenced via NSE), not typos.
utils::globalVariables(c("t", "x", "op", "avg_drift", "drift_se", "avg_diff", "diff_se",
                          "x1", "x2", "avg_drift1", "avg_drift2", "avg_diff11", "avg_diff22"))

#' @export
print.daddy <- function(x, ...) {
  if (isTRUE(x$vector)) {
    cat("<daddy> vector stochastic differential equation estimate\n")
    cat(sprintf("  %d observations (2D), t_int = %.4g\n", length(x$x1), x$t_int))
    cat(sprintf("  Dt = %d, dt = %d, bins = %d per axis (autocorrelation time ~ %d samples, from x1^2 + x2^2)\n",
                x$Dt, x$dt, x$bins, x$autocorr_time))
    for (nm in c("F1", "F2", "G11", "G22", "G12")) {
      if (!is.null(x$fits[[nm]])) {
        cat(sprintf("  %-4s fit: %s\n", nm, format(x$fits[[nm]])))
      } else {
        cat(sprintf("  %-4s fit:  (not fitted; see dd_fit(dd, \"%s\", degree = ...))\n", nm, nm))
      }
    }
    return(invisible(x))
  }

  cat("<daddy> stochastic differential equation estimate\n")
  cat(sprintf("  %d observations, t_int = %.4g\n", length(x$x), x$t_int))
  cat(sprintf("  Dt = %d, dt = %d, bins = %d (autocorrelation time ~ %d samples)\n",
              x$Dt, x$dt, x$bins, x$autocorr_time))
  if (!is.null(x$fits$drift)) {
    cat("  drift fit:     ", format(x$fits$drift), "\n")
  } else {
    cat("  drift fit:      (not fitted; see dd_fit(dd, \"drift\", degree = ...))\n")
  }
  if (!is.null(x$fits$diffusion)) {
    cat("  diffusion fit: ", format(x$fits$diffusion), "\n")
  } else {
    cat("  diffusion fit:  (not fitted; see dd_fit(dd, \"diffusion\", degree = ...))\n")
  }
  invisible(x)
}

#' @export
summary.daddy <- function(object, ...) {
  drift_df <- dd_drift(object)
  diff_df <- dd_diffusion(object)

  if (isTRUE(object$vector)) {
    cat("Summary of daddy object (vector)\n")
    cat("================================\n")
    print.daddy(object)
    cat("\nBinned drift/diffusion grid: ", nrow(drift_df), " (x1, x2) bins\n", sep = "")
    cat("  x1 range [", sprintf("%.4g, %.4g", min(drift_df$x1), max(drift_df$x1)), "]\n", sep = "")
    cat("  x2 range [", sprintf("%.4g, %.4g", min(drift_df$x2), max(drift_df$x2)), "]\n", sep = "")
    cat("\nFirst few rows of the binned drift estimate:\n")
    print(utils::head(drift_df, 5))
    cat("\nFirst few rows of the binned diffusion estimate:\n")
    print(utils::head(diff_df, 5))
    return(invisible(object))
  }

  cat("Summary of daddy object\n")
  cat("========================\n")
  print.daddy(object)
  cat("\nBinned drift (order parameter range ", sprintf("[%.4g, %.4g]", min(drift_df$op), max(drift_df$op)), "):\n", sep = "")
  print(utils::head(drift_df, 5))
  if (nrow(drift_df) > 5) cat("  ... (", nrow(drift_df) - 5, " more rows)\n", sep = "")
  cat("\nBinned diffusion:\n")
  print(utils::head(diff_df, 5))
  if (nrow(diff_df) > 5) cat("  ... (", nrow(diff_df) - 5, " more rows)\n", sep = "")
  invisible(object)
}

#' Plot a daddy object
#'
#' For a scalar object, produces a four-panel summary figure: the raw time
#' series, its histogram, the binned drift estimate (with the fitted
#' polynomial overlaid if present), and the binned diffusion estimate
#' (likewise). For a vector object, produces a phase portrait plus binned
#' heatmaps of the drift components (`F1`, `F2`) and diagonal diffusion
#' terms (`G11`, `G22`) -- fitted polynomials are not overlaid on the
#' heatmaps (there's no natural single line to draw for a 2D surface), but
#' see [dd_drift()]/[dd_diffusion()] for the underlying binned values and
#' `dd$fits` for the fitted `poly2d` objects.
#'
#' @param x A `"daddy"` object.
#' @param n_points Number of leading time-series points to show in the
#'   timeseries/phase-portrait panel (default 1000, to keep the plot legible
#'   for long series).
#' @param ... Unused; present for S3 method consistency.
#' @return A `patchwork` object (prints like a `ggplot`).
#' @export
plot.daddy <- function(x, n_points = 1000, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop("Plotting requires the 'ggplot2' and 'patchwork' packages.", call. = FALSE)
  }
  if (isTRUE(x$vector)) return(plot_daddy_vector(x, n_points = n_points))

  n_show <- min(n_points, length(x$x))
  ts_df <- data.frame(t = seq_len(n_show) * x$t_int, x = x$x[seq_len(n_show)])

  p_ts <- ggplot2::ggplot(ts_df, ggplot2::aes(x = t, y = x)) +
    ggplot2::geom_line(linewidth = 0.3) +
    ggplot2::labs(title = "Time series", x = "t", y = "x") +
    ggplot2::theme_minimal()

  p_hist <- ggplot2::ggplot(data.frame(x = x$x), ggplot2::aes(x = x)) +
    ggplot2::geom_histogram(bins = 40, fill = "grey60", color = "white") +
    ggplot2::labs(title = "Histogram", x = "x", y = "count") +
    ggplot2::theme_minimal()

  drift_df <- dd_drift(x)
  p_drift <- ggplot2::ggplot(drift_df, ggplot2::aes(x = op, y = avg_drift)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = avg_drift - drift_se,
                                        ymax = avg_drift + drift_se),
                            width = 0, color = "grey60") +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::labs(title = "Drift", x = "x", y = expression(f(x))) +
    ggplot2::theme_minimal()
  if (!is.null(x$fits$drift)) {
    xs <- seq(min(drift_df$op), max(drift_df$op), length.out = 200)
    p_drift <- p_drift + ggplot2::geom_line(
      data = data.frame(op = xs, avg_drift = stats::predict(x$fits$drift, xs)),
      color = "steelblue", linewidth = 0.8
    )
  }

  diff_df <- dd_diffusion(x)
  p_diff <- ggplot2::ggplot(diff_df, ggplot2::aes(x = op, y = avg_diff)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = avg_diff - diff_se,
                                        ymax = avg_diff + diff_se),
                            width = 0, color = "grey60") +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::labs(title = "Diffusion", x = "x", y = expression(g^2*(x))) +
    ggplot2::theme_minimal()
  if (!is.null(x$fits$diffusion)) {
    xs <- seq(min(diff_df$op), max(diff_df$op), length.out = 200)
    p_diff <- p_diff + ggplot2::geom_line(
      data = data.frame(op = xs, avg_diff = stats::predict(x$fits$diffusion, xs)),
      color = "firebrick", linewidth = 0.8
    )
  }

  (p_ts + p_hist) / (p_drift + p_diff)
}

#' @keywords internal
plot_daddy_vector <- function(x, n_points = 1000) {
  n_show <- min(n_points, length(x$x1))
  traj_df <- data.frame(x1 = x$x1[seq_len(n_show)], x2 = x$x2[seq_len(n_show)])

  p_traj <- ggplot2::ggplot(traj_df, ggplot2::aes(x = x1, y = x2)) +
    ggplot2::geom_path(linewidth = 0.2, alpha = 0.6) +
    ggplot2::labs(title = "Phase portrait", x = "x1", y = "x2") +
    ggplot2::theme_minimal()

  drift_df <- dd_drift(x)
  diff_df <- dd_diffusion(x)

  p_f1 <- ggplot2::ggplot(drift_df, ggplot2::aes(x = x1, y = x2, fill = avg_drift1)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2() +
    ggplot2::labs(title = "F1 (drift, x1 component)", x = "x1", y = "x2", fill = NULL) +
    ggplot2::theme_minimal()

  p_f2 <- ggplot2::ggplot(drift_df, ggplot2::aes(x = x1, y = x2, fill = avg_drift2)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2() +
    ggplot2::labs(title = "F2 (drift, x2 component)", x = "x1", y = "x2", fill = NULL) +
    ggplot2::theme_minimal()

  p_g11 <- ggplot2::ggplot(diff_df, ggplot2::aes(x = x1, y = x2, fill = avg_diff11)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(title = "G11 (diffusion, x1)", x = "x1", y = "x2", fill = NULL) +
    ggplot2::theme_minimal()

  p_g22 <- ggplot2::ggplot(diff_df, ggplot2::aes(x = x1, y = x2, fill = avg_diff22)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(title = "G22 (diffusion, x2)", x = "x1", y = "x2", fill = NULL) +
    ggplot2::theme_minimal()

  (p_traj + p_f1) / (p_f2 + p_g11) / (p_g22 + patchwork::plot_spacer())
}
