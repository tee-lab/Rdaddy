# daddyR demo: run this in RStudio (source the whole file, or step through it
# line by line to see each result as you go).

library(daddyR)

## 1. Load a bundled example dataset -----------------------------------------
# A simulated time series from a simple pairwise-interaction model (the kind
# of synthetic example used throughout the PyDaddy paper).
d <- dd_load_sample_data("scalar-pairwise")
cat("Loaded", length(d$x), "points, t_int =", d$t_int, "\n")
cat("Range of x:", paste(round(range(d$x), 3), collapse = " to "), "\n\n")

## 2. Estimate drift and diffusion --------------------------------------------
dd <- dd_analyse(d$x, t = d$t)
print(dd)
cat("\n")

## 3. Look at the binned estimates as plain data frames -----------------------
cat("First few rows of the binned drift estimate:\n")
print(head(dd_drift(dd)))
cat("\nFirst few rows of the binned diffusion estimate:\n")
print(head(dd_diffusion(dd)))
cat("\n")

## 4. Fit analytical (polynomial) expressions ---------------------------------
dd <- dd_fit(dd, "drift", degree = 3, threshold = 0.01)
dd <- dd_fit(dd, "diffusion", degree = 2, threshold = 0.01)
cat("Fitted drift function f(x)     =", format(dd$fits$drift), "\n")
cat("Fitted diffusion function g^2(x) =", format(dd$fits$diffusion), "\n\n")

## 5. Full summary -------------------------------------------------------------
summary(dd)

## 6. Plot it (opens in the RStudio Plots pane) -------------------------------
# Requires ggplot2 + patchwork -- if this errors, everything above this line
# still works fine without them.
print(plot(dd))

## 7. Diagnostics ---------------------------------------------------------------
cat("\nIs the extracted noise consistent with Gaussian white noise?\n")
print(dd_gaussianity_test(dd$estimate$noise_series))

cat("\nAutocorrelation time of the raw series:\n")
ac <- dd_autocorrelation(d$x)
print(ac$autocorr_time)

## 8. Export the drift/diffusion table ----------------------------------------
out <- dd_export_data(dd)
cat("\nExported drift/diffusion table (first few rows):\n")
print(head(out))
# Uncomment to actually write a CSV:
# dd_export_data(dd, file = "drift_diffusion_output.csv")

## 9. Try a second bundled dataset --------------------------------------------
cat("\n\n--- Same workflow on the 'scalar-ternary' dataset ---\n")
d2 <- dd_load_sample_data("scalar-ternary")
dd2 <- dd_analyse(d2$x, t = d2$t)
dd2 <- dd_fit(dd2, "drift", degree = 3, threshold = 0.01)
dd2 <- dd_fit(dd2, "diffusion", degree = 2, threshold = 0.01)
print(dd2)
print(plot(dd2))

## 10. Simulate a synthetic trajectory from the fitted scalar model ------------
cat("\n\n--- Simulating from the fitted scalar model (dd) ---\n")
sim <- dd_simulate(dd, t_int = d$t_int, timepoints = 2000, x0 = 0)
cat("Simulated", length(sim), "points; range:",
    paste(round(range(sim), 3), collapse = " to "), "\n")

## 11. Vector (2D) workflow -----------------------------------------------------
cat("\n\n--- Vector (2D) workflow, on the 'vector-pairwise' dataset ---\n")
dv <- dd_load_sample_data("vector-pairwise")
cat("Loaded", length(dv$x1), "points (2D), t_int =", dv$t_int, "\n")

ddv <- dd_analyse_vector(dv$x1, dv$x2, t = dv$t_int)
print(ddv)
cat("\n")

cat("First few rows of the binned drift estimate:\n")
print(head(dd_drift(ddv)))
cat("\nFirst few rows of the binned diffusion estimate:\n")
print(head(dd_diffusion(ddv)))

ddv <- dd_fit(ddv, "F1", degree = 3, threshold = 0.01)
ddv <- dd_fit(ddv, "F2", degree = 3, threshold = 0.01)
ddv <- dd_fit(ddv, "G11", degree = 2, threshold = 0.01)
ddv <- dd_fit(ddv, "G22", degree = 2, threshold = 0.01)
ddv <- dd_fit(ddv, "G12", degree = 2, threshold = 0.01)
cat("\nFitted F1(x1,x2)  =", format(ddv$fits$F1), "\n")
cat("Fitted F2(x1,x2)  =", format(ddv$fits$F2), "\n")
cat("Fitted G11(x1,x2) =", format(ddv$fits$G11), "\n")
cat("Fitted G22(x1,x2) =", format(ddv$fits$G22), "\n")
cat("Fitted G12(x1,x2) =", format(ddv$fits$G12), "\n\n")

summary(ddv)

cat("\nPhase-portrait + drift/diffusion heatmaps (RStudio Plots pane):\n")
print(plot(ddv))

cat("\nGaussianity test on each noise component:\n")
print(dd_gaussianity_test_vector(ddv))

cat("\nSimulating a synthetic 2D trajectory from the fitted vector model:\n")
simv <- dd_simulate(ddv, t_int = dv$t_int, timepoints = 2000, x0 = c(0, 0))
cat("Simulated", nrow(simv), "points; x1 range:",
    paste(round(range(simv[, "x1"]), 3), collapse = " to "),
    "  x2 range:", paste(round(range(simv[, "x2"]), 3), collapse = " to "), "\n")
