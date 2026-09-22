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
print(dd_gaussianity_test(dd$estimate$diff_series))

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
