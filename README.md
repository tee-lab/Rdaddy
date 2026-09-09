# daddyR

A native R port of [PyDaddy](https://github.com/tee-lab/PyDaddy), the
Python package accompanying Nabeel, Karichannavar, Palathingal, Jhawar,
Brückner, Danny Raj & Guttal (2025), "Discovering stochastic dynamical
equations from ecological time series data," *The American Naturalist*,
205: E100-E117 (<https://doi.org/10.1086/734083>).

daddyR estimates data-driven stochastic differential equations from time
series data: given a time series `x(t)`, it recovers the drift and
diffusion functions of

<p align="center"><code>dx/dt = f(x) + g(x)&middot;&eta;(t)</code></p>

via the Kramers-Moyal / Friedrich-Peinke binned-coefficient method, and fits
sparse polynomial expressions to `f` and `g` via sequentially-thresholded
ridge regression.

This is a from-scratch reimplementation, not a translation or a wrapper
around the Python package — see [`DESIGN_DECISIONS.md`](DESIGN_DECISIONS.md)
for a full account of where and why daddyR's behavior deviates from
PyDaddy's (several of which are bug fixes; PyDaddy's current release does
not run on numpy released since 2021 without patching, for instance).

**Current scope (v0.1.0):** scalar (1D) time series only. Vector (2D) data,
cross-diffusion, SDE simulation, and the full diagnostic suite are planned
for later releases — see the project roadmap.

## Installation

```r
# From a local checkout:
devtools::install(".")
```

## Quick start

```r
library(daddyR)

d <- dd_load_sample_data("scalar-pairwise")
dd <- dd_analyse(d$x, t = d$t)
print(dd)

dd_drift(dd)
dd_diffusion(dd)

dd <- dd_fit(dd, "drift", degree = 3, threshold = 0.01)
dd <- dd_fit(dd, "diffusion", degree = 2, threshold = 0.01)
print(dd)
plot(dd)

dd_gaussianity_test(dd$estimate$diff_series)
dd_autocorrelation(d$x)
```

## Development

```r
devtools::document()   # regenerate NAMESPACE / man pages from roxygen comments
devtools::test()       # run the test suite (includes checks against
                        # reference values from an actual PyDaddy run —
                        # see tests/testthat/test-reference-values.R)
devtools::check()      # full R CMD check
```

## License

GPL-3, matching the original PyDaddy package, of which this is a derivative
work.
