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

## Getting the code

If you were handed a folder (e.g. `daddyR/`) directly, skip to
[Installation](#installation) below.

If you're cloning from GitHub instead:

```bash
git clone https://github.com/<your-username>/daddyR.git
```

(or, if you don't use git, just download the repository as a ZIP from
GitHub — "Code" → "Download ZIP" — and unzip it). Either way you should end
up with a folder called `daddyR` containing `DESCRIPTION`, `NAMESPACE`, an
`R/` folder, and so on.

## Prerequisites

You need R (>= 4.1) and the **devtools** package. If you don't already have
devtools:

```r
install.packages("devtools")
```

Two more packages, **ggplot2** and **patchwork**, are needed only if you
want the `plot()` method to work; everything else in the package works
without them. If you don't have them yet:

```r
install.packages(c("ggplot2", "patchwork"))
```

## Installation

Open R or RStudio, and point `devtools::install()` at the folder you just
cloned/downloaded (use the actual path on your machine):

```r
devtools::install("/path/to/daddyR")
```

For example, if the folder is on your Desktop:

```r
devtools::install("~/Desktop/daddyR")
```

This compiles the documentation and installs daddyR like any other R
package. **This step is required** — simply calling `library(daddyR)`
without installing first will fail with `there is no package called
'daddyR'`. (If you're actively developing the package and just want to try
out changes without a full install, `devtools::load_all("/path/to/daddyR")`
loads it temporarily into your current R session instead — but for
following the demo below, a real `devtools::install()` is simplest.)

Once installed, you can load it in any R session the normal way:

```r
library(daddyR)
```

## Try the demo

The repository includes a self-contained demo script, `demo.R`, that walks
through the whole workflow — loading a bundled example dataset, estimating
drift and diffusion, fitting sparse polynomial equations, plotting, and
running diagnostics — on two different example time series.

After installing the package (previous step), open `demo.R` in RStudio and
either source the whole file:

```r
source("/path/to/demo.R")
```

or, better for a first look, open it in the RStudio editor and step through
it line by line (Cmd/Ctrl+Enter) so you can see each result — printed
summaries, fitted equations, and plots — as it appears.

You don't need any of your own data to try this: `demo.R` uses the example
datasets bundled inside the package itself (via `dd_load_sample_data()`).

## Quick start (minimal version)

If you'd rather type it yourself instead of running `demo.R`:

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

## Troubleshooting

- **`there is no package called 'daddyR'`** — you called `library(daddyR)`
  without installing it first. Run `devtools::install("/path/to/daddyR")`
  (see [Installation](#installation)) and try again.
- **Repeated prompts to install ggplot2/patchwork while running
  `devtools::document()`/`load_all()`/`check()`** — these two packages are
  optional (`Suggests`, not `Imports`): everything except `plot()` works
  without them. If you don't want the plots, you can answer "no"/skip the
  prompt and continue. If you do want plots, install them once with
  `install.packages(c("ggplot2", "patchwork"))` and the prompts should stop.
- **`plot()` errors with something like `object 'is_ggplot' is not exported
  by 'namespace:ggplot2'`** — this means your installed versions of
  ggplot2 and patchwork are out of sync with each other, not a daddyR bug.
  Reinstall both together: `install.packages(c("ggplot2", "patchwork"))`.

## Development

```r
devtools::document()   # regenerate NAMESPACE / man pages from roxygen comments
devtools::test()       # run the test suite (includes checks against
                        # reference values from an actual PyDaddy run --
                        # see tests/testthat/test-reference-values.R)
devtools::check()      # full R CMD check
```

See [`DESIGN_DECISIONS.md`](DESIGN_DECISIONS.md) for the rationale behind
every place daddyR's behavior intentionally differs from PyDaddy's.

## License

GPL-3, matching the original PyDaddy package, of which this is a derivative
work.
