# Provenance of inst/extdata/scalar_pairwise.csv and scalar_ternary.csv:
#
# These are copied unmodified from PyDaddy's own bundled sample data
# (pydaddy/data/model_data/scalar/pairwise.csv and ternary.csv in
# https://github.com/tee-lab/PyDaddy), which are simulated time series from
# simple pairwise- and ternary-interaction collective-motion models
# (Jhawar & Guttal, Phil. Trans. R. Soc. B, 2020), used as synthetic examples
# throughout the PyDaddy paper. Each file has two columns: the state
# variable and its timestamp, with no header row.
#
# No R data-preparation is currently needed for them; they're read directly
# by dd_load_sample_data() via system.file("extdata", ..., package = "daddyR").
# This script is kept as a record of provenance and as a place to add
# proper data-raw processing (e.g. converting to lazy-loaded .rda datasets)
# if/when that becomes worthwhile.
#
# Provenance of inst/extdata/vector_pairwise.csv and vector_ternary.csv:
#
# PyDaddy's own bundled vector sample data (pydaddy/data/model_data/vector/
# {pairwise,ternary}.csv) are ~600,000-row simulated trajectories from 2D
# pairwise- and ternary-interaction collective-motion models (Jhawar &
# Guttal, 2020), observed at a fixed t_int = 0.12 (see
# pydaddy.load_sample_dataset()'s source: vector datasets return a fixed
# `t = 0.12`, unlike the scalar datasets, which carry a timestamp column).
# At ~30 MB per file (18-digit-precision floats), these are far too large to
# bundle directly in an R package. inst/extdata/vector_{pairwise,ternary}.csv
# are a contiguous 20,000-point prefix of each, rounded to 6 decimal places,
# generated with:
#
#   for name in ["pairwise", "ternary"]:
#       src = f".../PyDaddy/pydaddy/data/model_data/vector/{name}.csv"
#       n = 20000
#       rows = []
#       with open(src) as f:
#           for i, line in enumerate(f):
#               if i >= n: break
#               x, y = (float(v) for v in line.strip().split())
#               rows.append((round(x, 6), round(y, 6)))
#       with open(f"inst/extdata/vector_{name}.csv", "w") as f:
#           for x, y in rows:
#               f.write(f"{x},{y}\n")
#
# A contiguous prefix (rather than a random subsample) preserves the
# series' temporal structure, which the drift/diffusion estimator and
# autocorrelation-time calculation both depend on. `t_int = 0.12` is
# hardcoded in `dd_load_sample_data()` for these datasets, matching
# PyDaddy's own convention, since there is no timestamp column to read it
# from.
