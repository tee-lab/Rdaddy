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
