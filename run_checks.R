# Run this from anywhere -- it sets the working directory itself.
# Requires: devtools (install.packages("devtools") if you don't have it).

pkg_path <- "~/Dropbox/Mac/Desktop/PyDaddy-R/daddyR"

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Please install.packages('devtools') first.")
}

cat("== devtools::document() (regenerates NAMESPACE/man from roxygen comments) ==\n")
devtools::document(pkg_path)

cat("\n== devtools::load_all() ==\n")
devtools::load_all(pkg_path)
  
cat("\n== devtools::test() ==\n")
print(devtools::test(pkg_path))

cat("\n== devtools::check() (full R CMD check; takes a minute or two) ==\n")
devtools::check(pkg_path)
