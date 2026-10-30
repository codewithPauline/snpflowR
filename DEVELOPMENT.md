# Development and verification

Run from a clean checkout. R >= 4.1 is declared; GitHub Actions checks current
release/devel and oldrel-1 on the listed platforms. The first local environment
uses R 4.3.3; this alone does not establish current-platform compatibility.

```r
install.packages(c("roxygen2", "testthat", "knitr", "rmarkdown"))
roxygen2::roxygenise()
testthat::test_local()
```

```sh
R CMD build .
R CMD check --as-cran snpflowR_0.0.0.9000.tar.gz
```

Build in the parent directory or move the tarball outside the checkout before
rebuilding. Package source excludes validation results and repository-only docs.
Pandoc is needed to build the HTML vignette, and TeX to check the PDF manual.
GitHub Actions follows the R-maintained action workflow conventions, runs examples,
tests, and the vignette, and treats warnings as failures. External submission
checks and manual review are still necessary before CRAN submission.

`R/` contains maintained roxygen documentation; `man/` and `NAMESPACE` are
generated and committed. All public functions have examples and help pages.
Use qualified package calls and avoid runtime downloads, global options,
implicit plotting, and uncontrolled parallel workers.

`validation/` contains reproducible experiments, including failure scenarios.
Run long simulations explicitly, outside ordinary `R CMD check`. Generated
intermediate counts are git-ignored; seeds, scripts, protocols, and every
replicate's outcome are retained. Do not replace negative results with a
handpicked successful example.

The optional [comparator audit](validation/COMPARATOR_AUDIT.md) executes pinned,
unmodified estimator sources after checksum verification. Install Rcpp and a C++
toolchain, then run `Rscript validation/run_comparator_audit.R` from the root.
Python 3 downloads and verifies the small external source cache. This is separate
from package dependencies and ordinary checks; the runtime package stays offline.
Do not refresh reference-output fixtures or widen tolerances merely to make a
failing comparison pass. Investigate and document the discrepancy first.

The prototype was developed with AI assistance. The human maintainer remains
responsible for checking authorship, methods, scientific claims, and compliance
with any journal's AI disclosure policy. No independent expert review is claimed.
