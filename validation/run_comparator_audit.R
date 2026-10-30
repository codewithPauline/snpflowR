# Run from the repository root after python3 validation/fetch_comparators.py.
# Executes unmodified upstream point-estimator sources in isolated environments.
# Requires an installed snpflowR, Rcpp, Python 3, and a C++ toolchain.
stopifnot(requireNamespace("snpflowR", quietly = TRUE),
          requireNamespace("Rcpp", quietly = TRUE))
python <- Sys.getenv("SNPFLOW_PYTHON", "python3")
if (system2(python, "validation/fetch_comparators.py") != 0L) {
  stop("Comparator source verification failed; no external source executed.")
}
source("validation/comparator_fixtures.R")
root <- "validation/generated/comparators"
dir.create("validation/generated/comparator-inputs", recursive = TRUE, showWarnings = FALSE)
div <- new.env(parent = globalenv())
cor <- new.env(parent = globalenv())
for (f in c("R/rgp.R", "R/divMigrate.R")) {
  sys.source(file.path(root, "diveRsity", f), envir = div)
}
for (f in c("src/myTab.cpp", "src/pwHt.cpp")) {
  Rcpp::sourceCpp(file.path(root, "diveRsity", f), env = div)
}
for (f in c("R/popcounts.R", "R/genetic_stats.R")) {
  sys.source(file.path(root, "corridoR", f), envir = cor)
}

capture <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = unique(warnings))
}
same <- function(x, y, tolerance = 1e-12) {
  (is.na(x) & is.na(y)) | (is.finite(x) & is.finite(y) & abs(x - y) <= tolerance) |
    (is.infinite(x) & is.infinite(y) & sign(x) == sign(y))
}
kind <- function(x) ifelse(is.nan(x), "not_a_number", ifelse(is.na(x), "missing",
                              ifelse(is.infinite(x), "infinite", "finite")))
records <- list()
warning_rows <- list()
inputs <- list()
fixtures <- comparator_fixtures()
for (id in names(fixtures)) {
  x <- fixtures[[id]]
  g <- fixture_genotypes(x)
  roundtrip <- snpflowR::genotype_counts(g$genotypes, g$populations)
  stopifnot(isTRUE(all.equal(unname(x$alt), unname(roundtrip$alt), tolerance = 0)),
            isTRUE(all.equal(unname(x$called), unname(roundtrip$called), tolerance = 0)))
  path <- file.path("validation/generated/comparator-inputs", paste0(id, ".gen"))
  write_fixture_genepop(g, path)
  fit <- snpflowR::estimate_flow(x)
  n <- nrow(x$alt)
  pops <- rownames(x$alt)
  dat <- div$rgp(path)
  # Verify the real upstream GENEPOP parser, including monomorphic and missing
  # loci, against the input frequencies before comparing statistical outputs.
  for (l in seq_len(ncol(x$alt))) {
    a <- dat$af[[l]]
    total <- colSums(a)
    alt <- if ("02" %in% rownames(a)) a["02", ] else rep(0, n)
    expected <- x$alt[, l] / x$called[, l]
    stopifnot(all(same(alt[total > 0], expected[total > 0])),
              identical(total == 0, unname(x$called[, l] == 0)))
    a[, total == 0] <- NA_real_
    dat$af[[l]] <- a
  }
  # divMigrate returns normalized scores, not raw GST. For raw GST only, call
  # its unmodified C++ heterozygosity helper and apply (mean(Ht)-mean(Hs))/mean(Ht).
  # This separate adapter is not described as an extra public divMigrate output.
  hh <- lapply(dat$af, div$pwHt, pw = utils::combn(n, 2) - 1L)
  mean_h <- function(field) apply(array(unlist(lapply(hh, `[[`, field)),
                                       c(n, n, length(hh))), c(1, 2), mean, na.rm = TRUE)
  ht <- mean_h("ht")
  div_gst <- (ht - mean_h("hs")) / ht
  diag(div_gst) <- 0
  d <- capture(div$divMigrate(path, stat = "gst", boots = 0, plot_network = FALSE))
  cc <- lapply(seq_len(ncol(x$alt)), function(l) {
    cbind(ref = x$called[, l] - x$alt[, l], alt = x$alt[, l])
  })
  pc <- cor$new_popcounts(cc, pops, as.integer(table(g$populations)[pops]))
  cor_result <- capture(cor$directional_gst(pc))
  outputs <- list(snpflowR = list(gst = fit$directional_gst, relative = fit$relative_migration),
                  diveRsity = list(gst = div_gst, relative = d$value$gRelMig),
                  corridoR = list(gst = cor_result$value$gst, relative = cor_result$value$relative))
  # Explicit directed entries avoid any matrix flattening/orientation ambiguity.
  for (software in names(outputs)) {
    for (i in seq_len(n)) for (j in setdiff(seq_len(n), i)) {
      v <- outputs[[software]]
      records[[length(records) + 1L]] <- data.frame(
        fixture = id, software = software, from = pops[i], to = pops[j],
        gst = v$gst[i, j], relative = v$relative[i, j],
        gst_kind = kind(v$gst[i, j]), relative_kind = kind(v$relative[i, j]),
        gst_abs_difference = abs(v$gst[i, j] - fit$directional_gst[i, j]),
        relative_abs_difference = abs(v$relative[i, j] - fit$relative_migration[i, j]),
        gst_agrees = same(v$gst[i, j], fit$directional_gst[i, j]),
        relative_agrees = same(v$relative[i, j], fit$relative_migration[i, j]),
        n_no_shared = fit$pairs$n_no_shared[match(paste(sort(c(pops[i], pops[j])), collapse = ":"),
          paste(fit$pairs$population_a, fit$pairs$population_b, sep = ":"))],
        stringsAsFactors = FALSE)
    }
  }
  for (software in c("diveRsity", "corridoR")) {
    w <- if (software == "diveRsity") d$warnings else cor_result$warnings
    if (length(w)) warning_rows[[length(warning_rows) + 1L]] <- data.frame(
      fixture = id, software = software, warning = w)
  }
  inputs[[id]] <- data.frame(fixture = id, populations = n, loci = ncol(x$alt),
                            seed = if (startsWith(id, "random_")) as.integer(sub("random_", "", id)) else NA,
                            genepop_md5 = unname(tools::md5sum(path)),
                            snpflow_scale_status = fit$scale_status)
}
results <- do.call(rbind, records)
warnings <- if (length(warning_rows)) do.call(rbind, warning_rows) else
  data.frame(fixture = character(), software = character(), warning = character())
write.csv(results, "validation/results/comparator-values.csv", row.names = FALSE, na = "NA")
write.csv(do.call(rbind, inputs), "validation/results/comparator-fixtures.csv", row.names = FALSE, na = "NA")
write.csv(warnings, "validation/results/comparator-warnings.csv", row.names = FALSE)
writeLines(sub("[[:blank:]]+$", "", capture.output(sessionInfo())),
           "validation/results/comparator-sessionInfo.txt")
summary <- do.call(rbind, lapply(split(results, results$software), function(z) {
  finite_g <- z$gst_abs_difference[is.finite(z$gst_abs_difference)]
  finite_r <- z$relative_abs_difference[is.finite(z$relative_abs_difference)]
  data.frame(software = z$software[1], fixtures = length(unique(z$fixture)),
             directed_entries = nrow(z), gst_disagreements = sum(!z$gst_agrees),
             relative_disagreements = sum(!z$relative_agrees),
             max_finite_gst_abs_difference = max(finite_g),
             max_finite_relative_abs_difference = max(finite_r))
}))
write.csv(summary, "validation/results/comparator-summary.csv", row.names = FALSE)
print(summary, row.names = FALSE)

# The initial run exposed both policy differences and floating-point boundary
# failures. See COMPARATOR_AUDIT.md. Preserve all outputs; require agreement in
# the ordinary domain, but do not require upstream roundoff bugs to reproduce
# identically under a different BLAS/compiler. Boundary tests independently
# protect snpflowR's declared behavior. The absolute tolerance stays 1e-12.
boundary <- results$fixture %in% c("identical_pair", "identical_all",
                                 "fixed_source_zero_diversity")
stopifnot(all(results$gst_agrees[!boundary]),
          all(results$relative_agrees[!boundary]),
          all(results$gst_agrees[results$software == "diveRsity"]))
different_gst <- results[!results$gst_agrees, ]
stopifnot(all(different_gst$software == "corridoR"),
          all(different_gst$fixture == "fixed_source_zero_diversity"),
          all(different_gst$from == "A"), all(different_gst$to == "B"))
# Even in boundary fixtures, all finite reference GST values must agree.
finite_diffs <- results$gst_abs_difference[is.finite(results$gst_abs_difference)]
stopifnot(all(finite_diffs <= 1e-12))
write.csv(results[!results$gst_agrees | !results$relative_agrees, ],
          "validation/results/comparator-discrepancies.csv", row.names = FALSE, na = "NA")
cat("PASS: ordinary-domain agreement; boundary differences recorded separately.\n")
