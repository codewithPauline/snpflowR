# Run after generate_msprime.py, from repository root with snpflowR installed.
library(snpflowR)
manifest <- read.csv("validation/results/msprime-manifest.csv", stringsAsFactors = FALSE)
rows <- vector("list", nrow(manifest))
for (k in seq_len(nrow(manifest))) {
  info <- manifest[k, ]
  d <- read.csv(info$file, stringsAsFactors = FALSE)
  a <- rbind(A = d$alt_a, B = d$alt_b)
  n <- rbind(A = d$called_a, B = d$called_b)
  colnames(a) <- colnames(n) <- d$locus
  x <- snp_counts(a, n, blocks = d$block)
  fit <- estimate_flow(x, n_boot = 999, resample = "block", seed = info$ancestry_seed + 4000000L)
  error <- max(abs(c(fit$pairs$gst_a_to_b - info$reference_gst_ab,
                     fit$pairs$gst_b_to_a - info$reference_gst_ba)))
  if (!is.finite(error) || error > 1e-12) stop("Independent arithmetic mismatch in ", info$file)
  p <- fit$pairs
  resolved <- p$direction != "unresolved"
  rows[[k]] <- data.frame(scenario = info$scenario, replicate = info$replicate,
    truth = info$truth, p, reference_max_error = error,
    correct = info$truth != "null" & p$direction == info$truth,
    wrong = info$truth != "null" & resolved & p$direction != info$truth,
    false_direction = info$truth == "null" & resolved, resolved = resolved,
    row.names = NULL)
}
raw <- do.call(rbind, rows)
write.csv(raw, "validation/results/msprime-replicates.csv", row.names = FALSE)
summary <- do.call(rbind, lapply(unique(raw$scenario), function(s) {
  d <- raw[raw$scenario == s, ]
  data.frame(scenario = s, n = nrow(d), supported = sum(d$resolved),
    unresolved = sum(!d$resolved), correct = sum(d$correct), wrong = sum(d$wrong),
    false_direction = sum(d$false_direction), max_reference_error = max(d$reference_max_error))
}))
write.csv(summary, "validation/results/msprime-summary.csv", row.names = FALSE)
capture.output(sessionInfo(), file = "validation/results/msprime-r-session-info.txt")
print(summary, row.names = FALSE)
