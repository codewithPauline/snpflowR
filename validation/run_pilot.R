# Run from the repository root: Rscript validation/run_pilot.R
# The default grid and seeds are fixed before inspecting results.
suppressPackageStartupMessages(library(snpflowR))
args <- commandArgs(trailingOnly = TRUE)
n_replicates <- if (length(args)) as.integer(args[1]) else 100L
stopifnot(length(n_replicates) == 1, is.finite(n_replicates), n_replicates > 0)
out <- "validation/results"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
scenarios <- data.frame(
  scenario = c("no_migration", "symmetric", "A_to_B", "B_to_A", "asymmetric_4to1",
               "null_unequal_sample", "null_unequal_Ne", "A_to_B_missing20", "founder_no_migration"),
  m_ab = c(0, .01, .02, 0, .02, 0, 0, .02, 0),
  m_ba = c(0, .01, 0, .02, .005, 0, 0, 0, 0),
  n_a = c(40, 40, 40, 40, 40, 10, 40, 40, 40),
  n_b = c(40, 40, 40, 40, 40, 80, 40, 40, 40),
  ne_a = c(500, 500, 500, 500, 500, 500, 100, 500, 500),
  ne_b = c(500, 500, 500, 500, 500, 500, 1000, 500, 500),
  missing = c(0, 0, 0, 0, 0, 0, 0, .2, 0), stringsAsFactors = FALSE)
write.csv(scenarios, file.path(out, "pilot-scenarios.csv"), row.names = FALSE)
rows <- vector("list", nrow(scenarios) * n_replicates)
index <- 0L
for (s in seq_len(nrow(scenarios))) {
  config <- scenarios[s, ]
  for (r in seq_len(n_replicates)) {
    seed <- 100000L + 10000L * s + r
    initial <- NULL
    if (config$scenario == "founder_no_migration") {
      set.seed(seed + 2000000L)
      ancestor <- runif(300, .05, .95)
      # B founded from 4 diploid individuals, followed by zero migration.
      initial <- rbind(ancestor, rbinom(300, 8, ancestor) / 8)
    }
    sim <- simulate_snp_flow(n_loci = 300, n_ind = c(config$n_a, config$n_b),
      effective_size = c(config$ne_a, config$ne_b), generations = 100,
      migration = matrix(c(0, config$m_ba, config$m_ab, 0), 2, 2),
      missing = config$missing, initial = initial, seed = seed)
    fit <- estimate_flow(sim$counts, n_boot = 999, resample = "locus", seed = seed + 1000000L)
    truth <- if (config$m_ab > config$m_ba) "A -> B" else if (config$m_ba > config$m_ab) "B -> A" else "null"
    p <- fit$pairs
    resolved <- p$direction != "unresolved"
    index <- index + 1L
    rows[[index]] <- data.frame(scenario = config$scenario, replicate = r,
      simulation_seed = seed, bootstrap_seed = seed + 1000000L, truth = truth,
      p, resolved = resolved,
      correct = truth != "null" & p$direction == truth,
      wrong = truth != "null" & resolved & p$direction != truth,
      false_direction = truth == "null" & resolved, row.names = NULL)
  }
  message(config$scenario, ": completed ", n_replicates, " replicates")
}
raw <- do.call(rbind, rows)
write.csv(raw, file.path(out, "pilot-replicates.csv"), row.names = FALSE, na = "NA")
wilson <- function(k, n) {
  z <- qnorm(.975); den <- 1 + z^2 / n
  center <- (k/n + z^2/(2*n)) / den
  half <- z * sqrt((k/n)*(1-k/n)/n + z^2/(4*n^2)) / den
  c(lower = center-half, upper = center+half)
}
summary <- do.call(rbind, lapply(scenarios$scenario, function(s) {
  d <- raw[raw$scenario == s, ]
  ci <- wilson(sum(d$resolved), nrow(d))
  data.frame(scenario = s, n = nrow(d), supported = sum(d$resolved),
    unresolved = sum(!d$resolved), correct = sum(d$correct), wrong = sum(d$wrong),
    false_direction = sum(d$false_direction), supported_rate = mean(d$resolved),
    supported_lower = ci[1], supported_upper = ci[2],
    median_delta = median(d$delta_gst, na.rm = TRUE), row.names = NULL)
}))
write.csv(summary, file.path(out, "pilot-summary.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "pilot-session-info.txt"))
writeLines(c("Pilot protocol v1; exploratory, not preregistered.",
  paste("Replicates per scenario:", n_replicates),
  "300 independent biallelic loci; 100 generations; 999 paired locus resamples.",
  "95% nominal intervals; default reporting safeguards; no MAF filtering.",
  "All outcomes, including failures to recover migration, must be reported.",
  "Wilson intervals describe Monte Carlo uncertainty in supported-call frequency.",
  "The bootstrap targets genetic contrasts, not known migration parameters.",
  "This pilot does not measure CI coverage for a population-genetic estimand.",
  paste("Run date:", Sys.Date())), file.path(out, "pilot-protocol.txt"))
print(summary, row.names = FALSE)
