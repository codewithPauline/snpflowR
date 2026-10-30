# Research helpers. Run from the repository root, never sourced by the package.
read_uncertainty_protocol <- function() {
  jsonlite::fromJSON("validation/uncertainty/protocol-v1.json")
}

uncertainty_truth <- function(p, q) {
  # Independent two-allele heterozygosity oracle, not the production algebra.
  g <- vapply(list(p, q), function(a) {
    r <- sqrt(p * q) / (sqrt(p * q) + sqrt((1 - p) * (1 - q)))
    hs <- 1 - (a^2 + (1-a)^2 + r^2 + (1-r)^2) / 2
    average <- (a + r) / 2
    ht <- 1 - average^2 - (1-average)^2
    sum(ht - hs) / sum(ht)
  }, numeric(1))
  unname(g[2] - g[1])
}

uncertainty_generate <- function(s, protocol, seed) {
  set.seed(seed)
  p <- protocol$p
  q <- switch(s$profile, complement = 1-p, identical = p,
              alternative = protocol$q_alternative, reversed = protocol$q_alternative)
  if (s$profile == "reversed") { old <- p; p <- q; q <- old }
  null <- s$profile %in% c("complement", "identical")
  ids <- if (s$panel == "fixed") rep(seq_along(p), each = protocol$blocks / length(p)) else
    sample.int(length(p), protocol$blocks, replace = TRUE)
  sample_population <- function(n, frequency) {
    shift <- if (s$row_shift > 0) sample(c(-1, 1), n, replace = TRUE) * s$row_shift else rep(0, n)
    probs <- outer(shift, frequency, "+")
    stopifnot(all(probs > 0 & probs < 1))
    g <- matrix(stats::rbinom(length(probs), 2, probs), nrow = n)
    if (s$missing > 0) g[matrix(stats::runif(length(g)) < s$missing, nrow = n)] <- NA_real_
    g
  }
  list(a = sample_population(s$n_a, p[ids]), b = sample_population(s$n_b, q[ids]),
       panel_truth = if (null) 0 else uncertainty_truth(p[ids], q[ids]),
       genomic_truth = if (null) 0 else uncertainty_truth(p, q), null = null,
       latent_p = p[ids], latent_q = q[ids], types = ids)
}

uncertainty_result <- function(fit, sim, s, protocol, replicate, simulation_seed, bootstrap_seed) {
  point <- fit$point
  truth <- if (s$panel == "fixed") sim$panel_truth else sim$genomic_truth
  rows <- lapply(colnames(fit$draws), function(method) {
    v <- fit$draws[, method]
    valid <- v[is.finite(v)]
    ci <- if (length(valid) >= 2) stats::quantile(valid, c(.025, .975), type = 7, names = FALSE) else c(NA_real_, NA_real_)
    eligible <- is.finite(point["delta"]) && point["gst_a"] > 0 && point["gst_b"] > 0 &&
      point["used"] >= 20 && point["used"] * s$copies >= 20 &&
      length(valid) >= ceiling(.95 * protocol$bootstrap_replicates) && length(valid) * .025 >= 5
    direction <- if (eligible && ci[1] > 0 && point["delta"] > 0) 1L else
      if (eligible && ci[2] < 0 && point["delta"] < 0) -1L else 0L
    defined <- all(is.finite(ci))
    data.frame(scenario = s$name, replicate = replicate, method = method,
      simulation_seed = simulation_seed, bootstrap_seed = bootstrap_seed,
      target = if (s$panel == "fixed") "panel" else "genomic",
      truth = truth, panel_truth = sim$panel_truth, estimate = unname(point["delta"]),
      lower = ci[1], upper = ci[2], valid_draws = length(valid), direction = direction,
      covered = defined && ci[1] <= truth && ci[2] >= truth,
      panel_covered = defined && ci[1] <= sim$panel_truth && ci[2] >= sim$panel_truth,
      interval_defined = defined, eligible = eligible,
      used_blocks = unname(point["used"]), no_shared_blocks = unname(point["no_shared"]),
      low_call_blocks = unname(point["low"]), shared_fixed_blocks = unname(point["shared_fixed"]),
      status = "ok", error = "", stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

uncertainty_wilson <- function(k, n) {
  z <- stats::qnorm(.975)
  center <- (k/n + z^2/(2*n)) / (1+z^2/n)
  half <- z * sqrt((k/n)*(1-k/n)/n + z^2/(4*n^2)) / (1+z^2/n)
  c(lower = center-half, upper = center+half)
}
