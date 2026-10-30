#' Estimate directional SNP differentiation and bootstrap signal contrasts
#'
#' Implements the geometric migrant-pool construction of Sundqvist et al. (2016).
#' This is an exploratory genetic summary, not a demographic migration estimator.
#'
#' @param x A `snpflow_counts` object.
#' @param n_boot Number of paired bootstrap replicates. Zero returns point
#'   estimates with all directions unresolved.
#' @param resample `"block"` resamples complete user-defined blocks. `"locus"`
#'   explicitly assumes loci are independent. Both condition on sampled individuals.
#' @param conf_level Nominal interval confidence level, strictly between 0 and 1.
#' @param adjustment `"bonferroni"` allocates the total nominal error across all
#'   unordered population pairs; `"none"` gives pointwise intervals. Neither
#'   choice guarantees empirical coverage for migration parameters.
#' @param min_called Minimum observed allele copies per population-locus.
#' @param min_loci Minimum usable loci for a supported signal label.
#' @param min_units Minimum contributing resampling units for a supported label.
#' @param min_effect Nonnegative minimum absolute GST contrast for a supported
#'   signal. Both the point estimate and interval must exceed this threshold.
#' @param seed Optional nonnegative integer seed. Restores the caller's RNG state.
#' @param keep_bootstrap Retain the replicate contrast matrix.
#' @return A `snpflow_fit` list with `directional_gst` (row is putative source,
#'   column is destination), symmetric `pairwise_gst`, `migration_transform`,
#'   `relative_migration`, pairwise diagnostics in `pairs`, `settings`, and
#'   `provenance`. The paired contrast is reverse GST minus forward GST;
#'   positive values favor the row population as a source under the model.
#' @details
#' Multilocus GST is the ratio of summed between-population diversity to summed
#' total diversity, with equally weighted populations at each locus. It is a
#' plug-in estimator, without finite-sample bias correction, and is not
#' Weir-Cockerham FST. SNPs fixed for opposite alleles have no migrant pool and
#' are excluded from both directional estimates with their count reported.
#' They remain in symmetric pairwise GST. This explicit exclusion differs from
#' the maximal-differentiation convention discussed in the original paper.
#'
#' The descriptive transformation is `(1 / GST - 1) / 4`. Zero GST yields Inf;
#' no arbitrary epsilon is substituted. Relative scaling is unavailable if any
#' off-diagonal transformation is infinite. Undefined comparisons remain NA.
#' Bootstrap intervals concern the GST contrast, avoiding this singularity.
#' A supported label requires at least 199 replicates, at least 95 percent finite
#' draws, at least five expected draws in each nominal tail, sufficient loci and
#' units, a nonsingular point estimate, and an interval beyond `min_effect`.
#' These are reporting safeguards, not established biological validity cutoffs.
#'
#' Resampling units are drawn with equal probability, with replacement. Every
#' locus in a selected block is retained, including repeated selections. The
#' same draws are used for both directions and all population pairs. Unequal
#' block lengths therefore give variable numbers of SNPs per replicate.
#' Intervals use percentile quantiles (R type 7). QC and sampled individuals are
#' fixed; individual sampling uncertainty and demographic confounding are not
#' resolved by this bootstrap. No p-values or posterior probabilities are returned.
#' A prospective sampling study found substantial interval miscalibration under
#' unequal sample sizes and other tested conditions. These are nominal intervals,
#' not calibrated confidence guarantees for the genetic contrast. The returned
#' provenance explicitly records `calibration_status = "not_established"`.
#' @references
#' Sundqvist L, Keenan K, Zackrisson M, Prodohl P, Kleinhans D (2016).
#' Directional genetic differentiation and relative migration.
#' Ecology and Evolution 6:3461-3475. \doi{10.1002/ece3.2096}.
#' @export
#' @examples
#' x <- simulate_snp_flow(n_loci = 100, seed = 12)$counts
#' fit <- estimate_flow(x, n_boot = 199, resample = "locus", seed = 13)
#' fit$pairs
estimate_flow <- function(x, n_boot = 0L, resample = c("block", "locus"),
                          conf_level = 0.95, adjustment = c("bonferroni", "none"),
                          min_called = 2L, min_loci = 20L, min_units = 20L,
                          min_effect = 0, seed = NULL, keep_bootstrap = FALSE) {
  .check_counts(x)
  .whole_number(n_boot, "n_boot")
  .whole_number(min_called, "min_called", 1)
  .whole_number(min_loci, "min_loci", 1)
  .whole_number(min_units, "min_units", 2)
  .probability(conf_level, "conf_level")
  if (conf_level <= 0 || conf_level >= 1) stop("conf_level must lie strictly between 0 and 1.", call. = FALSE)
  .probability(min_effect, "min_effect")
  .flag(keep_bootstrap, "keep_bootstrap")
  if (!is.null(seed)) .whole_number(seed, "seed")
  resample <- match.arg(resample)
  adjustment <- match.arg(adjustment)
  if (n_boot > 0 && resample == "block" && is.null(x$blocks)) {
    stop("Block bootstrap requires x$blocks. For independent SNPs, explicitly choose resample = 'locus'.", call. = FALSE)
  }
  pop <- rownames(x$alt)
  n_pop <- length(pop)
  n_loc <- ncol(x$alt)
  combos <- utils::combn(n_pop, 2)
  n_pair <- ncol(combos)
  unit_labels <- if (resample == "block" && !is.null(x$blocks)) x$blocks else seq_len(n_loc)
  unit_id <- match(unit_labels, unique(unit_labels))
  n_unit <- max(unit_id)
  f <- x$alt / x$called
  f[x$called == 0] <- NA_real_
  dg <- pg <- matrix(NA_real_, n_pop, n_pop, dimnames = list(pop, pop))
  diag(dg) <- diag(pg) <- 0
  components <- vector("list", n_pair)
  pair_rows <- vector("list", n_pair)
  for (k in seq_len(n_pair)) {
    i <- combos[1, k]
    j <- combos[2, k]
    z <- .pair_components(f[i, ], f[j, ], x$called[i, ] >= min_called &
                          x$called[j, ] >= min_called)
    totals <- colSums(z$components)
    g <- .component_gst(totals)
    dg[i, j] <- g[1]
    dg[j, i] <- g[2]
    pg[i, j] <- pg[j, i] <- z$pairwise
    components[[k]] <- rowsum(z$components, unit_id, reorder = FALSE)
    pair_rows[[k]] <- data.frame(population_a = pop[i], population_b = pop[j],
      gst_a_to_b = g[1], gst_b_to_a = g[2], delta_gst = g[2] - g[1],
      pairwise_gst = z$pairwise, n_loci = sum(z$used),
      n_units = length(unique(unit_id[z$used])),
      n_low_calls = sum(!z$observed), n_no_shared = sum(z$no_shared),
      n_shared_fixed = sum(z$shared_fixed), stringsAsFactors = FALSE)
  }
  pairs <- do.call(rbind, pair_rows)
  alpha <- 1 - conf_level
  if (adjustment == "bonferroni") alpha <- alpha / n_pair
  boot <- matrix(NA_real_, n_boot, n_pair)
  if (n_boot > 0) {
    boot <- .with_seed(seed, {
      draws <- matrix(NA_real_, n_boot, n_pair)
      for (b in seq_len(n_boot)) {
        weights <- tabulate(sample.int(n_unit, n_unit, replace = TRUE), nbins = n_unit)
        for (k in seq_len(n_pair)) {
          g <- .component_gst(as.numeric(crossprod(weights, components[[k]])))
          draws[b, k] <- g[2] - g[1]
        }
      }
      draws
    })
  }
  pairs$lower <- pairs$upper <- NA_real_
  pairs$n_boot_valid <- integer(n_pair)
  pairs$direction <- rep("unresolved", n_pair)
  pairs$reason <- rep("no_bootstrap", n_pair)
  for (k in seq_len(n_pair)) {
    draws <- boot[, k]
    valid <- draws[is.finite(draws)]
    nv <- length(valid)
    pairs$n_boot_valid[k] <- nv
    if (nv >= 2L) {
      ci <- stats::quantile(valid, c(alpha / 2, 1 - alpha / 2), names = FALSE, type = 7)
      pairs$lower[k] <- ci[1]
      pairs$upper[k] <- ci[2]
    }
    reason <- if (!is.finite(pairs$delta_gst[k])) {
      "undefined_differentiation"
    } else if (pairs$n_loci[k] < min_loci) {
      "insufficient_loci"
    } else if (min(pairs$gst_a_to_b[k], pairs$gst_b_to_a[k]) == 0) {
      "singular_migration_transform"
    } else if (n_boot == 0) {
      "no_bootstrap"
    } else if (pairs$n_units[k] < min_units) {
      "insufficient_resampling_units"
    } else if (n_boot < 199L) {
      "insufficient_bootstrap_replicates"
    } else if (nv < ceiling(0.95 * n_boot)) {
      "too_many_undefined_replicates"
    } else if (nv * alpha / 2 < 5 - 1e-12) {
      "insufficient_tail_resolution"
    } else if (pairs$lower[k] > min_effect && pairs$delta_gst[k] > min_effect) {
      pairs$direction[k] <- paste(pairs$population_a[k], "->", pairs$population_b[k])
      "bootstrap_supported_signal"
    } else if (pairs$upper[k] < -min_effect && pairs$delta_gst[k] < -min_effect) {
      pairs$direction[k] <- paste(pairs$population_b[k], "->", pairs$population_a[k])
      "bootstrap_supported_signal"
    } else {
      "interval_overlaps_threshold"
    }
    pairs$reason[k] <- reason
  }
  transform <- (1 / dg - 1) / 4
  diag(transform) <- NA_real_
  relative <- matrix(NA_real_, n_pop, n_pop, dimnames = dimnames(dg))
  finite <- is.finite(transform)
  scale_status <- "no_finite_positive_transform"
  if (any(is.infinite(transform))) {
    scale_status <- "unavailable_infinite_transform"
  } else if (any(finite) && max(transform[finite]) > 0) {
    relative[finite] <- transform[finite] / max(transform[finite])
    scale_status <- "scaled_by_largest_finite_transform"
  }
  structure(list(directional_gst = dg, pairwise_gst = pg,
    migration_transform = transform, relative_migration = relative,
    scale_status = scale_status, pairs = pairs,
    bootstrap = if (keep_bootstrap) boot else NULL,
    settings = list(n_boot = n_boot, resample = resample, conf_level = conf_level,
      adjustment = adjustment, interval_level = 1 - alpha, min_called = min_called,
      min_loci = min_loci, min_units = min_units, min_effect = min_effect, seed = seed),
    provenance = list(package_version = as.character(utils::packageVersion("snpflowR")),
      r_version = R.version.string, rng_kind = RNGkind(),
      input_loci = n_loc, population_order = pop, qc_history = x$history,
      method = "Sundqvist migrant pool; plug-in multilocus GST; no-shared loci excluded",
      interpretation = "Exploratory genetic signal; not a causal migration estimate",
      calibration_status = "not_established")),
    class = "snpflow_fit")
}

.pair_components <- function(p, q, observed) {
  observed <- observed & is.finite(p) & is.finite(q)
  p[!observed] <- q[!observed] <- 0
  shared_fixed <- observed & p == q & p %in% c(0, 1)
  s_alt <- sqrt(p * q)
  s_ref <- sqrt((1 - p) * (1 - q))
  denom <- s_alt + s_ref
  no_shared <- observed & denom == 0
  used <- observed & !shared_fixed & !no_shared
  r <- numeric(length(p))
  r[used] <- s_alt[used] / denom[used]
  # Ht - Hs = (p - r)^2 / 2 for biallelic frequencies, avoiding cancellation.
  t_a <- (p + r) * (1 - (p + r) / 2)
  t_b <- (q + r) * (1 - (q + r) / 2)
  comp <- cbind(n_a = (p - r)^2 / 2, t_a = t_a,
                 n_b = (q - r)^2 / 2, t_b = t_b)
  comp[!used, ] <- 0
  pair_denom <- sum(((p + q) * (1 - (p + q) / 2))[observed])
  pairwise <- if (pair_denom > 0) sum(((p - q)^2 / 2)[observed]) / pair_denom else NA_real_
  list(components = comp, observed = observed, used = used,
        no_shared = no_shared, shared_fixed = shared_fixed, pairwise = pairwise)
}

.component_gst <- function(totals) {
  c(if (totals[2] > 0) totals[1] / totals[2] else NA_real_,
    if (totals[4] > 0) totals[3] / totals[4] else NA_real_)
}

#' Print a directional signal analysis
#' @param x A `snpflow_fit` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.snpflow_fit <- function(x, ...) {
  cat("snpflowR: exploratory directional genetic signals\n")
  cat("Rows are putative sources; columns are destinations.\n")
  print(x$pairs[, c("population_a", "population_b", "delta_gst", "lower", "upper", "direction", "reason")], row.names = FALSE)
  cat("Signals do not establish contemporary migration or absolute rates.\n")
  cat("Intervals are nominal; calibration is not established.\n")
  invisible(x)
}

#' Plot paired directional contrasts
#'
#' @param x A `snpflow_fit` object.
#' @param ... Graphical parameters passed to [graphics::plot.default()].
#' @return `x`, invisibly. Creates a contrast interval plot; gray marks denote
#'   unresolved comparisons. Positive contrasts favor A to B for each label.
#' @export
#' @examples
#' x <- simulate_snp_flow(n_loci = 60, seed = 8)$counts
#' plot(estimate_flow(x))
plot.snpflow_fit <- function(x, ...) {
  p <- x$pairs
  ok <- is.finite(p$delta_gst)
  if (!any(ok)) stop("No finite contrasts to plot.", call. = FALSE)
  p <- p[ok, , drop = FALSE]
  y <- rev(seq_len(nrow(p)))
  span <- range(c(p$delta_gst, p$lower, p$upper, 0), finite = TRUE)
  if (diff(span) == 0) span <- span + c(-0.01, 0.01)
  col <- ifelse(p$direction == "unresolved", "gray50", "#126E82")
  graphics::plot(p$delta_gst, y, xlim = span, ylim = c(0.5, nrow(p) + 0.5),
    yaxt = "n", ylab = "", xlab = "Reverse GST - forward GST (positive favors A -> B)",
    pch = 19, col = col, ...)
  graphics::axis(2, at = y, labels = paste(p$population_a, "/", p$population_b), las = 1)
  graphics::abline(v = 0, lty = 2, col = "gray65")
  ci_ok <- is.finite(p$lower) & is.finite(p$upper)
  graphics::segments(p$lower[ci_ok], y[ci_ok], p$upper[ci_ok], y[ci_ok], col = col[ci_ok], lwd = 2)
  invisible(x)
}
