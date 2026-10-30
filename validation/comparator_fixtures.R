# Deterministic arithmetic fixtures; no demographic interpretation.
# Counts are chromosome copies, not read depth. All called totals are even so
# the same data can be represented as diploid hard calls and GENEPOP input.
comparator_fixtures <- function() {
  make <- function(alt, called = 40) {
    alt <- as.matrix(alt)
    if (length(called) == 1L) called <- matrix(called, nrow(alt), ncol(alt))
    dimnames(alt) <- dimnames(called) <- list(LETTERS[seq_len(nrow(alt))],
                                            paste0("snp", seq_len(ncol(alt))))
    snpflowR::snp_counts(alt, called)
  }
  base <- rbind(c(2, 6, 12, 22, 30, 36), c(12, 18, 20, 28, 32, 34),
                c(8, 26, 32, 14, 4, 24))
  missing_n <- matrix(40, 3, 6)
  missing_n[cbind(c(1, 2, 3), c(2, 3, 4))] <- 0
  missing_a <- base
  missing_a[missing_n == 0] <- 0
  out <- list(
    interior_three = make(base),
    unequal_samples = make(base * c(1, 4, 2), matrix(rep(c(40, 160, 80), 6), 3)),
    missing_pop_loci = make(missing_a, missing_n),
    shared_fixed_added = make(cbind(base, c(0, 0, 0), c(40, 40, 40))),
    opposite_fixed_added = make(cbind(base, c(0, 40, 0), c(40, 0, 40))),
    all_opposite_fixed = make(rbind(c(0, 40), c(40, 0))),
    all_shared_fixed = make(rbind(c(0, 40), c(0, 40), c(0, 40))),
    identical_pair = make(rbind(c(10, 20, 30, 10), c(10, 20, 30, 10), c(20, 30, 10, 30))),
    fixed_source_zero_diversity = make(rbind(c(0, 0, 0, 0), c(10, 20, 30, 10))),
    identical_all = make(rbind(c(10, 20, 30, 10), c(10, 20, 30, 10)))
  )
  # Independent, frozen arithmetic cases: 100 seeds, 3 populations, 50 loci.
  # These are randomized count fixtures, NOT 100 demographic simulations.
  for (seed in 810001:810100) {
    set.seed(seed)
    n <- matrix(2L * sample(10:40, 150, replace = TRUE), 3)
    a <- matrix(stats::rbinom(150, n, stats::runif(150, .05, .95)), 3)
    a[] <- pmax(1, pmin(n - 1, a)) # interior domain intentionally prespecified
    out[[paste0("random_", seed)]] <- make(a, n)
  }
  out
}

# Deterministic chromosome pairing for a parser round trip. Genotype structure
# is artificial and irrelevant to the point estimates under comparison.
fixture_genotypes <- function(x) {
  sizes <- apply(x$called, 1, max) / 2L
  pop <- rep(rownames(x$alt), sizes)
  g <- matrix(NA_integer_, sum(sizes), ncol(x$alt),
              dimnames = list(paste0("sample", seq_len(sum(sizes))), colnames(x$alt)))
  for (i in seq_along(sizes)) {
    rows <- which(pop == rownames(x$alt)[i])
    for (l in seq_len(ncol(g))) {
      n <- x$called[i, l] / 2L
      if (!n) next
      a <- x$alt[i, l]
      values <- pmin(2L, pmax(0L, a - 2L * (seq_len(n) - 1L)))
      g[rows[seq_len(n)], l] <- values
    }
  }
  list(genotypes = g, populations = pop)
}

write_fixture_genepop <- function(g, path) {
  lines <- c("snpflowR arithmetic comparator fixture", colnames(g$genotypes))
  for (pop in unique(g$populations)) {
    lines <- c(lines, "Pop")
    for (i in which(g$populations == pop)) {
      calls <- c("0101", "0102", "0202")[g$genotypes[i, ] + 1L]
      calls[is.na(calls)] <- "0000"
      lines <- c(lines, paste0(rownames(g$genotypes)[i], ", ", paste(calls, collapse = " ")))
    }
  }
  writeLines(lines, path)
}
