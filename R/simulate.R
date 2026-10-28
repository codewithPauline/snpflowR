#' Simulate neutral biallelic loci under a forward Wright-Fisher model
#'
#' A transparent, unlinked two-population simulator for development experiments.
#' It is not a fitted demographic model or a substitute for independent simulators.
#'
#' @param n_loci Number of independent loci.
#' @param n_ind Individuals sampled from populations A and B, length two.
#' @param effective_size Diploid effective population sizes, length two.
#' @param generations Number of drift-migration generations.
#' @param migration Two by two matrix. Entry `[i, j]` is the fraction of the
#'   destination j gene pool replaced by immigrants from source i each generation.
#'   Diagonal entries must be zero and incoming column sums cannot exceed one.
#'   Fractions are not numbers of migrants; unequal sizes change their meaning.
#' @param initial Optional two by `n_loci` matrix of initial frequencies.
#'   By default populations share ancestral frequencies sampled uniformly
#'   between 0.05 and 0.95, followed by independent population sampling.
#' @param mutation Symmetric per-generation biallelic mutation probability.
#' @param missing Independent probability that a sampled diploid genotype is missing.
#' @param seed Optional local seed, preserving the caller's RNG state.
#' @return List with population `counts`, sampled `genotypes`, `populations`,
#'   final `frequencies`, initial `frequencies_initial`, and simulation `truth`.
#'   Loci are retained regardless of whether they are polymorphic at sampling.
#' @details Each generation first mixes pre-generation frequencies according to
#'   incoming migration, then applies symmetric mutation, then binomial genetic
#'   drift in each population. All populations update synchronously. Genotypes
#'   are independent binomial samples with size two from final frequencies
#'   (Hardy-Weinberg sampling); missing calls are then masked. No selection,
#'   recombination, linkage, ascertainment, or population-size changes are modeled.
#' @export
#' @examples
#' m <- matrix(c(0, 0, 0.02, 0), 2, 2)
#' sim <- simulate_snp_flow(n_loci = 100, migration = m, seed = 10)
#' sim$truth$migration
#' sim$counts
simulate_snp_flow <- function(n_loci = 200L, n_ind = c(40L, 40L),
                              effective_size = c(500L, 500L), generations = 100L,
                              migration = matrix(c(0, 0, 0.02, 0), 2, 2),
                              initial = NULL, mutation = 0, missing = 0, seed = NULL) {
  .whole_number(n_loci, "n_loci", 1)
  .whole_number(generations, "generations")
  if (length(n_ind) != 2L || length(effective_size) != 2L) {
    stop("n_ind and effective_size must each have length two.", call. = FALSE)
  }
  for (i in 1:2) {
    .whole_number(n_ind[i], "n_ind", 1)
    .whole_number(effective_size[i], "effective_size", 1, floor(.Machine$integer.max / 2))
  }
  .probability(mutation, "mutation", 0.5)
  .probability(missing, "missing")
  if (!is.null(seed)) .whole_number(seed, "seed")
  if (!is.matrix(migration) || !is.numeric(migration) ||
      !identical(dim(migration), c(2L, 2L)) || anyNA(migration) ||
      any(!is.finite(migration)) || any(migration < 0) ||
      any(diag(migration) != 0) || any(colSums(migration) > 1)) {
    stop("migration must be a nonnegative finite 2 x 2 matrix, with zero diagonal and incoming column sums <= 1.", call. = FALSE)
  }
  if (!is.null(initial) && (!is.matrix(initial) || !is.numeric(initial) ||
      !identical(dim(initial), c(2L, as.integer(n_loci))) || anyNA(initial) ||
      any(!is.finite(initial)) || any(initial < 0 | initial > 1))) {
    stop("initial must be a 2 x n_loci frequency matrix in [0, 1].", call. = FALSE)
  }
  .with_seed(seed, {
    if (is.null(initial)) {
      ancestral <- stats::runif(n_loci, 0.05, 0.95)
      p <- rbind(stats::rbinom(n_loci, 2 * effective_size[1], ancestral) / (2 * effective_size[1]),
                 stats::rbinom(n_loci, 2 * effective_size[2], ancestral) / (2 * effective_size[2]))
    } else {
      p <- initial
    }
    p_start <- p
    mixing <- migration
    diag(mixing) <- 1 - colSums(migration)
    for (g in seq_len(generations)) {
      p <- t(mixing) %*% p
      p <- mutation + (1 - 2 * mutation) * p
      for (i in 1:2) p[i, ] <- stats::rbinom(n_loci, 2 * effective_size[i], p[i, ]) / (2 * effective_size[i])
    }
    genotypes <- do.call(rbind, lapply(1:2, function(i) {
      matrix(stats::rbinom(n_ind[i] * n_loci, 2, rep(p[i, ], each = n_ind[i])), n_ind[i], n_loci)
    }))
    if (missing > 0) genotypes[matrix(stats::runif(length(genotypes)) < missing, nrow(genotypes))] <- NA_real_
    rownames(genotypes) <- paste0("sample", seq_len(sum(n_ind)))
    colnames(genotypes) <- paste0("snp", seq_len(n_loci))
    populations <- rep(c("A", "B"), n_ind)
    dimnames(p) <- dimnames(p_start) <- list(c("A", "B"), colnames(genotypes))
    dimnames(migration) <- list(source = c("A", "B"), destination = c("A", "B"))
    list(counts = genotype_counts(genotypes, populations), genotypes = genotypes,
      populations = populations, frequencies = p, frequencies_initial = p_start,
      truth = list(model = "neutral unlinked diploid Wright-Fisher", migration = migration,
        n_loci = n_loci, n_ind = n_ind, effective_size = effective_size,
        generations = generations, mutation = mutation, missing = missing, seed = seed))
  })
}
