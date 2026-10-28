#' Construct population allele counts
#'
#' @param alt Numeric matrix, populations in rows and loci in columns, containing
#'   counts of the same alternate allele in every population at each locus.
#' @param called Numeric matrix of the same dimensions and dimnames, containing
#'   observed allele-copy counts. Zero denotes an entirely missing population-locus.
#'   These are chromosome counts, not sequencing read depths.
#' @param loci Optional data frame with one row per locus, in matrix column order.
#'   An existing `locus` column must match column names exactly.
#' @param blocks Optional vector of nonmissing block labels, one per locus.
#'   Blocks must be defined using genomic information, not observed direction.
#' @return An object of class `snpflow_counts` with `alt`, `called`, `loci`,
#'   `blocks`, and `history` fields. Counts are retained without imputation.
#' @export
#' @examples
#' a <- rbind(A = c(s1 = 2, s2 = 6), B = c(s1 = 8, s2 = 4))
#' x <- snp_counts(a, matrix(20, 2, 2, dimnames = dimnames(a)))
#' x
snp_counts <- function(alt, called, loci = NULL, blocks = NULL) {
  if (!is.matrix(alt) || !is.numeric(alt) || !is.matrix(called) ||
      !is.numeric(called) || !identical(dim(alt), dim(called)) ||
      nrow(alt) < 2L || ncol(alt) < 1L) {
    stop("alt and called must be numeric matrices with matching dimensions, at least two populations and one locus.",
         call. = FALSE)
  }
  if (anyNA(alt) || anyNA(called) || any(!is.finite(alt)) ||
      any(!is.finite(called)) || any(alt < 0) || any(called < 0) ||
      any(alt != floor(alt)) || any(called != floor(called)) ||
      any(alt > called) || any(called > 2^53)) {
    stop("Counts must be finite, exactly representable nonnegative integers with alt <= called; use 0/0 for missing data.",
         call. = FALSE)
  }
  if (!identical(dimnames(alt), dimnames(called))) {
    stop("alt and called dimnames must match exactly; align alleles and populations first.",
         call. = FALSE)
  }
  pops <- .labels(rownames(alt), nrow(alt), "Population names")
  if (is.null(colnames(alt))) colnames(alt) <- colnames(called) <- paste0("snp", seq_len(ncol(alt)))
  ids <- .labels(colnames(alt), ncol(alt), "Locus names")
  if (is.null(loci)) loci <- data.frame(locus = ids, stringsAsFactors = FALSE)
  if (!is.data.frame(loci) || nrow(loci) != ncol(alt)) {
    stop("loci must have one row per matrix column.", call. = FALSE)
  }
  if ("locus" %in% names(loci) && !identical(as.character(loci$locus), ids)) {
    stop("loci$locus must match matrix columns in order.", call. = FALSE)
  }
  loci$locus <- ids
  rownames(loci) <- NULL
  if (!is.null(blocks)) {
    if (length(blocks) != ncol(alt) || anyNA(blocks) ||
        any(!nzchar(trimws(as.character(blocks))))) {
      stop("blocks must contain one nonmissing, nonempty label per locus.", call. = FALSE)
    }
    blocks <- as.character(blocks)
  }
  structure(list(alt = alt, called = called, loci = loci, blocks = blocks,
                 history = list()), class = "snpflow_counts")
}

#' Convert diploid SNP genotypes to population allele counts
#'
#' @param genotypes Numeric matrix with individuals in rows, loci in columns,
#'   and hard-called alternate-allele dosages 0, 1, 2, or NA. Fractional dosages
#'   and genotype likelihoods are not supported.
#' @param populations Population label per row. If named, names must match the
#'   matrix row names exactly and in the same order.
#' @param loci,blocks See [snp_counts()].
#' @return A `snpflow_counts` object. Missing diploid genotypes contribute zero
#'   to both allele counts. Population order follows first appearance.
#' @export
#' @examples
#' g <- rbind(c(0, 1), c(1, NA), c(2, 1), c(1, 2))
#' genotype_counts(g, c("A", "A", "B", "B"))
genotype_counts <- function(genotypes, populations, loci = NULL, blocks = NULL) {
  if (!is.matrix(genotypes) || !is.numeric(genotypes) ||
      nrow(genotypes) < 2L || ncol(genotypes) < 1L ||
      any(!(genotypes[!is.na(genotypes)] %in% 0:2))) {
    stop("genotypes must be a numeric matrix of diploid hard calls 0, 1, 2, or NA.", call. = FALSE)
  }
  if (any(is.nan(genotypes))) stop("Use NA, not NaN, for missing genotypes.", call. = FALSE)
  if (length(populations) != nrow(genotypes) || anyNA(populations) ||
      any(!nzchar(trimws(as.character(populations))))) {
    stop("populations must contain a nonempty label for each individual.", call. = FALSE)
  }
  if (!is.null(rownames(genotypes))) .labels(rownames(genotypes), nrow(genotypes), "Sample names")
  if (!is.null(names(populations)) &&
      !identical(names(populations), rownames(genotypes))) {
    stop("Named populations must match genotype row names in order.", call. = FALSE)
  }
  populations <- as.character(populations)
  pops <- unique(populations)
  alt <- called <- matrix(0, length(pops), ncol(genotypes),
                          dimnames = list(pops, colnames(genotypes)))
  for (i in seq_along(pops)) {
    g <- genotypes[populations == pops[i], , drop = FALSE]
    alt[i, ] <- colSums(g, na.rm = TRUE)
    called[i, ] <- 2 * colSums(!is.na(g))
  }
  ans <- snp_counts(alt, called, loci, blocks)
  ans$sample_sizes <- table(factor(populations, levels = pops))
  ans
}

#' Inspect allele counts
#' @param x A `snpflow_counts` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.snpflow_counts <- function(x, ...) {
  cat("snpflowR allele counts:", nrow(x$alt), "populations;", ncol(x$alt), "loci\n")
  cat("Entirely missing population-loci:", sum(x$called == 0), "\n")
  if (!is.null(x$blocks)) cat("Resampling blocks:", length(unique(x$blocks)), "\n")
  invisible(x)
}

#' Filter loci with an audit trail
#'
#' @param x A `snpflow_counts` object.
#' @param min_called Minimum observed allele copies in every population.
#' @param min_maf Minimum pooled minor allele frequency, weighting each observed
#'   chromosome equally. Defaults to zero to avoid silently discarding rare alleles.
#' @param drop_monomorphic Remove loci monomorphic across all observed chromosomes.
#' @return A `snpflow_counts` object with a `qc` data frame for all input loci and
#'   the thresholds and audit appended to `history`. Stops if all loci fail.
#' @export
#' @examples
#' x <- simulate_snp_flow(n_loci = 50, seed = 4)$counts
#' filtered <- filter_snps(x, min_called = 10)
#' head(filtered$qc)
filter_snps <- function(x, min_called = 2L, min_maf = 0, drop_monomorphic = TRUE) {
  .check_counts(x)
  .whole_number(min_called, "min_called", lower = 1)
  .probability(min_maf, "min_maf", upper = 0.5)
  .flag(drop_monomorphic, "drop_monomorphic")
  total <- colSums(x$called)
  p <- colSums(x$alt) / total
  maf <- pmin(p, 1 - p)
  low <- apply(x$called, 2, min) < min_called
  mono <- is.finite(maf) & maf == 0
  rare <- is.finite(maf) & maf < min_maf
  keep <- !low & !is.na(maf) & !rare & !(drop_monomorphic & mono)
  audit <- data.frame(locus = colnames(x$alt), min_called = apply(x$called, 2, min),
                      pooled_maf = maf, low_calls = low, monomorphic = mono,
                      below_maf = rare, retained = keep, row.names = NULL)
  if (!any(keep)) stop("No loci pass QC; inspect counts and thresholds.", call. = FALSE)
  ans <- snp_counts(x$alt[, keep, drop = FALSE], x$called[, keep, drop = FALSE],
                    x$loci[keep, , drop = FALSE],
                    if (is.null(x$blocks)) NULL else x$blocks[keep])
  ans$qc <- audit
  ans$history <- c(x$history, list(list(operation = "filter_snps", min_called = min_called,
    min_maf = min_maf, drop_monomorphic = drop_monomorphic, audit = audit)))
  if (!is.null(x$sample_sizes)) ans$sample_sizes <- x$sample_sizes
  ans
}
