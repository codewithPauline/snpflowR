.whole_number <- function(x, name, lower = 0, upper = .Machine$integer.max) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < lower || x > upper) {
    stop(name, " must be one integer between ", lower, " and ", upper, ".",
         call. = FALSE)
  }
  invisible(x)
}

.probability <- function(x, name, upper = 1) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 0 || x > upper) {
    stop(name, " must be one number between 0 and ", upper, ".", call. = FALSE)
  }
  invisible(x)
}

.flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
}

.labels <- function(x, n, name) {
  if (length(x) != n || anyNA(x) || any(!nzchar(trimws(x))) || anyDuplicated(x)) {
    stop(name, " must contain ", n, " unique, nonempty labels.", call. = FALSE)
  }
  as.character(x)
}

# A supplied seed is local to this operation. NULL consumes the caller's RNG.
.with_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  .whole_number(seed, "seed")
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

.check_counts <- function(x) {
  if (!inherits(x, "snpflow_counts")) {
    stop("x must be created with snp_counts(), genotype_counts(), or read_snp_vcf().",
         call. = FALSE)
  }
  # Revalidate to catch accidental edits of the public list fields.
  snp_counts(x$alt, x$called, x$loci, x$blocks)
  invisible(x)
}
