make_counts <- function(p = c(.1, .3, .6, .9), q = c(.4, .5, .7, .8), copies = 1000, blocks = NULL) {
  a <- rbind(A = round(p * copies), B = round(q * copies))
  colnames(a) <- paste0("s", seq_along(p))
  snp_counts(a, matrix(copies, 2, length(p), dimnames = dimnames(a)), blocks = blocks)
}

# Independent, direct two-allele vector arithmetic (does not use package internals).
reference_gst <- function(p, q, directional = TRUE) {
  hs <- ht <- numeric(length(p))
  for (l in seq_along(p)) {
    a <- c(p[l], 1 - p[l])
    b <- c(q[l], 1 - q[l])
    if (directional) {
      b <- sqrt(a * b)
      b <- b / sum(b)
    }
    hs[l] <- ((1 - sum(a * a)) + (1 - sum(b * b))) / 2
    ht[l] <- 1 - sum(((a + b) / 2)^2)
  }
  sum(ht - hs) / sum(ht)
}
