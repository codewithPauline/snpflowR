test_that("genotype counts retain missingness and names", {
  g <- rbind(c(0, 1, NA), c(1, NA, NA), c(2, 1, 0), c(1, 2, 1))
  colnames(g) <- c("x", "y", "z")
  x <- genotype_counts(g, c("A", "A", "B", "B"))
  expect_equal(unname(x$alt), rbind(c(1, 1, 0), c(3, 3, 1)))
  expect_equal(unname(x$called), rbind(c(4, 2, 0), c(4, 4, 4)))
  expect_identical(x$loci$locus, colnames(g))
  expect_output(print(x), "Entirely missing population-loci: 1")
})

test_that("invalid input fails instead of silently realigning", {
  x <- make_counts()
  expect_error(snp_counts(x$alt + .1, x$called), "integers")
  expect_error(snp_counts(x$alt, x$called - 1000), "alt <= called")
  bad <- x$alt; bad[1, 1] <- NA
  expect_error(snp_counts(bad, x$called), "finite")
  bad <- x$called; rownames(bad) <- c("B", "A")
  expect_error(snp_counts(x$alt, bad), "dimnames")
  expect_error(snp_counts(x$alt, x$called, blocks = c("b", NA)), "blocks")
  expect_error(snp_counts(x$alt, x$called, loci = data.frame(locus = rev(colnames(x$alt)))), "match")
  expect_error(genotype_counts(matrix(.5, 2, 2), c("A", "B")), "hard calls")
  expect_error(genotype_counts(matrix(NaN, 2, 2), c("A", "B")), "NaN")
  expect_error(genotype_counts(matrix(0, 2, 2), c("A", NA)), "populations")
  expect_error(genotype_counts(matrix(0, 2, 2), c(a = "A", b = "B")), "row names")
  expect_error(genotype_counts(matrix(0, 2, 2), c("A", "A")), "two populations")
  bad <- x; bad$alt[1, 1] <- -1
  expect_error(estimate_flow(bad), "nonnegative")
  expect_error(estimate_flow(matrix(0, 2, 2)), "snp_counts")
})

test_that("QC is auditable and keeps rare alleles by default", {
  x <- make_counts(c(0, .001, .2, .3), c(0, 0, .3, .4))
  x$alt[1, 4] <- 0; x$called[1, 4] <- 0
  y <- filter_snps(x)
  expect_identical(y$qc$retained, c(FALSE, TRUE, TRUE, FALSE))
  expect_equal(ncol(y$alt), 2)
  expect_equal(length(y$history), 1)
  z <- filter_snps(y, min_maf = .01)
  expect_equal(ncol(z$alt), 1)
  expect_equal(length(z$history), 2)
  expect_error(filter_snps(x, min_maf = .6), "between")
  expect_error(filter_snps(x, min_called = 2000), "No loci")
  expect_equal(ncol(filter_snps(make_counts(0, 0), drop_monomorphic = FALSE)$alt), 1)
})
