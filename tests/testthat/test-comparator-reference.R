test_that("point estimates match pinned external reference outputs offline", {
  counts <- utils::read.csv(test_path("fixtures", "comparator-counts.csv"))
  expected <- utils::read.csv(test_path("fixtures", "comparator-expected.csv"))
  for (id in unique(counts$fixture)) {
    rows <- counts[counts$fixture == id, ]
    pops <- unique(rows$population)
    loci <- unique(rows$locus)
    alt <- called <- matrix(0, length(pops), length(loci), dimnames = list(pops, loci))
    idx <- cbind(match(rows$population, pops), match(rows$locus, loci))
    alt[idx] <- rows$alt
    called[idx] <- rows$called
    fit <- estimate_flow(snp_counts(alt, called))
    for (software in c("diveRsity", "corridoR")) {
      e <- expected[expected$fixture == id & expected$software == software, ]
      cells <- cbind(match(e$from, pops), match(e$to, pops))
      expect_length(e$gst, length(pops) * (length(pops) - 1L))
      expect_equal(unname(fit$directional_gst[cells]), e$gst, tolerance = 1e-12)
      expect_equal(unname(fit$relative_migration[cells]), e$relative, tolerance = 1e-12)
    }
  }
})

test_that("a singular pair makes global relative scaling unavailable", {
  a <- rbind(A = c(10, 20, 30, 10), B = c(10, 20, 30, 10), C = c(20, 30, 10, 30))
  n <- matrix(40, 3, 4, dimnames = dimnames(a))
  fit <- estimate_flow(snp_counts(a, n))
  expect_identical(unname(fit$directional_gst["A", "B"]), 0)
  expect_true(is.infinite(fit$migration_transform["A", "B"]))
  expect_true(is.finite(fit$migration_transform["A", "C"]))
  expect_identical(fit$scale_status, "unavailable_infinite_transform")
  expect_true(all(is.na(fit$relative_migration)))
  expect_true(all(fit$directional_gst >= 0))
})

test_that("a fixed source with a fixed pool gives undefined rather than zero GST", {
  a <- rbind(A = c(0, 0, 0, 0), B = c(10, 20, 30, 10))
  n <- matrix(40, 2, 4, dimnames = dimnames(a))
  fit <- estimate_flow(snp_counts(a, n), min_loci = 1)
  # A and the geometric pool are fixed for the same allele: Ht = Hs = 0.
  expect_true(is.na(fit$directional_gst["A", "B"]))
  expect_true(is.na(fit$migration_transform["A", "B"]))
  expect_true(is.na(fit$relative_migration["A", "B"]))
  expect_equal(unname(fit$directional_gst["B", "A"]), 15 / 41, tolerance = 1e-12)
  expect_identical(fit$pairs$direction, "unresolved")
  expect_identical(fit$pairs$reason, "undefined_differentiation")
})
