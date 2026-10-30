# Deterministic verification of the research engine; no evaluation seeds used.
library(snpflowR)
source("validation/uncertainty/helpers.R")
Rcpp::sourceCpp("validation/uncertainty/bootstrap.cpp")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
protocol <- read_uncertainty_protocol()
testthat::test_that("weighted genotype resampling equals literal row and column duplication", {
  for (seed in 50000001:50000012) {
    set.seed(seed)
    a <- matrix(sample(c(0:2, NA), 7*9, replace=TRUE), 7)
    b <- matrix(sample(c(0:2, NA), 11*9, replace=TRUE), 11)
    wa <- tabulate(sample.int(7, 7, TRUE), 7)
    wb <- tabulate(sample.int(11, 11, TRUE), 11)
    wk <- tabulate(sample.int(9, 9, TRUE), 9)
    z <- uncertainty_weighted_cpp(a, b, wa, wb, wk)
    ai <- a[rep(seq_len(7), wa), rep(seq_len(9), wk), drop=FALSE]
    bi <- b[rep(seq_len(11), wb), rep(seq_len(9), wk), drop=FALSE]
    x <- genotype_counts(rbind(ai, bi), c(rep("A", 7), rep("B", 11)))
    f <- estimate_flow(x)
    testthat::expect_equal(unname(z[c("gst_a", "gst_b", "delta")]),
      unname(unlist(f$pairs[c("gst_a_to_b", "gst_b_to_a", "delta_gst")])), tolerance=1e-12)
    testthat::expect_equal(unname(z[c("used", "low", "no_shared", "shared_fixed")]),
      unname(unlist(f$pairs[c("n_loci", "n_low_calls", "n_no_shared", "n_shared_fixed")])))
  }
})
testthat::test_that("truth and observed point estimates are independently checked", {
  for (i in seq_len(nrow(protocol$scenarios))) {
    s <- protocol$scenarios[i, ]
    sim <- uncertainty_generate(s, protocol, 50001000+i)
    z <- uncertainty_weighted_cpp(sim$a, sim$b, rep(1L,nrow(sim$a)), rep(1L,nrow(sim$b)), rep(1L,ncol(sim$a)))
    g <- rbind(sim$a, sim$b)[, rep(seq_len(ncol(sim$a)), each=s$copies)]
    x <- genotype_counts(g, c(rep("A", s$n_a),rep("B",s$n_b)), blocks=rep(seq_len(protocol$blocks), each=s$copies))
    f <- estimate_flow(x)
    testthat::expect_equal(unname(z["delta"]), f$pairs$delta_gst, tolerance=1e-12)
    # Known frequencies are multiples of 1/100: exact counts recover the oracle.
    alt <- rbind(A=round(sim$latent_p*100), B=round(sim$latent_q*100))
    latent <- estimate_flow(snp_counts(alt, matrix(100,2,ncol(alt),dimnames=dimnames(alt))))
    testthat::expect_equal(latent$pairs$delta_gst, sim$panel_truth, tolerance=1e-12)
  }
})
testthat::test_that("bootstrap pairing, identity and allele reversal are preserved", {
  sim <- uncertainty_generate(protocol$scenarios[11, ], protocol, 50002001)
  set.seed(50002002); f <- uncertainty_bootstrap_cpp(sim$a,sim$b,1,199)
  testthat::expect_identical(f$draws[,"locus"], f$draws[,"block"])
  set.seed(50002002); flip <- uncertainty_bootstrap_cpp(2-sim$a,2-sim$b,1,199)
  testthat::expect_equal(f$draws, flip$draws, tolerance=1e-12)
  set.seed(50002002); reverse <- uncertainty_bootstrap_cpp(sim$b,sim$a,1,199)
  # Identical row counts imply an identical RNG draw schedule; row weights
  # belong to swapped populations, so only the column-only draws negate exactly.
  testthat::expect_equal(f$draws[,"block"], -reverse$draws[,"block"], tolerance=1e-12)
  set.seed(50002002); copied <- uncertainty_bootstrap_cpp(sim$a,sim$b,4,199)
  testthat::expect_equal(f$point, copied$point)
  testthat::expect_true(any(copied$draws[,"locus"] != copied$draws[,"block"]))
  testthat::expect_error(uncertainty_weighted_cpp(sim$a,sim$b,rep(-1L,40),rep(1L,40),rep(1L,80)), "nonnegative")
})
testthat::test_that("support rule matches the package on identical conditional draws", {
  s <- protocol$scenarios[1, ]
  sim <- uncertainty_generate(s, protocol, 50003001)
  x <- genotype_counts(rbind(sim$a,sim$b), rep(c("A","B"),each=40))
  f <- estimate_flow(x,n_boot=999,resample="locus",seed=50003002,keep_bootstrap=TRUE)
  z <- uncertainty_weighted_cpp(sim$a,sim$b,rep(1L,40),rep(1L,40),rep(1L,80))
  draws <- matrix(rep(f$bootstrap,4), ncol=4, dimnames=list(NULL,protocol$methods))
  r <- uncertainty_result(list(point=z,draws=draws),sim,s,protocol,1,50003001,50003002)
  testthat::expect_equal(r$lower,rep(f$pairs$lower,4))
  testthat::expect_equal(r$upper,rep(f$pairs$upper,4))
  direction <- switch(f$pairs$direction,"A -> B"=1L,"B -> A"=-1L,0L)
  testthat::expect_identical(r$direction,rep(direction,4))
})
cat("Research-engine checks passed; evaluation seeds were not used.\n")
