# Post-evaluation computational audit only. It neither changes the protocol nor
# supplies additional evaluation outcomes. Replays every draw using base R.
library(snpflowR)
source("validation/uncertainty/helpers.R")
Rcpp::sourceCpp("validation/uncertainty/bootstrap.cpp")
RNGkind("Mersenne-Twister","Inversion","Rejection")
p <- read_uncertainty_protocol()
rows <- list()
literal <- function(a,b,wa,wb,wk) {
  a <- a[rep(seq_len(nrow(a)),wa),rep(seq_len(ncol(a)),wk),drop=FALSE]
  b <- b[rep(seq_len(nrow(b)),wb),rep(seq_len(ncol(b)),wk),drop=FALSE]
  estimate_flow(genotype_counts(rbind(a,b),c(rep("A",nrow(a)),rep("B",nrow(b)))))$pairs$delta_gst
}
for(id in c(2L,7L,11L,13L)) {
  s <- p$scenarios[p$scenarios$id==id,]
  seed <- p$evaluation_seed_base+p$scenario_seed_stride*id+1L
  bs <- p$bootstrap_seed_base+p$scenario_seed_stride*id+1L
  sim <- uncertainty_generate(s,p,seed)
  set.seed(bs)
  accelerated <- uncertainty_bootstrap_cpp(sim$a,sim$b,s$copies,p$bootstrap_replicates)
  set.seed(bs)
  ref <- matrix(NA_real_,p$bootstrap_replicates,4,dimnames=list(NULL,p$methods))
  weights <- function(n,draws=n) tabulate(1L+floor(stats::runif(draws)*n),nbins=n)
  for(r in seq_len(p$bootstrap_replicates)) {
    wa <- weights(s$n_a); wb <- weights(s$n_b); wk <- weights(p$blocks)
    wl <- if(s$copies==1) wk else weights(p$blocks,p$blocks*s$copies)
    ref[r,"locus"] <- literal(sim$a,sim$b,rep(1L,s$n_a),rep(1L,s$n_b),wl)
    ref[r,"block"] <- literal(sim$a,sim$b,rep(1L,s$n_a),rep(1L,s$n_b),wk)
    ref[r,"individual"] <- literal(sim$a,sim$b,wa,wb,rep(1L,p$blocks))
    ref[r,"crossed"] <- literal(sim$a,sim$b,wa,wb,wk)
  }
  stopifnot(identical(is.na(ref),is.na(accelerated$draws)),
            max(abs(ref-accelerated$draws),na.rm=TRUE)<1e-12)
  saved <- read.csv(file.path("validation/uncertainty/results",paste0(s$name,".csv")))
  saved <- saved[saved$replicate==1,]
  for(method in p$methods) {
    z <- saved[saved$method==method,]
    ci <- stats::quantile(ref[,method],c(.025,.975),type=7,names=FALSE,na.rm=TRUE)
    stopifnot(max(abs(ci-c(z$lower,z$upper)))<1e-12)
    rows[[length(rows)+1L]] <- data.frame(scenario=s$name,replicate=1,method=method,
      simulation_seed=seed,bootstrap_seed=bs,draws=nrow(ref),
      max_abs_difference=max(abs(ref[,method]-accelerated$draws[,method]),na.rm=TRUE),
      saved_ci_agrees=TRUE,point=accelerated$point["delta"],
      bootstrap_mean=mean(ref[,method],na.rm=TRUE),lower=ci[1],upper=ci[2],
      truth=if(s$panel=="fixed")sim$panel_truth else sim$genomic_truth,row.names=NULL)
  }
  message("Literal R replay verified: ",s$name)
}
write.csv(do.call(rbind,rows),"validation/uncertainty/results/literal-replay-check.csv",row.names=FALSE)
