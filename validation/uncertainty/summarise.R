# Rscript validation/uncertainty/summarise.R [evaluation|smoke]
source("validation/uncertainty/helpers.R")
args <- commandArgs(trailingOnly=TRUE)
mode <- if(length(args)) args[1] else "evaluation"
stopifnot(mode %in% c("evaluation","smoke"))
protocol <- read_uncertainty_protocol()
out <- if(mode=="evaluation") "validation/uncertainty/results" else "validation/generated/uncertainty-smoke"
failures <- list.files(out,pattern="failure.*json$",full.names=TRUE)
if(length(failures)) stop("Incomplete study: execution failure records exist.")
all_rows <- list(); summary_rows <- list(); comparisons <- list()
for (srow in seq_len(nrow(protocol$scenarios))) {
  s <- protocol$scenarios[srow,]
  n <- if(mode=="smoke") 2L else if(startsWith(s$name,"null_")) protocol$null_replicates else protocol$alternative_replicates
  x <- utils::read.csv(file.path(out,paste0(s$name,".csv")),stringsAsFactors=FALSE)
  stopifnot(nrow(x)==n*length(protocol$methods),all(x$status=="ok"),
            !anyDuplicated(x[c("replicate","method")]),
            setequal(x$method,protocol$methods),all(x$direction %in% -1:1))
  is_null <- startsWith(s$name,"null_")
  for (method in protocol$methods) {
    d <- x[x$method==method,]
    d <- d[order(d$replicate),]
    stopifnot(identical(d$replicate,seq_len(n)))
    cov <- sum(d$covered); support <- sum(d$direction!=0)
    cc <- uncertainty_wilson(cov,n); sc <- uncertainty_wilson(support,n)
    defined <- d$interval_defined
    finite <- is.finite(d$estimate)
    width <- d$upper[defined]-d$lower[defined]
    errors <- d$estimate[finite]-d$truth[finite]
    coverage_screen <- cc[1] >= protocol$screen_coverage_lower && cc[2] <= protocol$screen_coverage_upper
    false_screen <- if(is_null) sc[2] <= protocol$screen_false_support_upper else NA
    summary_rows[[length(summary_rows)+1L]] <- data.frame(
      scenario=s$name,method=method,n=n,null=is_null,panel=s$panel,target=d$target[1],
      covered=cov,coverage=cov/n,coverage_lower=cc[1],coverage_upper=cc[2],
      coverage_mcse=sqrt((cov/n)*(1-cov/n)/n),
      panel_coverage=mean(d$panel_covered),
      coverage_defined=if(any(defined)) mean(d$covered[defined]) else NA_real_,
      undefined_intervals=sum(!defined),undefined_estimates=sum(!finite),
      mean_width=mean(width),median_width=stats::median(width),
      bias=mean(errors),bias_mcse=stats::sd(errors)/sqrt(length(errors)),
      empirical_sd=stats::sd(d$estimate[finite]),
      supported=support,supported_rate=support/n,supported_lower=sc[1],supported_upper=sc[2],
      supported_mcse=sqrt((support/n)*(1-support/n)/n),
      false_direction=if(is_null) support else NA_integer_,
      correct=if(!is_null) sum(d$direction==sign(d$truth)) else NA_integer_,
      wrong=if(!is_null) sum(d$direction!=0 & d$direction!=sign(d$truth)) else NA_integer_,
      unresolved=sum(d$direction==0),
      sign_correct_rate=if(!is_null) mean(sign(d$estimate[finite])==sign(d$truth[finite])) else NA_real_,
      median_no_shared_blocks=stats::median(d$no_shared_blocks),
      coverage_screen=coverage_screen,false_support_screen=false_screen,row.names=NULL)
    if(method!="block") {
      b <- x[x$method=="block",]; b <- b[order(b$replicate),]
      dc <- as.numeric(d$covered)-as.numeric(b$covered)
      ds <- as.numeric(d$direction!=0)-as.numeric(b$direction!=0)
      comparisons[[length(comparisons)+1L]] <- data.frame(scenario=s$name,method=method,
        reference="block",n=n,coverage_difference=mean(dc),coverage_difference_mcse=stats::sd(dc)/sqrt(n),
        support_difference=mean(ds),support_difference_mcse=stats::sd(ds)/sqrt(n))
    }
  }
  all_rows[[s$name]] <- x
}
summary <- do.call(rbind,summary_rows)
utils::write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE,na="NA")
utils::write.csv(do.call(rbind,comparisons),file.path(out,"paired-comparisons.csv"),row.names=FALSE,na="NA")
screens <- do.call(rbind,lapply(protocol$methods,function(m) {
  x <- summary[summary$method==m,]
  data.frame(method=m,coverage_scenarios_passed=sum(x$coverage_screen),coverage_scenarios=nrow(x),
    null_scenarios_passed=sum(x$false_support_screen[x$null]),null_scenarios=sum(x$null),
    meets_all_screens=all(x$coverage_screen) && all(x$false_support_screen[x$null]))
}))
utils::write.csv(screens,file.path(out,"screening.csv"),row.names=FALSE)
datasets <- sum(vapply(all_rows,function(d)length(unique(d$replicate)),integer(1)))
jsonlite::write_json(list(protocol=protocol$protocol,mode=mode,scenarios=length(all_rows),datasets=datasets,
  method_results=sum(vapply(all_rows,nrow,integer(1))),complete=TRUE,execution_failures=0,
  note="Screening is scenario-specific; no public API or default is automatically promoted."),
  file.path(out,"completion.json"),pretty=TRUE,auto_unbox=TRUE)
message("Verified ",length(all_rows)," complete scenarios and ",datasets," datasets (",mode,").")
if(mode=="evaluation") print(screens,row.names=FALSE)
