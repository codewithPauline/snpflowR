# Rscript validation/uncertainty/run_study.R [evaluation|smoke] [scenario ids]
# Complete scenarios are atomic checkpoints. Evaluation seeds follow protocol v1.
args <- commandArgs(trailingOnly=TRUE)
mode <- if (length(args)) args[1] else "evaluation"
stopifnot(mode %in% c("evaluation","smoke"))
source("validation/uncertainty/helpers.R")
Rcpp::sourceCpp("validation/uncertainty/bootstrap.cpp")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
protocol <- read_uncertainty_protocol()
out <- if (mode == "evaluation") "validation/uncertainty/results" else "validation/generated/uncertainty-smoke"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
files <- file.path("validation/uncertainty",c("protocol-v1.json","PROTOCOL.md","bootstrap.cpp","helpers.R","run_study.R"))
fingerprint <- as.list(tools::md5sum(files))
stamp <- file.path(out,"source-fingerprint.json")
if (file.exists(stamp)) {
  previous <- jsonlite::read_json(stamp,simplifyVector=FALSE)
  stopifnot(identical(previous,fingerprint))
} else {
  temporary <- tempfile("fingerprint-",tmpdir=out)
  jsonlite::write_json(fingerprint,temporary,auto_unbox=TRUE,pretty=TRUE)
  stopifnot(file.rename(temporary,stamp))
}
scenarios <- protocol$scenarios
if (length(args)>1) scenarios <- scenarios[scenarios$id %in% as.integer(args[-1]),,drop=FALSE]
stopifnot(nrow(scenarios)>0)
for (srow in seq_len(nrow(scenarios))) {
  s <- scenarios[srow,]
  n <- if (mode=="smoke") 2L else if (startsWith(s$name,"null_")) protocol$null_replicates else protocol$alternative_replicates
  path <- file.path(out,paste0(s$name,".csv"))
  if (file.exists(path)) {
    done <- utils::read.csv(path)
    stopifnot(nrow(done)==n*4L,all(sort(unique(done$replicate))==seq_len(n)),all(done$status=="ok"))
    message("Verified checkpoint: ",s$name)
    next
  }
  started <- proc.time()[["elapsed"]]
  rows <- vector("list",n)
  for (r in seq_len(n)) {
    seed <- (if(mode=="smoke") protocol$development_seed_base else protocol$evaluation_seed_base) + protocol$scenario_seed_stride*s$id+r
    boot_seed <- (if(mode=="smoke") protocol$development_seed_base+10000000L else protocol$bootstrap_seed_base) + protocol$scenario_seed_stride*s$id+r
    # Smoke bootstrap range must not coincide with evaluation dataset seeds.
    if(mode=="smoke") boot_seed <- 55000000L+protocol$scenario_seed_stride*s$id+r
    rows[[r]] <- tryCatch({
      sim <- uncertainty_generate(s,protocol,seed)
      set.seed(boot_seed)
      fit <- uncertainty_bootstrap_cpp(sim$a,sim$b,s$copies,protocol$bootstrap_replicates)
      uncertainty_result(fit,sim,s,protocol,r,seed,boot_seed)
    },error=function(e) {
      if (r > 1) utils::write.csv(do.call(rbind,rows[seq_len(r-1)]),
        file.path(out,paste0(s$name,"-partial-before-failure.csv")),row.names=FALSE,na="NA")
      jsonlite::write_json(list(scenario=s$name,replicate=r,simulation_seed=seed,
        bootstrap_seed=boot_seed,error=conditionMessage(e)),
        file.path(out,paste0(s$name,"-failure-",r,".json")),auto_unbox=TRUE,pretty=TRUE)
      stop("Study incomplete: ",s$name," replicate ",r,": ",conditionMessage(e),call.=FALSE)
    })
    if(r %% 100L==0) message(s$name,": ",r,"/",n)
  }
  values <- do.call(rbind,rows)
  utils::write.csv(values,paste0(path,".tmp"),row.names=FALSE,na="NA")
  stopifnot(file.rename(paste0(path,".tmp"),path))
  elapsed <- proc.time()[["elapsed"]]-started
  jsonlite::write_json(list(scenario=s$name,replicates=n,method_rows=nrow(values),elapsed_seconds=elapsed,
    mode=mode,status="complete"),file.path(out,paste0(s$name,"-run.json")),auto_unbox=TRUE,pretty=TRUE)
  message(s$name," completed in ",round(elapsed,1)," seconds")
}
batch <- if(length(args)>1) paste(args[-1],collapse="-") else "all"
writeLines(sub("[[:blank:]]+$","",capture.output(sessionInfo())),file.path(out,paste0("sessionInfo-",batch,".txt")))
message("Requested scenario batch completed. Global completion is checked by summarise.R.")
