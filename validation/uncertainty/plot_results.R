# Static scientific figures; raw proportions and Monte Carlo intervals.
source("validation/uncertainty/helpers.R")
protocol <- read_uncertainty_protocol()
out <- "validation/uncertainty/results"
d <- utils::read.csv(file.path(out,"summary.csv"),stringsAsFactors=FALSE)
labels <- c("Fixed panel, 40:40", "Fixed panel, 10:80", "Random blocks, 40:40",
 "Random blocks, 10:80", "Identical frequencies, 40:40", "Identical frequencies, 10:80",
 "Linked SNP copies, 40:40", "Linked SNP copies, 10:80", "20% missing, 40:40",
 "20% missing, 10:80", "Individual structure, 40:40", "Individual structure, 10:80",
 "Alternative: fixed panel", "Alternative: random blocks", "Alternative: linked copies",
 "Alternative: reversed", "Alternative: individual structure")
names(labels) <- protocol$scenarios$name
colors <- c(locus="#0072B2",block="#009E73",individual="#D55E00",crossed="#CC79A7")
symbols <- c(locus=16,block=17,individual=15,crossed=18)
draw <- function(metric) {
  is_coverage <- metric=="coverage"
  ids <- protocol$scenarios$name[if(is_coverage) seq_len(17) else seq_len(12)]
  n <- length(ids)
  old <- par(mar=c(5,15,6,2),family="sans",las=1,cex=.95)
  on.exit(par(old))
  plot(NA,xlim=c(-.02,1.02),ylim=c(.5,n+.5),xaxs="i",yaxs="i",axes=FALSE,
       xlab=if(is_coverage) "Datasets whose interval contains the true contrast (%)" else
         "Datasets reporting a false genetic direction (%)",ylab="")
  axis(1,at=seq(0,1,.1),labels=seq(0,100,10),col="#B8C0C8",col.axis="#344253")
  axis(2,at=n:1,labels=labels[ids],tick=FALSE,cex.axis=.9)
  for(y in seq_len(n)) abline(h=y,col="#EDF0F3",lwd=1)
  abline(v=seq(0,1,.1),col="#EDF0F3",lwd=.7)
  abline(v=if(is_coverage) .95 else .05,col="#555555",lty=3,lwd=1.4)
  if(is_coverage) abline(h=5.5,col="#B8C0C8",lwd=1)
  shifts <- c(-.24,-.08,.08,.24)
  for(i in seq_along(protocol$methods)) {
    method <- protocol$methods[i]
    z <- d[d$method==method,]; z <- z[match(ids,z$scenario),]
    y <- n:1+shifts[i]
    value <- if(is_coverage) z$coverage else z$supported_rate
    lower <- if(is_coverage) z$coverage_lower else z$supported_lower
    upper <- if(is_coverage) z$coverage_upper else z$supported_upper
    segments(lower,y,upper,y,col=colors[method],lwd=1.6)
    points(value,y,pch=symbols[method],col=colors[method],cex=.9)
  }
  box(col="#B8C0C8")
  title(main=if(is_coverage) "Nominal 95% intervals: measured coverage" else
    "Directional support when the genetic contrast is zero",adj=0,line=4.1,cex.main=1.25)
  mtext(if(is_coverage) "Primary target: fixed genomic panel or genomic distribution, as prespecified" else
    "Twelve null scenarios, 1,000 datasets each; no demographic migration is simulated",
    side=3,line=2.8,adj=0,cex=.8,col="#526170")
  legend("top",inset=c(0,-.075),legend=c("Locus","Block","Individual","Crossed"),
    col=colors,pch=symbols,horiz=TRUE,bty="n",xpd=NA,cex=.9)
  mtext("Bars: Wilson 95% Monte Carlo intervals. Dotted line: nominal reference. No method is assumed calibrated.",
    side=1,line=4,adj=0,cex=.7,col="#526170")
}
for(metric in c("coverage","false-support")) {
  height <- if(metric=="coverage") 10.5 else 8.5
  grDevices::svg(file.path(out,paste0(metric,".svg")),width=12,height=height)
  draw(metric); grDevices::dev.off()
  # Preview for visual verification; exact figure data and SVG remain versioned.
  dir.create("validation/generated/uncertainty-figures",recursive=TRUE,showWarnings=FALSE)
  grDevices::png(file.path("validation/generated/uncertainty-figures",paste0(metric,".png")),
                width=1440,height=height*120,res=120,type="cairo")
  draw(metric); grDevices::dev.off()
}
