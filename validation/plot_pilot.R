# Run from repository root after run_pilot.R; needs only base R graphics.
d <- read.csv("validation/results/pilot-summary.csv", stringsAsFactors = FALSE)
draw_pilot <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old))
  par(mfrow = c(1, 2), mar = c(4.7, 12.5, 4.5, 1.2), oma = c(2.4, 0, 2.3, 0),
      family = "sans", las = 1, xaxs = "i", mgp = c(2.5, .6, 0))
  panels <- list(
    list(keys = c("no_migration", "symmetric", "null_unequal_sample", "null_unequal_Ne", "founder_no_migration"),
         labels = c("No migration", "Symmetric migration", "No flow; samples 10:80", "No flow; Ne 100:1000", "Founder history; no flow"),
         title = "Signals without asymmetric migration", col = "#B1493D"),
    list(keys = c("A_to_B", "B_to_A", "asymmetric_4to1", "A_to_B_missing20"),
         labels = c("A to B", "B to A", "Asymmetric 4:1", "A to B; 20% missing"),
         title = "Recovery when migration is asymmetric", col = "#176C78"))
  for (panel in panels) {
    z <- d[match(panel$keys, d$scenario), ]
    y <- rev(seq_len(nrow(z)))
    plot(z$supported_rate * 100, y, xlim = c(0, 105), ylim = c(.5, nrow(z)+.8),
         yaxt = "n", xlab = "Supported genetic contrasts (%)", ylab = "", type = "n",
         main = panel$title, cex.main = .93, axes = FALSE)
    axis(1, at = seq(0, 100, 20), cex.axis = .8)
    axis(2, at = y, labels = panel$labels, tick = FALSE, cex.axis = .85)
    abline(v = seq(0, 100, 20), col = "#E3E8EC", lwd = .8)
    segments(100*z$supported_lower, y, 100*z$supported_upper, y, col = panel$col, lwd = 2)
    points(100*z$supported_rate, y, pch = 19, col = panel$col, cex = 1.1)
    text(pmin(96, pmax(6, 100*z$supported_rate)), y+.23,
         paste0(z$supported, "/", z$n), cex = .78, col = panel$col)
    box(col = "#CBD3DA")
  }
  mtext("snpflowR | Exploratory forward-simulation pilot", outer = TRUE, side = 3,
        line = .5, adj = .02, font = 2, cex = 1.25)
  mtext("100 datasets per scenario; bars: 95% Wilson Monte Carlo intervals.\nA supported contrast is not proof of migration. Scenario-specific development results.",
        outer = TRUE, side = 1, line = .25, cex = .8, col = "#475569")
}
grDevices::svg("validation/results/pilot-diagnostics.svg", width = 14, height = 5.3)
draw_pilot()
invisible(grDevices::dev.off())
