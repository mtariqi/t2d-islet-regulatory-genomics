x <- read.delim(
  "results/influence/GSE221156/core285_lodo_assessment.tsv"
)

stopifnot(nrow(x) == 285)

metrics <- c(
  "Complete refits",
  "Direction stable",
  "High FDR retention"
)

counts <- c(
  sum(x$complete_refits),
  sum(x$stable_direction),
  sum(x$high_FDR_retention)
)

plot_figure <- function() {
  op <- par(mar = c(5, 10, 3, 3))
  on.exit(par(op))

  bp <- barplot(
    rev(counts),
    names.arg = rev(metrics),
    horiz = TRUE,
    las = 1,
    xlim = c(0, 310),
    col = "steelblue4",
    border = NA,
    xlab = "Number of genes",
    main = "Donor robustness of 285 T2D beta-cell candidates"
  )

  text(
    rev(counts) + 5,
    bp,
    labels = paste0(rev(counts), "/285"),
    pos = 4,
    cex = 1.1
  )

  abline(v = 285, lty = 3, col = "gray50")
}

pdf(
  "figures/phase7/core285_lodo_robustness.pdf",
  width = 10,
  height = 5
)
plot_figure()
dev.off()

png(
  "figures/phase7/core285_lodo_robustness.png",
  width = 2400,
  height = 1200,
  res = 240
)
plot_figure()
dev.off()
