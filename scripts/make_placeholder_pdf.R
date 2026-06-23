#!/usr/bin/env Rscript
# Write a minimal one-page placeholder PDF carrying a short message.
#
# Used as a fallback when a figure-producing step fails (e.g. the data is too
# large to render, or an R error / out-of-memory occurs) so the Snakemake rule
# still yields a valid, non-empty output instead of failing or leaving a
# truncated file.
#
# Usage: make_placeholder_pdf.R <out.pdf> [message]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("usage: make_placeholder_pdf.R <out.pdf> [message]")
}
out <- args[1]
msg <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "Figure not rendered"

pdf(out, width = 8, height = 6)
plot.new()
text(0.5, 0.55, "Figure not rendered", cex = 1.6, col = "grey25")
text(0.5, 0.42, msg, cex = 0.8, col = "grey55")
dev.off()
