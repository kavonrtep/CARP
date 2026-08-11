#!/usr/bin/env Rscript
# Equivalence test for the per-sequence tile grid in scripts/density_utils.R.
#
# density_track() no longer builds the whole-genome tileGenome() grid per track
# (94.26 M bins on a 94 Gbp assembly, x thousands of concurrent tasks -> the
# 568.5 GB peak measured on run-000156). It reproduces tileGenome's phase in
# closed form and materialises ONE SEQUENCE at a time. tileGenome derives the
# effective tile width from the WHOLE-genome total, so this arithmetic is the one
# genuinely new thing in that change: this test pins it against the real
# GenomicRanges::tileGenome over randomised inputs plus the boundary cases.
#
# PASS = seq_tile_ranges(params, sn) is identical to
#        unlist(tileGenome(seqlengths, tilewidth = step)) restricted to `sn`,
#        for every sequence.
suppressPackageStartupMessages({ library(GenomicRanges) })

ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  dirname(dirname(normalizePath(sub("^--file=", "", f[1]))))
})
source(file.path(ROOT, "scripts", "density_utils.R"))

failures <- 0L
checked  <- 0L

check_case <- function(label, sl, step) {
  ref    <- unlist(tileGenome(sl, tilewidth = step))
  params <- genome_tile_params(sl, step)
  ref_sn <- as.character(seqnames(ref))
  for (sn in names(sl)) {
    r <- ranges(ref[ref_sn == sn])
    o <- seq_tile_ranges(params, sn)
    checked <<- checked + 1L
    if (!identical(start(r), start(o)) || !identical(end(r), end(o))) {
      failures <<- failures + 1L
      cat(sprintf("FAIL %s / %s (step=%s): ref %d bins, got %d bins\n",
                  label, sn, format(step), length(r), length(o)))
      n <- min(length(r), length(o), 6L)
      if (n > 0) {
        cat("  ref:", paste(sprintf("%d-%d", start(r)[seq_len(n)], end(r)[seq_len(n)]),
                            collapse = " "), "\n")
        cat("  got:", paste(sprintf("%d-%d", start(o)[seq_len(n)], end(o)[seq_len(n)]),
                            collapse = " "), "\n")
      }
      return(invisible(FALSE))
    }
  }
  invisible(TRUE)
}

nm <- function(x) setNames(as.integer(x), paste0("seq", seq_along(x)))

# ── boundary cases ───────────────────────────────────────────────────────────
# tile breakpoints landing EXACTLY on chromosome breakpoints (tw divides both):
# the union in .superimpose_breakpoints deduplicates them.
check_case("aligned",            nm(c(1000, 1000, 1000)), 100)
# step == genome size (one tile), and step == 1 (one bin per base)
check_case("single-tile",        nm(c(300, 700)),         1000)
check_case("unit-step",          nm(c(37, 11, 52)),       1)
# a single sequence
check_case("one-sequence",       nm(2500),                333)
# steps larger than individual sequences (short scaffolds get one bin each)
check_case("short-scaffolds",    nm(c(10, 5, 100000, 7)), 1000)
# a sequence of length 0 (tileGenome emits nothing for it)
check_case("zero-length",        nm(c(500, 0, 500)),      64)
# genome > .Machine$integer.max: absolute coordinates must stay in double while
# the per-sequence coordinates stay integer. This is the run-000156 shape
# (5,096 sequences / 94.26 Gbp, longest 2,146,571,508 bp).
check_case("over-int-genome",    nm(c(2146571508, 2146417222, 1500000000, 999)), 1000000)
# the very last tile, whose unclamped end can round past genome_size
check_case("last-tile-rounding", nm(c(999983, 1000019)),  7919)

# ── randomised ───────────────────────────────────────────────────────────────
set.seed(20260811)
for (i in seq_len(200)) {
  n  <- sample.int(8, 1)
  sl <- nm(sample.int(5000, n, replace = TRUE))
  gs <- sum(as.numeric(sl))
  if (gs < 1) next
  step <- sample.int(as.integer(min(gs, 3000)), 1)
  check_case(sprintf("random[%d]", i), sl, step)
}
# randomised with many short scaffolds, the assembly shape that motivated this
for (i in seq_len(50)) {
  sl <- nm(c(sample.int(200000, 3, replace = TRUE),
             sample.int(2000, 40, replace = TRUE)))
  step <- sample(c(1000L, 10000L, 997L, 65536L), 1)
  check_case(sprintf("scaffolds[%d]", i), sl, step)
}

if (failures > 0L) {
  cat(sprintf("test_density_tiling: FAILED (%d of %d sequence grids differ)\n",
              failures, checked))
  quit(save = "no", status = 1)
}
cat(sprintf("test_density_tiling: PASSED (%d sequence grids identical to tileGenome)\n",
            checked))
