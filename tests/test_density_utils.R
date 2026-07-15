#!/usr/bin/env Rscript
# Equivalence test for the density-layer hang fix.
#
# rle_merge_granges (density_utils.R): the per-seqlevel
# `d[as.logical(seqnames(d) == sn)]` scan was O(n_tiles x n_seqnames) and hung at
# ~90M tiles x >10k contigs. Rewritten to a single pass over the sorted seqnames
# Rle runs. This sources the REAL new function and compares it to the old
# per-seqlevel implementation on random tiled inputs.
#
# (The occupied-only tileGenome optimization was investigated and REJECTED:
# tileGenome sets tile width from the whole-genome total, so tiling occupied-only
# changes the bins. calculate_density_batch.R keeps the whole-genome-then-filter
# grid.)
suppressPackageStartupMessages({ library(rtracklayer) })

ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  dirname(dirname(normalizePath(sub("^--file=", "", f[1]))))
})
source(file.path(ROOT, "scripts", "density_utils.R"))   # real new rle_merge_granges

# --- old rle_merge_granges (per-seqlevel scan), verbatim -----------------------
rle_merge_old <- function(d) {
  if (length(d) == 0) return(d)
  d <- sort(d)
  pieces <- list()
  for (sn in seqlevels(d)) {
    dd <- d[as.logical(seqnames(d) == sn)]
    if (length(dd) == 0) next
    r <- S4Vectors::Rle(score(dd))
    ends_idx   <- cumsum(runLength(r))
    starts_idx <- c(1L, head(ends_idx, -1L) + 1L)
    pieces[[sn]] <- GRanges(seqnames = sn,
      ranges = IRanges(start = start(dd)[starts_idx], end = end(dd)[ends_idx]),
      score  = runValue(r))
  }
  if (length(pieces) == 0) return(d[0])
  merged <- do.call(c, unname(pieces))
  seqlevels(merged) <- seqlevels(d)
  seqinfo(merged)   <- seqinfo(d)
  sort(merged)
}

sig <- function(gr) sprintf("%s:%d-%d=%g", as.character(seqnames(gr)),
                            start(gr), end(gr), score(gr))

make_tiles <- function(seed, chrs, step = 100L) {
  set.seed(seed)
  parts <- list()
  for (sn in chrs) {
    ntile <- sample(1:25, 1)
    starts <- seq(1L, by = step, length.out = ntile)
    # scores with repeats (incl. zero runs) so RLE merging has work to do
    sc <- sample(c(0, 0, 0, 0.25, 0.5, 1), ntile, replace = TRUE)
    parts[[sn]] <- GRanges(sn, IRanges(starts, width = step), score = sc)
  }
  g <- do.call(c, unname(parts))
  sl <- setNames(rep(100000L, length(chrs)), chrs)
  seqlengths(g) <- sl[seqlevels(g)]
  g[sample(length(g))]   # shuffle so sort() inside is exercised
}

test_rle_merge <- function() {
  chrsets <- list(c("chr1"), c("chr1","chr2","chr3"),
                  c("c10","c2","c1","chrX","chrM"))
  for (i in 1:40) {
    chrs <- chrsets[[(i %% length(chrsets)) + 1]]
    d <- make_tiles(3000 + i, chrs)
    a <- rle_merge_old(d)
    b <- rle_merge_granges(d)
    if (!identical(sig(a), sig(b)) || !identical(seqlevels(a), seqlevels(b))) {
      cat("MISMATCH trial", i, "\n"); print(sig(a)); print(sig(b))
      stop("rle_merge_granges: new != old")
    }
  }
  cat("  rle_merge_granges: new == old (40 trials, multi-seqname, zero/value runs)\n")
}

test_rle_merge()
cat("test_density_utils: PASSED\n")
