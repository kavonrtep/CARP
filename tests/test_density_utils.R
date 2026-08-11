#!/usr/bin/env Rscript
# Equivalence tests for the density layer (scripts/density_utils.R). Both compare
# the REAL shipped function against a verbatim copy of the implementation it
# replaced.
#
# 1. rle_merge_granges — the density-layer hang fix. The per-seqlevel
#    `d[as.logical(seqnames(d) == sn)]` scan was O(n_tiles x n_seqnames) and hung
#    at ~90M tiles x >10k contigs; rewritten to a single pass over the sorted
#    seqnames Rle runs.
#
# 2. density_track — the memory fix. The old code built the WHOLE-genome
#    tileGenome grid for every track and then discarded the bins the track did
#    not occupy (568.5 GB peak on the 94 Gbp run-000156). density_track()
#    reproduces tileGenome's phase in closed form and streams ONE SEQUENCE at a
#    time. An earlier "just tile the occupied sequences" attempt was correctly
#    rejected — tileGenome sets its tile width from the whole-genome total, so
#    naive subset tiling changes the bins — which is exactly why this equivalence
#    test exists. The grid arithmetic itself is pinned against the real
#    tileGenome by tests/test_density_tiling.R.
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

# ── 2. density_track vs the whole-genome-grid implementations it replaced ─────

# verbatim: calculate_density_batch.R::smooth_score2 (pre-change)
smooth_score2_old <- function(x, N_for_mean = 10){
  sc <- c(rep(0, N_for_mean-1), x, rep(0, N_for_mean-1))
  sc_smooth <- filter(sc, rep(1/N_for_mean, N_for_mean), sides=2)
  sc_smooth[(N_for_mean):(length(sc_smooth)-N_for_mean+1)]
}

# verbatim: calculate_density_batch.R::density_per_family (pre-change)
density_per_family_old <- function(g, chr_size_all, step, N_for_mean = 10){
  g <- reduce(g, ignore.strand = TRUE)
  occ <- seqlevels(g)
  not_used <- setdiff(names(chr_size_all), occ)
  chr_in_order <- chr_size_all[c(occ, not_used)]
  bins <- unlist(tileGenome(chr_in_order, tilewidth = step))
  bins <- bins[as.character(seqnames(bins)) %in% occ]
  bins <- keepSeqlevels(bins, occ, pruning.mode = "coarse")
  cvg <- coverage(g)
  d <- binnedAverage(bins, cvg, "score")
  s_part <- split(d$score, seqnames(d))
  d$score <- unlist(lapply(s_part, function(z) smooth_score2_old(z, N_for_mean)))
  d
}

# verbatim: calculate_density.R::get_density2 + its caller's seqlevel extension
# (pre-change). Keeps zero-score bins for scaffolds the input does not touch.
density_total_old <- function(g, chr_size_all, step, N_for_mean = 10){
  g <- reduce(g, ignore.strand = TRUE)
  chr_size <- chr_size_all[seqlevels(g)]
  not_used <- setdiff(names(chr_size_all), names(chr_size))
  chr_size_not_used <- chr_size_all[not_used]
  seqlevels(g) <- c(seqlevels(g), names(chr_size_not_used))
  chr_size_in_order <- chr_size_all[seqlevels(g)]
  cvg <- coverage(g)
  bins <- tileGenome(chr_size_in_order, tilewidth = step)
  d <- binnedAverage(unlist(bins), cvg, "score")
  s_part <- split(d$score, seqnames(d))
  d$score <- unlist(lapply(s_part, function(x) smooth_score2_old(x, N_for_mean)))
  d
}

# Random annotation-like input: a handful of sequences of very different sizes,
# features clustered so the density actually varies, and only SOME sequences
# occupied (the case the occupied-only filter used to handle).
make_annotation <- function(seed, sl, n_occ) {
  set.seed(seed)
  usable <- which(sl > 0L)                        # nothing can be annotated on a 0-bp seq
  occ <- names(sl)[sort(sample(usable, min(n_occ, length(usable))))]
  parts <- list()
  for (sn in occ) {
    n <- sample(1:60, 1)
    st <- sample.int(max(1L, sl[[sn]] - 500L), n)
    parts[[sn]] <- GRanges(sn, IRanges(st, width = sample(50:5000, n, replace = TRUE)),
                           strand = sample(c("+", "-", "*"), n, replace = TRUE))
  }
  g <- do.call(c, unname(parts))
  end(g) <- pmin(end(g), sl[as.character(seqnames(g))])
  g
}

test_density_track <- function() {
  slsets <- list(
    setNames(as.integer(c(120000, 43000, 250000, 9000, 71000)),
             c("chr1", "chr2", "chr3", "scaf1", "scaf2")),
    setNames(as.integer(c(1000, 999, 100000, 33333)), c("a", "b", "c", "d")),
    setNames(as.integer(c(500000)), "solo")
  )
  n <- 0L
  for (i in seq_len(24)) {
    sl   <- slsets[[(i %% length(slsets)) + 1]]
    step <- c(1000L, 10000L, 997L)[(i %% 3) + 1]
    n_occ <- max(1L, ((i - 1L) %% length(sl)) + 1L)
    g <- make_annotation(9000 + i, sl, n_occ)

    a <- rle_merge_granges(density_per_family_old(g, sl, step, 10))
    b <- density_track(g, sl, step, 10, keep_unoccupied = FALSE)
    if (!identical(sig(a), sig(b)) || !identical(seqlevels(a), seqlevels(b)) ||
        !identical(seqlengths(a), seqlengths(b)))
      stop(sprintf("density_track(keep_unoccupied=FALSE) != old, trial %d (step=%d, %d/%d occupied)",
                   i, step, n_occ, length(sl)))

    a2 <- rle_merge_granges(density_total_old(g, sl, step, 10))
    b2 <- density_track(g, sl, step, 10, keep_unoccupied = TRUE)
    if (!identical(sig(a2), sig(b2)) || !identical(seqlevels(a2), seqlevels(b2)) ||
        !identical(seqlengths(a2), seqlengths(b2)))
      stop(sprintf("density_track(keep_unoccupied=TRUE) != old, trial %d (step=%d, %d/%d occupied)",
                   i, step, n_occ, length(sl)))
    n <- n + 2L
  }
  cat(sprintf("  density_track: new == old (%d comparisons, partial/full occupancy, both modes)\n", n))
}

test_density_track()

# A zero-length sequence gets no bins from tileGenome, so it must contribute no
# interval and — the failure mode this guards — must not shift the sequences that
# follow it out of alignment (assigning NULL to a list element DELETES it in R).
# Not an equivalence check: the old whole-genome path CRASHED on this input.
# `smooth_score2(numeric(0))` indexes (N):(len-N+1) = 10:9, which counts DOWN and
# returns two spurious values, so `d$score <- unlist(...)` hit a length mismatch.
test_zero_length_seq <- function() {
  sl <- setNames(as.integer(c(80000, 0, 45000, 12000)), c("s1", "empty", "s2", "s3"))
  g  <- c(GRanges("s1", IRanges(c(1000, 40000), width = 3000)),
          GRanges("s3", IRanges(500, width = 2000)))
  for (ku in c(FALSE, TRUE)) {
    d <- density_track(g, sl, 1000L, 10, keep_unoccupied = ku)
    if ("empty" %in% as.character(seqnames(d)))
      stop("density_track emitted an interval on a zero-length sequence")
    for (sn in c("s1", "s3")) {
      dd <- d[as.character(seqnames(d)) == sn]
      if (length(dd) == 0 || min(start(dd)) != 1L || max(end(dd)) != sl[[sn]])
        stop(sprintf("density_track: %s not fully tiled (keep_unoccupied=%s)", sn, ku))
    }
  }
  cat("  density_track: zero-length sequence skipped without shifting the rest\n")
}
test_zero_length_seq()

cat("test_density_utils: PASSED\n")
