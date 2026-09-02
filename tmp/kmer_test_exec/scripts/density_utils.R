# density_utils.R — shared helpers for density BigWig generation.
#
# Contents:
#   genome_tile_params() / seq_tile_ranges()  — reproduce GenomicRanges::tileGenome's
#       bin grid ONE SEQUENCE AT A TIME, so a whole-genome grid is never allocated.
#   density_track()                           — the single density implementation
#       used by both calculate_density.R and calculate_density_batch.R.
#   rle_merge_granges()                       — lossless run-length collapse of the
#       tiled track (FR-1).
#
# Sourced by calculate_density.R and calculate_density_batch.R.
# See docs/archive/carp_feature_requests.md FR-1.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(S4Vectors)
})

# ── tileGenome, one sequence at a time ───────────────────────────────────────
#
# WHY: the density scripts used to build the WHOLE-genome grid
#   bins <- unlist(tileGenome(chr_in_order, tilewidth = step))
# for every track and then keep only the bins on the sequences that track
# actually occupies. On a 94.26 Gbp assembly at step=1000 that is 94.26 M bins
# (plus a 94.26 M-element character vector for the `%in%` filter) per task, and
# the per-family rule runs thousands of such tasks concurrently — measured
# 568.5 GB peak RSS on run-000156, 74% of a 768 GB host, of which 20.8x was
# grid that got discarded.
#
# An earlier attempt to tile "occupied sequences only" was rejected (and the
# rejection recorded in tests/test_density_utils.R) because tileGenome phases its
# bin boundaries by the WHOLE-genome total: tiling a subset changes the bins.
# That is true of naive re-tiling, but the phase is a closed form. From
# GenomicRanges::tileGenome (cut.last.tile.in.chrom = FALSE):
#
#   chrom_breakpoints <- cumsum(as.numeric(seqlengths))      # in the GIVEN order
#   genome_size       <- last(chrom_breakpoints)
#   ntile             <- as.integer(ceiling(genome_size / tilewidth))
#   tw                <- genome_size / ntile
#   tile_breakpoints  <- pmin(ceiling(tw * seq_len(ntile)), genome_size)
#   absolute_ends     <- sorted union of chrom_breakpoints and tile_breakpoints
#   per-sequence ends <- absolute_ends - offset(sequence)
#   per-sequence starts <- previous end + 1, restarting at 1 on each sequence
#
# So a single sequence needs only (genome_size, ntile, tw, its own offset, its
# own length), and because ceiling(tw * k) is strictly increasing (tw >= 1 for
# any step >= 1) the tile indices inside a sequence form a CONTIGUOUS range that
# can be generated directly. Peak grid memory becomes O(longest sequence / step)
# instead of O(genome / step), independent of genome size and of how much of the
# genome a given track occupies.
#
# Absolute positions are kept in double (a 94 Gbp genome overflows int); only
# the relative, per-sequence coordinates become IRanges integers.

genome_tile_params <- function(chr_size_in_order, step) {
  nms <- names(chr_size_in_order)
  if (is.null(nms) || anyDuplicated(nms) || any(is.na(nms) | !nzchar(nms)))
    stop("chr_size_in_order must be a uniquely, fully named vector of lengths")
  # tileGenome coerces to integer and rejects NA/negative; match it.
  sl <- suppressWarnings(setNames(as.integer(chr_size_in_order), nms))
  if (anyNA(sl) || any(sl < 0L))
    stop("sequence lengths must be non-negative and representable as integer")
  cb <- cumsum(as.numeric(sl))
  genome_size <- cb[[length(cb)]]
  if (!is.finite(step) || step < 1 || step > genome_size)
    stop("step must be >= 1 and <= genome size")
  ntile <- ceiling(genome_size / step)
  if (ntile > .Machine$integer.max) stop("this step is too small")
  ntile <- as.integer(ntile)
  list(seqlengths  = sl,
       offsets     = setNames(c(0, cb[-length(cb)]), nms),
       ends        = setNames(cb, nms),
       genome_size = genome_size,
       ntile       = ntile,
       tilewidth   = genome_size / ntile)
}

# Largest k in [0, ntile] with ceiling(tw * k) <= target. Seeded by the exact
# inverse and then nudged, so double rounding in either direction self-corrects
# rather than silently shifting a bin boundary.
.tile_index_at_or_below <- function(tw, target, ntile) {
  k <- floor(target / tw)
  if (!is.finite(k)) k <- 0
  k <- max(0, min(ntile, k))
  while (k > 0     && ceiling(tw * k)       >  target) k <- k - 1
  while (k < ntile && ceiling(tw * (k + 1)) <= target) k <- k + 1
  k
}

# The tile grid restricted to one sequence, as relative (1-based) IRanges —
# identical to subsetting unlist(tileGenome(...)) to that sequence.
seq_tile_ranges <- function(params, sn) {
  len <- params$seqlengths[[sn]]
  if (len <= 0L) return(IRanges())          # tileGenome emits nothing for these
  off     <- params$offsets[[sn]]
  end_abs <- params$ends[[sn]]
  tw      <- params$tilewidth
  k_lo <- .tile_index_at_or_below(tw, off,     params$ntile) + 1
  k_hi <- .tile_index_at_or_below(tw, end_abs, params$ntile)
  ends <- if (k_hi >= k_lo) ceiling(tw * seq.int(k_lo, k_hi)) else numeric(0)
  # Every chromosome boundary is also a bin end (tileGenome unions the chrom
  # breakpoints into the tile breakpoints, deduplicating). This also covers the
  # final tile, whose unclamped end can round past genome_size.
  if (length(ends) == 0L || ends[[length(ends)]] != end_abs) ends <- c(ends, end_abs)
  rel_end <- ends - off
  IRanges(start = c(1, rel_end[-length(rel_end)] + 1), end = rel_end)
}

# ── smoothing ────────────────────────────────────────────────────────────────
# Moving average over N_for_mean bins, zero-padded at both ends so the track
# keeps one value per bin. stats::filter is named explicitly — `filter` is a very
# commonly masked name.
smooth_score2 <- function(x, N_for_mean = 10) {
  sc <- c(rep(0, N_for_mean - 1), x, rep(0, N_for_mean - 1))
  sc_smooth <- stats::filter(sc, rep(1 / N_for_mean, N_for_mean), sides = 2)
  sc_smooth[(N_for_mean):(length(sc_smooth) - N_for_mean + 1)]
}

# ── run-length merge (FR-1) ──────────────────────────────────────────────────
# Collapse consecutive tiles carrying the *exactly equal* score into a single
# variable-width interval. BigWig natively stores variable-width intervals, so
# this is a lossless rewrite of the dense per-window track: non-zero values are
# unchanged at every position, and long runs of repeated values — most
# importantly all-zero runs — coalesce into one interval. Comparison is exact
# (no rounding).

# One sequence's worth: ranges `ir` with scores `sc`, already in position order.
# `si` (optional) is the seqinfo every piece should carry — giving all pieces the
# SAME seqinfo up front makes the final c() a plain concatenation instead of a
# pairwise seqlevel union, and silences the "no sequence levels in common"
# warning that would otherwise fire once per sequence (thousands of lines on a
# scaffold-rich assembly).
.rle_merge_one <- function(sn, ir, sc, si = NULL) {
  if (length(ir) == 0L)
    return(GRanges(seqnames = character(0), ranges = IRanges(), score = numeric(0),
                   seqinfo = if (is.null(si)) Seqinfo() else si))
  r <- S4Vectors::Rle(sc)
  ends_idx   <- cumsum(runLength(r))
  starts_idx <- c(1L, head(ends_idx, -1L) + 1L)
  rng <- IRanges(start = start(ir)[starts_idx], end = end(ir)[ends_idx])
  if (is.null(si)) GRanges(seqnames = sn, ranges = rng, score = runValue(r))
  else             GRanges(seqnames = sn, ranges = rng, score = runValue(r), seqinfo = si)
}

rle_merge_granges <- function(d) {
  # d: GRanges of contiguous, non-overlapping tiles with a numeric `score`
  # metadata column and seqlengths set (as produced by binnedAverage over
  # tileGenome). Returns a GRanges in which adjacent tiles with an exactly-equal
  # score, within the same seqname, are merged into one interval.
  if (length(d) == 0) return(d)
  d <- sort(d)
  # After sort(), each seqname occupies one contiguous run, so walk the runs of
  # the seqnames Rle in a single O(n) pass. The previous per-seqlevel
  # `d[as.logical(seqnames(d) == sn)]` scanned ALL of d once per seqlevel ->
  # O(n_tiles x n_seqnames), which hangs at ~90M tiles x >10k contigs.
  sn_rle    <- seqnames(d)
  run_len   <- runLength(sn_rle)
  run_end   <- cumsum(run_len)
  run_start <- run_end - run_len + 1L
  run_name  <- as.character(runValue(sn_rle))
  si <- seqinfo(d)
  pieces <- vector("list", length(run_name))
  for (k in seq_along(run_name)) {
    idx <- run_start[k]:run_end[k]
    pieces[[k]] <- .rle_merge_one(run_name[k], ranges(d)[idx], score(d)[idx], si)
  }
  if (length(pieces) == 0) return(d[0])
  sort(do.call(c, unname(pieces)))
}

# ── the density track ────────────────────────────────────────────────────────
# Smoothed union-coverage density over the tileGenome grid, run-length merged and
# ready for export(..., format="bigwig").
#
# g              annotation GRanges (any overlaps are collapsed to a union first,
#                so the score is a fraction in [0,1]: the unified annotation
#                deliberately tolerates overlap — L1 Simple_repeat / Low_complexity
#                over a TE, and all L2 nested children — and feeding that straight
#                to coverage() counts stacking depth and pushed the total track to
#                ~3.5x. ignore.strand=TRUE so '*'-strand features merge too.)
# chr_size_all   named vector of ALL sequence lengths (genome_seqlengths.rds)
# step           bin width in bp
# keep_unoccupied  TRUE  -> emit zero-score bins for sequences g does not touch
#                           (the whole-genome total tracks, calculate_density.R)
#                  FALSE -> emit only the occupied sequences; BigWig export drops
#                           dataless seqlevels anyway (per-family/per-class tracks,
#                           calculate_density_batch.R)
#
# Sequences are processed ONE AT A TIME. Every downstream step is already
# per-seqname — binnedAverage is per-bin independent, the smoothing splits by
# seqname, and the run-length merge never spans seqnames — so this is
# byte-identical to computing the whole grid at once, at a fraction of the peak
# memory. coverage() is still taken over the whole of g, so each sequence's Rle
# (and hence binnedAverage's trim()-based partial-bin scaling) is unchanged.
# ── BigWig item-width cap ────────────────────────────────────────────────────
# The UCSC/kent writer behind rtracklayer's BigWig export fails on very wide
# single items, and run-length merging is what creates them: a family with only a
# handful of features on a multi-gigabase chromosome collapses the gap between
# them into ONE zero-value interval spanning almost the whole sequence.
#
# Measured on run-000170 (94.3 Gbp, 2.1 Gb chromosomes): TRC_319 has 47 features
# across 22 sequences, which produced a single 2.07 Gb zero interval, and
# export(format="bigwig") died with "UCSC library operation failed" / "Internal
# error ucsc/bbiWrite.c 414" at step = 1000 while the same track wrote fine at
# step = 10000. Two of 1,599 families failed that way and took a 10-day run down
# with them at the very last rule, after every annotation output was complete.
#
# Splitting such an interval into consecutive pieces carrying the same value is
# representation-only: every base keeps exactly the value it had, so no track
# changes meaning and no consumer sees a difference. It is also a no-op wherever
# no run is that wide — i.e. on every genome small enough not to have the problem,
# whose outputs stay byte-identical.
#
# 100 Mb is far below where the writer misbehaves (2.07 Gb) and far above any
# genuine density run, so it costs a handful of extra intervals: a 2.07 Gb run
# becomes 21.
.BW_MAX_ITEM_WIDTH <- 100e6

split_wide_ranges <- function(gr, max_width = .BW_MAX_ITEM_WIDTH) {
  if (length(gr) == 0L) return(gr)
  w <- as.numeric(width(gr))
  if (all(w <= max_width)) return(gr)              # fast path: nothing to split
  n   <- as.integer(ceiling(w / max_width))        # pieces per input range
  idx <- rep(seq_along(gr), n)                     # parent of each piece
  k   <- sequence(n) - 1L                          # 0-based piece number
  # Arithmetic in double: start + k * max_width can exceed .Machine$integer.max
  # before pmin() brings it back under the sequence length.
  s <- as.numeric(start(gr))[idx] + k * max_width
  e <- pmin(s + max_width - 1, as.numeric(end(gr))[idx])
  out <- GRanges(seqnames = seqnames(gr)[idx],
                 ranges   = IRanges(start = as.integer(s), end = as.integer(e)),
                 seqinfo  = seqinfo(gr))
  mcols(out) <- mcols(gr)[idx, , drop = FALSE]
  out
}

density_track <- function(g, chr_size_all, step, N_for_mean = 10,
                          keep_unoccupied = FALSE) {
  g   <- reduce(g, ignore.strand = TRUE)
  occ <- seqlevels(g)
  not_used     <- setdiff(names(chr_size_all), occ)
  chr_in_order <- chr_size_all[c(occ, not_used)]   # occupied first, as the original
  params <- genome_tile_params(chr_in_order, step)

  cvg  <- coverage(g)
  keep <- if (keep_unoccupied) names(chr_in_order) else occ
  si   <- Seqinfo(seqnames = keep, seqlengths = unname(params$seqlengths[keep]))
  is_occ <- keep %in% occ

  pieces <- vector("list", length(keep))
  for (k in seq_along(keep)) {
    sn <- keep[[k]]
    ir <- seq_tile_ranges(params, sn)
    # A zero-length sequence gets no bins (tileGenome emits none either). Leave
    # the pre-allocated NULL in place and skip — assigning NULL to a list element
    # would DELETE it and shift every later index.
    if (length(ir) == 0L) next
    if (is_occ[[k]]) {
      bins <- GRanges(seqnames = sn, ranges = ir, seqinfo = si[sn])
      sc   <- smooth_score2(score(binnedAverage(bins, cvg[sn], "score")), N_for_mean)
    } else {
      # No coverage anywhere on this sequence: every bin averages to 0 and the
      # moving average of zeros is zero, so the whole sequence collapses to one
      # zero interval — exactly what the whole-genome path produced for it.
      sc <- numeric(length(ir))
    }
    pieces[[k]] <- .rle_merge_one(sn, ir, sc, si)
  }
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (length(pieces) == 0L)
    return(GRanges(seqnames = character(0), ranges = IRanges(),
                   score = numeric(0), seqinfo = si))
  # Cap item width before the caller exports: run-length merging can leave a
  # single interval spanning most of a multi-gigabase chromosome, which the
  # BigWig writer cannot serialise. Values are unchanged — see split_wide_ranges.
  split_wide_ranges(sort(do.call(c, unname(pieces))))
}
