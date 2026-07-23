# density_utils.R — shared helpers for density BigWig generation.
#
# rle_merge_granges(): collapse consecutive genome tiles that carry the
# *exactly equal* score into a single variable-width interval, per
# seqname. BigWig natively stores variable-width intervals, so this is a
# lossless rewrite of the dense per-window track: non-zero values are
# unchanged at every position, and long runs of repeated values — most
# importantly all-zero runs — coalesce into one interval. Sourced by
# calculate_density.R and calculate_density_batch.R.
#
# See docs/archive/carp_feature_requests.md FR-1.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(S4Vectors)
})

rle_merge_granges <- function(d) {
  # d: GRanges of contiguous, non-overlapping tiles with a numeric
  # `score` metadata column and seqlengths set (as produced by
  # binnedAverage over tileGenome). Returns a GRanges in which adjacent
  # tiles with an exactly-equal score, within the same seqname, are
  # merged into one interval [start_of_first, end_of_last]. Comparison
  # is exact (no rounding) so the track is lossless.
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
  pieces <- vector("list", length(run_name))
  for (k in seq_along(run_name)) {
    dd <- d[run_start[k]:run_end[k]]
    r <- S4Vectors::Rle(score(dd))
    ends_idx   <- cumsum(runLength(r))
    starts_idx <- c(1L, head(ends_idx, -1L) + 1L)
    pieces[[k]] <- GRanges(
      seqnames = run_name[k],
      ranges   = IRanges(start = start(dd)[starts_idx],
                         end   = end(dd)[ends_idx]),
      score    = runValue(r)
    )
  }
  if (length(pieces) == 0) return(d[0])
  merged <- do.call(c, unname(pieces))
  seqlevels(merged) <- seqlevels(d)
  seqinfo(merged)   <- seqinfo(d)
  sort(merged)
}
