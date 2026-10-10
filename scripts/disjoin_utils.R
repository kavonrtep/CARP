# disjoin_by_seq_groups() — disjoin() that stays fast on fragmented large genomes.
#
# Sourced by clean_rm_output.R and merge_repeat_annotations.R (both run
# gff_cleanup() over a whole-genome GRanges).
#
# Why: GenomicRanges' disjoin() splits by seqname and calls IRanges' disjoin on
# the CompressedIRangesList. That method has a fast path that packs every
# sequence onto ONE integer axis — but when the summed span of all sequences
# exceeds .Machine$integer.max (any genome > ~2.1 Gbp) it falls back to
# endoapply(), a per-sequence R loop costing ~5 ms per sequence. Harmless for a
# 4 Gbp genome in 10 chromosomes; ~25-30 min on a 22.4 Gbp assembly in 305,167
# sequences (merge_rm_and_dante took 4,804 s there). Measured with 400k ranges:
# 1,000 seqs 0.5 s, 10k 48 s, 40k 203 s, 160k 812 s (with.revmap = TRUE).
#
# Fix: disjoin consecutive groups of sequences whose summed span fits the fast
# path, then concatenate. The result is identical to disjoin(gr, ...):
#   * disjoin never combines ranges from different sequences;
#   * its output is seqname-major in seqlevel order (then strand, then
#     position), and the groups are consecutive in seqlevel order;
#   * revmap indices are mapped back from group-local to gr's indices;
#   * every group keeps gr's full Seqinfo, so c() merges identical Seqinfo.
disjoin_by_seq_groups <- function(gr, with.revmap = FALSE, ignore.strand = FALSE,
                                  max_span = 2e9) {
  if (length(gr) == 0)
    return(disjoin(gr, with.revmap = with.revmap, ignore.strand = ignore.strand))
  sn   <- as.integer(seqnames(gr))                       # seqlevel index per range
  nlev <- length(seqlevels(gr))
  # Per-seqlevel span as IRanges' fast path measures it: width(range) + 1.
  lo <- hi <- numeric(nlev)
  lo[] <- Inf; hi[] <- -Inf
  mn <- tapply(as.numeric(start(gr)), sn, min)
  mx <- tapply(as.numeric(end(gr)),   sn, max)
  used <- as.integer(names(mn))
  lo[used] <- mn; hi[used] <- mx
  span <- ifelse(is.finite(lo), hi - lo + 2, 0)
  if (sum(span) <= max_span)                              # already on the fast path
    return(disjoin(gr, with.revmap = with.revmap, ignore.strand = ignore.strand))
  # Greedy grouping of consecutive seqlevels (one sequence always fits: IRanges
  # coordinates are < 2^31).
  grp <- integer(nlev); g <- 1L; acc <- 0
  for (i in seq_len(nlev)) {
    if (acc > 0 && acc + span[i] > max_span) { g <- g + 1L; acc <- 0 }
    grp[i] <- g; acc <- acc + span[i]
  }
  range_grp <- grp[sn]
  parts <- lapply(sort(unique(range_grp)), function(k) {
    idx <- which(range_grp == k)
    d <- disjoin(gr[idx], with.revmap = with.revmap, ignore.strand = ignore.strand)
    if (with.revmap)
      mcols(d)$revmap <- relist(idx[unlist(mcols(d)$revmap, use.names = FALSE)],
                                mcols(d)$revmap)
    d
  })
  do.call(c, parts)
}
