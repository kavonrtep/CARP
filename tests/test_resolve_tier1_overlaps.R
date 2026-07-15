#!/usr/bin/env Rscript
# Equivalence test for the resolve_tier1_overlaps optimization in
# make_unified_annotation.R.
#
# The old resolver ran the greedy longest-first trim over EVERY tier-1 feature in
# a batch whenever a single overlapping pair existed, and grew `kept` one element
# per iteration -> O(N^2). The new one restricts the loop to the features that
# actually overlap (non-overlapping features pass through intact), which is exact
# because a non-overlapping feature cannot overlap an overlapping one. This test
# proves the two produce the IDENTICAL SET of features (order-independent; the
# real pipeline sorts the combined annotation afterwards) across many random
# overlap-dense inputs.
#
# trim_to_nonoverlap + the two resolvers are copied verbatim from
# make_unified_annotation.R (the script's main is unguarded, so it cannot be
# sourced in isolation). Keep them in sync with the source.
suppressPackageStartupMessages(library(GenomicRanges))

# ---- copied verbatim from make_unified_annotation.R -------------------------
trim_to_nonoverlap <- function(lower, higher, min_len = 50L) {
  if (length(lower) == 0 || length(higher) == 0) return(lower)
  higher_r <- reduce(higher, ignore.strand = TRUE)
  lower_overlaps_higher <- overlapsAny(lower, higher_r, ignore.strand = TRUE)
  if (!any(lower_overlaps_higher)) return(lower)
  n_lower      <- length(lower)
  lower_plain  <- granges(lower)
  higher_plain <- granges(higher_r)
  strand(lower_plain)  <- "*"
  strand(higher_plain) <- "*"
  combined <- c(lower_plain, higher_plain)
  dis    <- disjoin(combined, with.revmap = TRUE, ignore.strand = TRUE)
  revmap <- dis$revmap
  has_higher <- any(revmap > n_lower)
  has_lower  <- any(revmap <= n_lower)
  base_mask  <- has_lower & !has_higher
  if (!any(base_mask)) return(lower[integer(0)])
  kept     <- dis[base_mask]
  kept_rev <- revmap[base_mask]
  lower_idx <- as.integer(min(kept_rev))
  source_was_trimmed <- lower_overlaps_higher[lower_idx]
  piece_widths       <- width(kept)
  keep_by_len        <- !source_was_trimmed | piece_widths >= min_len
  kept               <- kept[keep_by_len]
  lower_idx          <- lower_idx[keep_by_len]
  source_was_trimmed <- source_was_trimmed[keep_by_len]
  if (length(kept) == 0) return(lower[integer(0)])
  mcols(kept) <- mcols(lower)[lower_idx, , drop = FALSE]
  new_strand <- as.character(strand(lower))[lower_idx]
  new_strand[source_was_trimmed] <- "*"
  strand(kept) <- new_strand
  kept$revmap  <- NULL
  kept
}

# OLD resolver (pre-optimization)
resolve_old <- function(t1, min_len) {
  if (length(t1) <= 1) return(t1)
  h <- suppressWarnings(findOverlaps(t1, ignore.strand = TRUE,
                                     drop.self = TRUE, drop.redundant = TRUE))
  if (length(h) == 0) return(t1)
  t1s  <- t1[order(width(t1), decreasing = TRUE)]
  kept <- t1s[1]
  for (i in 2:length(t1s)) {
    piece <- trim_to_nonoverlap(t1s[i], kept, min_len)
    if (length(piece) > 0) kept <- suppressWarnings(c(kept, piece))
  }
  kept
}

# NEW resolver (copied from make_unified_annotation.R)
resolve_new <- function(t1, min_len) {
  if (length(t1) <= 1) return(t1)
  h <- suppressWarnings(findOverlaps(t1, ignore.strand = TRUE,
                                     drop.self = TRUE, drop.redundant = TRUE))
  if (length(h) == 0) return(t1)
  involved    <- sort(unique(c(queryHits(h), subjectHits(h))))
  passthrough <- t1[-involved]
  t1s  <- t1[involved][order(width(t1[involved]), decreasing = TRUE)]
  pieces <- vector("list", length(t1s))
  pieces[[1]] <- t1s[1]
  kept <- t1s[1]
  for (i in 2:length(t1s)) {
    piece <- trim_to_nonoverlap(t1s[i], kept, min_len)
    if (length(piece) > 0) {
      pieces[[i]] <- piece
      kept <- suppressWarnings(c(kept, piece))
    }
  }
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  suppressWarnings(c(passthrough, do.call(c, pieces)))
}

# ---- canonical, order-independent signature of a GRanges --------------------
sig <- function(gr) {
  if (length(gr) == 0) return(character(0))
  sort(sprintf("%s:%d-%d:%s|%s|%s|%s|%s",
               as.character(seqnames(gr)), start(gr), end(gr), as.character(strand(gr)),
               gr$classification, gr$Name, gr$source_tier, gr$source_tool))
}

rand_t1 <- function(rng_i, n) {
  set.seed(1000 + rng_i)
  sn  <- sample(c("chr1", "chr2", "c3"), n, replace = TRUE)
  st  <- sample.int(4000L, n, replace = TRUE)
  wd  <- sample(50:900, n, replace = TRUE)     # overlap-dense
  cls <- sample(c("Class_I/LTR/Ty1_copia/Ale", "Class_II/Subclass_1/TIR/hAT",
                  "Class_I/LINE"), n, replace = TRUE)
  gr <- GRanges(sn, IRanges(st, st + wd),
                strand = sample(c("+", "-"), n, replace = TRUE))
  mcols(gr) <- DataFrame(classification = cls, Name = cls,
                         source_tier = 1L, source_tool = "DANTE")
  gr
}

main <- function() {
  min_len <- 50L
  ntrials <- 60
  for (i in seq_len(ntrials)) {
    n <- sample(2:40, 1)
    t1 <- rand_t1(i, n)
    a <- resolve_old(t1, min_len)
    b <- resolve_new(t1, min_len)
    if (!identical(sig(a), sig(b))) {
      cat("MISMATCH at trial", i, "n=", n, "\n")
      cat("old:\n"); print(sig(a))
      cat("new:\n"); print(sig(b))
      stop("resolve_tier1_overlaps: new != old")
    }
  }
  cat(sprintf("  resolve_tier1_overlaps: new == old (%d random overlap-dense trials)\n",
              ntrials))
  cat("test_resolve_tier1_overlaps: PASSED\n")
}

main()
