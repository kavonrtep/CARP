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
# Always decomposes `lower`'s internal overlaps (disjoin lower alone when nothing
# overlaps higher) so the result is independent of the batch's feature mix. Keep
# in sync with the source.
trim_to_nonoverlap <- function(lower, higher, min_len = 50L) {
  if (length(lower) == 0) return(lower)
  n_lower     <- length(lower)
  lower_plain <- granges(lower)
  strand(lower_plain) <- "*"
  higher_r <- if (length(higher) > 0) reduce(higher, ignore.strand = TRUE) else GRanges()
  lower_overlaps_higher <- if (length(higher_r) > 0)
    overlapsAny(lower, higher_r, ignore.strand = TRUE) else logical(n_lower)
  if (any(lower_overlaps_higher)) {
    higher_plain <- granges(higher_r)
    strand(higher_plain) <- "*"
    combined <- c(lower_plain, higher_plain)
  } else {
    combined <- lower_plain
  }
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

# Deterministic longest-first order with coordinate tie-break (matches the
# production resolver). Ordering equal-width features by (seqname,start,end,
# strand,classification,source_tool) makes the greedy trim independent of input
# order -> batch/thread independent. The last two keys break the residual tie
# between genuinely identical intervals from different tools/classes.
det_order <- function(gr)
  order(-width(gr), as.character(seqnames(gr)), start(gr), end(gr),
        as.character(strand(gr)), as.character(gr$classification),
        as.character(gr$source_tool))

# OLD resolver shape (pre-optimization: grows `kept` incrementally), same
# deterministic ordering.
resolve_old <- function(t1, min_len) {
  if (length(t1) <= 1) return(t1)
  h <- suppressWarnings(findOverlaps(t1, ignore.strand = TRUE,
                                     drop.self = TRUE, drop.redundant = TRUE))
  if (length(h) == 0) return(t1)
  t1s  <- t1[det_order(t1)]
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
  inv  <- t1[involved]
  t1s  <- inv[det_order(inv)]
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
  ntrials <- 30
  for (i in seq_len(ntrials)) {
    n <- sample(2:22, 1)   # capped: the O(k^2) greedy has heavy GRanges overhead
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

  # Order-invariance (determinism): the greedy trim RESULT must not depend on the
  # order features arrive in — that order is batch/thread-dependent. Use inputs
  # dense with WIDTH TIES (where a naive width-only sort tie-breaks by position)
  # and check the resolved SET is identical under input permutation.
  n_inv <- 12L
  for (i in seq_len(n_inv)) {
    set.seed(5000 + i)
    n  <- sample(4:7, 1)
    # spread over a wider window so only some features overlap (keeps the O(k^2)
    # greedy cheap) while WIDTH TIES still exercise the tie-break path
    st <- sample.int(2600, n, replace = TRUE)
    wd <- sample(c(400L, 400L, 600L), n, replace = TRUE)  # ties
    gr <- GRanges("c1", IRanges(st, st + wd),
                  strand = sample(c("+", "-"), n, replace = TRUE))
    mcols(gr) <- DataFrame(classification = "Class_I/LTR", Name = "x",
                           source_tier = 1L, source_tool = "DANTE")
    ref  <- sig(resolve_new(gr, min_len))
    perm <- sig(resolve_new(gr[sample(seq_len(n))], min_len))  # SAME features, shuffled
    if (!identical(ref, perm)) {
      cat("ORDER-DEPENDENT at trial", i, "\n"); print(ref); print(perm)
      stop("resolve_tier1_overlaps: result depends on input order")
    }
  }
  cat(sprintf("  resolve_tier1_overlaps: order-invariant (%d tie permutation trials)\n", n_inv))

  # ── Regression: te_sat trim must run AFTER resolve, never before ──────────
  # process_batch's TE-derived-satellite pre-pass trims t1 against te_sat_r with
  # trim_to_nonoverlap, which disjoins c(lower, higher). If t1 still carries its
  # INTERNAL overlaps at that point, they are silently decomposed into disjoint
  # pieces (both overlappers kept, split at the overlap) — defeating the greedy
  # longest-first resolution. Because te_sat is non-empty only in batches that
  # contain a TE-derived-TRC array, whether a given sequence's t1 got pre-disjoined
  # depended on the thread-count batching -> non-deterministic on large genomes
  # (measured: 816 tier-1 features flipped, run116 threads 1 vs N). The fix runs
  # resolve_tier1_overlaps BEFORE the te_sat trim. This test locks that ordering.
  mk1 <- function(seq, s, e, strand, cls = "Class_I/LTR/Ty3_gypsy/chromovirus/CRM",
                  tool = "DANTE_LTR") {
    g <- GRanges(seq, IRanges(s, e), strand = strand)
    mcols(g) <- DataFrame(classification = cls, Name = cls,
                          source_tier = 1L, source_tool = tool)
    g
  }
  keep_c1 <- function(g) sig(g[as.character(seqnames(g)) == "c1"])
  # On c1: two overlapping tier-1 elements (the CFRUME023 shape) — a 7031 bp and a
  # wider 9631 bp CRM element overlapping in 5557-7405. These are FAR from any
  # te_sat and never overlap one.
  pair <- suppressWarnings(c(mk1("c1", 375, 7405, "+"), mk1("c1", 5557, 15187, "+")))
  # On c2: a structural TE that DOES sit under a te_sat array. Its presence makes
  # trim_to_nonoverlap's global `any(lower_overlaps_higher)` TRUE, so the disjoin
  # path runs over the WHOLE batch — decomposing c1's untouched overlapping pair
  # if the trim runs before the greedy. This is exactly how a te_sat elsewhere in
  # a single (thread=1) batch corrupted far-away sequences' tier-1 resolution.
  trigger  <- mk1("c2", 50100, 50500, "+")
  t1       <- suppressWarnings(c(pair, trigger))
  te_sat_r <- reduce(granges(mk1("c2", 50000, 51000, "+")), ignore.strand = TRUE)

  # BUGGY order (trim then resolve): the global disjoin pre-splits c1's pair -> 3.
  buggy <- resolve_new(trim_to_nonoverlap(t1, te_sat_r, min_len), min_len)
  # FIXED order (resolve then trim): greedy longest-first on c1 -> 2 pieces.
  fixed <- trim_to_nonoverlap(resolve_new(t1, min_len), te_sat_r, min_len)
  # Baseline: no te_sat present at all (a batch with no TE-derived array).
  none  <- resolve_new(pair, min_len)

  if (length(fixed[as.character(seqnames(fixed)) == "c1"]) != 2)
    stop(sprintf("regression: fixed order must yield 2 greedy c1 pieces, got %d",
                 length(fixed[as.character(seqnames(fixed)) == "c1"])))
  # The whole point of the fix: c1's tier-1 result is INVARIANT to whether the
  # batch also contains a te_sat-covered TE elsewhere (i.e. to batching).
  if (!identical(keep_c1(fixed), keep_c1(none)))
    stop("regression: c1 tier-1 result changed when a te_sat exists elsewhere in the batch")
  # And it demonstrably differs from the buggy pre-disjoin (guards against a
  # future refactor silently reintroducing the trim-before-resolve ordering).
  if (identical(keep_c1(buggy), keep_c1(fixed)))
    stop("regression: buggy pre-disjoin should over-fragment c1 vs the fixed order")
  cat(sprintf("  te_sat-trim-after-resolve: fixed c1=%d pieces (batch-invariant), buggy c1=%d\n",
              length(fixed[as.character(seqnames(fixed)) == "c1"]),
              length(buggy[as.character(seqnames(buggy)) == "c1"])))

  cat("test_resolve_tier1_overlaps: PASSED\n")
}

main()
