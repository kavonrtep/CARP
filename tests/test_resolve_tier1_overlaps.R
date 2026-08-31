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

# NEW resolver: loaded from the REAL script, not copied.
#
# This used to be a verbatim copy with a "keep in sync" comment. That silently
# stopped testing the shipped code — the candidate-bounded rewrite of
# resolve_tier1_overlaps (the fix for the 94 Gbp run-000156 blowup, where one
# batch spent 4.7 h in this function) left the copy untouched and the test still
# passed. Extract the real definition instead: parse the script and eval only the
# top-level assignments we need, so the script's unguarded main never runs.
load_from_script <- function(names_wanted) {
  here <- tryCatch(dirname(normalizePath(sub("^--file=", "",
             grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
           error = function(e) ".")
  src <- file.path(here, "..", "scripts", "make_unified_annotation.R")
  if (!file.exists(src)) stop("cannot locate make_unified_annotation.R at ", src)
  env <- new.env(parent = globalenv())
  found <- character(0)
  for (e in parse(src)) {
    if (is.call(e) && length(e) >= 3 &&
        as.character(e[[1]])[1] %in% c("<-", "=") &&
        is.name(e[[2]]) && as.character(e[[2]]) %in% names_wanted) {
      eval(e, envir = env)
      found <- c(found, as.character(e[[2]]))
    }
  }
  missing <- setdiff(names_wanted, found)
  if (length(missing) > 0)
    stop("not found in make_unified_annotation.R: ", paste(missing, collapse = ", "))
  env
}

.real <- load_from_script(c("trim_to_nonoverlap", "resolve_tier1_overlaps"))
resolve_new <- .real$resolve_tier1_overlaps
# The naive reference above must trim with the SAME primitive as the real
# resolver, or a difference in trim_to_nonoverlap would masquerade as a resolver
# difference. Point the local copy at the real one.
trim_to_nonoverlap <- .real$trim_to_nonoverlap

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

# ---- observed vs inferred extent (claim_width) -------------------------------
# DANTE_LINE and the DANTE_TIR fallback INFER most of their extent by aligning
# the flanks of other copies; DANTE_LTR and primary DANTE_TIR observe theirs
# structurally. The greedy is longest-first, so without claim_width an inferred
# flank outranks an observed element purely on length it never observed.
mk_cw <- function(s, e, tool, cw) {
  g <- GRanges("chr1", IRanges(s, e))
  g$Name <- tool; g$classification <- tool; g$source_tier <- 1L
  g$source_tool <- tool; g$claim_width <- as.integer(cw); g
}
span_of <- function(gr, tool) sum(width(gr[gr$source_tool == tool]))

check_claim_width <- function(min_len) {
  # 1. a chimeric LINE (22 kb span, 2 kb observed core) must not displace a
  #    structurally delimited 10 kb DANTE_LTR element
  r <- resolve_new(c(mk_cw(1, 22000, "DANTE_LINE", 2000),
                     mk_cw(15000, 25000, "DANTE_LTR", 10001)), min_len)
  if (span_of(r, "DANTE_LTR") != 10001L)
    stop("inferred LINE flank displaced an observed DANTE_LTR element")

  # 2. but an OBSERVED LINE core still beats a smaller competitor
  r <- resolve_new(c(mk_cw(1, 5000, "DANTE_LINE", 4000),
                     mk_cw(4500, 5200, "DANTE_LTR", 701)), min_len)
  if (span_of(r, "DANTE_LINE") < 4000L)
    stop("an observed LINE core was carved by a smaller competitor")

  # 3. primary DANTE_TIR (observed) over a fallback element (inferred flank)
  r <- resolve_new(c(mk_cw(1, 20000, "DANTE_TIR", 20000),
                     mk_cw(15000, 30000, "DANTE_TIR", 1200)), min_len)
  if (max(width(r)) != 20000L)
    stop("an inferred TIR fallback flank displaced a primary DANTE_TIR element")

  # 4. absent claim_width must reproduce the width-only behaviour exactly --
  #    otherwise this whole check could pass vacuously on inputs that never set it
  a <- c(mk_cw(1, 22000, "DANTE_LINE", 2000), mk_cw(15000, 25000, "DANTE_LTR", 10001))
  b <- a; b$claim_width <- NULL
  if (identical(sig(resolve_new(a, min_len)), sig(resolve_new(b, min_len))))
    stop("claim_width had no effect - the ordering key is not being read")
  if (!identical(sig(resolve_new(b, min_len)), sig(resolve_old(b, min_len))))
    stop("without claim_width the resolver diverged from the width-only reference")

  # 5. still order-invariant when claim_width drives the ordering
  base <- c(mk_cw(1, 22000, "DANTE_LINE", 2000),
            mk_cw(15000, 25000, "DANTE_LTR", 10001),
            mk_cw(21000, 24000, "DANTE_TIR", 3001))
  ref <- sig(resolve_new(base, min_len))
  for (perm in list(c(2, 1, 3), c(3, 2, 1), c(2, 3, 1))) {
    if (!identical(sig(resolve_new(base[perm], min_len)), ref))
      stop("claim_width ordering is not order-invariant")
  }
  cat("  claim_width: inferred flank yields to observed extent; observed core kept;\n")
  cat("               absent claim_width == width-only reference; order-invariant\n")
}

# ---- subset_seqs drops list-valued columns BEFORE the row subset -------------
# A GFF3 attribute containing a comma (GFF3's multi-value separator) imports as a
# CompressedCharacterList mcols column. Subsetting a DataFrame that holds one dies
# on the Bioconductor 3.14 stack the container pins via r-base 4.1
# (S4Vectors 0.32.4) -- "'end' must be <= 'length(x)'" -- but not on 3.18, so a
# newer dev environment cannot see it. subset_seqs must therefore drop non-standard
# columns BEFORE it subsets rows, never after.
check_subset_seqs_drops_first <- function() {
  # load the function AND the column whitelist it closes over, into one env
  env <- load_from_script(c("subset_seqs", ".META_COLS"))
  subset_seqs <- env$subset_seqs
  meta_cols   <- env$.META_COLS

  gr <- GRanges(c("chr1", "chr1", "chr2"), IRanges(c(1, 100, 1), c(50, 150, 50)))
  gr$Name <- c("a", "b", "c"); gr$classification <- gr$Name
  gr$source_tier <- 1L; gr$source_tool <- "DANTE_LINE"
  # the shape a comma-valued GFF3 attribute takes after rtracklayer import
  gr$Extension_capped <- IRanges::CharacterList(c("5prime", "3prime"),
                                                character(0), "5prime")
  out <- subset_seqs(gr, "chr1")
  if (length(out) != 2L) stop("subset_seqs returned the wrong number of features")
  if ("Extension_capped" %in% colnames(mcols(out)))
    stop("subset_seqs kept a non-standard, list-valued column")
  if (!all(colnames(mcols(out)) %in% meta_cols))
    stop("subset_seqs kept a column outside .META_COLS")
  # The drop must happen BEFORE the row subset. On Bioconductor 3.18 both orders
  # work, so assert the ordering in the source rather than by behaviour.
  body_txt <- paste(deparse(subset_seqs), collapse = "\n")
  drop_at <- regexpr("keep_cols", body_txt)
  subs_at <- regexpr("seqnames\\(gr\\) %in% seqs", body_txt)
  if (drop_at < 0 || subs_at < 0 || drop_at > subs_at)
    stop("subset_seqs must drop non-standard mcols BEFORE subsetting rows")
  cat("  subset_seqs: list-valued columns dropped before the row subset\n")
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

  check_claim_width(min_len)
  check_subset_seqs_drops_first()

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
