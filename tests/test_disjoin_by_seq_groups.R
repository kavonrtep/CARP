#!/usr/bin/env Rscript
# disjoin_by_seq_groups() (scripts/disjoin_utils.R) must return EXACTLY what
# disjoin() returns — ranges, order, strand, Seqinfo and revmap.
#
# The helper exists because IRanges' disjoin falls back to a per-sequence R loop
# once the summed span of all sequences exceeds .Machine$integer.max (merge_rm_and_dante
# spent ~4,800 s there on a 22.4 Gbp, 305k-sequence assembly). It disjoins
# consecutive groups of sequences that each fit the fast path and concatenates.
#
# Here the inputs stay small so the reference disjoin() is on its fast path, and
# a tiny `max_span` forces the helper to split into many groups. Covers mixed
# strands, ignore.strand, unused seqlevels, a seqlevel order unlike the data,
# a sequence larger than max_span (its own group), and empty input.
suppressPackageStartupMessages(library(GenomicRanges))
source(file.path("scripts", "disjoin_utils.R"))

set.seed(42)
fails <- 0L
check <- function(label, a, b) {
  if (!identical(a, b)) {
    fails <<- fails + 1L
    cat("FAIL", label, "\n"); print(all.equal(a, b))
  }
}

make <- function(ns, n, span, strands = c("+", "-", "*")) {
  lv <- sprintf("JBDLLP01%07d.1", sample(1e6, ns))
  sq <- sample(lv[-(1:3)], n, replace = TRUE)                # 3 unused seqlevels
  GRanges(sq, IRanges(sample(span, n, TRUE), width = sample(1:400, n, TRUE)),
          strand = sample(strands, n, TRUE), seqinfo = Seqinfo(sample(lv)))
}

cases <- list(
  mixed   = make(300, 6000, 20000),
  star    = make(300, 6000, 20000, "*"),
  one_big = suppressWarnings(c(make(50, 2000, 5000),
              GRanges("big", IRanges(sample(3e5, 500, TRUE), width = 300))))
)
for (nm in names(cases)) {
  g <- cases[[nm]]
  for (rev in c(FALSE, TRUE)) for (ign in c(FALSE, TRUE)) {
    lbl <- sprintf("%s revmap=%s ignore.strand=%s", nm, rev, ign)
    check(lbl, disjoin(g, with.revmap = rev, ignore.strand = ign),
          disjoin_by_seq_groups(g, with.revmap = rev, ignore.strand = ign,
                                max_span = 1e5))
  }
}

# The grouped path must actually be taken (else the comparison is vacuous):
# with max_span = 1e5 the 300 x ~20 kb sequences need many groups.
g <- cases$mixed
spans <- tapply(end(g), as.character(seqnames(g)), max) -
         tapply(start(g), as.character(seqnames(g)), min) + 2
if (sum(spans) <= 1e5) { fails <- fails + 1L; cat("FAIL fixture does not force grouping\n") }

# Default max_span on input that already fits -> plain disjoin path.
check("fits fast path", disjoin(g, with.revmap = TRUE),
      disjoin_by_seq_groups(g, with.revmap = TRUE))

# Empty input.
e <- GRanges(seqinfo = Seqinfo(c("a", "b")))
check("empty", disjoin(e, with.revmap = TRUE), disjoin_by_seq_groups(e, with.revmap = TRUE))

if (fails > 0) stop(sprintf("test_disjoin_by_seq_groups: %d FAILED", fails))
cat("test_disjoin_by_seq_groups: PASSED\n")
