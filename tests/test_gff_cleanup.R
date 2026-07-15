#!/usr/bin/env Rscript
# Equivalence test for the gff_cleanup rewrite in clean_rm_output.R and
# merge_repeat_annotations.R (identical function in both).
#
# The old version did as.list(gff_disjoin$revmap) — tens of millions of tiny R
# vectors — inside an mclapply that COW-forked the whole GRanges into every
# worker (296 GB max_rss on a 3.9 Gb genome). The new version collapses per-range
# Names/strands with a single vectorized extractList + unstrsplit that stays in
# the compressed S4Vectors list domain. This test proves the two produce an
# identical GRanges (ranges + all mcols) across random overlap-dense inputs whose
# disjoint ranges mix several classifications (exercising the resolve_name LCA).
#
# Functions are copied from the source (scripts run main at top level, so they
# cannot be sourced in isolation). Keep resolve_name / the two gff_cleanup bodies
# in sync with the source.
suppressPackageStartupMessages({ library(rtracklayer); library(parallel) })

resolve_name <- function(x){
  if (length(x)==1){ return(x) } else {
    y <- sapply(x, strsplit, split="/", fixed = TRUE)
    ny <- table(unlist(sapply(y, function(x)paste(seq_along(x), x))))
    if (max(ny)<length(x)){ return("Unknown") } else {
      k <- which(ny==length(x))
      r <- max(as.numeric((gsub(" .+", "", names(k)))))
      out <- paste(y[[1]][1:r], collapse="/"); return(out)
    }
  }
}

gff_cleanup_old <- function(gff, mc.cores = 1L){
  gff_disjoin <- disjoin(gff, with.revmap=TRUE)
  gff_annot <- mclapply(as.list(gff_disjoin$revmap),
                        FUN = function(x) list(name = gff$Name[x],
                                               strand = strand(gff[x])),
                        mc.cores = mc.cores)
  gff_names <- lapply(gff_annot, `[[`, "name")
  gff_strands <- lapply(gff_annot, `[[`, "strand")
  new_annot <- sapply(sapply(gff_names, unique), paste, collapse="|")
  new_annot_uniq <- unique(new_annot)
  lca_annot <- sapply(strsplit(new_annot_uniq, "|", fixed = TRUE), resolve_name)
  names(lca_annot) <- new_annot_uniq
  new_annot_lca <- lca_annot[new_annot]
  strand_attribute <- sapply(sapply(gff_strands, unique), paste, collapse="|")
  gff_disjoin$source <- "RM"; gff_disjoin$type <- "repeat"
  gff_disjoin$score <- NA;    gff_disjoin$phase <- NA
  gff_disjoin$Name <- new_annot_lca
  gff_disjoin$Original_names <- new_annot
  gff_disjoin$strands <- strand_attribute
  gff_disjoin$revmap <- NULL
  gff_disjoin
}

gff_cleanup_new <- function(gff){
  gff_disjoin <- disjoin(gff, with.revmap=TRUE)
  names_cl  <- extractList(gff$Name, gff_disjoin$revmap)
  new_annot <- unstrsplit(unique(names_cl), sep = "|")
  strand_cl <- extractList(as.character(strand(gff)), gff_disjoin$revmap)
  strand_attribute <- unstrsplit(unique(strand_cl), sep = "|")
  new_annot_uniq <- unique(new_annot)
  lca_annot <- sapply(strsplit(new_annot_uniq, "|", fixed = TRUE), resolve_name)
  names(lca_annot) <- new_annot_uniq
  new_annot_lca <- lca_annot[new_annot]
  gff_disjoin$source <- "RM"; gff_disjoin$type <- "repeat"
  gff_disjoin$score <- NA;    gff_disjoin$phase <- NA
  gff_disjoin$Name <- new_annot_lca
  gff_disjoin$Original_names <- new_annot
  gff_disjoin$strands <- strand_attribute
  gff_disjoin$revmap <- NULL
  gff_disjoin
}

sig <- function(gr) {
  sprintf("%s:%d-%d|%s|%s|%s|%s|%s",
          as.character(seqnames(gr)), start(gr), end(gr),
          gr$Name, gr$Original_names, gr$strands, gr$source, gr$type)
}

NAMES <- c("Class_I/LTR/Ty1_copia/Ale", "Class_I/LTR/Ty1_copia/Angela",
           "Class_I/LTR/Ty3_gypsy/Tekay", "Class_II/Subclass_1/TIR/hAT",
           "Simple_repeat", "Class_I/LTR/Ty1_copia")

main <- function() {
  for (i in seq_len(40)) {
    set.seed(2000 + i)
    n  <- sample(3:40, 1)
    sn <- sample(c("chr1", "chr2", "c3"), n, replace = TRUE)
    st <- sample.int(300L, n, replace = TRUE)          # dense -> many-name overlaps
    wd <- sample(10:120, n, replace = TRUE)
    gr <- GRanges(sn, IRanges(st, st + wd))            # strand "*" (as in the scripts)
    gr$Name <- sample(NAMES, n, replace = TRUE)
    # scripts pre-merge same-Name neighbours before gff_cleanup:
    gr <- unlist(reduce(split(gr, gr$Name))); gr$Name <- names(gr)
    a <- gff_cleanup_old(gr, mc.cores = 1L)
    b <- gff_cleanup_new(gr)
    if (!identical(sig(a), sig(b))) {
      cat("MISMATCH trial", i, "\n"); print(sig(a)); print(sig(b))
      stop("gff_cleanup: new != old")
    }
  }
  cat("  gff_cleanup: new == old (40 random overlap-dense trials, LCA exercised)\n")
  cat("test_gff_cleanup: PASSED\n")
}
main()
