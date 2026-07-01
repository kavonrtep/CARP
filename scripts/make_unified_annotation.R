#!/usr/bin/env Rscript
# make_unified_annotation.R
# Produces a unified, tier-prioritised repeat annotation GFF3 from all pipeline layers.
# See annotation_rules.md and unified_annotation_implementation_plan.md for design rationale.

suppressPackageStartupMessages({
  library(optparse)
  library(rtracklayer)
  library(GenomicRanges)
  library(parallel)
})

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "classification.R"))

# ── Progress / timing helpers ─────────────────────────────────────────────────
# Simple wall-clock timers. All messages go to stderr so they interleave
# correctly with mclapply worker output.
.script_start <- proc.time()[3]
elapsed <- function() sprintf("[+%6.1fs]", proc.time()[3] - .script_start)

log_msg <- function(...) message(elapsed(), " ", ...)

# Time a single expression: x <- timed("label", expr); x is the result.
timed <- function(label, expr) {
  t0 <- proc.time()[3]
  res <- force(expr)
  dt <- proc.time()[3] - t0
  log_msg(sprintf("%6.2fs  %s", dt, label))
  res
}

# Start/stop pair for multi-line operations:
#   t0 <- tic("label"); ...; toc(t0)
tic <- function(label) {
  log_msg("▶ ", label)
  list(label = label, t0 = proc.time()[3])
}
toc <- function(tic_state, detail = NULL) {
  dt <- proc.time()[3] - tic_state$t0
  if (is.null(detail))
    log_msg(sprintf("◀ %s  (%.2fs)", tic_state$label, dt))
  else
    log_msg(sprintf("◀ %s  (%.2fs) — %s", tic_state$label, dt, detail))
  invisible(dt)
}

# ── 0. CLI arguments ─────────────────────────────────────────────────────────

option_list <- list(
  make_option("--ltr",              type="character", help="DANTE_LTR GFF3"),
  make_option("--ltr_tandems",      type="character", default=NULL,
              help="DANTE_LTR_tandems.gff3 (tandem LTR-RT containers; optional)"),
  make_option("--tir",              type="character", help="DANTE_TIR GFF3"),
  make_option("--line",             type="character", help="DANTE_LINE GFF3"),
  make_option("--dante",            type="character", help="DANTE filtered domains GFF3"),
  make_option("--tc_default",       type="character", help="TideCluster default GFF3"),
  make_option("--tc_short",         type="character", help="TideCluster short monomer GFF3"),
  make_option("--tc_rm",            type="character", help="RM on TideCluster library GFF3"),
  make_option("--tc_rdna_default",  type="character", default=NULL,
              help="TideCluster default rDNA TSV (TRC/rDNA_type/coverage; authoritative rDNA-TRC calls; optional)"),
  make_option("--tc_rdna_short",    type="character", default=NULL,
              help="TideCluster short-monomer rDNA TSV (authoritative rDNA-TRC calls; optional)"),
  make_option("--rm",               type="character", help="RepeatMasker+DANTE merged GFF3"),
  make_option("--th_default",       type="character", help="TideHunter default residuals GFF3"),
  make_option("--th_short",         type="character", help="TideHunter short residuals GFF3"),
  make_option("--th_raw_default",   type="character", default="",
              help="Raw TideHunter default GFF3 (independent tandem evidence for the RM-on-TideCluster tandem gate; optional)"),
  make_option("--th_raw_short",     type="character", default="",
              help="Raw TideHunter short GFF3 (independent tandem evidence; optional)"),
  make_option("--rm_tc_tandem_gate", type="character", default="TRUE",
              help="When TRUE, a Tier-4 RM-on-TideCluster satellite may override a Tier-5 TE call only where it has independent tandem evidence (raw TideHunter); unsupported RM_TC is demoted below the TE. [TRUE]"),
  make_option("--fai",              type="character", help="Genome FAI file"),
  make_option("--output",           type="character", help="Output unified GFF3"),
  make_option("--threads",          type="integer",   default=4L,
              help="Number of parallel threads [4]"),
  make_option("--batch_size",       type="double",    default=200e6,
              help="Target genome bp per processing batch [200000000]"),
  make_option("--chunk_threshold",  type="double",    default=500e6,
              help="Genome bp above which chunking is activated [500000000]"),
  make_option("--min_feature_length", type="integer", default=50L,
              help="Min bp for a trimmed feature fragment to be retained [50]")
)
opt <- parse_args(OptionParser(option_list=option_list))

for (arg in c("ltr","tir","line","dante","tc_default","tc_short","tc_rm",
              "rm","th_default","th_short","fai","output")) {
  if (is.null(opt[[arg]])) stop("Missing required argument: --", arg)
}

# ── 1. Utility functions ──────────────────────────────────────────────────────
# Classification handling lives in scripts/classification.R — canonicalise()
# replaces the former fix_sep / convert_tir_cls helpers (see CLAUDE.md).

# Attach standard metadata columns to a GRanges
set_meta <- function(gr, name_vec, cls_vec, tier, tool) {
  gr$Name           <- name_vec
  gr$classification <- cls_vec
  gr$source_tier    <- tier
  gr$source_tool    <- tool
  gr
}

# Import GFF3 safely; return empty GRanges for empty/missing/error files.
# Does NOT pre-scan lines: DANTE_filtered.gff3 has ~800 comment lines at the
# top, so any fixed-N scan would incorrectly classify it as header-only.
safe_import <- function(path) {
  if (is.null(path) || !file.exists(path) || file.size(path) == 0)
    return(GRanges())
  tryCatch(
    import.gff3(path),
    error = function(e) {
      message("  Warning: could not import ", basename(path), ": ", conditionMessage(e))
      GRanges()
    }
  )
}

# Standard column names that every GRanges in process_batch carries.
# All other tool-specific columns are dropped before combining to avoid
# CharacterList vs character schema incompatibilities that cause NSBS errors.
.META_COLS <- c("ID", "Name", "classification", "source_tier",
                "source_tool", "element_type", "TE_origin", "TE_origin_structure",
                "structure", "copy_number")

# Subset a GRanges to a set of sequence names.
# Also standardizes seqlevels and mcols schema so that c() / reduce()
# across GRanges from different tools never triggers seqinfo/NSBS errors.
# Ensures all retained columns are plain atomic vectors (not CharacterList).
subset_seqs <- function(gr, seqs) {
  if (length(gr) == 0) return(gr)
  gr <- gr[as.character(seqnames(gr)) %in% seqs]
  if (length(gr) == 0) return(gr)
  # Rebuild mcols with only the standardised columns, coerced to plain types
  md <- list()
  for (col in .META_COLS) {
    if (col %in% colnames(mcols(gr))) {
      val <- mcols(gr)[[col]]
      md[[col]] <- if (col == "source_tier") as.integer(val) else as.character(val)
    }
  }
  # Rebuild with the standardised columns; a track with none (e.g. the th_raw
  # tandem-evidence intervals) gets its mcols cleared — assigning an empty
  # DataFrame() to a non-zero-length object errors, so use NULL.
  if (length(md) > 0) {
    mcols(gr) <- do.call(S4Vectors::DataFrame, md)
  } else {
    mcols(gr) <- NULL
  }
  # Standardize seqlevels to the batch set
  suppressWarnings(seqlevels(gr, pruning.mode = "coarse") <- seqs)
  gr
}

# ── 2. Per-tool loading functions ─────────────────────────────────────────────

# DANTE_LTR.gff3 is read untouched (all individual elements). Tandem LTR-RT
# (LTR_RT_TR) arrays come from the small companion file written by
# resolve_ltr_tandems.py: one container per array carrying `members=` (the member
# element IDs). We split the individual elements into members (IDs listed in a
# container) vs standalone; containers + standalone are Level-1 tier-1, and the
# member copies become Level-2 children nested under their container (parent
# resolved by overlap in finalise_output).
load_tier1_ltr <- function(path, tandems_path = NULL) {
  message("Loading DANTE_LTR: ", path)
  raw <- safe_import(path)
  if (length(raw) == 0)
    return(list(top = GRanges(), children = GRanges(), members = GRanges()))

  all_te   <- raw[raw$type == "transposable_element"]
  children <- raw[raw$type %in% c("protein_domain", "long_terminal_repeat",
                                   "target_site_duplication", "primer_binding_site")]

  # ── tandem containers (small side file) ────────────────────────────────────
  containers <- GRanges()
  member_ids <- character(0)
  if (!is.null(tandems_path)) {
    tand <- safe_import(tandems_path)
    if (length(tand) > 0) {
      tand <- tand[tand$type == "transposable_element"]
      # rtracklayer parses the comma-separated members= attribute into a
      # CharacterList; flatten to the set of member element IDs either way.
      m <- tand$members
      member_ids <- if (is.null(m)) character(0)
                    else if (is(m, "CharacterList")) unique(unlist(m))
                    else unique(unlist(strsplit(as.character(m), ",", fixed = TRUE)))
      ccls <- canonicalise(tand$Final_Classification, source = "DANTE_LTR")
      cn   <- if (!is.null(tand$copy_number)) as.character(tand$copy_number)
              else rep(NA_character_, length(tand))
      containers <- set_meta(tand, ccls, ccls, 1L, "DANTE_LTR")
      containers$members      <- NULL   # CharacterList attr — used above, drop before c()
      containers$element_type <- NA_character_
      containers$structure    <- "LTR_RT_TR"
      containers$copy_number  <- cn
    }
  }

  is_member  <- as.character(all_te$ID) %in% member_ids
  members    <- all_te[is_member]
  standalone <- all_te[!is_member]

  meta_el <- function(gr) {   # individual elements (standalone or tandem members)
    if (length(gr) == 0) return(gr)
    cls <- canonicalise(gr$Final_Classification, source = "DANTE_LTR")
    gr  <- set_meta(gr, cls, cls, 1L, "DANTE_LTR")
    gr$element_type <- ifelse(grepl("^TE_partial_", gr$ID, perl = TRUE), "partial", "complete")
    gr$structure    <- NA_character_
    gr$copy_number  <- NA_character_
    gr
  }
  standalone <- meta_el(standalone)
  members    <- meta_el(members)
  top        <- suppressWarnings(c(standalone, containers))

  # Strip verbose alignment attributes to reduce memory
  verbose_attrs <- c("DB_Seq", "Region_Seq", "Query_Seq", "Best_Hit_DB_Pos",
                     "cls_prefilter", "neighbors_count",
                     "upstream_domain", "downstream_domain", "domain_order")
  for (a in verbose_attrs) mcols(children)[[a]] <- NULL

  message("  ", length(standalone), " standalone + ", length(containers),
          " tandem LTR_RT_TR container(s); ", length(members),
          " member copies -> Level 2")
  list(top = top, children = children, members = members)
}

load_tier1_tir <- function(path) {
  message("Loading DANTE_TIR: ", path)
  raw <- safe_import(path)
  if (length(raw) == 0) return(GRanges())

  # DANTE_TIR_combined.gff3 bundles the primary elements plus fallback
  # elements and their protein_domain children. We only want the top-level
  # TIR annotations (sequence_feature rows); protein_domain children carry
  # Final_Classification, not Classification, and would yield NAs otherwise.
  raw <- raw[raw$type == "sequence_feature"]
  if (length(raw) == 0) return(GRanges())

  raw$type <- "transposable_element"
  cls <- canonicalise(raw$Classification, source = "DANTE_TIR")
  raw <- set_meta(raw, cls, cls, 1L, "DANTE_TIR")
  message("  ", length(raw), " TIR elements")
  raw
}

load_tier1_line <- function(path) {
  message("Loading DANTE_LINE: ", path)
  raw <- safe_import(path)
  if (length(raw) == 0) return(list(top = GRanges(), children = GRanges()))

  top      <- raw[raw$type == "LINE_element"]
  children <- raw[raw$type == "protein_domain"]
  top$type <- "transposable_element"
  cls <- rep("Class_I/LINE", length(top))
  top <- set_meta(top, cls, cls, 1L, "DANTE_LINE")
  message("  ", length(top), " LINE elements")
  list(top = top, children = children)
}

load_tier2_dante <- function(path) {
  message("Loading DANTE filtered: ", path)
  raw <- safe_import(path)
  if (length(raw) == 0) return(GRanges())

  raw <- raw[raw$type == "protein_domain"]
  cls <- canonicalise(raw$Final_Classification, source = "DANTE")
  raw <- set_meta(raw, cls, cls, 2L, "DANTE")
  message("  ", length(raw), " domain features")
  raw
}

# Authoritative TRC -> rDNA-class map from TideCluster's `<prefix>_rdna.tsv`
# (columns TRC / rDNA_type / coverage), written by identify_rdna() in tc_utils.py.
# This is TideCluster's definitive per-TRC rDNA call and is what we trust. It is
# strictly more complete than the clustering GFF3's `rDNA_type` attribute: a TRC
# can be called rDNA from its consensus dimer yet carry no rDNA_type-tagged
# feature in the clustering GFF3, and the RM-on-TideCluster GFF3 (tier 4) carries
# no `rDNA_type` at all — so without this map every tier-4 rDNA array would be
# mislabelled `Satellite/TideCluster/<TRC>`. Returns a named character vector
# TRC_<n> -> "rDNA/45S_rDNA"/"rDNA/5S_rDNA"; empty when the file is absent/unreadable
# (e.g. `--no_rdna`, detection failed, or older TideCluster), in which case
# normalise_tc_satellite falls back to the per-feature `rDNA_type` attribute.
load_rdna_map <- function(path) {
  if (is.null(path) || !file.exists(path) || file.size(path) == 0)
    return(character(0))
  tab <- tryCatch(
    read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
               colClasses = "character", quote = "", comment.char = ""),
    error = function(e) {
      message("  Warning: could not read rDNA TSV ", basename(path), ": ",
              conditionMessage(e))
      NULL
    })
  if (is.null(tab) || nrow(tab) == 0 ||
      !all(c("TRC", "rDNA_type") %in% names(tab)))
    return(character(0))
  rtype <- toupper(trimws(tab$rDNA_type))
  cls   <- ifelse(rtype == "45S", "rDNA/45S_rDNA",
           ifelse(rtype == "5S",  "rDNA/5S_rDNA", NA_character_))
  keep  <- !is.na(cls)
  m <- setNames(cls[keep], trimws(tab$TRC[keep]))
  m[!duplicated(names(m))]
}

# Map a TideCluster clustering / RM-on-TC feature to (Name, classification).
# rDNA arrays (45S/5S) are surfaced as array-level rDNA/45S_rDNA / rDNA/5S_rDNA (no internal
# 18S/ITS/5.8S/IGS/25S substructure, matching TideCluster's design), which routes
# them to the rDNA class downstream (calculate_statistics_and_make_groups.R keys
# on `^rDNA`) and keeps them out of the Tandem_repeats aggregation. The rDNA call
# comes from the authoritative `rdna_map` (TideCluster's `<prefix>_rdna.tsv`),
# applied identically to tier-3 clustering and tier-4 RM-on-TC features by TRC id
# (RM-on-TC reannotates with the *default* dimer library, so it shares the default
# run's TRC namespace). The per-feature `rDNA_type` attribute is a fallback for
# when the TSV is unavailable. All other TRCs keep Satellite/TideCluster/<TRC>.
normalise_tc_satellite <- function(gr, tier, tool, rdna_map = character(0)) {
  if (length(gr) == 0) return(gr)
  gr$type <- "repeat_region"
  trc   <- as.character(gr$Name)
  cls   <- paste0("Satellite/TideCluster/", trc)
  # Authoritative rDNA call, keyed by TRC id.
  if (length(rdna_map)) {
    hit <- trc %in% names(rdna_map)
    cls[hit] <- unname(rdna_map[trc[hit]])
  }
  # Fallback to the per-feature rDNA_type attribute only where the map didn't
  # already resolve the TRC (TSV absent / older TideCluster / --no_rdna).
  rtype <- if (!is.null(gr$rDNA_type)) as.character(gr$rDNA_type) else rep(NA_character_, length(gr))
  unresolved <- !grepl("^rDNA_", cls)
  cls[unresolved & !is.na(rtype) & rtype == "45S"] <- "rDNA/45S_rDNA"
  cls[unresolved & !is.na(rtype) & rtype == "5S"]  <- "rDNA/5S_rDNA"
  # Name ALWAYS keeps the bare TRC id (TRC_<n>) — downstream apps and
  # split_gff_by_name.R (--name-prefix TRC_) key on it, so it must stay byte-stable
  # across versions for every satellite (rDNA included). The rDNA / TE-derived
  # distinction lives in `classification` (rDNA/45S_rDNA|rDNA/5S_rDNA) and the TE_origin
  # attribute, never in Name. Plain-satellite records are therefore identical to
  # the previous release; rDNA only changes `classification`. The stats splitter
  # routes rDNA by classification (it carries Name=TRC_<n>).
  gr <- set_meta(gr, trc, cls, tier, tool)
  gr
}

load_tier3_tidecluster <- function(path_default, path_short,
                                   rdna_default = character(0),
                                   rdna_short   = character(0)) {
  message("Loading TideCluster default: ", path_default)
  def <- safe_import(path_default)
  message("Loading TideCluster short monomer: ", path_short)
  sho <- safe_import(path_short)
  d <- normalise_tc_satellite(def, 3L, "TideCluster_default", rdna_default)
  s <- normalise_tc_satellite(sho, 3L, "TideCluster_short",   rdna_short)
  message("  ", if (length(d)) sum(grepl("^rDNA_", d$classification)) else 0,
          " default + ",
          if (length(s)) sum(grepl("^rDNA_", s$classification)) else 0,
          " short rDNA-labelled arrays")
  list(default = d, short = s)
}

load_tier4_tc_rm <- function(path, rdna_default = character(0)) {
  message("Loading RM on TideCluster library: ", path)
  raw <- safe_import(path)
  if (length(raw) == 0) return(GRanges())

  # RM-on-TC reannotates with the *default* dimer library, so its TRC namespace
  # is the default clustering run's — use the default rDNA map here.
  raw <- normalise_tc_satellite(raw, 4L, "TideCluster_RM", rdna_default)
  message("  ", length(raw), " RM-on-TideCluster features")
  raw
}

load_tier5_rm <- function(path) {
  message("Loading RepeatMasker+DANTE: ", path)
  raw <- safe_import(path)
  if (length(raw) == 0) return(list(te = GRanges(), simple = GRanges()))

  # Name already uses / separator (processed by clean_rm_output.R)
  raw <- set_meta(raw, as.character(raw$Name), as.character(raw$Name), 5L, "RepeatMasker")
  is_simple <- grepl("^Simple_repeat|^Low_complexity", raw$Name)
  message("  ", sum(!is_simple), " TE features, ", sum(is_simple), " Simple/Low features")
  list(te = raw[!is_simple], simple = raw[is_simple])
}

load_tier6_tidehunter <- function(path_default, path_short) {
  message("Loading TideHunter default: ", path_default)
  def <- safe_import(path_default)
  message("Loading TideHunter short: ", path_short)
  sho <- safe_import(path_short)

  normalise_th <- function(gr) {
    if (length(gr) == 0) return(gr)
    gr$type <- "repeat_region"
    set_meta(gr,
             rep("Satellite/Unknown", length(gr)),
             rep("Satellite/Unknown", length(gr)),
             6L, "TideHunter")
  }
  combined <- c(normalise_th(def), normalise_th(sho))
  if (length(combined) == 0) return(GRanges())

  # Merge overlapping TideHunter intervals from the two runs
  combined_r <- reduce(combined)
  combined_r$Name           <- "Satellite/Unknown"
  combined_r$classification <- "Satellite/Unknown"
  combined_r$source_tier    <- 6L
  combined_r$source_tool    <- "TideHunter"
  message("  ", length(combined_r), " TideHunter residuals (after merge)")
  combined_r
}

# ── 3. Batching ───────────────────────────────────────────────────────────────

read_fai <- function(path) {
  fai <- read.table(path, header = FALSE, sep = "\t",
                    col.names = c("name", "length", "offset",
                                  "bases_per_line", "bytes_per_line"))
  setNames(as.integer(fai$length), fai$name)
}

# Greedy bin-packing: sequences > batch_target get their own batch;
# smaller sequences are accumulated until the batch would overflow.
make_batches <- function(seqlengths, batch_target_bp = 200e6) {
  sl      <- sort(seqlengths, decreasing = TRUE)
  batches <- list()
  current <- character(0)
  cur_sz  <- 0L

  for (i in seq_along(sl)) {
    nm  <- names(sl)[i]
    len <- sl[i]
    if (len >= batch_target_bp) {
      batches <- c(batches, list(nm))
    } else if (cur_sz + len > batch_target_bp && length(current) > 0) {
      batches <- c(batches, list(current))
      current <- nm
      cur_sz  <- len
    } else {
      current <- c(current, nm)
      cur_sz  <- cur_sz + len
    }
  }
  if (length(current) > 0) batches <- c(batches, list(current))
  batches
}

# ── 4. Resolution helpers ─────────────────────────────────────────────────────

# Trim each feature in `lower` to its non-overlapping portion against `higher`.
# Metadata is propagated from the first overlapping source feature to each
# resulting fragment (within-tier overlaps are resolved later by
# resolve_within_tier via disjoin+LCA).
#
# Implementation: single vectorized disjoin() over the union of lower+higher
# with a revmap, then mask-filter to keep pieces that come from a lower feature
# but NOT from any higher feature. Avoids the per-feature lapply + keepSeqlevels
# + setdiff that dominated the previous runtime.
#
# min_len semantics (matching the legacy per-feature implementation): the
# filter only applies to pieces that came from a feature that actually got
# trimmed. Features that don't overlap any higher region pass through
# unchanged — even if they are shorter than min_len.
trim_to_nonoverlap <- function(lower, higher, min_len = 50L) {
  if (length(lower) == 0 || length(higher) == 0) return(lower)

  higher_r <- reduce(higher, ignore.strand = TRUE)

  # Which lower features actually overlap something in higher? Pieces derived
  # from features NOT in this set are "intact" — skip the min_len filter.
  lower_overlaps_higher <- overlapsAny(lower, higher_r, ignore.strand = TRUE)

  # Fast exit: no lower feature overlaps any higher region → return as-is.
  if (!any(lower_overlaps_higher)) return(lower)

  n_lower      <- length(lower)
  lower_plain  <- granges(lower)
  higher_plain <- granges(higher_r)
  strand(lower_plain)  <- "*"
  strand(higher_plain) <- "*"
  combined <- c(lower_plain, higher_plain)

  # Single disjoin with revmap: indices 1..n_lower are lower features,
  # indices n_lower+1..end are higher regions.
  dis    <- disjoin(combined, with.revmap = TRUE, ignore.strand = TRUE)
  revmap <- dis$revmap  # IntegerList

  # Vectorized list-ops: which pieces come from lower only, none from higher?
  has_higher <- any(revmap > n_lower)
  has_lower  <- any(revmap <= n_lower)
  base_mask  <- has_lower & !has_higher

  if (!any(base_mask)) return(lower[integer(0)])

  kept     <- dis[base_mask]
  kept_rev <- revmap[base_mask]
  # Kept pieces have only lower-source indices → min() gives the first lower parent.
  lower_idx <- as.integer(min(kept_rev))

  # Apply min_len only to pieces whose source lower feature was actually
  # trimmed (i.e. overlapped higher_r). Intact features pass through.
  source_was_trimmed <- lower_overlaps_higher[lower_idx]
  piece_widths       <- width(kept)
  keep_by_len        <- !source_was_trimmed | piece_widths >= min_len
  kept               <- kept[keep_by_len]
  lower_idx          <- lower_idx[keep_by_len]
  source_was_trimmed <- source_was_trimmed[keep_by_len]

  if (length(kept) == 0) return(lower[integer(0)])

  mcols(kept) <- mcols(lower)[lower_idx, , drop = FALSE]
  # Strand handling (matches the legacy per-feature implementation):
  #   - intact pieces (from features not overlapping higher) keep their
  #     original strand from the source lower feature;
  #   - trimmed pieces get "*" because the legacy setdiff-based path
  #     returned "*" for them under ignore.strand = TRUE.
  new_strand <- as.character(strand(lower))[lower_idx]
  new_strand[source_was_trimmed] <- "*"
  strand(kept) <- new_strand
  kept$revmap  <- NULL  # defensive
  kept
}

# For each feature in `children`, return the ID of the Level 1 parent
# (the overlapping feature in `parents` with the maximum intersection width).
# Fully vectorized — no R-level loop. Called once per pool type in
# finalise_output(); do NOT call this once per child.
get_parent_id <- function(children, parents) {
  n <- length(children)
  if (n == 0 || length(parents) == 0)
    return(rep(NA_character_, n))

  hits <- findOverlaps(children, parents, ignore.strand = TRUE)
  if (length(hits) == 0)
    return(rep(NA_character_, n))

  qi <- queryHits(hits)
  si <- subjectHits(hits)

  # Vectorized intersection width using start/end position arrays
  ch_start <- start(children)[qi]
  ch_end   <- end(children)[qi]
  pa_start <- start(parents)[si]
  pa_end   <- end(parents)[si]
  ov_w     <- pmin(ch_end, pa_end) - pmax(ch_start, pa_start) + 1L

  # Sort (qi ASC, ov_w DESC) → first entry per qi group is the best overlap.
  ord        <- order(qi, -ov_w)
  qi_ord     <- qi[ord]
  si_ord     <- si[ord]
  first_seen <- !duplicated(qi_ord)

  result <- rep(NA_character_, n)
  result[qi_ord[first_seen]] <- parents$ID[si_ord[first_seen]]
  result
}

# Lowest-Common-Ancestor of a set of slash-delimited classification strings.
# E.g. c("Class_I/LTR/Ty1/copia/Ale", "Class_I/LTR/Ty1/copia/TAR") → "Class_I/LTR/Ty1/copia"
lca_classification <- function(classifications) {
  cls <- unique(as.character(classifications))
  cls <- cls[!is.na(cls) & nchar(cls) > 0]
  if (length(cls) == 0) return("Unknown")
  if (length(cls) == 1) return(cls)
  parts <- strsplit(cls, "/", fixed = TRUE)
  min_len <- min(lengths(parts))
  if (min_len == 0) return("Unknown")
  common <- character(0)
  for (i in seq_len(min_len)) {
    vals <- vapply(parts, `[[`, character(1), i)
    if (length(unique(vals)) == 1) common <- c(common, vals[1]) else break
  }
  if (length(common) == 0) return("Unknown")
  paste(common, collapse = "/")
}

# ── TR-from-structural-TE resolution (Item 3) ────────────────────────────────
# A length-qualified TideCluster clustering TRC (Tier 3) that tandemly stacks
# multiple structural TEs (Tier 1: DANTE_LTR/TIR/LINE) of one family is itself a
# TE-derived satellite (e.g. the OZ408687.1 Ale array → TRC_13). Policy: the
# satellite WINS the region — it is tagged TE_origin=<LCA class of the covered
# structural TEs> and the structural TEs underneath are removed from the unified
# output (they remain in the DANTE_* GFF3s). Judged against STRUCTURAL
# annotation only (Tier 1), never RepeatMasker, and only for clustering TRCs,
# which already carry TideCluster's -m/-M length filters; TideHunter residuals
# (Tier 6) never reach this tier.
TE_ORIGIN_MIN_ELEMENTS  <- 2L    # ≥ this many distinct structural TEs under the TRC
TE_ORIGIN_MIN_LCA_DEPTH <- 3L    # shared lineage ≥ this many '/'-levels (e.g. Class_I/LTR/Ty1_copia)
# ≥ this fraction of the TRC's array bp covered by those structural TEs. DANTE
# annotates only *complete* elements, so a tandem array of mostly-degraded
# monomers is only partially covered: calibrated on OZ408687.1 TRC_13 (a true
# Ale-derived array at 0.35 coverage) vs the genome's two genuine large
# satellites with incidental TE insertions (TRC_1 11 Mb / TRC_3 1.9 Mb, both
# ~0.00) — a wide, clean gap. 0.25 sits in it with margin on both sides.
TE_ORIGIN_MIN_COVERAGE  <- 0.25

.cls_depth <- function(x)
  if (is.na(x) || !nzchar(x)) 0L else length(strsplit(x, "/", fixed = TRUE)[[1]])

# Identify TE-derived TRCs among Tier-3 clustering satellites `t3` against the
# structural TE tier `t1`. Returns a named character vector: TRC-Name ->
# TE_origin class (one entry per qualifying TRC). Grouping is by satellite Name
# (the TRC id); run separately on default vs short to avoid the TRC_<n> name
# collision between the two clustering runs. rDNA arrays share a Name but never
# qualify (they do not overlap structural TEs).
identify_te_derived_trcs <- function(t3, t1) {
  if (length(t3) == 0 || length(t1) == 0) return(character(0))
  t3p <- granges(t3); strand(t3p) <- "*"
  t1p <- granges(t1); strand(t1p) <- "*"
  nm  <- as.character(t3$Name)
  out <- character(0)
  for (name in unique(nm)) {
    arrays <- reduce(t3p[nm == name])
    h <- unique(subjectHits(findOverlaps(arrays, t1p)))
    if (length(h) < TE_ORIGIN_MIN_ELEMENTS) next
    lca <- lca_classification(t1$classification[h])
    if (.cls_depth(lca) < TE_ORIGIN_MIN_LCA_DEPTH) next
    trc_bp <- sum(as.numeric(width(arrays)))
    if (trc_bp <= 0) next
    cov_bp <- sum(as.numeric(width(intersect(arrays, t1p[h]))))
    if (cov_bp / trc_bp < TE_ORIGIN_MIN_COVERAGE) next
    out[[name]] <- lca
  }
  out
}

# Generic tier-1 structural-overlap fallback. DANTE_LTR (tandem containers +
# standalone) + DANTE_TIR + DANTE_LINE features must not overlap each other in the
# unified output. Greedy longest-first: keep the longest element complete, trim
# shorter overlapping ones to their non-overlapping remainder (>= min_len). This
# catches structural overlaps that are NOT clean tandem LTR-RT arrays (those were
# already collapsed upstream by resolve_ltr_tandems.py) — cross-lineage,
# cross-tool, or partial overlaps. Fast no-op when no tier-1 features overlap.
resolve_tier1_overlaps <- function(t1, min_len) {
  if (length(t1) <= 1) return(t1)
  h <- suppressWarnings(findOverlaps(t1, ignore.strand = TRUE,
                                     drop.self = TRUE, drop.redundant = TRUE))
  if (length(h) == 0) return(t1)
  t1s  <- t1[order(width(t1), decreasing = TRUE)]   # longest first stays complete
  kept <- t1s[1]
  for (i in 2:length(t1s)) {
    piece <- trim_to_nonoverlap(t1s[i], kept, min_len)
    if (length(piece) > 0) kept <- suppressWarnings(c(kept, piece))
  }
  kept
}

# After trimming lower-tier features against higher tiers, some overlaps can be
# re-introduced (two fragments of different features now cover the same region).
# disjoin + LCA resolves these, producing strictly non-overlapping output.
resolve_within_tier <- function(gr) {
  if (length(gr) == 0) return(gr)
  h <- suppressWarnings(findOverlaps(gr, gr, ignore.strand = TRUE))
  h <- h[queryHits(h) < subjectHits(h)]
  if (length(h) == 0) return(gr)
  dis <- suppressWarnings(disjoin(gr, with.revmap = TRUE, ignore.strand = TRUE))
  revmap <- as.list(dis$revmap)
  dis$classification <- vapply(revmap, function(idx) lca_classification(gr$classification[idx]), character(1))
  dis$Name           <- dis$classification
  dis$source_tier    <- vapply(revmap, function(idx) min(gr$source_tier[idx]),  integer(1))
  dis$source_tool    <- vapply(revmap, function(idx) gr$source_tool[idx[1]],    character(1))
  dis$revmap         <- NULL
  dis
}

# ── 5. Per-batch resolution ───────────────────────────────────────────────────

process_batch <- function(seqs, data, min_len) {
  batch_label <- sprintf("batch[%s%s]",
                         paste(head(seqs, 2), collapse=","),
                         if (length(seqs) > 2) sprintf(",...+%d", length(seqs)-2) else "")
  t_batch <- tic(sprintf("%s  (%d seqs)", batch_label, length(seqs)))

  sub     <- timed(sprintf("  %s  subset+normalise all tiers", batch_label),
                   lapply(data, subset_seqs, seqs = seqs))
  t1      <- sub$t1
  t1_ltr  <- if (length(t1) > 0) t1[t1$source_tool == "DANTE_LTR"] else GRanges()
  t1_members <- sub$t1_members           # tandem LTR-RT member copies -> Level 2
  t2      <- sub$t2
  t3_def  <- sub$t3_def
  t3_sho  <- sub$t3_sho
  t4      <- sub$t4
  t5_te   <- sub$t5_te
  t5_sc   <- sub$t5_sc
  t6      <- sub$t6

  log_msg(sprintf("  %s  input sizes: t1=%d t2=%d t3_def=%d t3_sho=%d t4=%d t5_te=%d t5_sc=%d t6=%d",
                  batch_label, length(t1), length(t2), length(t3_def), length(t3_sho),
                  length(t4), length(t5_te), length(t5_sc), length(t6)))

  level1 <- GRanges()
  level2 <- GRanges()

  # ── Pre-pass: TR-from-structural-TE (Item 3) ─────────────────────────────
  # Flag length-qualified Tier-3 clustering TRCs that tandemly stack ≥2 same-
  # family structural TEs. Those satellites WIN: tag them TE_origin and trim the
  # covered structural Tier-1/2 features out (they remain in the DANTE_* GFF3s).
  # Run default vs short separately (TRC_<n> names collide across the two runs).
  te_sat <- GRanges()
  trc_origin_def <- if (length(t3_def) > 0 && length(t1) > 0) identify_te_derived_trcs(t3_def, t1) else character(0)
  trc_origin_sho <- if (length(t3_sho) > 0 && length(t1) > 0) identify_te_derived_trcs(t3_sho, t1) else character(0)
  if (length(trc_origin_def) + length(trc_origin_sho) > 0) {
    # A TE-derived satellite that overlaps a tandem LTR-RT container is the
    # "full LTR-RT in tandems" sub-type (complete LTR-RTs sharing LTRs); one over
    # scattered complete TEs is the general (often degraded) TE-derived case.
    containers_t1 <- if (!is.null(t1$structure))
      t1[!is.na(t1$structure) & t1$structure == "LTR_RT_TR"] else GRanges()
    tag_te <- function(gr, omap) {
      if (length(gr) == 0) return(gr)
      gr$TE_origin <- unname(omap[as.character(gr$Name)])   # NA where not TE-derived
      sv <- rep(NA_character_, length(gr))
      if (length(containers_t1) > 0) {
        on_cont <- !is.na(gr$TE_origin) &
          overlapsAny(gr, containers_t1, ignore.strand = TRUE)
        sv[on_cont] <- "tandem_LTR_RT"
      }
      gr$TE_origin_structure <- sv
      gr
    }
    t3_def <- tag_te(t3_def, trc_origin_def)
    t3_sho <- tag_te(t3_sho, trc_origin_sho)
    def_te <- if (length(t3_def) > 0) !is.na(t3_def$TE_origin) else logical(0)
    sho_te <- if (length(t3_sho) > 0) !is.na(t3_sho$TE_origin) else logical(0)
    te_sat <- suppressWarnings(c(t3_def[def_te], t3_sho[sho_te]))
    t3_def <- t3_def[!def_te]
    t3_sho <- t3_sho[!sho_te]
    # Structural tiers — AND tandem member copies — lose their portion under the
    # TE-derived satellites (the satellite wins; TE_origin records the family).
    # Dropping members here is essential: otherwise a member whose container is
    # trimmed away orphans (its parent is gone), producing a malformed feature.
    te_sat_r <- reduce(granges(te_sat), ignore.strand = TRUE)
    t1 <- trim_to_nonoverlap(t1, te_sat_r, min_len)
    t2 <- trim_to_nonoverlap(t2, te_sat_r, min_len)
    if (length(t1_members) > 0)
      t1_members <- t1_members[!overlapsAny(t1_members, te_sat_r, ignore.strand = TRUE)]
    t1_ltr <- if (length(t1) > 0) t1[t1$source_tool == "DANTE_LTR"] else GRanges()
    log_msg(sprintf("  %s  pre-pass: %d TE-derived TRC(s) → te_sat=%d feats; t1→%d t2→%d t1_members→%d",
                    batch_label, length(trc_origin_def) + length(trc_origin_sho),
                    length(te_sat), length(t1), length(t2), length(t1_members)))
  }

  # Generic tier-1 structural-overlap fallback (always): no two tier-1 elements
  # overlap in the output. Clean tandem LTR-RT arrays were already collapsed
  # upstream (resolve_ltr_tandems.py); this handles any residual cross-lineage /
  # cross-tool / partial overlap. No-op when nothing overlaps (the common case).
  n_t1_before <- length(t1)
  t1 <- resolve_tier1_overlaps(t1, min_len)
  t1_ltr <- if (length(t1) > 0) t1[t1$source_tool == "DANTE_LTR"] else GRanges()
  if (length(t1) != n_t1_before)
    log_msg(sprintf("  %s  tier1 overlap resolve: %d → %d features",
                    batch_label, n_t1_before, length(t1)))

  # ── Step 1: Tier 1 — structure-based elements ────────────────────────────
  # TE-derived satellites win their region: seed Level 1 with them, then add the
  # (te_sat-trimmed) structural Tier-1 elements.
  if (length(te_sat) > 0) level1 <- suppressWarnings(c(level1, te_sat))
  if (length(t1) > 0)     level1 <- suppressWarnings(c(level1, t1))
  # Tandem LTR-RT member copies → Level 2, nested under their container (the
  # container is a Tier-1 DANTE_LTR feature now in level1; finalise_output
  # resolves each member's Parent by overlap via the DANTE_LTR parent lookup).
  if (length(t1_members) > 0) {
    t1_members$temp_parent_tool <- "DANTE_LTR"
    level2 <- suppressWarnings(c(level2, t1_members))
  }
  log_msg(sprintf("  %s  step1 done: level1=%d (te_sat=%d + tier1=%d), +%d LTR members → L2",
                  batch_label, length(level1), length(te_sat), length(t1),
                  length(t1_members)))

  # ── Step 2: Tier 2 — DANTE domains, trim against Tier 1 ──────────────────
  t2_trimmed <- timed(sprintf("  %s  step2: trim Tier2 (%d) vs Tier1 (%d)",
                               batch_label, length(t2), length(t1)),
                       resolve_within_tier(trim_to_nonoverlap(t2, t1, min_len)))
  if (length(t2_trimmed) > 0) level1 <- suppressWarnings(c(level1, t2_trimmed))
  log_msg(sprintf("  %s  step2 done: +%d → level1=%d",
                  batch_label, length(t2_trimmed), length(level1)))

  # ── Step 3: Tier 3 default — TideCluster satellite clusters ──────────────
  # Rare case: cluster entirely within an LTR element → Level 2 nested
  t3_def_l1 <- t3_def
  if (length(t3_def) > 0 && length(t1_ltr) > 0) {
    in_hits <- suppressWarnings(findOverlaps(t3_def, t1_ltr, type = "within", ignore.strand = TRUE))
    in_idx  <- as.integer(unique(queryHits(in_hits)))
    if (length(in_idx) > 0) {
      t3_in            <- t3_def[in_idx]
      t3_in$temp_parent_tool <- "DANTE_LTR"  # resolved to ID in finalise_output
      level2   <- suppressWarnings(c(level2, t3_in))
      t3_def_l1 <- t3_def[setdiff(seq_along(t3_def), in_idx)]
    }
  }
  # Standard tier priority for the non-TE-derived remainder: trim against
  # everything already placed (te_sat + Tier 1 + Tier 2). A satellite that merely
  # abuts or contains a single structural element yields that contested span to
  # it (the inserted/abutting TE wins), keeping the rest. TE-derived satellites
  # already won their whole region in the pre-pass; the within-LTR clusters above
  # stay nested. This is what removes the residual satellite-vs-structural-TE
  # overlaps (the non-TE-derived counterpart of the TR-from-TE fix).
  n_nested <- length(t3_def) - length(t3_def_l1)
  if (length(t3_def_l1) > 0) {
    higher_3  <- if (length(level1) > 0) suppressWarnings(reduce(level1, ignore.strand = TRUE)) else GRanges()
    t3_def_l1 <- timed(sprintf("  %s  step3: trim Tier3 default (%d) vs placed",
                               batch_label, length(t3_def_l1)),
                       trim_to_nonoverlap(t3_def_l1, higher_3, min_len))
  }
  if (length(t3_def_l1) > 0) level1 <- suppressWarnings(c(level1, t3_def_l1))
  log_msg(sprintf("  %s  step3 done: +%d (L1, trimmed vs higher), +%d (L2 nested in LTR) → level1=%d",
                  batch_label, length(t3_def_l1), n_nested, length(level1)))

  # ── Step 4: Tier 3 short monomer — trim against everything placed so far ──
  higher_4 <- if (length(level1) > 0) suppressWarnings(reduce(level1, ignore.strand = TRUE)) else GRanges()
  t3s_trimmed <- timed(sprintf("  %s  step4: trim Tier3 short (%d) vs placed (%d regions)",
                                batch_label, length(t3_sho), length(higher_4)),
                        trim_to_nonoverlap(t3_sho, higher_4, min_len))
  if (length(t3s_trimmed) > 0) level1 <- suppressWarnings(c(level1, t3s_trimmed))
  log_msg(sprintf("  %s  step4 done: +%d → level1=%d",
                  batch_label, length(t3s_trimmed), length(level1)))

  # ── Step 5: Tier 4 — RM on TideCluster library, trim against Tiers 1–3 ──
  # Tandem gate: a Tier-4 satellite may override a Tier-5 TE only where it has
  # independent tandem evidence (raw TideHunter, `th_raw`). Supported arrays keep
  # Tier-4 priority here; UNsupported arrays are held and placed only AFTER
  # Tier-5 TE (Step 6b), so a TE wins any contested span while a genuine
  # satellite that TideHunter merely missed still survives in open sequence. When
  # the gate is off / no evidence, th_raw is empty → all t4 stays supported (no-op).
  th_raw_b <- sub$th_raw
  if (length(t4) > 0 && length(th_raw_b) > 0) {
    has_td      <- overlapsAny(t4, th_raw_b, ignore.strand = TRUE)
    t4_tandem   <- t4[has_td]
    t4_notandem <- t4[!has_td]
  } else {
    t4_tandem   <- t4
    t4_notandem <- t4[integer(0)]
  }
  higher_1_3  <- timed(sprintf("  %s  step5: reduce(level1) [%d]", batch_label, length(level1)),
                        if (length(level1) > 0) suppressWarnings(reduce(level1, ignore.strand = TRUE)) else GRanges())
  t4_trimmed  <- timed(sprintf("  %s  step5: trim Tier4 tandem (%d/%d) vs Tiers1-3 (%d regions)",
                                batch_label, length(t4_tandem), length(t4), length(higher_1_3)),
                        resolve_within_tier(trim_to_nonoverlap(t4_tandem, higher_1_3, min_len)))
  if (length(t4_trimmed) > 0) level1 <- suppressWarnings(c(level1, t4_trimmed))
  log_msg(sprintf("  %s  step5 done: +%d (tandem-supported; %d held for step6b) → level1=%d",
                  batch_label, length(t4_trimmed), length(t4_notandem), length(level1)))

  # ── Step 6: Tier 5 TE hits — trim against Tiers 1–4 ─────────────────────
  higher_1_4 <- timed(sprintf("  %s  step6: reduce(level1) [%d]", batch_label, length(level1)),
                       if (length(level1) > 0) suppressWarnings(reduce(level1, ignore.strand = TRUE)) else GRanges())
  t5_trimmed <- timed(sprintf("  %s  step6: trim Tier5 TE (%d) vs Tiers1-4 (%d regions)",
                               batch_label, length(t5_te), length(higher_1_4)),
                       resolve_within_tier(trim_to_nonoverlap(t5_te, higher_1_4, min_len)))
  if (length(t5_trimmed) > 0) level1 <- suppressWarnings(c(level1, t5_trimmed))
  log_msg(sprintf("  %s  step6 done: +%d → level1=%d",
                  batch_label, length(t5_trimmed), length(level1)))

  # ── Step 6b: RM-on-TideCluster satellites WITHOUT tandem evidence ────────
  # Demoted below Tier-5 TE by the tandem gate: trim against everything placed so
  # far (Tiers 1–4-tandem AND Tier-5 TE), so a TE wins any contested span. What
  # remains is satellite in sequence no TE claims (genuine array TideHunter
  # missed). Empty when the gate is off. These keep their Tier-4 satellite tag.
  if (length(t4_notandem) > 0) {
    higher_6b   <- if (length(level1) > 0) suppressWarnings(reduce(level1, ignore.strand = TRUE)) else GRanges()
    t4nt_trimmed <- timed(sprintf("  %s  step6b: trim ungated Tier4 (%d) vs Tiers1-5 (%d regions)",
                                  batch_label, length(t4_notandem), length(higher_6b)),
                          resolve_within_tier(trim_to_nonoverlap(t4_notandem, higher_6b, min_len)))
    if (length(t4nt_trimmed) > 0) level1 <- suppressWarnings(c(level1, t4nt_trimmed))
    log_msg(sprintf("  %s  step6b done: +%d (ungated RM_TC kept where no TE) → level1=%d",
                    batch_label, length(t4nt_trimmed), length(level1)))
  }

  # ── Step 7: Tier 5 Simple/Low complexity ─────────────────────────────────
  # Overlapping an existing Level 1 feature → Level 2 nested
  # Not overlapping → Level 1
  if (length(t5_sc) > 0) {
    cur_l1 <- timed(sprintf("  %s  step7: reduce(level1) [%d]", batch_label, length(level1)),
                     if (length(level1) > 0) suppressWarnings(reduce(level1, ignore.strand = TRUE)) else GRanges())
    if (length(cur_l1) > 0) {
      sc_hits <- timed(sprintf("  %s  step7: findOverlaps Simple/Low (%d) vs L1 (%d)",
                                batch_label, length(t5_sc), length(cur_l1)),
                        suppressWarnings(findOverlaps(t5_sc, cur_l1, ignore.strand = TRUE)))
      sc_in   <- as.integer(unique(queryHits(sc_hits)))
      if (length(sc_in) > 0) {
        sc_l2            <- t5_sc[sc_in]
        sc_l2$temp_parent_tool <- "any"
        level2           <- suppressWarnings(c(level2, sc_l2))
        sc_l1            <- t5_sc[setdiff(seq_along(t5_sc), sc_in)]
      } else {
        sc_l1 <- t5_sc
      }
    } else {
      sc_l1 <- t5_sc
    }
    if (length(sc_l1) > 0) level1 <- suppressWarnings(c(level1, sc_l1))
    log_msg(sprintf("  %s  step7 done: +%d (L1), +%d (L2 nested) → level1=%d, level2=%d",
                    batch_label, length(sc_l1), length(t5_sc) - length(sc_l1),
                    length(level1), length(level2)))
  }

  # ── Step 8: Tier 6 — TideHunter residuals ────────────────────────────────
  # Within LTR body → Level 2 nested
  # Not covered by any higher tier → Level 1
  if (length(t6) > 0) {
    all_higher <- timed(sprintf("  %s  step8: reduce(level1) [%d]", batch_label, length(level1)),
                         if (length(level1) > 0) suppressWarnings(reduce(level1, ignore.strand = TRUE)) else GRanges())

    t6_rest <- t6
    if (length(t1_ltr) > 0) {
      ltr_hits <- suppressWarnings(findOverlaps(t6_rest, t1_ltr, type = "within", ignore.strand = TRUE))
      ltr_idx  <- as.integer(unique(queryHits(ltr_hits)))
      if (length(ltr_idx) > 0) {
        t6_l2                    <- t6_rest[ltr_idx]
        t6_l2$temp_parent_tool   <- "DANTE_LTR"
        level2                   <- suppressWarnings(c(level2, t6_l2))
        t6_rest                  <- t6_rest[setdiff(seq_along(t6_rest), ltr_idx)]
      }
    }
    # Keep only TideHunter not covered by anything in Level 1
    if (length(t6_rest) > 0 && length(all_higher) > 0) {
      cov_hits <- suppressWarnings(findOverlaps(t6_rest, all_higher, ignore.strand = TRUE))
      cov_idx  <- as.integer(unique(queryHits(cov_hits)))
      if (length(cov_idx) > 0)
        t6_rest <- t6_rest[setdiff(seq_along(t6_rest), cov_idx)]
    }
    if (length(t6_rest) > 0) level1 <- suppressWarnings(c(level1, t6_rest))
    log_msg(sprintf("  %s  step8 done: +%d → level1=%d, level2=%d",
                    batch_label, length(t6_rest), length(level1), length(level2)))
  }

  toc(t_batch, sprintf("final: level1=%d, level2=%d", length(level1), length(level2)))
  list(level1 = level1, level2 = level2)
}

# ── 6. Output assembly ────────────────────────────────────────────────────────

finalise_output <- function(level1, level2, seqlengths_vec, output_path) {
  t_fin <- tic(sprintf("finalise_output: L1=%d, L2=%d", length(level1), length(level2)))

  # ── Sort Level 1 ──────────────────────────────────────────────────────────
  all_seqs <- names(seqlengths_vec)
  seqlevels(level1, pruning.mode = "coarse") <- all_seqs
  seqlengths(level1) <- seqlengths_vec
  level1 <- timed(sprintf("sort Level 1 (%d features)", length(level1)),
                  sort(sortSeqlevels(level1)))

  # ── Assign globally unique IDs to Level 1 ────────────────────────────────
  level1$ID <- paste0("UA_L1_", formatC(seq_along(level1), width = 8, flag = "0"))

  # ── Assign feature types ──────────────────────────────────────────────────
  level1$type <- ifelse(
    grepl("^Satellite|^Simple_repeat|^Low_complexity|^rDNA|^Unknown",
          level1$classification),
    "repeat_region",
    "transposable_element"
  )
  level1$source <- level1$source_tool

  # ── Assign Parent IDs to Level 2 features ────────────────────────────────
  if (length(level2) > 0) {
    seqlevels(level2, pruning.mode = "coarse") <- all_seqs
    seqlengths(level2) <- seqlengths_vec
    level2 <- timed("sort Level 2", sort(sortSeqlevels(level2)))
    level2$ID <- paste0("UA_L2_", formatC(seq_along(level2), width = 8, flag = "0"))

    # Split Level 2 by pool type and call get_parent_id in batch for each
    # pool. Previous per-feature loop was O(n_l2 × n_l1) findOverlaps calls;
    # this is 1–2 batched findOverlaps calls total.
    ltr_l1 <- level1[level1$source_tool == "DANTE_LTR"]

    parent_tool <- if (!is.null(level2$temp_parent_tool))
      as.character(level2$temp_parent_tool) else rep(NA_character_, length(level2))
    parent_tool[is.na(parent_tool)] <- "any"

    is_ltr     <- parent_tool == "DANTE_LTR"
    parent_ids <- rep(NA_character_, length(level2))
    if (any(is_ltr)) {
      parent_ids[is_ltr] <- timed(
        sprintf("Level 2 → LTR parent lookup (%d children, %d parents)",
                sum(is_ltr), length(ltr_l1)),
        get_parent_id(level2[is_ltr], ltr_l1))
    }
    if (any(!is_ltr)) {
      parent_ids[!is_ltr] <- timed(
        sprintf("Level 2 → ANY parent lookup (%d children, %d parents)",
                sum(!is_ltr), length(level1)),
        get_parent_id(level2[!is_ltr], level1))
    }
    level2$Parent <- parent_ids
    level2$temp_parent_tool <- NULL
    level2$type   <- ifelse(
      grepl("^Satellite|^Simple_repeat|^Low_complexity|^rDNA|^Unknown",
            level2$classification),
      "repeat_region", "transposable_element"
    )
    level2$source <- level2$source_tool
    # Promote any Level-2 feature whose parent wasn't found to Level 1 (it can't
    # nest). type/source are already set above and we re-ID to a fresh UA_L1_, so
    # the promoted record is a valid Level-1 feature — not a parentless UA_L2_
    # with type/source='.' (the malformed-orphan bug). With members dropped under
    # TE-derived satellites this path is now rare, but stays correct if hit.
    orphans <- is.na(level2$Parent)
    if (any(orphans)) {
      message("  Promoting ", sum(orphans), " Level 2 orphans to Level 1")
      orphan_feats <- level2[orphans]
      orphan_feats$Parent <- NULL
      orphan_feats$ID <- paste0("UA_L1_",
                                formatC(length(level1) + seq_along(orphan_feats),
                                        width = 8, flag = "0"))
      level1 <- suppressWarnings(c(level1, orphan_feats))
      level2 <- level2[!orphans]
    }
  }

  # ── Clean up temp columns from Level 1 ───────────────────────────────────
  level1$temp_parent_tool <- NULL

  # ── Combine, final sort, export ───────────────────────────────────────────
  all_feats <- if (length(level2) > 0) suppressWarnings(c(level1, level2)) else level1
  seqlevels(all_feats, pruning.mode = "coarse") <- all_seqs
  seqlengths(all_feats) <- seqlengths_vec
  all_feats <- timed(sprintf("final sort (L1+L2 = %d features)", length(all_feats)),
                     sort(sortSeqlevels(all_feats)))

  timed(sprintf("export GFF3 to %s", output_path),
        export(all_feats, output_path, format = "gff3"))

  # Prepend provenance header lines so a consumer reading just the
  # GFF3 (not the whole output_dir) still knows which pipeline
  # version + run produced these calls. Falls back silently if
  # run_provenance.json is missing — the rule keeps its existing
  # behaviour without it.
  prov_path <- file.path(dirname(output_path), "run_provenance.json")
  if (file.exists(prov_path) &&
      requireNamespace("jsonlite", quietly = TRUE)) {
    prov <- tryCatch(jsonlite::fromJSON(prov_path, simplifyVector = TRUE),
                      error = function(e) NULL)
    if (!is.null(prov)) {
      .or_unknown <- function(x) if (is.null(x) || length(x) == 0) "unknown" else as.character(x)
      header <- c(
        sprintf("##pipeline-version %s", .or_unknown(prov$pipeline_version)),
        sprintf("##git-sha %s", .or_unknown(prov$git_sha)),
        sprintf("##run-started %s", .or_unknown(prov$run_started))
      )
      gff_lines <- readLines(output_path)
      # Insert after the existing ##gff-version line (always first).
      writeLines(c(gff_lines[1], header, gff_lines[-1]), output_path)
    }
  }

  toc(t_fin, sprintf("%d total features written", length(all_feats)))
  invisible(all_feats)
}

# ── 7. Sanity checks ──────────────────────────────────────────────────────────

sanity_check <- function(gr, seqlengths_vec, overlaps_tsv = NULL) {
  log_msg("── Sanity checks ──────────────────────────────────────────────")
  # gr is the in-memory GRanges returned by finalise_output — no GFF3 re-parse
  # required (saves ~30s on a large file and avoids CharacterList coercions).
  t0 <- tic("sanity check")

  parent_col <- if ("Parent" %in% colnames(mcols(gr))) gr$Parent else NULL
  if (is.null(parent_col)) {
    has_parent <- rep(FALSE, length(gr))
  } else if (is(parent_col, "CharacterList")) {
    has_parent <- lengths(parent_col) > 0
  } else {
    has_parent <- !is.na(parent_col)
  }
  l1 <- gr[!has_parent]
  l2 <- gr[has_parent]

  # ── Residual-overlap report (Item 3b) ────────────────────────────────────
  # Every Level-1 vs Level-1 overlap where NEITHER partner is
  # Simple_repeat / Low_complexity. Non-fatal. After the TR-from-TE fix this
  # should contain no satellite/rDNA-vs-structural-TE rows (the bug class);
  # same-tool partial-TE overlaps may legitimately remain and are surfaced for
  # review. Always (re)writes the file, header-only when there are no overlaps.
  if (!is.null(overlaps_tsv)) {
    hdr     <- "seqid\tstart\tend\tname_a\tclass_a\ttier_a\tname_b\tclass_b\ttier_b\toverlap_bp"
    allowed <- grepl("^Simple_repeat|^Low_complexity", as.character(l1$classification))
    cand    <- l1[!allowed]
    n_ov    <- 0L
    if (length(cand) >= 2) {
      ch <- findOverlaps(cand, cand, ignore.strand = TRUE)
      ch <- ch[queryHits(ch) < subjectHits(ch)]
      if (length(ch) > 0) {
        qi <- queryHits(ch); si <- subjectHits(ch)
        df <- data.frame(
          seqid      = as.character(seqnames(cand))[qi],
          start      = pmax(start(cand)[qi], start(cand)[si]),
          end        = pmin(end(cand)[qi],   end(cand)[si]),
          name_a     = as.character(cand$Name)[qi],
          class_a    = as.character(cand$classification)[qi],
          tier_a     = as.character(cand$source_tier)[qi],
          name_b     = as.character(cand$Name)[si],
          class_b    = as.character(cand$classification)[si],
          tier_b     = as.character(cand$source_tier)[si],
          overlap_bp = pmin(end(cand)[qi], end(cand)[si]) -
                       pmax(start(cand)[qi], start(cand)[si]) + 1L,
          stringsAsFactors = FALSE)
        write.table(df, overlaps_tsv, sep = "\t", quote = FALSE,
                    row.names = FALSE, col.names = TRUE)
        n_ov <- nrow(df)
      }
    }
    if (n_ov == 0L) writeLines(hdr, overlaps_tsv)
    message(sprintf("Overlap report: %d non-simple/low Level-1 overlap pair(s) → %s",
                    n_ov, overlaps_tsv))
    if (n_ov > 0L)
      warning(n_ov, " residual non-simple/low overlap(s) in unified annotation; see ",
              basename(overlaps_tsv))
  }

  # 1. TE features (tier 1–5 TE) should be non-overlapping; Simple/Low can overlap.
  is_simple <- grepl("^Simple_repeat|^Low_complexity|^Satellite|^Unknown",
                     as.character(l1$classification))
  l1_te     <- l1[!is_simple]
  te_hits   <- findOverlaps(l1_te, l1_te, ignore.strand = TRUE)
  te_hits   <- te_hits[queryHits(te_hits) < subjectHits(te_hits)]
  if (length(te_hits) > 0) {
    tier_q <- as.character(l1_te$source_tier[queryHits(te_hits)])
    tier_s <- as.character(l1_te$source_tier[subjectHits(te_hits)])
    pair_tbl <- table(paste0("T", tier_q, "-T", tier_s))
    message("Note: ", length(te_hits), " overlapping TE pairs remain (mostly same-tool partial elements):")
    print(pair_tbl)
  } else {
    message("OK: Level 1 TE features are non-overlapping (", length(l1_te), " features)")
  }
  message("Level 1 total: ", length(l1), " features; ",
          sum(is_simple), " Simple/Low/Satellite (may overlap each other)")

  # 2. Total Level 1 coverage ≤ genome size
  genome_bp <- sum(as.numeric(seqlengths_vec))
  l1_bp     <- sum(as.numeric(width(reduce(l1, ignore.strand = TRUE))))
  pct       <- round(l1_bp / genome_bp * 100, 2)
  message("OK: Level 1 coverage = ", format(l1_bp, big.mark=","), " bp (",
          pct, "% of genome)")

  # 3. Level 2 Parent references valid Level 1 IDs
  if (length(l2) > 0) {
    l2_parents <- if (is(parent_col[has_parent], "CharacterList"))
      unlist(parent_col[has_parent]) else parent_col[has_parent]
    bad_parents <- sum(!l2_parents %in% l1$ID)
    if (bad_parents > 0) {
      warning(bad_parents, " Level 2 features reference non-existent Parent IDs")
    } else {
      message("OK: All ", length(l2), " Level 2 features have valid Parent IDs")
    }
  }

  # 4. Summary — use the source column (GFF3 col 2) which is always a plain vector
  message("\nLevel 1 breakdown by source (tool):")
  tool_col <- if (!is.null(l1$source)) l1$source else l1$source_tool
  print(sort(table(as.character(tool_col)), decreasing = TRUE))
  message("\nLevel 1 breakdown by source_tier attribute:")
  tier_col <- suppressWarnings(l1$source_tier)
  if (!is.null(tier_col)) print(sort(table(as.character(tier_col))))
  toc(t0)
  message("──────────────────────────────────────────────────────────────\n")
}

# ── 8. Main ───────────────────────────────────────────────────────────────────

log_msg("=== make_unified_annotation.R ===")

# Load all tiers
log_msg("── Loading inputs ─────────────────────────────────────────────")
t_load <- tic("loading all tiers")
ltr_data  <- timed("load DANTE_LTR",   load_tier1_ltr(opt$ltr, opt$ltr_tandems))
tir_data  <- timed("load DANTE_TIR",   load_tier1_tir(opt$tir))
line_data <- timed("load DANTE_LINE",  load_tier1_line(opt$line))

t1 <- suppressWarnings(c(ltr_data$top, tir_data, line_data$top))
log_msg("Tier 1 total: ", length(t1), " top-level features")

t2 <- timed("load DANTE filtered (Tier 2)", load_tier2_dante(opt$dante))
log_msg("Tier 2 total: ", length(t2), " features")

# Authoritative rDNA-TRC calls from TideCluster's <prefix>_rdna.tsv (default and
# short runs). Applied to tier-3 and tier-4 alike; empty maps fall back to the
# per-feature rDNA_type attribute inside normalise_tc_satellite().
rdna_default <- load_rdna_map(opt$tc_rdna_default)
rdna_short   <- load_rdna_map(opt$tc_rdna_short)
log_msg("Authoritative rDNA TRCs: ", length(rdna_default), " default + ",
        length(rdna_short), " short")

tc_data <- timed("load TideCluster (Tier 3)",
                 load_tier3_tidecluster(opt$tc_default, opt$tc_short,
                                        rdna_default, rdna_short))
log_msg("Tier 3 total: ", length(tc_data$default), " default + ",
        length(tc_data$short), " short monomer features")

t4 <- timed("load RM-on-TC (Tier 4)", load_tier4_tc_rm(opt$tc_rm, rdna_default))
log_msg("Tier 4 total: ", length(t4), " features")

rm_data <- timed("load RepeatMasker+DANTE (Tier 5)", load_tier5_rm(opt$rm))
log_msg("Tier 5 total: ", length(rm_data$te), " TE + ",
        length(rm_data$simple), " Simple/Low features")

t6 <- timed("load TideHunter (Tier 6)",
            load_tier6_tidehunter(opt$th_default, opt$th_short))
log_msg("Tier 6 total: ", length(t6), " features")

# Independent tandem evidence for the RM-on-TideCluster (Tier 4) tandem gate.
# The RM-on-TC dimer library can tile a *non*-tandem, AT-rich TE (e.g. a Tekay
# LTR-RT) as a long apparent "array"; when such a satellite overrides the TE it
# is spurious over-masking. Raw TideHunter (a de-novo tandem detector, run
# independently of the dimer library) marks where real tandem structure exists.
# We use it to gate whether a Tier-4 satellite may override a Tier-5 TE. Reduced
# to plain intervals; empty when the gate is off or the files are absent.
rm_tc_tandem_gate <- toupper(opt$rm_tc_tandem_gate) %in% c("TRUE","T","YES","1")
th_raw <- GRanges()
if (rm_tc_tandem_gate) {
  raw_paths <- Filter(function(p) !is.null(p) && nzchar(p) && file.exists(p),
                      c(opt$th_raw_default, opt$th_raw_short))
  if (length(raw_paths) > 0) {
    th_raw <- timed("load raw TideHunter (tandem evidence)",
                    suppressWarnings(reduce(granges(do.call(c, lapply(raw_paths,
                      function(p) { g <- safe_import(p); mcols(g) <- NULL; g }))),
                      ignore.strand = TRUE)))
    log_msg("RM_TC tandem gate: ON, ", length(th_raw), " raw-TideHunter evidence regions")
  } else {
    log_msg("RM_TC tandem gate: requested but no raw-TideHunter files given — gate inert")
  }
} else {
  log_msg("RM_TC tandem gate: OFF")
}
toc(t_load)

# Read genome FAI
log_msg("── Reading FAI ────────────────────────────────────────────────")
seqlengths_vec <- read_fai(opt$fai)
genome_bp      <- sum(as.numeric(seqlengths_vec))
log_msg("Genome: ", length(seqlengths_vec), " sequences, ",
        format(genome_bp, big.mark=","), " bp total")

# Decide whether to chunk
data <- list(
  t1         = t1,
  t1_members = ltr_data$members,   # tandem LTR-RT member copies -> Level 2
  t2         = t2,
  t3_def     = tc_data$default,
  t3_sho     = tc_data$short,
  t4         = t4,
  t5_te      = rm_data$te,
  t5_sc      = rm_data$simple,
  t6         = t6,
  th_raw     = th_raw               # tandem evidence for the Tier-4 gate (may be empty)
)

# Decide batching. Even when the genome is below the "chunked processing"
# threshold we still want to split the work across threads — otherwise
# --threads > 1 is wasted on small-to-medium genomes. Pick a batch target
# that creates roughly one batch per thread, capped by opt$batch_size.
if (opt$threads <= 1 || length(seqlengths_vec) == 1) {
  # Nothing to parallelise — one batch, one worker.
  batches <- list(names(seqlengths_vec))
  log_msg("── Single batch (threads=1 or single-sequence genome) ─────────")
} else {
  effective_target <- min(
    opt$batch_size,
    max(1e6, ceiling(genome_bp / opt$threads))
  )
  batches <- make_batches(seqlengths_vec, effective_target)
  log_msg(sprintf(
    "── Batching: %d batch(es), target ~%.1f Mb/batch, %d threads ─────────",
    length(batches), effective_target / 1e6, opt$threads))
  # Optional: report the largest batch to help diagnose skew.
  batch_bp <- vapply(batches, function(b) sum(as.numeric(seqlengths_vec[b])), numeric(1))
  log_msg(sprintf(
    "   batch sizes: min=%.1f Mb  median=%.1f Mb  max=%.1f Mb  (features cluster on largest)",
    min(batch_bp)/1e6, stats::median(batch_bp)/1e6, max(batch_bp)/1e6))
}

# Run resolution — always via mclapply so the code path is uniform.
log_msg(sprintf("── Running tier resolution (%d threads, %d batch(es)) ─────",
                opt$threads, length(batches)))
t_resolve <- tic("tier resolution (all batches)")
results <- mclapply(batches,
                    function(seqs) process_batch(seqs, data, opt$min_feature_length),
                    mc.cores = min(opt$threads, length(batches)))
toc(t_resolve)

# Check for worker errors (mclapply returns error objects on failure)
err_idx <- which(vapply(results, function(r) inherits(r, "try-error"), logical(1)))
if (length(err_idx) > 0) {
  for (i in err_idx) message("Batch ", i, " failed: ", as.character(results[[i]]))
  stop("One or more batches failed — see messages above")
}

# Combine batch results
log_msg("── Combining batch results ─────────────────────────────────────")
level1_all <- timed("combine level1 across batches",
                    suppressWarnings(do.call(c, lapply(results, `[[`, "level1"))))
level2_all <- timed("combine level2 across batches",
                    suppressWarnings(do.call(c, lapply(results, `[[`, "level2"))))
log_msg("Combined: ", length(level1_all), " Level 1, ",
        length(level2_all), " Level 2 features")

# Finalise and export
log_msg("── Finalising output ──────────────────────────────────────────")
all_feats <- finalise_output(level1_all, level2_all, seqlengths_vec, opt$output)

# Sanity checks on the in-memory GRanges (avoids re-parsing the GFF3)
sanity_check(all_feats, seqlengths_vec,
             overlaps_tsv = sub("\\.gff3$", ".overlaps.tsv", opt$output))

log_msg("Done. Output: ", opt$output,
        sprintf("  (total wall time: %.1fs)", proc.time()[3] - .script_start))
