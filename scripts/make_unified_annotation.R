#!/usr/bin/env Rscript
# make_unified_annotation.R
# Produces a unified, tier-prioritised repeat annotation GFF3 from all pipeline layers.
# See docs/archive/annotation_rules.md and docs/archive/unified_annotation_implementation_plan.md for design rationale.

suppressPackageStartupMessages({
  library(optparse)
  library(rtracklayer)
  library(GenomicRanges)
  library(parallel)
})

.script_dir <- dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.script_dir, "classification.R"))
# Memory instrumentation + cgroup-aware budget detection (mem_rss_mb, mem_hwm_mb,
# mem_str, mem_budget_mb, ...). Shared with the density scripts, which size their
# own worker pools the same way.
source(file.path(.script_dir, "mem_utils.R"))

# ── Progress / timing helpers ─────────────────────────────────────────────────
# Simple wall-clock timers. All messages go to stderr so they interleave
# correctly with mclapply worker output.
.script_start <- proc.time()[3]
elapsed <- function() sprintf("[+%6.1fs]", proc.time()[3] - .script_start)

log_msg <- function(...) message(elapsed(), " ", ...)

# ── memory instrumentation ───────────────────────────────────────────────────
# mem_rss_mb / mem_hwm_mb / mem_str and the cgroup-aware budget detection
# (mem_budget_mb) live in mem_utils.R, sourced above and shared with the density
# scripts. Only the log formatting is local, so the figures carry this script's
# elapsed-time prefix.
# label is padded so the columns line up when grepping '[mem]' out of the log.
log_mem <- function(label) {
  s <- mem_str()
  if (!is.na(s)) log_msg(sprintf("[mem] %-44s %s", label, s))
}

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
  make_option("--tc_trc_table_default", type="character", default=NULL,
              help="TideCluster default report trc_table.tsv (authoritative per-TRC period monomer_tarean->monomer_kite; monomer column + domain-rhythm gate; optional)"),
  make_option("--tc_trc_table_short",   type="character", default=NULL,
              help="TideCluster short-monomer report trc_table.tsv (authoritative per-TRC period; optional)"),
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
  # --threads sets BATCH COMPOSITION (batch target = genome_bp / threads); the two
  # options below cap only the number of CONCURRENT WORKERS, so tuning them can
  # never change the result — see the memory-gating block near the mclapply call.
  make_option("--max_workers",      type="integer",   default=8L,
              help=paste("Hard ceiling on concurrent mclapply workers, independent of",
                         "--threads. Each forked worker's peak RSS converges on the",
                         "PARENT heap (measured: 48.3 GB workers vs a 48.4 GB parent,",
                         "identical across batches spanning 143 Mb to 2.15 Gb), so",
                         "workers are expensive while tier resolution is only ~13% of",
                         "this rule's wall time. 0 = no ceiling. [8]")),
  make_option("--mem_budget_gb",    type="double",    default=0,
              help=paste("Memory budget (GB) for sizing the worker pool. 0 = detect:",
                         "the tightest of the cgroup limit (walking up the hierarchy,",
                         "so a PBS/Slurm job scope limit is honoured) and /proc/meminfo",
                         "MemAvailable. [0]")),
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
  # `%in%` on the seqnames Rle matches on run values (cheap); the old
  # `as.character(seqnames(gr))` expanded the whole tier's seqnames to a plain
  # character vector on every batch call -> O(n_batches x tier_features).
  gr <- gr[seqnames(gr) %in% seqs]
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
  # Preserve the original DANTE domain identity (GAG/PROT/INT/RT/RH/TPase/...)
  # before set_meta overwrites Name with the classification. Kept only on the
  # top-level t2 object (dropped by subset_seqs' .META_COLS standardisation, which
  # the per-batch resolver uses) — write_te_derived_trc_table reads it globally.
  dom <- if (!is.null(raw$Name)) as.character(raw$Name) else rep(NA_character_, length(raw))
  raw <- set_meta(raw, cls, cls, 2L, "DANTE")
  raw$domain <- dom
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

# Clip each feature in `lower` to its portion not covered by `higher`, AND
# decompose `lower`'s own internal overlaps into disjoint pieces at the same
# time. Metadata (and strand) for each resulting piece come from its lowest-index
# source feature in `lower` (an order-stable, deterministic choice); pieces from
# a feature that actually overlapped `higher` get strand "*".
#
# DETERMINISM — why the two former `return(lower)` short-circuits are gone: they
# returned `lower` with its INTERNAL overlaps INTACT (when `higher` was empty, or
# when no lower feature overlapped higher), whereas the main path disjoined those
# overlaps. Which path ran depended on the batch's feature mix — `higher` is
# `reduce(level1)`, batch-composition-dependent, and batches are split by thread
# count. So the SAME `lower` could come out disjoined-with-min-index-metadata in
# one run but left internally-overlapping (and then disjoined by resolve_within_tier
# with LCA + strand "*") in another — a strand/classification flip on same-tier
# features that overlap only each other (surfaced by the unified_multibatch
# fixture; see docs/archive/tier1_resolution_determinism_audit.md). We now ALWAYS
# decompose `lower`, so `resolve_within_tier(trim_to_nonoverlap(...))` is
# batch-invariant and `lower` is always returned internally non-overlapping
# (matching what a single-batch threads=1 run already produced). When nothing
# overlaps `higher` we disjoin `lower` ALONE (never against a large `higher`), so
# the fast path stays cheap.
#
# min_len semantics (unchanged): the filter applies only to pieces whose source
# feature actually overlapped higher; features not overlapping higher pass
# through even if shorter than min_len.
trim_to_nonoverlap <- function(lower, higher, min_len = 50L) {
  if (length(lower) == 0) return(lower)

  n_lower     <- length(lower)
  lower_plain <- granges(lower)
  strand(lower_plain) <- "*"

  higher_r <- if (length(higher) > 0) reduce(higher, ignore.strand = TRUE) else GRanges()
  # Which lower features actually overlap something in higher? Pieces derived
  # from features NOT in this set are "intact" — skip the min_len filter and keep
  # their original strand.
  lower_overlaps_higher <- if (length(higher_r) > 0)
    overlapsAny(lower, higher_r, ignore.strand = TRUE) else logical(n_lower)

  if (any(lower_overlaps_higher)) {
    higher_plain <- granges(higher_r)
    strand(higher_plain) <- "*"
    combined <- c(lower_plain, higher_plain)   # trim against higher + decompose lower
  } else {
    # Nothing to trim: still decompose lower's internal overlaps so the result is
    # identical whether or not this batch placed any higher-tier feature next to
    # `lower`. Disjoining lower alone (not lower+higher) gives the same pieces —
    # no piece is higher-covered — at a fraction of the cost.
    combined <- lower_plain
  }

  # Single disjoin with revmap: indices 1..n_lower are lower features,
  # indices n_lower+1..end (if any) are higher regions.
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
TE_ORIGIN_LCA_MIN_SHARE <- 0.10   # a covered family must be >= this share of the
                                  # covered bp to count toward the LCA (ignore strays)

.cls_depth <- function(x)
  if (is.na(x) || !nzchar(x)) 0L else length(strsplit(x, "/", fixed = TRUE)[[1]])

# Identify TE-derived TRCs among Tier-3 clustering satellites `t3` against the
# structural TE tier `t1`. Returns a named character vector: TRC-Name ->
# TE_origin class (one entry per qualifying TRC). Grouping is by satellite Name
# (the TRC id); run separately on default vs short to avoid the TRC_<n> name
# collision between the two clustering runs. rDNA arrays share a Name but never
# qualify (they do not overlap structural TEs).
#
# A candidate must pass BOTH:
#   (a) structural coverage — >= TE_ORIGIN_MIN_ELEMENTS same-family (LCA depth
#       >= TE_ORIGIN_MIN_LCA_DEPTH) structural TEs covering >= TE_ORIGIN_MIN_COVERAGE
#       of the array bp; and
#   (b) domain rhythm — the TE's DANTE domains recur through the tandem (occupancy
#       >= TE_RHYTHM_MIN_OCC of P-windows, in >= TE_RHYTHM_MIN_FRAC of the arrays),
#       tiling at the authoritative period `period_map[[name]]` (trc_table
#       monomer_tarean->monomer_kite). This is what separates a genuine TE-derived
#       tandem (TE in ~every monomer) from a plain satellite merely INTERRUPTED by
#       a few TE insertions (TE clumped in a few blocks); the latter is NOT tagged
#       and falls through to normal tier resolution, which splits the satellite
#       around the interrupting TEs. Coverage alone cannot tell these apart —
#       a degraded genuine tandem also has low coverage — but domain rhythm is
#       degradation-robust. No period available (rare; TRC absent from trc_table)
#       -> the rhythm test is skipped and the coverage decision stands (we do not
#       demote what we cannot measure).
identify_te_derived_trcs <- function(t3, t1, t2 = GRanges(), period_map = integer(0)) {
  if (length(t3) == 0 || length(t1) == 0) return(character(0))
  t3p <- granges(t3); strand(t3p) <- "*"
  t1p <- granges(t1); strand(t1p) <- "*"
  t2p <- if (length(t2) > 0) { g <- granges(t2); strand(g) <- "*"; g } else GRanges()
  # Restrict domains to those overlapping any clustering array (a tiny fraction of
  # the genome) up front, so the per-TRC te_domain_rhythm countOverlaps stays cheap
  # when this runs globally over a multi-Gbp genome's full t2 domain set.
  if (length(t2p) > 0)
    t2p <- t2p[overlapsAny(t2p, t3p, ignore.strand = TRUE)]
  nm  <- as.character(t3$Name)
  out <- character(0)
  t1cls <- as.character(t1$classification)
  for (name in unique(nm)) {
    arrays <- reduce(t3p[nm == name])
    h <- unique(subjectHits(findOverlaps(arrays, t1p)))
    if (length(h) < TE_ORIGIN_MIN_ELEMENTS) next
    trc_bp <- sum(as.numeric(width(arrays)))
    if (trc_bp <= 0) next
    cov_bp <- sum(as.numeric(width(intersect(arrays, t1p[h]))))
    if (cov_bp / trc_bp < TE_ORIGIN_MIN_COVERAGE) next
    # Robust LCA over the DOMINANT covered families. Strict LCA over *every*
    # overlapping structural TE is fragile: a single stray inserted element of an
    # unrelated family (e.g. 1 Ty1_copia/SIRE among 73 Ty3_gypsy/CRM in a CRM
    # array) collapses the LCA to Class_I/LTR (depth 2) and would drop a genuine
    # TE-derived TRC — and, now that the decision is global over all of a TRC's
    # arrays, one such stray anywhere in the family kills it. So compute the LCA
    # only over classes that each cover >= TE_ORIGIN_LCA_MIN_SHARE of the covered
    # bp; rare strays are ignored. Genuine multi-family stacks (all classes above
    # the share) still take their true shared lineage.
    cls_h <- t1cls[h]
    ucls  <- unique(cls_h)
    cov_c <- vapply(ucls, function(cc)
      sum(as.numeric(width(intersect(arrays, t1p[h[cls_h == cc]])))), numeric(1))
    keep_cls <- ucls[cov_c >= TE_ORIGIN_LCA_MIN_SHARE * sum(cov_c)]
    if (length(keep_cls) == 0) keep_cls <- ucls
    lca <- lca_classification(keep_cls)
    if (.cls_depth(lca) < TE_ORIGIN_MIN_LCA_DEPTH) next
    # Domain-rhythm gate (skipped only when no period is available for this TRC).
    P <- if (length(period_map) > 0) period_map[[name]] else NULL
    if (!is.null(P) && !is.na(P) && P > 0) {
      rhy <- te_domain_rhythm(arrays, t2p, P)
      if (!is.na(rhy$occ) &&
          (rhy$occ < TE_RHYTHM_MIN_OCC || rhy$frac < TE_RHYTHM_MIN_FRAC)) next
    }
    out[[name]] <- lca
  }
  out
}

# ── TE-derived tandem-repeat summary table ───────────────────────────────────
# One row per (TRC id, run) for every TideCluster clustering satellite the
# TR-from-structural-TE pre-pass tagged with TE_origin (identify_te_derived_trcs
# above). Aggregates genome-wide over that TRC's tagged arrays. Written next to
# the unified GFF3 as Repeat_Annotation_Unified.te_derived_trc.csv and rendered
# in the Tandem-repeats section of the HTML report. All data is already in
# memory: the final combined level1 (the TE_origin arrays), the untrimmed global
# structural tier `t1` (covered bp + complete-element count), the global DANTE
# domain tier `t2` (domain set via t2$domain), and the two kite monomer maps.

# Canonical protein-domain order (REXdb/DANTE vocabulary). Domains not listed
# sort after these, alphabetically — so the emitted field is stable run-to-run.
.DOMAIN_ORDER <- c("GAG","PROT","AP","INT","RT","RH","aRH",
                   "CHD","CHDCR","CHDII","TPase","ENDO","EN","HEL1","HEL2")
.order_domains <- function(d) {
  d <- unique(d[!is.na(d) & nzchar(d)])
  if (length(d) == 0) return(character(0))
  rnk <- match(d, .DOMAIN_ORDER)
  d[order(ifelse(is.na(rnk), length(.DOMAIN_ORDER) + 1L, rnk), d)]
}

# Authoritative per-TRC tandem monomer period (bp) from TideCluster's report
# table `trc_table.tsv`. Per TRC: `monomer_tarean` (TAREAN family consensus) when
# present, else `monomer_kite` (most-frequent KITE *founder* period), else
# `prevalent_founder`. Returns a named integer vector TRC_ID -> period; empty when
# the file is absent/unusable (optional — `--no_rdna`/older TideCluster/purged).
#
# Why NOT the kite `monomer_size` CSV: that column is the top k-mer *peak*, which
# can lock onto a short SSR sub-period (measured: 79 bp reported for a genuine
# 13134 bp TIR-derived monomer). Tiling the domain-rhythm occupancy test at 79 bp
# wrongly reads a real TIR tandem as sparse; the TAREAN/founder period is correct.
# trc_table.tsv also survives `cleanup_intermediates: maximal` (the kite tree does not).
read_trc_periods <- function(trc_table_tsv) {
  empty <- setNames(integer(0), character(0))
  if (is.null(trc_table_tsv) || !nzchar(trc_table_tsv) || !file.exists(trc_table_tsv))
    return(empty)
  tab <- tryCatch(read.table(trc_table_tsv, header = TRUE, sep = "\t", check.names = FALSE,
                             stringsAsFactors = FALSE, quote = "", comment.char = ""),
                  error = function(e) NULL)
  if (is.null(tab) || nrow(tab) == 0 || !("TRC_ID" %in% names(tab))) return(empty)
  pick <- function(row_i) {
    for (col in c("monomer_tarean", "monomer_kite", "prevalent_founder")) {
      if (col %in% names(tab)) {
        v <- suppressWarnings(as.integer(as.character(tab[[col]][row_i])))
        if (!is.na(v) && v > 0) return(v)
      }
    }
    NA_integer_
  }
  vals <- vapply(seq_len(nrow(tab)), pick, integer(1))
  ids  <- as.character(tab$TRC_ID)
  keep <- !is.na(vals) & nzchar(ids)
  setNames(vals[keep], ids[keep])
}

# Domain-rhythm of a TRC's arrays at tandem period P: is the TE signal spread
# through the tandem (TE-derived) or clumped in a few blocks (satellite merely
# interrupted by TE insertions)? Tile each array into P-bp windows; a window is
# "occupied" if a DANTE protein domain overlaps it. Uses domains (not full
# elements) because they survive element decay — a degraded TE-tandem still
# carries its RT/INT/GAG/TPase domains in ~every monomer. Returns:
#   occ  = pooled occupied/total windows across the TRC's arrays
#   frac = fraction of the TRC's arrays whose own occupancy >= 0.5
# Calibrated on 3 genomes (100 TE_origin TRCs): derived 0.67-1.0, all frac=1.0.
te_domain_rhythm <- function(arr, doms, P) {
  if (is.na(P) || P <= 0 || length(arr) == 0)
    return(list(occ = NA_real_, frac = NA_real_))
  strand(arr) <- "*"
  totN <- 0L; totW <- 0L; per <- numeric(length(arr))
  for (i in seq_along(arr)) {
    wins <- suppressWarnings(unlist(tile(arr[i], width = P)))  # P-bp windows
    strand(wins) <- "*"
    occw <- suppressWarnings(countOverlaps(wins, doms, ignore.strand = TRUE)) > 0
    n <- length(wins)
    totN <- totN + n; totW <- totW + sum(occw)
    per[i] <- if (n > 0) mean(occw) else 0
  }
  list(occ = if (totN > 0) totW / totN else 0,
       frac = if (length(per) > 0) mean(per >= 0.5) else 0)
}
TE_RHYTHM_MIN_OCC  <- 0.5   # >= this pooled per-monomer domain occupancy
TE_RHYTHM_MIN_FRAC <- 0.5   # AND >= this fraction of the TRC's arrays in-rhythm

write_te_derived_trc_table <- function(level1, t1, t1_members, t2, period_default, period_short, out_csv) {
  hdr <- c("trc_id","run","n_arrays","total_array_bp","monomer_length_bp",
           "te_classification","te_origin_structure","protein_domains",
           "n_complete_elements","n_expected_monomers","complete_bp_fraction",
           "domain_occupancy","frac_arrays_in_rhythm")
  write_empty <- function() { writeLines(paste(hdr, collapse = ","), out_csv)
                              log_msg("No TE-derived TRCs; wrote header-only ", basename(out_csv)) }

  if (length(level1) == 0 || is.null(level1$TE_origin)) return(write_empty())
  te <- level1[!is.na(level1$TE_origin)]
  if (length(te) == 0) return(write_empty())

  run_of <- ifelse(as.character(te$source_tool) == "TideCluster_short", "short", "default")
  key    <- paste(run_of, as.character(te$Name), sep = "\t")

  # Complete structural element copies as plain intervals, for the covered-bp
  # fraction and the element count: standalone DANTE_LTR/TIR/LINE elements PLUS
  # the LTR-RT member copies of any LTR_RT_TR container — but NOT the containers
  # themselves. A container wraps a whole tandem array; counting it as one
  # element would tally an entire array of copies as a single element (and its
  # member copies live in t1_members, Level 2, not in t1).
  elem0 <- if (length(t1) > 0) {
    st <- if (is.null(t1$structure)) rep(NA_character_, length(t1)) else as.character(t1$structure)
    granges(t1[is.na(st) | st != "LTR_RT_TR"])
  } else GRanges()
  elem <- suppressWarnings(c(elem0, if (length(t1_members) > 0) granges(t1_members) else GRanges()))
  if (length(elem) > 0) strand(elem) <- "*"
  t2p <- if (length(t2) > 0) { g <- granges(t2); strand(g) <- "*"; g } else GRanges()
  t2_dom <- if (length(t2) > 0 && !is.null(t2$domain)) as.character(t2$domain) else character(0)
  mode1 <- function(v) { v <- v[!is.na(v) & nzchar(v)]
                         if (length(v) == 0) "" else names(sort(table(v), decreasing = TRUE))[1] }

  rows <- lapply(sort(unique(key)), function(k) {
    parts <- strsplit(k, "\t", fixed = TRUE)[[1]]
    run <- parts[1]; trc <- parts[2]
    sel <- key == k
    arr <- reduce(granges(te[sel]), ignore.strand = TRUE)   # union of tagged arrays
    strand(arr) <- "*"
    array_bp <- sum(as.numeric(width(arr)))

    te_cls    <- mode1(as.character(te$TE_origin[sel]))
    te_struct <- mode1(as.character(te$TE_origin_structure[sel]))

    # Complete structural element copies overlapping the arrays -> count + covered bp.
    n_complete <- 0L; cov_bp <- 0
    if (length(elem) > 0 && length(arr) > 0) {
      h <- unique(subjectHits(findOverlaps(arr, elem)))
      n_complete <- length(h)
      if (n_complete > 0) cov_bp <- sum(as.numeric(width(intersect(arr, elem[h]))))
    }
    frac <- if (array_bp > 0) cov_bp / array_bp else NA_real_

    # Distinct DANTE protein domains present within the arrays.
    dom_str <- ""
    if (length(t2p) > 0 && length(t2_dom) > 0 && length(arr) > 0) {
      dh <- unique(queryHits(findOverlaps(t2p, arr)))
      if (length(dh) > 0) dom_str <- paste(.order_domains(t2_dom[dh]), collapse = "|")
    }

    # Authoritative tandem period (trc_table monomer_tarean->monomer_kite).
    per_map  <- if (run == "short") period_short else period_default
    mono_n   <- if (length(per_map) > 0 && !is.null(per_map[[trc]])) per_map[[trc]] else NA_integer_
    mono_n   <- suppressWarnings(as.integer(mono_n))
    n_exp    <- if (!is.na(mono_n) && mono_n > 0) as.integer(round(array_bp / mono_n)) else NA_integer_

    # Domain-rhythm metrics at that period (the gate criterion; NA if no period).
    rhy <- if (!is.na(mono_n) && mono_n > 0) te_domain_rhythm(arr, t2p, mono_n)
            else list(occ = NA_real_, frac = NA_real_)

    data.frame(
      trc_id                = trc,
      run                   = run,
      n_arrays              = sum(sel),
      total_array_bp        = as.integer(array_bp),
      monomer_length_bp     = if (is.na(mono_n)) "" else format(mono_n, trim = TRUE, scientific = FALSE),
      te_classification     = te_cls,
      te_origin_structure   = te_struct,
      protein_domains       = dom_str,
      n_complete_elements   = n_complete,
      n_expected_monomers   = if (is.na(n_exp)) "" else as.character(n_exp),
      complete_bp_fraction  = if (is.na(frac)) "" else sprintf("%.4f", frac),
      domain_occupancy      = if (is.na(rhy$occ))  "" else sprintf("%.4f", rhy$occ),
      frac_arrays_in_rhythm = if (is.na(rhy$frac)) "" else sprintf("%.4f", rhy$frac),
      stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  # Deterministic order: run (default before short), then numeric TRC id.
  trc_num <- suppressWarnings(as.integer(sub("^TRC_", "", df$trc_id)))
  df <- df[order(df$run, trc_num, df$trc_id), , drop = FALSE]
  write.table(df, out_csv, sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
  log_msg(sprintf("Wrote %d TE-derived TRC row(s) -> %s", nrow(df), basename(out_csv)))
  invisible()
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
  # Only features that overlap another need the greedy trim; the rest pass
  # through intact. A non-overlapping feature cannot overlap an overlapping one
  # (or it would appear in `h`), so a trim never has to consider the pass-through
  # set — restricting the loop to the overlapping subset is exact. This turns the
  # old O(N^2) (grow `kept` by one element per feature AND trim each feature
  # against the whole growing set) into O(k^2), k = overlapping tier-1 features in
  # the batch (usually tiny; a single overlapping pair used to drag the whole
  # batch's tier-1 set through the loop).
  involved    <- sort(unique(c(queryHits(h), subjectHits(h))))
  passthrough <- t1[-involved]
  # Process longest-first, with a FULLY deterministic tie-break on coordinates.
  # The greedy trim RESULT depends on processing order (an earlier-kept feature
  # trims later overlappers), so equal-width features must not be ordered by their
  # incoming position — that position is batch-composition-dependent (batches are
  # split by thread count), which made the trimmed output non-deterministic
  # run-to-run on multi-sequence genomes. Ordering ties by (seqname,start,end,
  # strand) — then classification,source_tool for genuinely identical intervals
  # (e.g. a DANTE_LTR complete + a DANTE_TIR call at the same span) — makes the
  # result independent of input order and thus of batching.
  inv  <- t1[involved]
  ord  <- order(-width(inv), as.character(seqnames(inv)), start(inv), end(inv),
                as.character(strand(inv)), as.character(inv$classification),
                as.character(inv$source_tool))
  t1s  <- inv[ord]
  # NOTE (run-000156, 94 Gbp): this loop was suspected of an O(k^2) blowup because
  # each batch spent ~17,045 s of its ~17,050 s between the "input sizes" and
  # "tier1 overlap resolve" log lines. Measured against the real tier-1 features of
  # that run's largest chromosome (28,022 features, 139 overlap pairs, involved=259
  # i.e. 0.9%, max degree 2) it runs in ~21 s — the observed 17,045 s was the host
  # thrashing under 55 concurrent workers, not this function. A candidate-bounded
  # variant was written and measured at 18 s on the same data (12% faster) and
  # reverted as not worth the added branching in a determinism-critical path. Do
  # not "optimise" this loop without first measuring `involved` on real data.
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

process_batch <- function(seqs, data, min_len, trc_origin_def = character(0), trc_origin_sho = character(0)) {
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

  # ── Tier-1 structural-overlap resolution (greedy longest-first) ──────────
  # MUST run BEFORE any trim_to_nonoverlap(t1, ...) call below. trim_to_nonoverlap
  # disjoins `c(lower, higher)`, so if `lower` (t1) still carries INTERNAL overlaps
  # it is silently decomposed into disjoint pieces — both overlapping elements kept
  # and split at the overlap, with the contested span assigned by input index
  # (min(revmap)) — which defeats the greedy AND is batch-composition-dependent.
  # The te_sat pre-pass below trims t1 against te_sat_r, and te_sat is non-empty
  # only in batches that happen to contain a TE-derived-TRC array; whether a given
  # sequence's t1 got pre-disjoined therefore depended on the thread-count batching
  # (single batch always has te_sat; a per-thread batch may not) — the source of
  # the residual tier-1 non-determinism on large multi-sequence genomes. Resolving
  # here first makes t1 internally non-overlapping, so the te_sat trim only carves
  # (never disjoins) and the greedy result is batch-invariant. No-op when nothing
  # overlaps (the common case). Clean tandem LTR-RT arrays were already collapsed
  # upstream (resolve_ltr_tandems.py); this handles residual cross-lineage /
  # cross-tool / partial overlaps.
  n_t1_before <- length(t1)
  t1 <- resolve_tier1_overlaps(t1, min_len)
  if (length(t1) != n_t1_before)
    log_msg(sprintf("  %s  tier1 overlap resolve: %d → %d features",
                    batch_label, n_t1_before, length(t1)))

  # ── Pre-pass: TR-from-structural-TE (Item 3) ─────────────────────────────
  # Flag length-qualified Tier-3 clustering TRCs that tandemly stack ≥2 same-
  # family structural TEs. Those satellites WIN: tag them TE_origin and trim the
  # covered structural Tier-1/2 features out (they remain in the DANTE_* GFF3s).
  # The TRC->TE_origin decision (`trc_origin_def`/`trc_origin_sho`) is computed
  # ONCE at the top level over ALL of each TRC's arrays genome-wide, then passed
  # in — NOT recomputed per batch. A TRC's arrays can span sequences that land in
  # different batches, and the batch split itself is thread-count-dependent
  # (threads=1 -> one batch, else N); a per-batch decision (coverage AND the
  # domain-rhythm gate judged over only the arrays in one batch) would make
  # TE_origin thread-dependent. Global judging keeps it deterministic. Here we
  # only APPLY the decision to this batch's arrays (tagging + trimming).
  te_sat <- GRanges()
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
    # t1 is already internally non-overlapping (resolved above), so this trim only
    # carves te_sat out — it does NOT disjoin t1's structure. t2 is deliberately
    # NOT trimmed here: trim_to_nonoverlap would disjoin t2's INTERNAL overlaps,
    # and te_sat is present only in batches that contain a TE-derived array, so a
    # pre-disjoin here would make the Step-2 resolve_within_tier a no-op in some
    # batches but not others (batch-dependent strand/geometry — the tier-2 twin of
    # the tier-1 bug). Instead te_sat is folded into Step 2's `higher` set (it is
    # in level1 by then), so t2 is trimmed-and-resolved exactly once, batch-invariantly.
    t1 <- trim_to_nonoverlap(t1, te_sat_r, min_len)
    if (length(t1_members) > 0)
      t1_members <- t1_members[!overlapsAny(t1_members, te_sat_r, ignore.strand = TRUE)]
    t1_ltr <- if (length(t1) > 0) t1[t1$source_tool == "DANTE_LTR"] else GRanges()
    log_msg(sprintf("  %s  pre-pass: %d TE-derived TRC(s) → te_sat=%d feats; t1→%d t1_members→%d",
                    batch_label, length(trc_origin_def) + length(trc_origin_sho),
                    length(te_sat), length(t1), length(t1_members)))
  }

  # t1 is now fully resolved (greedy longest-first, above) and te_sat-carved
  # (pre-pass, above); both preserve internal non-overlap, so no second
  # resolve pass is needed. Refresh the DANTE_LTR view used for within-LTR
  # nesting in Steps 3 and 8.
  t1_ltr <- if (length(t1) > 0) t1[t1$source_tool == "DANTE_LTR"] else GRanges()

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

  # ── Step 2: Tier 2 — DANTE domains, trim against placed Tier 1 (+ te_sat) ─
  # Trim against reduce(level1) — i.e. te_sat ∪ Tier-1 — NOT bare t1. Folding the
  # te_sat carve in here (rather than a separate pre-pass trim) means t2 is
  # disjoined exactly once, by resolve_within_tier, so its geometry/strand is
  # batch-invariant. When there are no TE-derived TRCs, level1 == t1 here, so this
  # is identical to the previous `trim_to_nonoverlap(t2, t1)`.
  higher_2 <- if (length(level1) > 0) suppressWarnings(reduce(level1, ignore.strand = TRUE)) else GRanges()
  t2_trimmed <- timed(sprintf("  %s  step2: trim Tier2 (%d) vs placed (%d regions)",
                               batch_label, length(t2), length(higher_2)),
                       resolve_within_tier(trim_to_nonoverlap(t2, higher_2, min_len)))
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

  # Report this worker's peak RSS back to the parent. In a forked child VmHWM
  # covers the inherited (copy-on-write) parent pages as well as the child's own
  # allocations, so it is the per-worker footprint to divide a host budget by.
  peak <- mem_hwm_mb()
  toc(t_batch, sprintf("final: level1=%d, level2=%d%s",
                       length(level1), length(level2),
                       if (is.na(peak)) "" else sprintf(", peak_rss=%.1fG", peak/1024)))
  list(level1 = level1, level2 = level2, peak_rss_mb = peak)
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
  # LTR_RT_TR membership tags live only on Level-2 members (set below). Level-1
  # features carry them as NA so level1 / level2 keep one mcols schema for the
  # c() combine (orphan promotion + final) and the GFF3 export.
  level1$in_structure <- NA_character_
  level1$member_of    <- NA_character_

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
    # ── LTR_RT_TR membership tags ────────────────────────────────────────────
    # The ONLY Level-2 DANTE_LTR features are tandem LTR-RT member copies
    # (standalone LTR-RTs are Level 1; the LTR protein-domain/LTR children are
    # not emitted to the unified GFF3), so source_tool=="DANTE_LTR" at Level 2
    # uniquely identifies a tandem member. Tag it so a consumer can spot a member
    # directly (in_structure) and jump to its container (member_of == the Parent
    # just assigned) without having to look the parent up.
    level2$in_structure <- NA_character_
    level2$member_of    <- NA_character_
    is_member <- as.character(level2$source_tool) == "DANTE_LTR" & !is.na(level2$Parent)
    level2$in_structure[is_member] <- "LTR_RT_TR"
    level2$member_of[is_member]    <- as.character(level2$Parent)[is_member]

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

  # Canonicalise the mcols columns before export so the GFF3 is BYTE-identical
  # across thread counts, not just data-identical. Two batch-dependent effects:
  #   1. Column ORDER — rtracklayer emits column-9 attributes in mcols order, which
  #      was inherited from whichever batch combined first (e.g. `ID` landed first
  #      under multi-batch but after source_tool under a single batch).
  #   2. Column SET — a c()-combine across batches can carry an entirely-NA column
  #      that a single batch does not (all its features lacked that field). An
  #      all-NA column emits no key=value, but if it sorts AFTER the last real
  #      attribute rtracklayer still writes a trailing ';' for it (verified), so the
  #      raw line differs by a stray ';' (seen only on LTR_RT_TR containers, where
  #      copy_number is last). Dropping all-NA columns removes that ambiguity and
  #      is loss-free — those columns contribute nothing to the output.
  {
    keep <- vapply(mcols(all_feats), function(v) !all(is.na(v)), logical(1))
    keep <- keep | names(keep) %in% c("type", "source")  # never drop GFF3 reserved cols
    mcols(all_feats) <- mcols(all_feats)[, keep, drop = FALSE]
    cur  <- colnames(mcols(all_feats))
    pref <- c("source", "type", "score", "phase", "ID", "Parent", "Name",
              "classification", "source_tier", "source_tool", "element_type",
              "in_structure", "member_of", "TE_origin", "TE_origin_structure",
              "structure", "copy_number")
    ord  <- c(intersect(pref, cur), sort(setdiff(cur, pref)))
    mcols(all_feats) <- mcols(all_feats)[, ord, drop = FALSE]
  }

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
      # Insert after the existing ##gff-version line (always first) by streaming
      # into a sibling temp file — readLines(output_path) used to pull the entire
      # multi-GB unified GFF3 into memory just to prepend 3 lines.
      tmp_out <- tempfile(tmpdir = dirname(output_path))
      con_in  <- file(output_path, "r")
      con_out <- file(tmp_out, "w")
      # Safety net for the error path (a mid-stream readLines/writeLines throw
      # leaves the connections open). On the success path the explicit close()
      # below has already run, so these close() calls hit an already-closed
      # connection — `close()` errors on that ("invalid connection"), so wrap in
      # try(). NB: do NOT guard with isOpen() here: isOpen() itself *raises*
      # "invalid connection" on a closed/destroyed connection (it does not return
      # FALSE), which is exactly what crashed make_unified after a clean write.
      on.exit({ try(close(con_in),  silent = TRUE)
                try(close(con_out), silent = TRUE) }, add = TRUE)
      first <- readLines(con_in, n = 1L)
      writeLines(c(first, header), con_out)
      repeat {
        chunk <- readLines(con_in, n = 100000L)
        if (length(chunk) == 0L) break
        writeLines(chunk, con_out)
      }
      close(con_in); close(con_out)
      file.rename(tmp_out, output_path)
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
log_mem("after loading all tiers")

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

# Authoritative per-TRC tandem periods (trc_table monomer_tarean->monomer_kite),
# for the domain-rhythm gate and the TE-derived TRC table's monomer column.
# Default and short runs have independent TRC_<n> spaces.
period_def   <- read_trc_periods(opt$tc_trc_table_default)
period_short <- read_trc_periods(opt$tc_trc_table_short)
log_msg("TRC periods (trc_table): ", length(period_def), " default + ",
        length(period_short), " short")

# TE-derived TRC decision — computed GLOBALLY, once, over every array of each TRC
# genome-wide (NOT per batch). process_batch only applies the decision to its
# arrays. This is essential for determinism: batch composition is thread-count-
# dependent (threads=1 -> one batch), and a per-batch coverage+rhythm judgement
# over a partial set of a TRC's arrays would flip TE_origin between thread counts.
trc_origin_def <- if (length(tc_data$default) > 0 && length(t1) > 0)
  identify_te_derived_trcs(tc_data$default, t1, t2, period_def) else character(0)
trc_origin_sho <- if (length(tc_data$short) > 0 && length(t1) > 0)
  identify_te_derived_trcs(tc_data$short, t1, t2, period_short) else character(0)
log_msg("TE-derived TRCs (global): ", length(trc_origin_def), " default + ",
        length(trc_origin_sho), " short")

# Run resolution — always via mclapply so the code path is uniform.
log_msg(sprintf("── Running tier resolution (%d threads, %d batch(es)) ─────",
                opt$threads, length(batches)))

# ── Worker-pool sizing ───────────────────────────────────────────────────────
# NOTE: this caps CONCURRENCY ONLY. Batch composition is fixed above from
# --threads/--batch_size, so nothing here can change the result — only how many
# batches run at once.
#
# Why a cap is needed at all (measured on run-000156, 94.26 Gbp, 768 GB host):
# mclapply forks, and R's GC writes to object headers in the child, so the
# inherited parent heap gets privately copied. Each worker's peak RSS therefore
# converges on the PARENT's, regardless of how much data its own batch holds —
# measured 48.3 GB for every one of 55 batches (min = median = max) against a
# 48.4 GB parent, for batches spanning 143 Mb to 2.15 Gb of sequence. At 55
# workers that demands ~2.66 TB on a 768 GB box: the run thrashed (388 M minor
# faults at only 4 workers ≈ 1.5 TB of page copying), 15 workers were OOM-killed
# and the rule died after 5h50m. The same work at 4 workers took 885 s.
n_workers <- min(opt$threads, length(batches))
parent_rss <- mem_rss_mb()
log_mem(sprintf("at fork (before gating: %d worker(s))", n_workers))

# 1. Hard ceiling — this rule gains little from wide parallelism (tier resolution
#    was 885 s of a 6,949 s run; loading and GFF3 export are serial), while each
#    worker costs a parent-heap copy. 0 disables.
if (opt$max_workers > 0 && n_workers > opt$max_workers) {
  log_msg(sprintf("[mem] capping workers %d -> %d (--max_workers)",
                  n_workers, opt$max_workers))
  n_workers <- opt$max_workers
}

# 2. Memory gate — budget / per-worker cost, using the parent RSS as the
#    per-worker estimate (justified by the measurement above). Leave 20% headroom
#    for the parent's own growth during combine/finalise, which peaked well above
#    the fork-point RSS (121.9 GB vs 48.4 GB on run-000156).
{
  bud <- mem_budget_mb(opt$mem_budget_gb)
  if (!is.na(bud$mb) && !is.na(parent_rss) && parent_rss > 0) {
    cap <- max(1L, as.integer(floor((bud$mb * 0.8) / parent_rss)))
    log_msg(sprintf("[mem] budget %.1fG (%s), parent rss %.1fG -> room for %d worker(s)",
                    bud$mb/1024, bud$src, parent_rss/1024, cap))
    if (cap < n_workers) {
      log_msg(sprintf("[mem] capping workers %d -> %d (memory gate)", n_workers, cap))
      n_workers <- cap
    }
  } else {
    log_msg(sprintf("[mem] no memory gate applied (budget=%s, parent rss=%s)",
                    ifelse(is.na(bud$mb), "unknown", sprintf("%.1fG", bud$mb/1024)),
                    ifelse(is.na(parent_rss), "unknown", sprintf("%.1fG", parent_rss/1024))))
  }
}
log_msg(sprintf("── Worker pool: %d concurrent worker(s) for %d batch(es) ─────",
                n_workers, length(batches)))
t_resolve <- tic("tier resolution (all batches)")
results <- mclapply(batches,
                    function(seqs) process_batch(seqs, data, opt$min_feature_length,
                                                 trc_origin_def, trc_origin_sho),
                    mc.cores = n_workers)
toc(t_resolve)
log_mem("after tier resolution (parent)")

# Per-worker peaks reported back by the children (NA where /proc is unreadable).
{
  peaks <- suppressWarnings(as.numeric(vapply(results, function(r)
    if (is.list(r) && !is.null(r$peak_rss_mb)) r$peak_rss_mb else NA_real_, numeric(1))))
  peaks <- peaks[is.finite(peaks)]
  if (length(peaks) > 0)
    log_msg(sprintf("[mem] worker peak rss: min=%.1fG median=%.1fG max=%.1fG over %d batch(es)",
                    min(peaks)/1024, median(peaks)/1024, max(peaks)/1024, length(peaks)))
}

# Validate every batch result POSITIVELY.
#
# An mclapply child that dies on a signal (OOM killer, segfault) does NOT come
# back as a "try-error" — mclapply only warns "scheduled cores ... did not
# deliver results" and leaves a non-conforming element in its place. The old
# try-error-only check therefore passed with 15 of 55 workers dead on the 94 Gbp
# run-000156; execution continued, do.call(c, ...) degraded to a plain list
# ("Combined: 40 Level 1" — batches, not features), and the run died much later
# in finalise_output with a misleading S4 dispatch error. Had the survivors
# happened to concatenate cleanly it would have written an annotation missing
# 27% of the genome and exited 0. So: require the expected count, and require
# every element to actually carry GRanges level1/level2.
if (length(results) != length(batches))
  stop(sprintf("mclapply returned %d result(s) for %d batch(es) — worker loss",
               length(results), length(batches)))
bad_idx <- which(!vapply(results, function(r) {
  is.list(r) && !is.null(r$level1) && !is.null(r$level2) &&
    is(r$level1, "GRanges") && is(r$level2, "GRanges")
}, logical(1)))
if (length(bad_idx) > 0) {
  for (i in bad_idx) {
    detail <- if (inherits(results[[i]], "try-error")) as.character(results[[i]])
              else sprintf("no result returned (worker died — check for OOM); got %s",
                           paste(class(results[[i]]), collapse = "/"))
    seqs_i <- batches[[i]]
    label  <- paste0(paste(head(seqs_i, 3), collapse = ","),
                     if (length(seqs_i) > 3) sprintf(",...+%d", length(seqs_i) - 3) else "")
    message("Batch ", i, " [", label, "] failed: ", detail)
  }
  stop(sprintf("%d of %d batch(es) did not return a usable result — refusing to ",
               length(bad_idx), length(results)),
       "write a truncated annotation. See the per-batch messages above.")
}

# Combine batch results
log_msg("── Combining batch results ─────────────────────────────────────")
level1_all <- timed("combine level1 across batches",
                    suppressWarnings(do.call(c, lapply(results, `[[`, "level1"))))
level2_all <- timed("combine level2 across batches",
                    suppressWarnings(do.call(c, lapply(results, `[[`, "level2"))))
log_msg("Combined: ", length(level1_all), " Level 1, ",
        length(level2_all), " Level 2 features")

# TE-derived tandem-repeat summary (TRCs the pre-pass tagged with TE_origin).
# Computed from the intact combined level1 + the untrimmed global t1/t2 tiers.
timed("write TE-derived TRC table",
      write_te_derived_trc_table(
        level1_all, t1, ltr_data$members, t2,
        period_def, period_short,
        sub("\\.gff3$", ".te_derived_trc.csv", opt$output)))

# Finalise and export
log_mem("after combining batch results")

log_msg("── Finalising output ──────────────────────────────────────────")
all_feats <- finalise_output(level1_all, level2_all, seqlengths_vec, opt$output)

# Sanity checks on the in-memory GRanges (avoids re-parsing the GFF3)
sanity_check(all_feats, seqlengths_vec,
             overlaps_tsv = sub("\\.gff3$", ".overlaps.tsv", opt$output))

log_mem("final (whole-run peak)")
log_msg("Done. Output: ", opt$output,
        sprintf("  (total wall time: %.1fs)", proc.time()[3] - .script_start))
