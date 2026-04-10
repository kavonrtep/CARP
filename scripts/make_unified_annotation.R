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

# ── 0. CLI arguments ─────────────────────────────────────────────────────────

option_list <- list(
  make_option("--ltr",              type="character", help="DANTE_LTR GFF3"),
  make_option("--tir",              type="character", help="DANTE_TIR GFF3"),
  make_option("--line",             type="character", help="DANTE_LINE GFF3"),
  make_option("--dante",            type="character", help="DANTE filtered domains GFF3"),
  make_option("--tc_default",       type="character", help="TideCluster default GFF3"),
  make_option("--tc_short",         type="character", help="TideCluster short monomer GFF3"),
  make_option("--tc_rm",            type="character", help="RM on TideCluster library GFF3"),
  make_option("--rm",               type="character", help="RepeatMasker+DANTE merged GFF3"),
  make_option("--th_default",       type="character", help="TideHunter default residuals GFF3"),
  make_option("--th_short",         type="character", help="TideHunter short residuals GFF3"),
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

# Convert | separator to /  in classification paths (DANTE_LTR / DANTE_LINE)
fix_sep <- function(x) gsub("|", "/", as.character(x), fixed = TRUE)

# Convert DANTE_TIR underscore-encoded hierarchy to slash-separated
# e.g. "Class_II_Subclass_1_TIR_hAT" → "Class_II/Subclass_1/TIR/hAT"
convert_tir_cls <- function(x) {
  x <- as.character(x)
  x <- sub("^Class_II_Subclass_1_", "Class_II/Subclass_1/", x)
  x <- sub("^Class_II_Subclass_2_", "Class_II/Subclass_2/", x)
  x
}

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
                "source_tool", "element_type")

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
  mcols(gr) <- if (length(md) > 0) do.call(S4Vectors::DataFrame, md) else S4Vectors::DataFrame()
  # Standardize seqlevels to the batch set
  suppressWarnings(seqlevels(gr, pruning.mode = "coarse") <- seqs)
  gr
}

# ── 2. Per-tool loading functions ─────────────────────────────────────────────

load_tier1_ltr <- function(path) {
  message("Loading DANTE_LTR: ", path)
  raw <- safe_import(path)
  if (length(raw) == 0) return(list(top = GRanges(), children = GRanges()))

  top      <- raw[raw$type == "transposable_element"]
  children <- raw[raw$type %in% c("protein_domain", "long_terminal_repeat",
                                   "target_site_duplication", "primer_binding_site")]
  cls <- fix_sep(top$Final_Classification)
  top <- set_meta(top, cls, cls, 1L, "DANTE_LTR")
  top$element_type <- ifelse(grepl("^TE_partial_", top$ID, perl = TRUE),
                             "partial", "complete")

  # Strip verbose alignment attributes to reduce memory
  verbose_attrs <- c("DB_Seq", "Region_Seq", "Query_Seq", "Best_Hit_DB_Pos",
                     "cls_prefilter", "neighbors_count",
                     "upstream_domain", "downstream_domain", "domain_order")
  for (a in verbose_attrs) mcols(children)[[a]] <- NULL

  message("  ", length(top), " top-level LTR elements (",
          sum(top$element_type == "complete"), " complete, ",
          sum(top$element_type == "partial"), " partial)")
  list(top = top, children = children)
}

load_tier1_tir <- function(path) {
  message("Loading DANTE_TIR: ", path)
  raw <- safe_import(path)
  if (length(raw) == 0) return(GRanges())

  raw$type <- "transposable_element"
  cls <- convert_tir_cls(raw$Classification)
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
  cls <- fix_sep(raw$Final_Classification)
  raw <- set_meta(raw, cls, cls, 2L, "DANTE")
  message("  ", length(raw), " domain features")
  raw
}

load_tier3_tidecluster <- function(path_default, path_short) {
  message("Loading TideCluster default: ", path_default)
  def <- safe_import(path_default)
  message("Loading TideCluster short monomer: ", path_short)
  sho <- safe_import(path_short)

  normalise_tc <- function(gr, tool) {
    if (length(gr) == 0) return(gr)
    gr$type <- "repeat_region"
    cls <- paste0("Satellite/TideCluster/", as.character(gr$Name))
    set_meta(gr, as.character(gr$Name), cls, 3L, tool)
  }
  list(
    default = normalise_tc(def, "TideCluster_default"),
    short   = normalise_tc(sho, "TideCluster_short")
  )
}

load_tier4_tc_rm <- function(path) {
  message("Loading RM on TideCluster library: ", path)
  raw <- safe_import(path)
  if (length(raw) == 0) return(GRanges())

  raw$type <- "repeat_region"
  cls <- paste0("Satellite/TideCluster/", as.character(raw$Name))
  raw <- set_meta(raw, as.character(raw$Name), cls, 4L, "TideCluster_RM")
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
# Fragments shorter than min_len bp are discarded.
# Metadata is propagated from the source feature to all resulting fragments.
trim_to_nonoverlap <- function(lower, higher, min_len = 50L) {
  if (length(lower) == 0 || length(higher) == 0) return(lower)

  higher_r <- reduce(higher, ignore.strand = TRUE)
  hits     <- findOverlaps(lower, higher_r, ignore.strand = TRUE)
  if (length(hits) == 0) return(lower)

  ov_idx  <- as.integer(unique(queryHits(hits)))
  intact  <- lower[setdiff(seq_along(lower), ov_idx)]

  # Columns to propagate from the source feature to trimmed fragments
  keep_cols <- intersect(c("Name", "classification", "source_tier", "source_tool",
                            "element_type"),
                         colnames(mcols(lower)))

  trimmed_list <- lapply(ov_idx, function(i) {
    feat    <- lower[i]
    feat_sl <- as.character(seqnames(feat))
    # Restrict higher_r to the seqname of this feature so setdiff seqlevels match
    hr_sub  <- keepSeqlevels(higher_r,
                              intersect(seqlevels(higher_r), feat_sl),
                              pruning.mode = "coarse")
    seqlevels(feat, pruning.mode = "coarse") <- feat_sl
    remaining <- GenomicRanges::setdiff(feat, hr_sub, ignore.strand = TRUE)
    if (length(remaining) == 0) return(NULL)
    remaining <- remaining[width(remaining) >= min_len]
    if (length(remaining) == 0) return(NULL)
    for (col in keep_cols) {
      val <- mcols(lower)[[col]][i]
      mcols(remaining)[[col]] <- rep(val, length(remaining))
    }
    remaining
  })

  trimmed_list <- Filter(Negate(is.null), trimmed_list)
  if (length(trimmed_list) == 0) return(intact)
  suppressWarnings(c(intact, do.call(c, trimmed_list)))
}

# For each feature in `children`, return the ID of the Level 1 parent
# (the overlapping feature in `parents` with the maximum intersection width).
get_parent_id <- function(children, parents) {
  n <- length(children)
  if (n == 0 || length(parents) == 0)
    return(rep(NA_character_, n))

  hits <- findOverlaps(children, parents)
  if (length(hits) == 0)
    return(rep(NA_character_, n))

  # Manual intersection width (avoids pintersect version incompatibilities)
  ch_r   <- children[queryHits(hits)]
  par_r  <- parents[subjectHits(hits)]
  ov_w   <- pmax(0L, pmin(end(ch_r), end(par_r)) - pmax(start(ch_r), start(par_r)) + 1L)
  qi     <- queryHits(hits)
  result <- rep(NA_character_, n)
  for (ci in unique(qi)) {
    idx  <- which(qi == ci)
    best <- subjectHits(hits)[idx[which.max(ov_w[idx])]]
    result[ci] <- parents$ID[best]
  }
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
  message("  Processing batch: ", paste(head(seqs, 3), collapse=", "),
          if (length(seqs) > 3) paste0(" ... (", length(seqs), " seqs)") else "")

  sub     <- lapply(data, subset_seqs, seqs = seqs)
  t1      <- sub$t1
  t1_ltr  <- if (length(t1) > 0) t1[t1$source_tool == "DANTE_LTR"] else GRanges()
  t2      <- sub$t2
  t3_def  <- sub$t3_def
  t3_sho  <- sub$t3_sho
  t4      <- sub$t4
  t5_te   <- sub$t5_te
  t5_sc   <- sub$t5_sc
  t6      <- sub$t6

  level1 <- GRanges()
  level2 <- GRanges()

  # ── Step 1: Tier 1 — structure-based elements (no trimming) ──────────────
  if (length(t1) > 0) level1 <- t1

  # ── Step 2: Tier 2 — DANTE domains, trim against Tier 1 ──────────────────
  t2_trimmed <- resolve_within_tier(trim_to_nonoverlap(t2, t1, min_len))
  if (length(t2_trimmed) > 0) level1 <- suppressWarnings(c(level1, t2_trimmed))

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
  if (length(t3_def_l1) > 0) level1 <- suppressWarnings(c(level1, t3_def_l1))

  # ── Step 4: Tier 3 short monomer — trim against Tier 3 default + Tier 1 ──
  t3s_trimmed <- trim_to_nonoverlap(t3_sho, suppressWarnings(c(t3_def_l1, t1)), min_len)
  if (length(t3s_trimmed) > 0) level1 <- suppressWarnings(c(level1, t3s_trimmed))

  # ── Step 5: Tier 4 — RM on TideCluster library, trim against Tiers 1–3 ──
  higher_1_3  <- if (length(level1) > 0) suppressWarnings(reduce(level1, ignore.strand = TRUE)) else GRanges()
  t4_trimmed  <- resolve_within_tier(trim_to_nonoverlap(t4, higher_1_3, min_len))
  if (length(t4_trimmed) > 0) level1 <- suppressWarnings(c(level1, t4_trimmed))

  # ── Step 6: Tier 5 TE hits — trim against Tiers 1–4 ─────────────────────
  higher_1_4 <- if (length(level1) > 0) suppressWarnings(reduce(level1, ignore.strand = TRUE)) else GRanges()
  t5_trimmed <- resolve_within_tier(trim_to_nonoverlap(t5_te, higher_1_4, min_len))
  if (length(t5_trimmed) > 0) level1 <- suppressWarnings(c(level1, t5_trimmed))

  # ── Step 7: Tier 5 Simple/Low complexity ─────────────────────────────────
  # Overlapping an existing Level 1 feature → Level 2 nested
  # Not overlapping → Level 1
  if (length(t5_sc) > 0) {
    cur_l1 <- if (length(level1) > 0) suppressWarnings(reduce(level1, ignore.strand = TRUE)) else GRanges()
    if (length(cur_l1) > 0) {
      sc_hits <- suppressWarnings(findOverlaps(t5_sc, cur_l1, ignore.strand = TRUE))
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
  }

  # ── Step 8: Tier 6 — TideHunter residuals ────────────────────────────────
  # Within LTR body → Level 2 nested
  # Not covered by any higher tier → Level 1
  if (length(t6) > 0) {
    all_higher <- if (length(level1) > 0) suppressWarnings(reduce(level1, ignore.strand = TRUE)) else GRanges()

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
  }

  list(level1 = level1, level2 = level2)
}

# ── 6. Output assembly ────────────────────────────────────────────────────────

finalise_output <- function(level1, level2, seqlengths_vec, output_path) {
  message("Finalising: ", length(level1), " Level 1, ",
          length(level2), " Level 2 features")

  # ── Sort Level 1 ──────────────────────────────────────────────────────────
  all_seqs <- names(seqlengths_vec)
  seqlevels(level1, pruning.mode = "coarse") <- all_seqs
  seqlengths(level1) <- seqlengths_vec
  level1 <- sort(sortSeqlevels(level1))

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
    level2 <- sort(sortSeqlevels(level2))
    level2$ID <- paste0("UA_L2_", formatC(seq_along(level2), width = 8, flag = "0"))

    # For Level 2 features with temp_parent_tool == "DANTE_LTR": parent from LTR tier1
    ltr_l1   <- level1[level1$source_tool == "DANTE_LTR"]
    any_l1   <- level1

    parent_ids <- character(length(level2))
    for (i in seq_along(level2)) {
      ptype <- if (!is.null(level2$temp_parent_tool)) level2$temp_parent_tool[i] else "any"
      ptype <- if (is.na(ptype)) "any" else ptype
      pool  <- if (ptype == "DANTE_LTR") ltr_l1 else any_l1
      pid   <- get_parent_id(level2[i], pool)
      parent_ids[i] <- if (!is.na(pid[1])) pid[1] else NA_character_
    }
    level2$Parent <- parent_ids
    # Remove Level 2 features with no parent found (can't be nested)
    orphans <- is.na(level2$Parent)
    if (any(orphans)) {
      message("  Promoting ", sum(orphans), " Level 2 orphans to Level 1")
      orphan_feats <- level2[orphans]
      orphan_feats$Parent <- NULL
      orphan_feats$temp_parent_tool <- NULL
      level1 <- c(level1, orphan_feats)
      level2 <- level2[!orphans]
    }
    level2$temp_parent_tool <- NULL
    level2$type   <- ifelse(
      grepl("^Satellite|^Simple_repeat|^Low_complexity|^rDNA|^Unknown",
            level2$classification),
      "repeat_region", "transposable_element"
    )
    level2$source <- level2$source_tool
  }

  # ── Clean up temp columns from Level 1 ───────────────────────────────────
  level1$temp_parent_tool <- NULL

  # ── Combine, final sort, export ───────────────────────────────────────────
  all_feats <- if (length(level2) > 0) c(level1, level2) else level1
  seqlevels(all_feats, pruning.mode = "coarse") <- all_seqs
  seqlengths(all_feats) <- seqlengths_vec
  all_feats <- sort(sortSeqlevels(all_feats))

  message("Exporting to: ", output_path)
  export(all_feats, output_path, format = "gff3")

  invisible(all_feats)
}

# ── 7. Sanity checks ──────────────────────────────────────────────────────────

sanity_check <- function(output_path, seqlengths_vec) {
  message("\n── Sanity checks ──────────────────────────────────────────────")
  gr <- import.gff3(output_path)

  # After GFF3 round-trip, multi-valued attributes come back as CharacterList;
  # use the standard source column (col 2) and unlist for safe access.
  parent_col <- suppressWarnings(gr$Parent)
  if (is(parent_col, "CharacterList")) {
    has_parent <- lengths(parent_col) > 0
  } else {
    has_parent <- !is.na(parent_col)
  }
  l1 <- gr[!has_parent]
  l2 <- gr[has_parent]

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
  print(sort(table(as.character(l1$source)), decreasing = TRUE))
  message("\nLevel 1 breakdown by source_tier attribute:")
  tier_col <- suppressWarnings(l1$source_tier)
  if (!is.null(tier_col)) print(sort(table(as.character(tier_col))))
  message("──────────────────────────────────────────────────────────────\n")
}

# ── 8. Main ───────────────────────────────────────────────────────────────────

message("=== make_unified_annotation.R ===")

# Load all tiers
message("\n── Loading inputs ─────────────────────────────────────────────")
ltr_data  <- load_tier1_ltr(opt$ltr)
tir_data  <- load_tier1_tir(opt$tir)
line_data <- load_tier1_line(opt$line)

t1 <- c(ltr_data$top, tir_data, line_data$top)
message("Tier 1 total: ", length(t1), " top-level features")

t2 <- load_tier2_dante(opt$dante)
message("Tier 2 total: ", length(t2), " features")

tc_data <- load_tier3_tidecluster(opt$tc_default, opt$tc_short)
message("Tier 3 total: ", length(tc_data$default), " default + ",
        length(tc_data$short), " short monomer features")

t4 <- load_tier4_tc_rm(opt$tc_rm)
message("Tier 4 total: ", length(t4), " features")

rm_data <- load_tier5_rm(opt$rm)
message("Tier 5 total: ", length(rm_data$te), " TE + ",
        length(rm_data$simple), " Simple/Low features")

t6 <- load_tier6_tidehunter(opt$th_default, opt$th_short)
message("Tier 6 total: ", length(t6), " features")

# Read genome FAI
message("\n── Reading FAI ────────────────────────────────────────────────")
seqlengths_vec <- read_fai(opt$fai)
genome_bp      <- sum(as.numeric(seqlengths_vec))
message("Genome: ", length(seqlengths_vec), " sequences, ",
        format(genome_bp, big.mark=","), " bp total")

# Decide whether to chunk
data <- list(
  t1     = t1,
  t2     = t2,
  t3_def = tc_data$default,
  t3_sho = tc_data$short,
  t4     = t4,
  t5_te  = rm_data$te,
  t5_sc  = rm_data$simple,
  t6     = t6
)

if (genome_bp <= opt$chunk_threshold) {
  message("\n── Processing genome as single batch (", format(genome_bp, big.mark=","),
          " bp ≤ threshold) ─────────")
  batches <- list(names(seqlengths_vec))
} else {
  message("\n── Activating chunked processing (genome ", format(genome_bp, big.mark=","),
          " bp > threshold) ───────")
  batches <- make_batches(seqlengths_vec, opt$batch_size)
  message("Created ", length(batches), " batches (target ", opt$batch_size/1e6, " Mb each)")
}

# Run resolution — parallel over batches
message("\n── Running tier resolution (", opt$threads, " threads, ",
        length(batches), " batch(es)) ─────────────────────")

if (length(batches) == 1) {
  results <- list(process_batch(batches[[1]], data, opt$min_feature_length))
} else {
  results <- mclapply(batches,
                      function(seqs) process_batch(seqs, data, opt$min_feature_length),
                      mc.cores = opt$threads)
}

# Combine batch results
message("\n── Combining results ──────────────────────────────────────────")
level1_all <- do.call(c, lapply(results, `[[`, "level1"))
level2_all <- do.call(c, lapply(results, `[[`, "level2"))
message("Combined: ", length(level1_all), " Level 1, ",
        length(level2_all), " Level 2 features")

# Finalise and export
message("\n── Finalising output ──────────────────────────────────────────")
finalise_output(level1_all, level2_all, seqlengths_vec, opt$output)

# Sanity checks
sanity_check(opt$output, seqlengths_vec)

message("Done. Output: ", opt$output)
