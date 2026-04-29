#!/usr/bin/env Rscript
# make_repeat_report.R
# Generates a standalone interactive HTML repeat annotation report.
# Usage: make_repeat_report.R --output_dir <dir> [options]

suppressPackageStartupMessages({
  library(optparse)
  library(rtracklayer)
  library(IRanges)
  library(GenomicRanges)
  library(jsonlite)
})

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "classification.R"))

# ═══════════════════════════════════════════════════════════════════════════
# A. ARGUMENT PARSING
# ═══════════════════════════════════════════════════════════════════════════

option_list <- list(
  make_option("--output_dir",       type = "character", help = "Pipeline output directory [required]"),
  make_option("--bin_width",        type = "integer",   default = 100000L,
              help = "Density track bin width in bp [default %default]"),
  make_option("--min_len_chart",    type = "integer",   default = 500000L,
              help = "Min seq length for own bar in composition chart [default %default]"),
  make_option("--min_len_tracks",   type = "integer",   default = 1000000L,
              help = "Min seq length for density track panel [default %default]"),
  make_option("--max_tracks",       type = "integer",   default = 50L,
              help = "Max sequences in density track dropdown [default %default]"),
  make_option("--top_sat_clusters", type = "integer",   default = 10L,
              help = "Top N satellite clusters shown in section 5 [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$output_dir)) stop("--output_dir is required")
outdir <- normalizePath(opt$output_dir, mustWork = TRUE)

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1])) a else b

# ═══════════════════════════════════════════════════════════════════════════
# B. DATA LOADING
# ═══════════════════════════════════════════════════════════════════════════

# ── B1. Genome sequence info from FAI ──────────────────────────────────────
load_genome_info <- function(outdir) {
  fai <- file.path(outdir, "genome_cleaned.fasta.fai")
  if (!file.exists(fai)) stop("FAI file not found: ", fai)
  df <- read.table(fai, header = FALSE, sep = "\t",
                   col.names = c("seqname", "length", "offset", "linebases", "linewidth"))
  df[, c("seqname", "length")]
}

# ── B2. Repeat composition from summary_statistics.csv ────────────────────
load_composition <- function(outdir) {
  f <- file.path(outdir, "summary_statistics.csv")
  df <- read.table(f, header = TRUE, sep = "\t", quote = '"', check.names = FALSE)
  colnames(df) <- c("type", "bp", "pct")
  df$bp  <- as.numeric(df$bp)
  df$pct <- as.numeric(df$pct)
  df <- df[!is.na(df$bp), ]
  df
}

# ── B3. DANTE_LTR stats from GFF3 ─────────────────────────────────────────
# Count only *complete* LTR-RTs for the "Complete TEs" column of the
# classification table. DANTE_LTR labels the element's completeness via
# its Rank attribute:
#   DL, DLT, DLP, DLTP  → complete (domains + LTRs + optional PBS / TSD)
#   D                   → partial — an LTR-RT region with protein
#                         domains identified but no reconstructed LTR
#                         sequences; the element boundaries are therefore
#                         unknown.
# Including Rank=D over-counts by ~3x on a typical plant genome (most
# of the D-rank calls are fragmented / solo-LTR-less copies without
# reliable coordinates).
DANTE_LTR_COMPLETE_RANKS <- c("DL", "DLT", "DLP", "DLTP")

load_dante_ltr_stats <- function(outdir) {
  gff <- file.path(outdir, "DANTE_LTR", "DANTE_LTR.gff3")
  if (!file.exists(gff) || file.size(gff) == 0) return(NULL)
  gr <- tryCatch(import(gff), error = function(e) NULL)
  if (is.null(gr) || length(gr) == 0) return(NULL)
  te <- gr[gr$type == "transposable_element"]
  if (length(te) == 0) return(NULL)
  if (!is.null(te$Rank)) {
    te <- te[as.character(te$Rank) %in% DANTE_LTR_COMPLETE_RANKS]
  }
  if (length(te) == 0) return(NULL)
  cls <- te$Final_Classification
  if (is.null(cls)) cls <- te$Name
  cls_path <- canonicalise(cls, source = "DANTE_LTR", validate = FALSE)
  tmp_df <- data.frame(path = cls_path, count = 1L, bp = as.numeric(width(te)))
  agg <- aggregate(cbind(count, bp) ~ path, data = tmp_df, FUN = sum)
  agg
}

# ── B4. DANTE_TIR stats from summary txt ──────────────────────────────────
load_dante_tir_stats <- function(outdir) {
  f <- file.path(outdir, "DANTE_TIR", "TIR_classification_summary.txt")
  if (!file.exists(f) || file.size(f) == 0) return(NULL)
  df <- read.table(f, header = TRUE, sep = "\t", check.names = FALSE)
  colnames(df) <- c("raw_name", "count")
  df$path <- canonicalise(df$raw_name, source = "DANTE_TIR", validate = FALSE)
  df[, c("path", "count")]
}

# ── B5. DANTE_LINE stats from GFF3 ────────────────────────────────────────
load_dante_line_stats <- function(outdir) {
  gff <- file.path(outdir, "DANTE_LINE", "DANTE_LINE.gff3")
  if (!file.exists(gff) || file.size(gff) == 0) return(list(regions = 0, bp = 0))
  gr <- tryCatch(import(gff), error = function(e) NULL)
  if (is.null(gr) || length(gr) == 0) return(list(regions = 0, bp = 0))
  le <- gr[gr$type == "LINE_element"]
  if (length(le) == 0) return(list(regions = 0, bp = 0))
  list(regions = length(le), bp = sum(as.numeric(width(le))))
}

# ── B6. Discover BigWig files ──────────────────────────────────────────────
discover_bw_files <- function(outdir, bin_width) {
  suffix <- if (bin_width == 100000) "_100k.bw" else "_10k.bw"
  bw_dir_rm <- file.path(outdir, "Repeat_Annotation_NoSat_split_by_class_bigwig",
                          if (bin_width == 100000) "100k" else "10k")
  bw_dir_tc <- file.path(outdir, "TideCluster", "default",
                          "TideCluster_clustering_split_files_bigwig",
                          if (bin_width == 100000) "100k" else "10k")
  tc_agg    <- file.path(outdir, "TideCluster", "default",
                          paste0("TideCluster_clustering", suffix))

  find_bw <- function(dir, pattern) {
    if (!dir.exists(dir)) return(character(0))
    list.files(dir, pattern = pattern, full.names = TRUE)
  }

  rm_files <- find_bw(bw_dir_rm, "[.]bw$")
  names(rm_files) <- sub(paste0(suffix, "$"), "", basename(rm_files))

  tc_cluster_files <- find_bw(bw_dir_tc, "[.]bw$")
  names(tc_cluster_files) <- sub(paste0(suffix, "$"), "", basename(tc_cluster_files))

  list(
    rm         = rm_files,
    tc_agg     = if (file.exists(tc_agg)) tc_agg else NULL,
    tc_cluster = tc_cluster_files
  )
}

# ═══════════════════════════════════════════════════════════════════════════
# C. DATA PROCESSING
# ═══════════════════════════════════════════════════════════════════════════

# ── C1. Apply sequence length thresholds ──────────────────────────────────
apply_seq_thresholds <- function(genome_info, min_len_chart, min_len_tracks, max_tracks) {
  gi <- genome_info[order(-genome_info$length), ]
  list(
    chart_seqs  = gi$seqname[gi$length >= min_len_chart],
    track_seqs  = head(gi$seqname[gi$length >= min_len_tracks], max_tracks),
    n_other     = sum(gi$length < min_len_chart),
    other_bp    = sum(gi$length[gi$length < min_len_chart])
  )
}

# ── C2. Build sunburst hierarchy from composition CSV ─────────────────────
build_sunburst_data <- function(comp) {
  pal <- c(
    "Class_I"       = "#4575b4",
    "Class_II"      = "#d73027",
    "Tandem_repeats" = "#1a9850",
    "Simple_repeat" = "#878787",
    "Low_complexity"= "#bdbdbd",
    "rDNA"          = "#762a83",
    "rDNA_45S"      = "#762a83",
    "rDNA_5S"       = "#9e6faf",
    "Unknown"       = "#d9d9d9"
  )
  get_color <- function(id) {
    top <- strsplit(id, "/")[[1]][1]
    col <- pal[top]
    if (is.na(col)) "#aaaaaa" else col
  }

  csv_ids    <- comp$type
  csv_bp     <- comp$bp
  csv_pct    <- comp$pct

  # Compute the closure: all parent paths needed
  all_prefixes <- unique(unlist(lapply(csv_ids, function(id) {
    parts <- strsplit(id, "/")[[1]]
    if (length(parts) <= 1) return(character(0))
    sapply(seq_len(length(parts) - 1), function(n) paste(parts[1:n], collapse = "/"))
  })))
  synthetic_ids <- setdiff(all_prefixes, csv_ids)

  # Combine CSV nodes + synthetic nodes
  all_ids <- c(csv_ids, synthetic_ids)
  all_bp  <- c(csv_bp,  rep(0, length(synthetic_ids)))
  all_pct <- c(csv_pct, rep(0, length(synthetic_ids)))

  labels  <- sapply(strsplit(all_ids, "/"), function(x) x[length(x)])
  parents <- sapply(all_ids, function(id) {
    p <- strsplit(id, "/")[[1]]
    if (length(p) == 1) "root" else paste(p[-length(p)], collapse = "/")
  })
  colors  <- sapply(all_ids, get_color)

  list(
    ids     = c("root", all_ids),
    labels  = c("All repeats", labels),
    parents = c("", parents),
    values  = c(0, all_bp),
    text    = c("", ifelse(all_pct > 0, sprintf("%.3f%%", all_pct), "")),
    colors  = c("#cccccc", colors)
  )
}

# ── C3. Build composition tree for hierarchical table ─────────────────────
build_comp_tree <- function(comp, ltr_stats, tir_stats, line_stats = NULL,
                            genome_size = NULL) {
  # Compute the closure of all node paths (CSV rows + synthetic parents)
  csv_ids  <- comp$type
  csv_bp   <- setNames(comp$bp,  comp$type)
  csv_pct  <- setNames(comp$pct, comp$type)

  all_prefixes <- unique(unlist(lapply(csv_ids, function(id) {
    parts <- strsplit(id, "/")[[1]]
    if (length(parts) <= 1) return(character(0))
    sapply(seq_len(length(parts) - 1), function(n) paste(parts[1:n], collapse = "/"))
  })))
  synthetic_ids <- setdiff(all_prefixes, csv_ids)
  all_ids <- c(csv_ids, synthetic_ids)

  # Parent map
  parent_of <- sapply(all_ids, function(id) {
    p <- strsplit(id, "/")[[1]]
    if (length(p) == 1) NA_character_ else paste(p[-length(p)], collapse = "/")
  })

  # Children map
  children_of <- lapply(setNames(all_ids, all_ids), function(id) {
    all_ids[!is.na(parent_of) & parent_of == id]
  })

  # Subtree bp sum (recursive)
  subtree_bp <- function(id) {
    own <- unname(csv_bp[id] %||% 0)  # unname() prevents sapply name-collision
    kids <- children_of[[id]]
    if (length(kids) == 0) return(own)
    own + sum(sapply(kids, subtree_bp))
  }
  subtree_bp_cache <- sapply(all_ids, subtree_bp)

  if (is.null(genome_size)) genome_size <- sum(as.numeric(comp$bp))  # fallback
  # genome_size is the total assembly size, NOT the repeat content

  # Build DANTE lookup: path → count
  dante_counts <- setNames(integer(0), character(0))
  if (!is.null(ltr_stats) && nrow(ltr_stats) > 0) {
    for (i in seq_len(nrow(ltr_stats))) {
      p <- ltr_stats$path[i]
      dante_counts[p] <- (dante_counts[p] %||% 0L) + ltr_stats$count[i]
    }
  }
  if (!is.null(tir_stats) && nrow(tir_stats) > 0) {
    for (i in seq_len(nrow(tir_stats))) {
      p <- tir_stats$path[i]
      dante_counts[p] <- (dante_counts[p] %||% 0L) + tir_stats$count[i]
    }
  }
  # LINE elements: count at the Class_I/LINE node if it exists
  if (!is.null(line_stats) && !is.null(line_stats$regions) && line_stats$regions > 0) {
    line_path <- grep("LINE$", all_ids, value = TRUE)
    if (length(line_path) > 0)
      dante_counts[line_path[1]] <- as.integer(line_stats$regions)
  }

  # DFS pre-order traversal
  rows <- list()
  dfs <- function(id, depth) {
    own_bp   <- unname(csv_bp[id]  %||% 0)
    tot_bp   <- unname(subtree_bp_cache[id])
    kids     <- children_of[[id]]
    has_kids <- length(kids) > 0
    label    <- strsplit(id, "/")[[1]]
    label    <- label[length(label)]
    dc_raw   <- dante_counts[id]
    dc       <- if (length(dc_raw) > 0 && !is.na(dc_raw)) unname(dc_raw) else NA_integer_

    if (has_kids) {
      # Total row
      rows[[length(rows) + 1]] <<- list(
        path        = id,
        label       = label,
        row_type    = "total",
        depth       = depth,
        bp          = tot_bp,
        pct         = tot_bp / genome_size * 100,
        dante_count = NA_integer_
      )
      # Unspecified row (own CSV value) only if non-zero
      if (own_bp > 0) {
        rows[[length(rows) + 1]] <<- list(
          path        = id,
          label       = label,
          row_type    = "unspecified",
          depth       = depth + 1L,
          bp          = own_bp,
          pct         = unname(csv_pct[id] %||% 0),
          dante_count = NA_integer_
        )
      }
      # Recurse into children sorted by subtree bp descending
      kids_sorted <- kids[order(-subtree_bp_cache[kids])]
      for (k in kids_sorted) dfs(k, depth + 1L)
    } else {
      # Leaf row — show DANTE count if available
      rows[[length(rows) + 1]] <<- list(
        path        = id,
        label       = label,
        row_type    = "leaf",
        depth       = depth,
        bp          = own_bp,
        pct         = unname(csv_pct[id] %||% 0),
        dante_count = dc
      )
    }
  }

  # Find top-level nodes (no parent in all_ids)
  top_nodes <- all_ids[is.na(parent_of)]
  top_nodes <- top_nodes[order(-subtree_bp_cache[top_nodes])]
  for (id in top_nodes) dfs(id, 0L)

  do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
}

# ── C4. Bin a BigWig file for a set of sequences ──────────────────────────
bin_bw <- function(bw_path, seqnames, seq_lengths, bin_width) {
  if (is.null(bw_path) || !file.exists(bw_path)) return(NULL)
  bw <- tryCatch(import(bw_path, as = "RleList"), error = function(e) NULL)
  if (is.null(bw)) return(NULL)

  result <- vector("list", length(seqnames))
  names(result) <- seqnames
  for (s in seqnames) {
    if (!s %in% names(bw)) {
      result[[s]] <- rep(0, ceiling(seq_lengths[s] / bin_width))
      next
    }
    rle_v  <- bw[[s]]
    len    <- seq_lengths[s]
    starts <- seq(1L, len, by = bin_width)
    ends   <- pmin(starts + bin_width - 1L, len)
    v      <- Views(rle_v, start = starts, end = ends)
    result[[s]] <- pmax(0, viewMeans(v, na.rm = TRUE))
  }
  result
}

# Sum multiple binned BW arrays element-wise
sum_bw_arrays <- function(arrays_list, seqnames) {
  arrays_list <- Filter(Negate(is.null), arrays_list)
  if (length(arrays_list) == 0) return(NULL)
  Reduce(function(a, b) {
    mapply(function(x, y) pmin(1, x + y), a, b, SIMPLIFY = FALSE)
  }, arrays_list)
}

# ── Helpers for concatenated density track plots ───────────────────────────

# Centred rolling mean — NA at both ends (creates gaps at sequence boundaries)
smooth_vec <- function(x, N) {
  if (N <= 1L || length(x) < N) return(as.numeric(x))
  padded <- c(rep(0, N - 1L), x, rep(0, N - 1L))
  sm     <- as.numeric(stats::filter(padded, rep(1/N, N), sides = 2))
  sm     <- sm[N:(N + length(x) - 1L)]
  sm[1:(N - 1L)] <- NA
  sm[(length(sm) - N + 2L):length(sm)] <- NA
  sm
}

# Ramer-Douglas-Peucker simplification for a monotonic-x time series.
# NA values in y are always preserved (they mark visual segment breaks).
# Returns integer indices of points to keep.
rdp_simplify <- function(x, y, epsilon) {
  n <- length(x)
  if (n <= 2L || epsilon <= 0) return(seq_len(n))
  keep <- logical(n)

  rdp_seg <- function(a, b) {
    keep[a] <<- TRUE
    keep[b] <<- TRUE
    if (b <= a + 1L) return(invisible(NULL))
    xi    <- x[(a + 1L):(b - 1L)]
    yi    <- y[(a + 1L):(b - 1L)]
    frac  <- (xi - x[a]) / (x[b] - x[a])
    dists <- abs(yi - (y[a] + frac * (y[b] - y[a])))
    dists[is.na(dists)] <- Inf   # NA interior points are always splitting candidates
    j <- which.max(dists)
    if (dists[j] > epsilon) {
      mid <- a + j
      rdp_seg(a, mid)
      rdp_seg(mid, b)
    }
  }

  is_na <- is.na(y)
  keep[is_na] <- TRUE   # preserve NA gap markers
  non_na <- which(!is_na)
  if (length(non_na) > 0) {
    # identify start/end of each run of non-NA values
    brk <- c(TRUE, diff(non_na) > 1L)
    seg_s <- non_na[brk]
    seg_e <- non_na[c(brk[-1L], TRUE)]
    for (k in seq_along(seg_s))
      rdp_seg(seg_s[k], seg_e[k])
  }
  sort(which(keep))
}

# Find lineage-level BigWig files in a repeat-annotation directory
discover_lineage_bw_files <- function(rm_dir, bin_width) {
  suffix <- if (bin_width == 100000L) "_100k.bw" else "_10k.bw"
  known  <- c("Ale","Alesia","Angela","Bianca","Bryco","Gymco-I","Gymco-II","Gymco-III",
              "Gymco-IV","Ikeros","Ivana","Lyco","Osser","SIRE","TAR","Tork","Chlamyvir",
              "chromo-unclass","CRM","Galadriel","Reina","Tcn1","Tekay","Athila","Ogre",
              "Retand","TatI","TatII","TatIII","Phygy","Selgy")
  if (!dir.exists(rm_dir)) return(list())
  pat    <- paste0("LTR\\.Ty.*", gsub("\\.", "\\\\.", suffix), "$")
  fnames <- list.files(rm_dir, pattern = pat, full.names = FALSE)
  if (length(fnames) == 0) return(list())
  lnames <- gsub("_100k\\.bw$|_10k\\.bw$", "", fnames)
  lnames <- gsub(".*\\.", "", lnames)
  keep   <- lnames %in% known
  fnames <- fnames[keep]; lnames <- lnames[keep]
  if (length(fnames) == 0) return(list())
  ord    <- order(match(lnames, known))
  setNames(as.list(file.path(rm_dir, fnames[ord])), lnames[ord])
}

# Find TRC satellite BigWig files and label them with monomer sizes
discover_trc_bw_files <- function(outdir, bin_width) {
  res     <- if (bin_width == 100000L) "100k" else "10k"
  suffix  <- paste0("_", res, ".bw")
  trc_dir <- file.path(outdir, "TideCluster", "default",
                        "TideCluster_clustering_split_files_bigwig", res)
  mon_f   <- file.path(outdir, "TideCluster", "default", "TideCluster_kite",
                        "monomer_size_best_estimate_stat.csv")
  if (!dir.exists(trc_dir)) return(list())
  bw_f   <- list.files(trc_dir, pattern = "\\.bw$", full.names = FALSE)
  tnames <- gsub(paste0(gsub("\\.", "\\\\.", suffix), "$"), "", bw_f)
  tidx   <- suppressWarnings(as.numeric(sub("TRC_", "", tnames)))
  valid  <- !is.na(tidx)
  bw_f   <- bw_f[valid]; tnames <- tnames[valid]; tidx <- tidx[valid]
  if (length(bw_f) == 0) return(list())
  ord    <- order(tidx)
  N      <- min(20L, length(bw_f))
  bw_f   <- bw_f[ord][seq_len(N)]; tnames <- tnames[ord][seq_len(N)]
  mon_sizes <- setNames(rep("?", N), tnames)
  if (file.exists(mon_f)) {
    ms <- tryCatch(read.table(mon_f, header = TRUE, sep = "\t", check.names = FALSE),
                   error = function(e) NULL)
    if (!is.null(ms) && all(c("TRC_ID", "position") %in% names(ms))) {
      for (tn in tnames) {
        sub <- ms[ms$TRC_ID == tn, ]
        if (nrow(sub) > 0)
          mon_sizes[tn] <- names(sort(table(sub$position), decreasing = TRUE)[1])
      }
    }
  }
  setNames(as.list(file.path(trc_dir, bw_f)),
           paste0(tnames, " (", mon_sizes, "bp)"))
}

# ── C5. Build density track arrays for selected sequences ─────────────────
build_density_arrays <- function(bw_files, seqnames, seq_lengths_vec, bin_width) {
  result <- list()
  for (label in names(bw_files)) {
    paths <- bw_files[[label]]
    if (is.null(paths)) next
    if (length(paths) == 1) {
      arr <- bin_bw(paths, seqnames, seq_lengths_vec, bin_width)
    } else {
      arr_list <- lapply(paths, bin_bw, seqnames = seqnames,
                         seq_lengths = seq_lengths_vec, bin_width = bin_width)
      arr <- sum_bw_arrays(arr_list, seqnames)
    }
    if (!is.null(arr)) result[[label]] <- arr
  }
  result
}

# ── C6. Per-sequence coverage fractions for stacked bar ───────────────────
compute_per_seq_coverage <- function(bw_files_named, seqnames, seq_lengths_vec) {
  mat <- matrix(0, nrow = length(seqnames), ncol = length(bw_files_named),
                dimnames = list(seqnames, names(bw_files_named)))

  for (cat in names(bw_files_named)) {
    paths <- bw_files_named[[cat]]
    paths <- paths[!is.null(paths) & file.exists(paths)]
    if (length(paths) == 0) next
    combined <- rep(0, length(seqnames))
    names(combined) <- seqnames
    for (p in paths) {
      gr <- tryCatch(import(p), error = function(e) NULL)
      if (is.null(gr) || length(gr) == 0) next
      by_seq <- tapply(score(gr) * width(gr), as.character(seqnames(gr)), sum, na.rm = TRUE)
      for (s in intersect(names(by_seq), seqnames)) {
        combined[s] <- combined[s] + by_seq[s]
      }
    }
    mat[, cat] <- pmin(1, combined / seq_lengths_vec[seqnames])
  }
  mat
}

# ── C7. Satellite cluster coverage ────────────────────────────────────────
build_satellite_coverage <- function(tc_cluster_files, seqnames, seq_lengths_vec,
                                      bin_width, top_n) {
  if (length(tc_cluster_files) == 0) return(list(totals = NULL, arrays = NULL))

  totals <- vapply(names(tc_cluster_files), function(cl) {
    gr <- tryCatch(import(tc_cluster_files[[cl]]), error = function(e) NULL)
    if (is.null(gr) || length(gr) == 0) return(0)
    sum(score(gr) * width(gr), na.rm = TRUE)
  }, numeric(1))

  top_clusters <- names(sort(totals, decreasing = TRUE))[seq_len(min(top_n, length(totals)))]

  per_seq_mat <- matrix(0, nrow = length(seqnames), ncol = length(top_clusters),
                        dimnames = list(seqnames, top_clusters))
  for (cl in top_clusters) {
    gr <- tryCatch(import(tc_cluster_files[[cl]]), error = function(e) NULL)
    if (is.null(gr)) next
    by_seq <- tapply(score(gr) * width(gr), as.character(seqnames(gr)), sum, na.rm = TRUE)
    for (s in intersect(names(by_seq), seqnames)) {
      per_seq_mat[s, cl] <- by_seq[s] / seq_lengths_vec[s]
    }
  }

  density_arrays <- build_density_arrays(
    setNames(as.list(tc_cluster_files[top_clusters]), top_clusters),
    seqnames, seq_lengths_vec, bin_width
  )

  list(totals = totals[top_clusters], per_seq = per_seq_mat, arrays = density_arrays)
}

# ═══════════════════════════════════════════════════════════════════════════
# D. JSON / PLOTLY GENERATION
# ═══════════════════════════════════════════════════════════════════════════

J <- function(x) toJSON(x, auto_unbox = TRUE, null = "null")

TRACK_COLORS <- c(
  "Tandem_repeats"  = "#1a9850",
  "Ty1/copia"       = "#4575b4",
  "Ty3/gypsy"       = "#313695",
  "LINE"            = "#74add1",
  "TIR (Class II)"  = "#d73027",
  "Pararetrovirus"  = "#f46d43",
  "Simple repeats"  = "#878787",
  "Low complexity"  = "#bdbdbd",
  "rDNA"            = "#762a83",
  "Unknown"         = "#d9d9d9",
  "Other"           = "#ffffbf"
)

# ── D1. Sunburst JSON ──────────────────────────────────────────────────────
json_sunburst <- function(sb) {
  trace <- list(
    type          = "sunburst",
    ids           = sb$ids,
    labels        = sb$labels,
    parents       = sb$parents,
    values        = sb$values,
    branchvalues  = "remainder",
    hovertemplate = paste0(
      "<b>%{label}</b><br>",
      "%{value:,.0f} bp<br>",
      "%{customdata}<br>",
      "<extra></extra>"
    ),
    customdata    = sb$text,
    textinfo      = "label+percent entry",
    marker        = list(colors = sb$colors),
    maxdepth      = 3
  )
  layout <- list(
    margin = list(t = 10, l = 10, r = 10, b = 10),
    height = 520
  )
  list(traces = list(trace), layout = layout)
}

# ── D2. Per-sequence stacked composition bar ──────────────────────────────
json_composition_bar <- function(cov_mat, seq_info, n_other, other_total_bp,
                                  genome_avg_frac) {
  cats      <- colnames(cov_mat)
  seq_names <- rownames(cov_mat)
  total_rep <- rowSums(cov_mat)
  ord       <- order(-total_rep)
  cov_mat   <- cov_mat[ord, , drop = FALSE]
  seq_names <- seq_names[ord]

  seq_len_mb <- seq_info$length[match(seq_names, seq_info$seqname)] / 1e6
  y_labels   <- sprintf("%s (%.1f Mb)", seq_names, seq_len_mb)

  max_x <- max(total_rep, genome_avg_frac, 0.05) * 1.08

  traces <- lapply(seq_along(cats), function(i) {
    cat <- cats[i]
    list(
      type        = "bar",
      orientation = "h",
      name        = cat,
      x           = cov_mat[, cat],
      y           = y_labels,
      marker      = list(color = TRACK_COLORS[[cat]] %||% "#aaaaaa"),
      hovertemplate = paste0(cat, ": %{x:.3f}<extra></extra>")
    )
  })

  if (n_other > 0) {
    traces[[length(traces) + 1]] <- list(
      type        = "bar",
      orientation = "h",
      name        = sprintf("Other sequences (%d)", n_other),
      x           = list(0),
      y           = list(sprintf("Other (%d seqs, %.0f Mb)",
                                  n_other, other_total_bp / 1e6)),
      marker      = list(color = "#eeeeee"),
      hovertemplate = "Sequences below length threshold<extra></extra>"
    )
  }

  # Genome-average vertical line as a shape
  shapes <- list(list(
    type = "line",
    x0 = genome_avg_frac, x1 = genome_avg_frac,
    y0 = 0, y1 = 1,
    yref = "paper",
    line = list(color = "#333333", width = 1.5, dash = "dash")
  ))

  annotations <- list(list(
    x = genome_avg_frac, y = 1,
    xref = "x", yref = "paper",
    text = sprintf("genome avg %.1f%%", genome_avg_frac * 100),
    showarrow = FALSE,
    xanchor = "left",
    font = list(size = 10, color = "#333333")
  ))

  layout <- list(
    barmode     = "stack",
    height      = max(300, length(seq_names) * 22 + 80),
    margin      = list(l = 200, r = 20, t = 30, b = 60),
    xaxis       = list(title = "Fraction of sequence", range = c(0, max_x)),
    yaxis       = list(title = "", automargin = TRUE),
    legend      = list(orientation = "h", y = -0.15),
    showlegend  = TRUE,
    shapes      = shapes,
    annotations = annotations
  )
  list(traces = traces, layout = layout)
}

# ── D3. Concatenated density track plot (matches summary_plots.pdf design) ─
# All sequences concatenated on a single x-axis; tracks stacked vertically
# with a SHARED y-scale (global ymax across all tracks).
# Two lines per track: thin gray (raw) + thicker red (smoothed).
# Track labels on the right; chromosome names on x-axis.
json_concat_density_plot <- function(bw_map, seq_info, bin_width, n_smooth = 10) {
  bw_map <- Filter(function(p) !is.null(p) && file.exists(unlist(p)), bw_map)
  if (length(bw_map) == 0) return(NULL)

  track_names <- names(bw_map)
  n_tracks    <- length(track_names)

  # All sequences sorted by length descending (matches PDF order)
  seq_info <- seq_info[order(-seq_info$length), ]
  seqnames <- seq_info$seqname
  seq_len  <- setNames(seq_info$length, seqnames)

  # Cumulative x offsets (Mb): cum_mb[i] = start of sequence i
  cum_mb  <- c(0, cumsum(seq_len / 1e6))
  x_ticks <- unname((cum_mb[-length(cum_mb)] + cum_mb[-1]) / 2)  # midpoints
  x_sep   <- unname(cum_mb[-c(1, length(cum_mb))])         # internal boundaries
  total_mb <- cum_mb[length(cum_mb)]

  # Build raw + smoothed arrays per track, concatenated across sequences
  x_all       <- numeric(0)
  track_y_raw <- list()
  track_y_sm  <- list()
  built_x     <- FALSE

  for (lbl in track_names) {
    p   <- unlist(bw_map[[lbl]])
    bw  <- tryCatch(import(p, as = "RleList"), error = function(e) NULL)
    raw_cat <- numeric(0)
    sm_cat  <- numeric(0)

    for (si in seq_along(seqnames)) {
      s      <- seqnames[si]
      len    <- seq_len[s]
      starts <- seq(1L, len, by = bin_width)
      ends   <- pmin(starts + bin_width - 1L, len)
      n_bins <- length(starts)

      if (!is.null(bw) && s %in% names(bw)) {
        v   <- Views(bw[[s]], start = starts, end = ends)
        raw <- pmax(0, viewMeans(v, na.rm = TRUE))
      } else {
        raw <- rep(0, n_bins)
      }
      sm      <- smooth_vec(raw, n_smooth)
      raw_cat <- c(raw_cat, raw)
      sm_cat  <- c(sm_cat,  sm)
      if (!built_x)
        x_all <- c(x_all, round((starts + ends) / 2 / 1e6 + cum_mb[si], 3))
    }
    built_x            <- TRUE
    track_y_raw[[lbl]] <- raw_cat
    track_y_sm[[lbl]]  <- sm_cat
  }

  # Global ymax shared across all tracks
  ymax <- max(sapply(track_y_raw, max, na.rm = TRUE), na.rm = TRUE)
  if (!is.finite(ymax) || ymax == 0) ymax <- 1
  y_step <- 1.1 * ymax

  # RDP epsilon relative to the global y range
  eps_raw <- ymax * 0.01    # 1 % of range — aggressive (raw is decorative)
  eps_sm  <- ymax * 0.002   # 0.2 % of range — conservative (smooth is informative)

  # Build traces: for each track, raw (gray thin) + smooth (red thick)
  traces <- list()
  tidx   <- 0L
  for (ti in seq_along(track_names)) {
    lbl   <- track_names[ti]
    y_off <- (ti - 1) * y_step

    raw_y  <- track_y_raw[[lbl]] + y_off
    idx_rw <- rdp_simplify(x_all, raw_y, eps_raw)
    # Raw gray line
    tidx <- tidx + 1L
    traces[[tidx]] <- list(
      type = "scattergl", mode = "lines",
      x = x_all[idx_rw], y = raw_y[idx_rw],
      showlegend = FALSE, hoverinfo = "skip",
      line = list(width = 0.5, color = "#00000030"),
      name = paste0(lbl, "_raw")
    )

    sm_y   <- track_y_sm[[lbl]] + y_off
    idx_sm <- rdp_simplify(x_all, sm_y, eps_sm)
    # Smoothed red line (NA values create gaps at sequence boundaries)
    tidx <- tidx + 1L
    traces[[tidx]] <- list(
      type = "scattergl", mode = "lines",
      x = x_all[idx_sm], y = sm_y[idx_sm],
      showlegend = FALSE, connectgaps = FALSE,
      name = lbl,
      line = list(width = 1.5, color = "#cc0000bb"),
      hovertemplate = paste0(lbl, ": %{y:.4f} at %{x:.2f} Mb<extra></extra>")
    )
  }

  # Shapes drawn BELOW all traces (layer="below") so they don't obscure density
  # lines even on short contigs.
  # Horizontal grid lines at every track boundary
  h_shapes <- lapply(seq(0L, n_tracks), function(ti)
    list(type = "line", layer = "below",
         x0 = 0, x1 = 1, xref = "paper",
         y0 = ti * y_step, y1 = ti * y_step, yref = "y",
         line = list(color = "#bbbbbb", width = 0.8)))

  # Vertical chromosome boundary lines — clearly visible grid
  v_shapes <- lapply(x_sep, function(xb)
    list(type = "line", layer = "below",
         x0 = xb, x1 = xb, xref = "x",
         y0 = 0, y1 = 1, yref = "paper",
         line = list(color = "#999999", width = 1.0)))

  # Track label annotations (right side, in data y coordinates)
  annots <- lapply(seq_along(track_names), function(ti) {
    list(
      x = 1.01, xref = "paper", xanchor = "left",
      y = (ti - 1) * y_step + ymax * 0.5, yref = "y",
      text = track_names[ti], showarrow = FALSE,
      font = list(size = 9, color = "#333333"),
      bgcolor = "rgba(255,255,255,0.85)", borderpad = 2
    )
  })

  layout <- list(
    height   = max(200L, n_tracks * 55L + 90L),
    margin   = list(l = 40, r = 140, t = 10, b = 95),
    xaxis    = list(
      range     = c(0, total_mb),
      tickvals  = x_ticks,
      ticktext  = seqnames,
      tickangle = -90,
      tickfont  = list(size = 7),
      title     = "Genome position (Mb)",
      showgrid  = FALSE
    ),
    yaxis    = list(
      range          = c(-0.02 * ymax, n_tracks * y_step + 0.05 * ymax),
      showticklabels = FALSE,
      zeroline = FALSE, showgrid = FALSE, title = ""
    ),
    showlegend  = FALSE,
    annotations = annots,
    shapes      = unname(c(h_shapes, v_shapes))
  )

  list(traces = traces, layout = layout)
}

# ── D4. Satellite cluster stacked bar ─────────────────────────────────────
json_satellite_bar <- function(sat_data, seq_info) {
  if (is.null(sat_data$per_seq) || nrow(sat_data$per_seq) == 0) return(NULL)
  mat       <- sat_data$per_seq
  seq_names <- rownames(mat)
  clusters  <- colnames(mat)
  ord       <- order(-rowSums(mat))
  mat       <- mat[ord, , drop = FALSE]
  seq_names <- seq_names[ord]

  seq_len_mb <- seq_info$length[match(seq_names, seq_info$seqname)] / 1e6
  y_labels   <- sprintf("%s (%.1f Mb)", seq_names, seq_len_mb)

  cluster_pal <- colorRampPalette(c("#1a9850","#66bd63","#a6d96a","#d9ef8b",
                                     "#fee08b","#fdae61","#f46d43","#d73027"))(length(clusters))

  traces <- lapply(seq_along(clusters), function(i) {
    list(
      type        = "bar",
      orientation = "h",
      name        = clusters[i],
      x           = mat[, clusters[i]],
      y           = y_labels,
      marker      = list(color = cluster_pal[i]),
      hovertemplate = paste0(clusters[i], ": %{x:.4f}<extra></extra>")
    )
  })
  layout <- list(
    barmode   = "stack",
    height    = max(280, length(seq_names) * 22 + 80),
    margin    = list(l = 200, r = 20, t = 30, b = 60),
    xaxis     = list(title = "Fraction of sequence"),
    yaxis     = list(title = "", automargin = TRUE),
    legend    = list(orientation = "h", y = -0.18),
    showlegend = TRUE
  )
  list(traces = traces, layout = layout)
}

# ═══════════════════════════════════════════════════════════════════════════
# E. HTML ASSEMBLY
# ═══════════════════════════════════════════════════════════════════════════

# ── E1. Load Plotly.js ─────────────────────────────────────────────────────
load_plotly_js <- function() {
  script_dir   <- tryCatch(normalizePath(dirname(sys.frame(1)$ofile), mustWork = FALSE),
                            error = function(e) ".")
  pipeline_dir <- dirname(script_dir)
  candidates   <- c(
    file.path(pipeline_dir, "data", "plotly.min.js"),
    file.path("/opt/pipeline/data", "plotly.min.js"),
    file.path(getwd(), "data", "plotly.min.js")
  )
  local_js <- Find(file.exists, candidates)
  if (!is.null(local_js)) {
    message("Using local Plotly.js: ", local_js)
    return(paste(readLines(local_js, warn = FALSE), collapse = "\n"))
  }
  url  <- "https://cdn.plot.ly/plotly-2.35.2.min.js"
  dest <- tempfile(fileext = ".js")
  ok   <- tryCatch({ download.file(url, dest, quiet = TRUE, method = "auto"); TRUE },
                    error = function(e) FALSE)
  if (ok && file.exists(dest)) {
    message("Downloaded Plotly.js from CDN (consider caching in data/plotly.min.js)")
    return(paste(readLines(dest, warn = FALSE), collapse = "\n"))
  }
  warning("Plotly.js not found locally and could not be downloaded.")
  sprintf('<script src="%s"></script>', url)
}

# ── E2. Plotly div + init script ──────────────────────────────────────────
plotly_div <- function(div_id, chart) {
  traces_json <- J(chart$traces)
  layout_json <- J(chart$layout)
  sprintf(
    '<div id="%s" style="width:100%%;"></div>\n<script>\nPlotly.newPlot("%s", %s, %s, {responsive:true, displayModeBar:true});\n</script>',
    div_id, div_id, traces_json, layout_json
  )
}

# ── E3. Summary cards ──────────────────────────────────────────────────────
html_cards <- function(genome_info, comp, ltr_stats, tir_stats, line_stats,
                        seq_thresholds) {
  genome_size_mb    <- sum(as.numeric(genome_info$length)) / 1e6
  total_repeat_pct  <- sum(comp$bp) / sum(as.numeric(genome_info$length)) * 100
  n_seqs  <- nrow(genome_info)
  n_shown <- length(seq_thresholds$chart_seqs)

  ltr_count  <- if (!is.null(ltr_stats) && nrow(ltr_stats) > 0) sum(ltr_stats$count) else "N/A"
  tir_count  <- if (!is.null(tir_stats) && nrow(tir_stats) > 0) sum(tir_stats$count) else "N/A"
  line_count <- line_stats$regions %||% "N/A"

  sprintf('
<div class="cards">
  <div class="card">
    <div class="card-label">Genome size</div>
    <div class="card-value">%.0f Mb</div>
    <div class="card-sub">%d sequences</div>
  </div>
  <div class="card">
    <div class="card-label">Repeat content</div>
    <div class="card-value">%.1f%%</div>
    <div class="card-sub">similarity-based</div>
  </div>
  <div class="card">
    <div class="card-label">LTR retrotransposons</div>
    <div class="card-value">%s</div>
    <div class="card-sub">structure-based (DANTE_LTR)</div>
  </div>
  <div class="card">
    <div class="card-label">TIR elements</div>
    <div class="card-value">%s</div>
    <div class="card-sub">structure-based (DANTE_TIR)</div>
  </div>
  <div class="card">
    <div class="card-label">LINE elements</div>
    <div class="card-value">%s</div>
    <div class="card-sub">structure-based (DANTE_LINE)</div>
  </div>
  <div class="card">
    <div class="card-label">Sequences in charts</div>
    <div class="card-value">%d / %d</div>
    <div class="card-sub">&ge; %.0f kb shown individually</div>
  </div>
</div>',
    genome_size_mb, n_seqs,
    total_repeat_pct,
    ltr_count, tir_count, line_count,
    n_shown, n_seqs, opt$min_len_chart / 1000
  )
}

# ── E4. Hierarchical composition table ────────────────────────────────────
html_comp_table <- function(tree_df) {
  fmt_bp <- function(x) format(round(as.numeric(x)), big.mark = ",", scientific = FALSE)

  rows_html <- apply(tree_df, 1, function(r) {
    depth   <- as.integer(r["depth"])
    indent  <- depth * 20
    bp      <- as.numeric(r["bp"])
    pct     <- as.numeric(r["pct"])
    rtype   <- r["row_type"]
    label   <- r["label"]
    dc      <- r["dante_count"]

    # Display label
    if (rtype == "total") {
      disp <- sprintf('<b>%s</b> <span class="row-tag row-total">Total</span>', label)
    } else if (rtype == "unspecified") {
      disp <- sprintf('<i>%s</i> <span class="row-tag row-unspec">Unspecified</span>', label)
    } else {
      disp <- label
    }

    dante_cell <- if (!is.na(dc) && dc != "NA") {
      sprintf('<td style="text-align:right">%s</td>', format(as.integer(dc), big.mark = ","))
    } else {
      '<td></td>'
    }

    row_class <- switch(rtype, total = "tr-total", unspecified = "tr-unspec", "")
    sprintf('<tr class="%s"><td style="padding-left:%dpx">%s</td><td style="text-align:right">%s</td><td style="text-align:right">%.3f%%</td>%s</tr>',
            row_class, indent, disp, fmt_bp(bp), pct, dante_cell)
  })

  sprintf('
<div style="overflow-x:auto">
<table id="comp-table" class="data-table">
<thead>
<tr>
  <th>Classification</th>
  <th style="text-align:right">Total bp</th>
  <th style="text-align:right">%% genome</th>
  <th style="text-align:right">Complete TEs</th>
</tr>
</thead>
<tbody>
%s
</tbody>
</table>
</div>
<p class="caption">Complete TEs: structure-based count from DANTE_LTR / DANTE_TIR where classification matches.</p>',
    paste(rows_html, collapse = "\n"))
}

# ── E5. DANTE structure-based summary table ───────────────────────────────
html_dante_summary <- function(ltr_stats, tir_stats, line_stats, outdir) {
  ltr_count <- if (!is.null(ltr_stats) && nrow(ltr_stats) > 0) sum(ltr_stats$count) else 0
  ltr_mb    <- if (!is.null(ltr_stats) && nrow(ltr_stats) > 0)
                 sprintf("%.1f Mb", sum(ltr_stats$bp) / 1e6) else "—"
  tir_count <- if (!is.null(tir_stats) && nrow(tir_stats) > 0) sum(tir_stats$count) else 0
  tir_mb    <- "—"  # bp not available from TIR summary txt
  line_count <- line_stats$regions %||% 0
  line_mb    <- if (!is.null(line_stats$bp) && line_stats$bp > 0)
                  sprintf("%.1f Mb", line_stats$bp / 1e6) else "—"

  # Build sub-report links for tool column
  ltr_link <- if (file.exists(file.path(outdir, "DANTE_LTR", "DANTE_LTR_summary.html")))
    'DANTE_LTR (<a href="DANTE_LTR/DANTE_LTR_summary.html">report</a>)' else "DANTE_LTR"
  tc_link  <- if (file.exists(file.path(outdir, "TideCluster", "default", "TideCluster_index.html")))
    'TideCluster (<a href="TideCluster/default/TideCluster_index.html">report</a>)' else "TideCluster"
  tir_link <- if (file.exists(file.path(outdir, "DANTE_TIR", "report.html")))
    'DANTE_TIR (<a href="DANTE_TIR/report.html">report</a>)' else "DANTE_TIR"

  sprintf('
<p class="annot-note">Structure-based annotation identifies complete element boundaries from
protein domain homology. Element counts below reflect intact, full-length elements only.
Per-lineage counts appear in the <a href="#s2">Full Classification Table</a> (Complete&nbsp;TEs column).</p>
<table class="data-table" style="max-width:600px">
<thead><tr><th>Tool</th><th style="text-align:right">Elements</th><th style="text-align:right">Coverage</th><th>Notes</th></tr></thead>
<tbody>
<tr><td>%s</td><td style="text-align:right">%s</td><td style="text-align:right">%s</td><td>Intact LTR retrotransposons</td></tr>
<tr><td>%s</td><td style="text-align:right">%s</td><td style="text-align:right">%s</td><td>Intact TIR transposons</td></tr>
<tr><td>%s</td><td style="text-align:right">%s</td><td style="text-align:right">%s</td><td>RT protein-domain loci (boundaries not delimited)</td></tr>
</tbody>
</table>',
    ltr_link,  format(ltr_count,  big.mark = ","), ltr_mb,
    tir_link,  format(tir_count,  big.mark = ","), tir_mb,
    "DANTE_LINE", format(line_count, big.mark = ","), line_mb
  )
}

# ── E6. Satellite cluster table ────────────────────────────────────────────
html_sat_table <- function(sat_data, genome_size) {
  if (is.null(sat_data$totals)) return("<p>No satellite data available.</p>")
  totals <- sort(sat_data$totals, decreasing = TRUE)
  rows <- mapply(function(cl, bp) {
    sprintf('<tr><td>%s</td><td style="text-align:right">%.2f</td><td style="text-align:right">%.3f%%</td></tr>',
            cl, bp / 1e6, bp / genome_size * 100)
  }, names(totals), totals, SIMPLIFY = FALSE)
  sprintf('
<div style="overflow-x:auto">
<table class="data-table">
<thead><tr><th>Cluster</th><th style="text-align:right">Total (Mb)</th><th style="text-align:right">%% genome</th></tr></thead>
<tbody>%s</tbody>
</table>
</div>', paste(rows, collapse = "\n"))
}

# ── E7. Assemble full HTML ─────────────────────────────────────────────────
assemble_html <- function(plotly_js, cards_html, sunburst_div, comp_table_html,
                           dante_summary_html, comp_bar_div,
                           density_top_div, density_lineage_div, density_trc_div,
                           sat_table_html, sat_bar_div,
                           outdir, bin_width, genome_avg_pct) {

  # Sub-report links — use real paths, not symlinks, so relative resources work
  sub_links <- ""
  real_paths <- list(
    list(file = file.path("DANTE_LTR", "DANTE_LTR_summary.html"), label = "DANTE_LTR"),
    list(file = file.path("TideCluster", "default", "TideCluster_index.html"),
         label = "TideCluster"),
    list(file = file.path("DANTE_TIR", "report.html"), label = "DANTE_TIR")
  )
  for (rp in real_paths) {
    if (file.exists(file.path(outdir, rp$file))) {
      href <- gsub("\\\\", "/", rp$file)  # normalise on Windows too
      sub_links <- paste0(sub_links,
        sprintf('<a class="report-link" href="%s">%s</a>\n', href, rp$label))
    }
  }

  sort_js <- '
function sortTable(id, col) {
  var tbl = document.getElementById(id).tBodies[0];
  var rows = Array.from(tbl.rows);
  var asc = tbl.dataset.sortCol == col && tbl.dataset.sortDir == "asc";
  rows.sort(function(a, b) {
    var va = a.cells[col].innerText.replace(/[,%]/g,"");
    var vb = b.cells[col].innerText.replace(/[,%]/g,"");
    var na = parseFloat(va), nb = parseFloat(vb);
    if (!isNaN(na) && !isNaN(nb)) return asc ? na-nb : nb-na;
    return asc ? va.localeCompare(vb) : vb.localeCompare(va);
  });
  rows.forEach(function(r){ tbl.appendChild(r); });
  tbl.dataset.sortCol = col;
  tbl.dataset.sortDir = asc ? "desc" : "asc";
}'

  sprintf('<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Repeat Annotation Report</title>
<style>
*{box-sizing:border-box;margin:0;padding:0}
body{font-family:"Segoe UI",Arial,sans-serif;background:#f5f5f5;color:#333;font-size:14px}
.topbar{background:#2c3e50;color:white;padding:14px 24px;display:flex;align-items:center;gap:24px;flex-wrap:wrap;position:sticky;top:0;z-index:100}
.topbar h1{font-size:1.15em;font-weight:600}
.topbar nav a{color:#ccc;text-decoration:none;font-size:0.88em;padding:4px 8px;border-radius:3px}
.topbar nav a:hover{background:#3d5166;color:white}
.container{max-width:1400px;margin:0 auto;padding:24px}
section{background:white;border-radius:8px;padding:24px;margin-bottom:24px;box-shadow:0 1px 4px rgba(0,0,0,.08)}
h2{color:#2c3e50;font-size:1.1em;margin-bottom:16px;padding-bottom:8px;border-bottom:2px solid #eee}
h3{color:#34495e;font-size:0.95em;margin:14px 0 8px}
.annot-note{background:#eaf2ff;border-left:4px solid #4575b4;padding:10px 14px;border-radius:0 4px 4px 0;margin-bottom:16px;font-size:0.88em;color:#2c3e50}
.cards{display:flex;gap:14px;flex-wrap:wrap;margin-bottom:4px}
.card{background:#f8f9fa;border-radius:6px;padding:14px 18px;min-width:160px;flex:1;border:1px solid #e9ecef}
.card-label{font-size:0.75em;color:#6c757d;text-transform:uppercase;letter-spacing:.04em;margin-bottom:4px}
.card-value{font-size:1.5em;font-weight:700;color:#2c3e50;line-height:1.1}
.card-sub{font-size:0.75em;color:#888;margin-top:3px}
.two-col{display:grid;grid-template-columns:1fr 1fr;gap:24px}
@media(max-width:800px){.two-col{grid-template-columns:1fr}}
.report-link{display:inline-block;background:#2c3e50;color:white;text-decoration:none;padding:7px 14px;border-radius:4px;margin:4px 6px 4px 0;font-size:0.88em}
.report-link:hover{background:#3d5166}
.data-table{width:100%%;border-collapse:collapse;font-size:0.85em}
.data-table th{background:#2c3e50;color:white;padding:7px 10px;text-align:left;white-space:nowrap}
.data-table td{padding:5px 10px;border-bottom:1px solid #f0f0f0}
.data-table tbody tr:hover{background:#f8f9fa}
.tr-total td{background:#f0f4f8}
.tr-total td:first-child{font-weight:600}
.tr-unspec td{color:#666}
.row-tag{font-size:0.72em;padding:1px 5px;border-radius:3px;vertical-align:middle;margin-left:4px}
.row-total{background:#d0e4f7;color:#1a4a72}
.row-unspec{background:#f0f0f0;color:#666}
.caption{font-size:0.8em;color:#888;margin-top:6px;font-style:italic}
</style>
<script>%s</script>
<script>%s</script>
</head>
<body>
<div class="topbar">
  <h1>Repeat Annotation Report</h1>
  <nav>
    <a href="#s1">Overview</a>
    <a href="#s2">Composition</a>
    <a href="#s3">Structure-based</a>
    <a href="#s4">Genomic distribution</a>
    <a href="#s5">Tandem repeats</a>
    <a href="#s6">Reports</a>
  </nav>
</div>
<div class="container">

<!-- SECTION 1: OVERVIEW -->
<section id="s1">
<h2>Overview</h2>
<p class="annot-note">
  Two complementary annotation strategies are used.
  <b>Similarity-based</b> (RepeatMasker): identifies repeats by sequence similarity
  to a constructed library — reflected in composition charts and density tracks.
  <b>Structure-based</b> (DANTE family): identifies complete elements by protein domain homology —
  reported in Section 3 and in the Complete&nbsp;TEs column of the classification table.
  The same genomic region may be covered by both.
</p>
%s
</section>

<!-- SECTION 2: REPEAT COMPOSITION -->
<section id="s2">
<h2>Repeat Composition — Similarity-based</h2>
<div class="two-col">
  <div>
    <h3>Hierarchy (click to explore)</h3>
    %s
    <p class="caption">Each arc represents elements classified only at that level
    (not further refined). Click a node to zoom into its subtree.</p>
  </div>
  <div>
    <h3>Full classification table</h3>
    <p class="caption" style="margin-bottom:8px">Bold rows show subtree totals.
    <i>Unspecified</i> rows are elements not classified to a finer level.
    Complete&nbsp;TEs column: intact elements from DANTE_LTR / DANTE_TIR.</p>
    %s
  </div>
</div>
</section>

<!-- SECTION 3: STRUCTURE-BASED ANNOTATION -->
<section id="s3">
<h2>Structure-based Annotation (DANTE)</h2>
%s
</section>

<!-- SECTION 4: GENOMIC DISTRIBUTION -->
<section id="s4">
<h2>Genomic Distribution</h2>
<h3>Repeat content per sequence</h3>
%s
<p class="caption">Fraction of each sequence covered by each repeat class (similarity-based).
Dashed line = genome-wide average (%.1f%%). Sequences &lt; %.0f kb aggregated into Other.
Per-sequence values can differ from the genome average because satellite arrays
are concentrated on specific chromosomes.</p>
<h3 style="margin-top:24px">Density — major repeat categories</h3>
<p class="caption">Gray: raw %s-kb bins. Red: smoothed. All sequences concatenated; vertical lines separate chromosomes.</p>
%s
<h3 style="margin-top:24px">Density — LTR lineages</h3>
%s
<h3 style="margin-top:24px">Density — satellite clusters (TideCluster)</h3>
%s
</section>

<!-- SECTION 5: SATELLITES -->
<section id="s5">
<h2>Tandem repeats (TideCluster)</h2>
<div class="two-col">
  <div>
    <h3>Top satellite clusters</h3>
    %s
  </div>
  <div>
    <h3>Cluster distribution per sequence</h3>
    %s
  </div>
</div>
</section>

<!-- SECTION 6: SUB-REPORTS -->
<section id="s6">
<h2>Detailed Sub-reports</h2>
<p>Open individual tool reports for element-level details:</p>
<div style="margin-top:12px">%s</div>
</section>

</div><!-- /container -->
</body>
</html>',
    plotly_js, sort_js,
    cards_html,
    sunburst_div,
    comp_table_html,
    dante_summary_html,
    comp_bar_div,
    genome_avg_pct,
    opt$min_len_chart / 1000,
    bin_width / 1000,
    density_top_div     %||% "<p>No top-level density data available.</p>",
    density_lineage_div %||% "<p>No lineage density data available.</p>",
    density_trc_div     %||% "<p>No TRC density data available.</p>",
    sat_table_html,
    sat_bar_div %||% "<p>No satellite density data available.</p>",
    sub_links
  )
}

# ═══════════════════════════════════════════════════════════════════════════
# F. MAIN
# ═══════════════════════════════════════════════════════════════════════════

main <- function() {
  message("=== Repeat Annotation Report Generator ===")
  message("Output directory: ", outdir)

  message("Loading genome info...")
  genome_info  <- load_genome_info(outdir)
  genome_size  <- sum(as.numeric(genome_info$length))

  message("Loading repeat composition...")
  comp         <- load_composition(outdir)

  message("Loading DANTE_LTR stats...")
  ltr_stats    <- load_dante_ltr_stats(outdir)

  message("Loading DANTE_TIR stats...")
  tir_stats    <- load_dante_tir_stats(outdir)

  message("Loading DANTE_LINE stats...")
  line_stats   <- load_dante_line_stats(outdir)

  message("Discovering BigWig files...")
  bw_files     <- discover_bw_files(outdir, opt$bin_width)

  seq_thresh   <- apply_seq_thresholds(genome_info, opt$min_len_chart,
                                        opt$min_len_tracks, opt$max_tracks)
  chart_seqs   <- seq_thresh$chart_seqs
  track_seqs   <- seq_thresh$track_seqs
  message(sprintf("Chart seqs: %d / Track seqs: %d", length(chart_seqs), length(track_seqs)))

  seq_len_named <- setNames(genome_info$length, genome_info$seqname)

  # ── Per-sequence composition ───────────────────────────────────────────
  message("Computing per-sequence coverage...")
  rm_bw <- bw_files$rm

  find_bw <- function(name) {
    m <- rm_bw[name]
    if (length(m) > 0 && !is.na(m) && file.exists(m)) m else NULL
  }

  tir_files <- rm_bw[grepl("^Class_II", names(rm_bw))]
  tir_files <- tir_files[file.exists(tir_files)]

  bar_bw_map <- list(
    "Tandem_repeats"  = bw_files$tc_agg,
    "Ty1/copia"       = find_bw("All_Ty1_Copia"),
    "Ty3/gypsy"       = find_bw("All_Ty3_Gypsy"),
    "LINE"            = find_bw("Class_I.LINE"),
    "TIR (Class II)"  = if (length(tir_files) > 0) tir_files else NULL,
    "Pararetrovirus"  = find_bw("Class_I.pararetrovirus"),
    "Simple repeats"  = find_bw("Simple_repeat") %||% find_bw("Simple_repeats"),
    "Low complexity"  = find_bw("Low_complexity"),
    "rDNA"            = find_bw("rDNA"),
    "Unknown"         = find_bw("Unknown")
  )
  bar_bw_map <- Filter(Negate(is.null), bar_bw_map)

  cov_mat <- if (length(bar_bw_map) > 0 && length(chart_seqs) > 0)
    compute_per_seq_coverage(bar_bw_map, chart_seqs, seq_len_named)
  else
    matrix(0, 0, 0)

  genome_avg_frac <- sum(comp$bp) / genome_size

  # ── Density track plots (3 panels matching summary_plots.pdf) ─────────────
  # Always use 100k resolution for the concatenated density plots.
  message("Building density track data...")
  d_bin   <- 100000L
  d_res   <- "100k"
  rm_100k <- file.path(outdir, "Repeat_Annotation_NoSat_split_by_class_bigwig", d_res)
  tc_agg_100k <- file.path(outdir, "TideCluster", "default",
                             paste0("TideCluster_clustering_", d_res, ".bw"))
  find_100k <- function(stem) {
    p <- file.path(rm_100k, paste0(stem, "_", d_res, ".bw"))
    if (file.exists(p)) p else NULL
  }

  # Panel 1 — major categories (matches page 1 of summary_plots.pdf)
  tir_100k <- {
    tf <- list.files(rm_100k, pattern = paste0("^Class_II.*_", d_res, "\\.bw$"),
                     full.names = TRUE)
    tf <- tf[file.exists(tf)]
    if (length(tf) > 0) tf[1] else NULL   # use first; concat handled inside
  }
  top_bw_map <- Filter(Negate(is.null), list(
    "Mobile elements" = find_100k("Mobile_elements"),
    "Low complexity"  = find_100k("Low_complexity"),
    "Simple repeats"  = find_100k("Simple_repeat") %||% find_100k("Simple_repeats"),
    "rDNA"            = find_100k("rDNA"),
    "Ty1/copia"       = find_100k("All_Ty1_Copia"),
    "Ty3/gypsy"       = find_100k("All_Ty3_Gypsy"),
    "Tandem_repeats"  = if (file.exists(tc_agg_100k)) tc_agg_100k else NULL
  ))

  # Panel 2 — LTR lineages (matches page 2)
  lineage_bw_map <- discover_lineage_bw_files(rm_100k, d_bin)

  # Panel 3 — TRC satellite clusters (matches page 3)
  trc_bw_map <- discover_trc_bw_files(outdir, d_bin)

  # ── Satellite cluster data (for Section 5 bar + table) ────────────────
  message("Building satellite cluster data...")
  sat_data <- build_satellite_coverage(bw_files$tc_cluster, chart_seqs,
                                        seq_len_named, opt$bin_width,
                                        opt$top_sat_clusters)

  # ── Build Plotly charts ────────────────────────────────────────────────
  message("Generating Plotly JSON...")
  sb_data        <- build_sunburst_data(comp)
  tree_df        <- build_comp_tree(comp, ltr_stats, tir_stats, line_stats,
                                     genome_size = genome_size)
  sunburst_chart <- json_sunburst(sb_data)

  comp_bar_chart <- if (nrow(cov_mat) > 0)
    json_composition_bar(cov_mat, genome_info, seq_thresh$n_other,
                          seq_thresh$other_bp, genome_avg_frac)
  else NULL

  message("Building density panel 1 (top-level categories)...")
  density_top_chart     <- json_concat_density_plot(top_bw_map,     genome_info, d_bin)
  message("Building density panel 2 (LTR lineages)...")
  density_lineage_chart <- json_concat_density_plot(lineage_bw_map, genome_info, d_bin)
  message("Building density panel 3 (TRC clusters)...")
  # Reverse so TRC_1 ends up at the top (matches summary_plots.pdf: rev(trc_bw))
  density_trc_chart     <- json_concat_density_plot(rev(trc_bw_map), genome_info, d_bin)

  sat_bar_chart <- if (!is.null(sat_data$per_seq) && nrow(sat_data$per_seq) > 0)
    json_satellite_bar(sat_data, genome_info)
  else NULL

  # ── Render HTML ────────────────────────────────────────────────────────
  message("Loading Plotly.js...")
  plotly_js <- load_plotly_js()

  cards_html       <- html_cards(genome_info, comp, ltr_stats, tir_stats,
                                   line_stats, seq_thresh)
  comp_table_html  <- html_comp_table(tree_df)
  dante_summ_html  <- html_dante_summary(ltr_stats, tir_stats, line_stats, outdir)
  sat_table_html   <- html_sat_table(sat_data, genome_size)

  sunburst_div      <- plotly_div("sunburst-plot", sunburst_chart)
  comp_bar_div      <- if (!is.null(comp_bar_chart))
    plotly_div("comp-bar-plot", comp_bar_chart) else "<p>No composition bar data.</p>"
  density_top_div     <- if (!is.null(density_top_chart))
    plotly_div("density-top-plot",     density_top_chart)     else NULL
  density_lineage_div <- if (!is.null(density_lineage_chart))
    plotly_div("density-lineage-plot", density_lineage_chart) else NULL
  density_trc_div     <- if (!is.null(density_trc_chart))
    plotly_div("density-trc-plot",     density_trc_chart)     else NULL
  sat_bar_div       <- if (!is.null(sat_bar_chart))
    plotly_div("sat-bar-plot", sat_bar_chart) else NULL

  html_out <- assemble_html(
    plotly_js, cards_html, sunburst_div, comp_table_html,
    dante_summ_html,
    comp_bar_div,
    density_top_div, density_lineage_div, density_trc_div,
    sat_table_html, sat_bar_div,
    outdir, d_bin, genome_avg_frac * 100
  )

  out_file <- file.path(outdir, "repeat_annotation_report.html")
  writeLines(html_out, out_file)
  message("Report written to: ", out_file)
}

main()
