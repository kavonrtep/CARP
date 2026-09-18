#!/usr/bin/env Rscript
# Regression test for the "Complete TEs" column of the report's composition
# table (build_comp_tree in make_repeat_report.R).
#
# The column used to render BLANK on every internal node: dfs() attached
# dante_count only to `leaf` rows, and passed NA_integer_ on both the `total`
# and `unspecified` rows. So a classification that carries its own bp AND has
# annotated children showed its bp but no element count anywhere. It was also
# data-dependent -- with no annotated child the same node became a leaf and the
# count DID appear, so identical biology rendered two different ways.
#
# This matters most under DANTE_LTR core mode, where an element's call routinely
# stops at an internal node (Ty3_gypsy, chromovirus): on a Draparnaldia assembly
# 1,320 of 1,530 complete elements sit on such nodes. It is a pre-existing bug
# though -- DANTE tier-2 domains and RepeatMasker already produce internal-node
# rows today, just tiny ones (0.0007-0.13% of the genome in current runs).
#
# Now: the `unspecified` row carries the node's OWN count (parallel to its own
# bp) and the `total` row carries the subtree sum (parallel to its subtree bp).
#
# build_comp_tree() and its two globals are copied VERBATIM from
# scripts/make_repeat_report.R (its main is unguarded, so it cannot be sourced).
# Keep in sync with the source.

# ---- copied verbatim from scripts/make_repeat_report.R ----------------------
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1])) a else b

CATEGORY_ORDER <- c("Class_I", "Class_II", "rDNA", "rDNA_45S", "rDNA_5S",
                    "Tandem_repeats", "Simple_repeat", "Low_complexity", "Unknown")

build_comp_tree <- function(comp, ltr_stats, tir_stats, line_stats = NULL,
                            genome_size = NULL, ltr_tr_stats = NULL) {
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

  # LTR_RT_TR member counts: path → number of complete copies inside tandem arrays
  ltr_tr_counts <- setNames(integer(0), character(0))
  if (!is.null(ltr_tr_stats) && nrow(ltr_tr_stats) > 0) {
    for (i in seq_len(nrow(ltr_tr_stats))) {
      p <- ltr_tr_stats$path[i]
      ltr_tr_counts[p] <- (ltr_tr_counts[p] %||% 0L) + ltr_tr_stats$tr_members[i]
    }
  }

  # Subtree sums for the count columns, on the same basis as `subtree_bp`: a
  # "Total" row aggregates its whole subtree, so the element counts must
  # aggregate too. They used to render blank on every internal node, which reads
  # as "no elements" on a node whose children carry them all. Returns NA when
  # nothing in the subtree has a count, so the cell stays empty rather than
  # printing 0.
  subtree_count <- function(id, tbl) {
    own  <- as.numeric(unname(tbl[id]))
    kids <- children_of[[id]]
    vals <- c(own, if (length(kids) > 0)
                     vapply(kids, subtree_count, numeric(1), tbl = tbl)
                   else numeric(0))
    if (all(is.na(vals))) return(NA_real_)
    sum(vals, na.rm = TRUE)
  }
  as_count <- function(x) if (is.na(x)) NA_integer_ else as.integer(x)
  subtree_dc_cache  <- sapply(all_ids, subtree_count, tbl = dante_counts)
  subtree_trc_cache <- sapply(all_ids, subtree_count, tbl = ltr_tr_counts)

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
    trc_raw  <- ltr_tr_counts[id]
    trc      <- if (length(trc_raw) > 0 && !is.na(trc_raw)) unname(trc_raw) else NA_integer_
    tot_dc   <- as_count(unname(subtree_dc_cache[id]))
    tot_trc  <- as_count(unname(subtree_trc_cache[id]))

    if (has_kids) {
      # Total row
      rows[[length(rows) + 1]] <<- list(
        path         = id,
        label        = label,
        row_type     = "total",
        depth        = depth,
        bp           = tot_bp,
        pct          = tot_bp / genome_size * 100,
        dante_count  = tot_dc,
        ltr_tr_count = tot_trc
      )
      # Unspecified row (own CSV value) only if non-zero
      if (own_bp > 0) {
        rows[[length(rows) + 1]] <<- list(
          path         = id,
          label        = label,
          row_type     = "unspecified",
          depth        = depth + 1L,
          bp           = own_bp,
          pct          = unname(csv_pct[id] %||% 0),
          dante_count  = dc,
          ltr_tr_count = trc
        )
      }
      # Recurse into children sorted by subtree bp descending
      kids_sorted <- kids[order(-subtree_bp_cache[kids])]
      for (k in kids_sorted) dfs(k, depth + 1L)
    } else {
      # Leaf row — show DANTE count (and LTR_RT_TR-member count) if available
      rows[[length(rows) + 1]] <<- list(
        path         = id,
        label        = label,
        row_type     = "leaf",
        depth        = depth,
        bp           = own_bp,
        pct          = unname(csv_pct[id] %||% 0),
        dante_count  = dc,
        ltr_tr_count = trc
      )
    }
  }

  # Find top-level nodes (no parent in all_ids). Fixed category order
  # (CATEGORY_ORDER); anything unlisted is appended, ordered by size.
  top_nodes <- all_ids[is.na(parent_of)]
  top_nodes <- top_nodes[order(match(top_nodes, CATEGORY_ORDER),
                               -subtree_bp_cache[top_nodes])]
  for (id in top_nodes) dfs(id, 0L)

  do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
}

# -----------------------------------------------------------------------------

ok <- 0L
check <- function(what, got, want) {
  if (!isTRUE(all.equal(got, want))) {
    stop(sprintf("FAIL %s: got %s, want %s", what,
                 paste(format(got), collapse = ","),
                 paste(format(want), collapse = ",")), call. = FALSE)
  }
  ok <<- ok + 1L
}
row_of <- function(df, path, type) {
  r <- df[df$path == path & df$row_type == type, , drop = FALSE]
  if (nrow(r) != 1L) stop(sprintf("expected exactly 1 %s row for %s, got %d",
                                  type, path, nrow(r)), call. = FALSE)
  r
}

# A core-mode-shaped tree: the superfamily and the chromovirus node each carry
# their own bp (elements whose lineage could not be resolved) and both have an
# annotated child. Class_II/.../hAT is the control: children, but no counts
# anywhere in the subtree.
comp <- data.frame(
  type = c("Class_I/LTR/Ty3_gypsy",
           "Class_I/LTR/Ty3_gypsy/chromovirus",
           "Class_I/LTR/Ty3_gypsy/chromovirus/Chlamyvir",
           "Class_II/Subclass_1/TIR/hAT/child"),
  bp   = c(2000000, 9000000, 2000000, 50000),
  pct  = c(0.5, 2.25, 0.5, 0.0125),
  stringsAsFactors = FALSE
)
ltr_stats <- data.frame(
  path  = c("Class_I/LTR/Ty3_gypsy",
            "Class_I/LTR/Ty3_gypsy/chromovirus",
            "Class_I/LTR/Ty3_gypsy/chromovirus/Chlamyvir"),
  count = c(225L, 985L, 205L),
  bp    = c(2000000, 9000000, 2000000),
  stringsAsFactors = FALSE
)
ltr_tr_stats <- data.frame(
  path       = "Class_I/LTR/Ty3_gypsy/chromovirus/Chlamyvir",
  tr_members = 12L, tr_arrays = 3L,
  stringsAsFactors = FALSE
)
tree <- build_comp_tree(comp, ltr_stats, tir_stats = NULL, line_stats = NULL,
                        genome_size = 400000000, ltr_tr_stats = ltr_tr_stats)

# 1. The node's own elements appear on its "unspecified" row (was NA).
check("chromovirus unspecified count",
      row_of(tree, "Class_I/LTR/Ty3_gypsy/chromovirus", "unspecified")$dante_count, 985L)
check("Ty3_gypsy unspecified count",
      row_of(tree, "Class_I/LTR/Ty3_gypsy", "unspecified")$dante_count, 225L)

# 2. "Total" rows aggregate the subtree, like their bp column does.
check("chromovirus total count",
      row_of(tree, "Class_I/LTR/Ty3_gypsy/chromovirus", "total")$dante_count, 985L + 205L)
check("Ty3_gypsy total count",
      row_of(tree, "Class_I/LTR/Ty3_gypsy", "total")$dante_count, 225L + 985L + 205L)
check("synthetic ancestor total count",
      row_of(tree, "Class_I/LTR", "total")$dante_count, 1415L)

# 3. Leaf rows are untouched.
check("Chlamyvir leaf count",
      row_of(tree, "Class_I/LTR/Ty3_gypsy/chromovirus/Chlamyvir", "leaf")$dante_count, 205L)

# 4. No count anywhere in a subtree stays NA, so the cell renders empty, not 0.
check("no-count subtree stays NA",
      is.na(row_of(tree, "Class_II/Subclass_1/TIR/hAT", "total")$dante_count), TRUE)

# 5. LTR_RT_TR member counts roll up the same way.
check("LTR_RT_TR rolls up to the superfamily total",
      row_of(tree, "Class_I/LTR/Ty3_gypsy", "total")$ltr_tr_count, 12L)
check("LTR_RT_TR absent on an unrelated unspecified row",
      is.na(row_of(tree, "Class_I/LTR/Ty3_gypsy", "unspecified")$ltr_tr_count), TRUE)

# 6. The bp columns must not have moved (regression guard on existing output).
check("chromovirus total bp",
      row_of(tree, "Class_I/LTR/Ty3_gypsy/chromovirus", "total")$bp, 11000000)
check("chromovirus unspecified bp",
      row_of(tree, "Class_I/LTR/Ty3_gypsy/chromovirus", "unspecified")$bp, 9000000)

cat(sprintf("test_comp_tree_counts.R: %d checks passed\n", ok))
