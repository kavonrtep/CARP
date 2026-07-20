#!/usr/bin/env Rscript
# Regression test for the TRC monomer-size label lookup in make_repeat_report.R
# (mirrored in make_summary_plots.R).
#
# The satellite density label used to render as `TRC_1 (?bp)` because both report
# scripts read the OLD TideCluster kite file `monomer_size_best_estimate_stat.csv`
# / column `position`; current TideCluster (>=1.15) writes
# `monomer_size_top3_estimats.csv` / column `monomer_size`. The file never existed
# on a modern run, so the `?` default stuck. read_trc_monomer_sizes() now reads
# the current file and returns per-TRC modes; the caller drops the `(bp)` suffix
# when a size is unknown instead of printing `?`.
#
# read_trc_monomer_sizes() is copied VERBATIM from make_repeat_report.R (its main
# is unguarded, so it cannot be sourced). Keep in sync with the source.

# ---- copied verbatim from scripts/make_repeat_report.R ----------------------
read_trc_monomer_sizes <- function(kite_csv) {
  empty <- setNames(character(0), character(0))
  if (is.null(kite_csv) || !file.exists(kite_csv)) return(empty)
  ms <- tryCatch(read.table(kite_csv, header = TRUE, sep = "\t", check.names = FALSE,
                            stringsAsFactors = FALSE, quote = "", comment.char = ""),
                 error = function(e) NULL)
  if (is.null(ms) || nrow(ms) == 0 || !all(c("TRC_ID", "monomer_size") %in% names(ms)))
    return(empty)
  ids <- unique(as.character(ms$TRC_ID))
  vals <- vapply(ids, function(id) {
    v <- as.character(ms$monomer_size[as.character(ms$TRC_ID) == id])
    v <- v[!is.na(v) & nzchar(v)]
    if (length(v) == 0) return(NA_character_)
    names(sort(table(v), decreasing = TRUE))[1]        # mode across the TRC's arrays
  }, character(1))
  setNames(vals, ids)
}
# -----------------------------------------------------------------------------

# label builder mirrors discover_trc_bw_files(): drop (bp) when size unknown.
label_for <- function(tn, size_map) {
  sz <- if (tn %in% names(size_map)) size_map[[tn]] else NA_character_
  if (is.na(sz) || !nzchar(sz)) tn else paste0(tn, " (", sz, "bp)")
}

# current-schema kite CSV: TRC_1 two arrays both 174 (mode 174); TRC_2 174/174/350
# (mode 174, not the outlier); TRC_3 a single 42.
tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
writeLines(c(
  "TRC_ID\tseqid\tarray_length\tmonomer_size\tconsensus_period_bp",
  "TRC_1\tchr1\t46500\t174\t174",
  "TRC_1\tchr1\t9380\t174\t174",
  "TRC_2\tchr2\t5000\t174\t174",
  "TRC_2\tchr2\t5000\t174\t174",
  "TRC_2\tchr2\t5000\t350\t350",
  "TRC_3\tchr3\t2000\t42\t42"), tmp)

m <- read_trc_monomer_sizes(tmp)
stopifnot(identical(unname(m["TRC_1"]), "174"))
stopifnot(identical(unname(m["TRC_2"]), "174"))   # mode, not the 350 outlier
stopifnot(identical(unname(m["TRC_3"]), "42"))
stopifnot(identical(label_for("TRC_1", m), "TRC_1 (174bp)"))
stopifnot(identical(label_for("TRC_3", m), "TRC_3 (42bp)"))
stopifnot(identical(label_for("TRC_9", m), "TRC_9"))   # absent TRC -> no (bp), no ?
cat("  current-schema kite CSV: per-TRC monomer-size mode + label OK\n")

# missing file -> empty map -> label drops (bp), never prints "?"
m0 <- read_trc_monomer_sizes(file.path(tempdir(), "does_not_exist.csv"))
stopifnot(length(m0) == 0, identical(label_for("TRC_1", m0), "TRC_1"))
cat("  missing kite CSV: empty map, no (bp), no '?'\n")

# stale-schema file (old `position` column, no `monomer_size`) -> empty map
stale <- tempfile(fileext = ".csv"); on.exit(unlink(stale), add = TRUE)
writeLines(c("TRC_ID\tposition", "TRC_1\t174"), stale)
ms <- read_trc_monomer_sizes(stale)
stopifnot(length(ms) == 0, identical(label_for("TRC_1", ms), "TRC_1"))
cat("  stale-schema CSV (position col): ignored, no crash\n")

cat("test_trc_monomer_label: PASSED\n")
