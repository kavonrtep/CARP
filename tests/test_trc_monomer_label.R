#!/usr/bin/env Rscript
# Regression test for the TRC monomer-size label on the satellite density panel
# (make_repeat_report.R discover_trc_bw_files() and make_summary_plots.R).
#
# The label used to come from the kite `monomer_size_top3_estimats.csv` column
# `monomer_size` (per-array top k-mer peak, mode across arrays). That peak locks
# onto short internal sub-repeats: run-000076 labelled five 45S rDNA TRCs
# (~10.8 kb unit) as "22bp". Both reports now use read_trc_periods()
# (scripts/trc_periods.R) — the same trc_table.tsv source and fallback
# (prevalent_founder -> monomer_tarean -> monomer_kite) as make_unified_annotation.R.
source(file.path("scripts", "trc_periods.R"))

# label builder mirrors discover_trc_bw_files(): drop (bp) when no estimate.
label_for <- function(tn, size_map) {
  if (tn %in% names(size_map)) paste0(tn, " (", size_map[[tn]], "bp)") else tn
}

tmp <- tempfile(fileext = ".tsv"); on.exit(unlink(tmp), add = TRUE)
writeLines(c(
  "TRC_ID\trepeat_type\tmonomer_kite\tmonomer_tarean\tprevalent_founder",
  "TRC_1\tTR\t170\t20\t170",        # founder 170, not the 20 bp TAREAN sub-repeat
  "TRC_2\tTR\t22\t7810\t10808",     # 45S rDNA: founder ~10.8 kb, never the 22 bp kite peak
  "TRC_3\tTR\t351\t351\t",          # no founder -> tarean
  "TRC_4\tTR\t2784\t\t",            # only kite
  "TRC_5\tTR\t\t\t"), tmp)          # no estimate at all

m <- read_trc_periods(tmp)
stopifnot(identical(label_for("TRC_1", m), "TRC_1 (170bp)"))
stopifnot(identical(label_for("TRC_2", m), "TRC_2 (10808bp)"))
stopifnot(identical(label_for("TRC_3", m), "TRC_3 (351bp)"))
stopifnot(identical(label_for("TRC_4", m), "TRC_4 (2784bp)"))
stopifnot(identical(label_for("TRC_5", m), "TRC_5"))   # no estimate -> no (bp), no ?
stopifnot(identical(label_for("TRC_9", m), "TRC_9"))   # absent TRC
cat("  trc_table fallback prevalent_founder -> tarean -> kite: labels OK\n")

# missing file (older TideCluster / failed report) -> bare labels, never "?"
m0 <- read_trc_periods(file.path(tempdir(), "does_not_exist.tsv"))
stopifnot(length(m0) == 0, identical(label_for("TRC_1", m0), "TRC_1"))
cat("  missing trc_table.tsv: empty map, no (bp), no '?'\n")

cat("test_trc_monomer_label: PASSED\n")
