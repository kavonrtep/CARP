# Per-TRC tandem monomer period, shared by make_unified_annotation.R (TE_origin
# domain-rhythm gate), make_repeat_report.R and make_summary_plots.R (the
# "TRC_n (<bp>bp)" density-track labels). One definition, so the annotation and
# both reports always quote the same monomer for a TRC.

# Authoritative per-TRC tandem monomer period (bp) from TideCluster's report
# table `trc_table.tsv`. Per TRC: `prevalent_founder` (the founder period most
# arrays of the TRC agree on) when present, else `monomer_tarean` (TAREAN family
# consensus), else `monomer_kite` (most-frequent KITE founder period). Founder
# first because TAREAN can also lock onto a sub-repeat (run-000076: TRC_1 tarean
# 20 bp vs founder 170 bp; 45S rDNA tarean ~7.8 kb vs founder ~10.8 kb). Returns a named integer vector TRC_ID -> period; empty when
# the file is absent/unusable (optional — `--no_rdna`/older TideCluster/purged).
#
# Why NOT the kite `monomer_size` CSV: that column is the top k-mer *peak*, which
# can lock onto a short SSR sub-period (measured: 79 bp reported for a genuine
# 13134 bp TIR-derived monomer; 22 bp for the ~10.8 kb 45S rDNA unit). Tiling the
# domain-rhythm occupancy test at 79 bp wrongly reads a real TIR tandem as sparse;
# the founder/TAREAN period is correct.
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
    for (col in c("prevalent_founder", "monomer_tarean", "monomer_kite")) {
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
