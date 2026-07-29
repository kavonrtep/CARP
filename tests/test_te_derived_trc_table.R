#!/usr/bin/env Rscript
# Unit test for the TE-derived TRC logic in make_unified_annotation.R:
#   .order_domains          — canonical protein-domain ordering
#   read_trc_periods        — per-TRC tandem period from trc_table.tsv
#                             (monomer_tarean -> monomer_kite -> prevalent_founder)
#   te_domain_rhythm        — per-monomer domain occupancy (the gate criterion)
#   write_te_derived_trc_table — the CSV summary (element counts + occupancy cols)
#
# Behaviours covered:
#   1. Domain ordering: canonical, unknown domains alphabetical after.
#   2. Period source: TAREAN consensus wins; KITE founder / prevalent_founder
#      fall back; all-empty rows are dropped. (The kite `monomer_size` SSR
#      sub-period trap is why the source is trc_table, not the kite CSV.)
#   3. Domain rhythm: a TE-derived tandem has domains in ~every monomer window
#      (occupancy -> 1); a satellite interrupted by a single TE block has low
#      occupancy; per-array `frac` needs MOST arrays in-rhythm.
#   4. Container case: the LTR_RT_TR container is NOT counted as one element;
#      its member copies are (the run116 TRC_9 undercount fix).
#   5. CSV carries the monomer (from the period map) and the occupancy columns.
#
# The definitions below are copied VERBATIM from make_unified_annotation.R (its
# main is unguarded, so the script cannot be sourced in isolation). Keep in sync.
suppressPackageStartupMessages(library(GenomicRanges))
log_msg <- function(...) invisible(NULL)   # stub

# ---- copied verbatim from make_unified_annotation.R -------------------------
.DOMAIN_ORDER <- c("GAG","PROT","AP","INT","RT","RH","aRH",
                   "CHD","CHDCR","CHDII","TPase","ENDO","EN","HEL1","HEL2")
.order_domains <- function(d) {
  d <- unique(d[!is.na(d) & nzchar(d)])
  if (length(d) == 0) return(character(0))
  rnk <- match(d, .DOMAIN_ORDER)
  d[order(ifelse(is.na(rnk), length(.DOMAIN_ORDER) + 1L, rnk), d)]
}

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
    arr <- reduce(granges(te[sel]), ignore.strand = TRUE)
    strand(arr) <- "*"
    array_bp <- sum(as.numeric(width(arr)))

    te_cls    <- mode1(as.character(te$TE_origin[sel]))
    te_struct <- mode1(as.character(te$TE_origin_structure[sel]))

    n_complete <- 0L; cov_bp <- 0
    if (length(elem) > 0 && length(arr) > 0) {
      h <- unique(subjectHits(findOverlaps(arr, elem)))
      n_complete <- length(h)
      if (n_complete > 0) cov_bp <- sum(as.numeric(width(intersect(arr, elem[h]))))
    }
    frac <- if (array_bp > 0) cov_bp / array_bp else NA_real_

    dom_str <- ""
    if (length(t2p) > 0 && length(t2_dom) > 0 && length(arr) > 0) {
      dh <- unique(queryHits(findOverlaps(t2p, arr)))
      if (length(dh) > 0) dom_str <- paste(.order_domains(t2_dom[dh]), collapse = "|")
    }

    per_map  <- if (run == "short") period_short else period_default
    mono_n   <- if (length(per_map) > 0 && !is.null(per_map[[trc]])) per_map[[trc]] else NA_integer_
    mono_n   <- suppressWarnings(as.integer(mono_n))
    n_exp    <- if (!is.na(mono_n) && mono_n > 0) as.integer(round(array_bp / mono_n)) else NA_integer_

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
  trc_num <- suppressWarnings(as.integer(sub("^TRC_", "", df$trc_id)))
  df <- df[order(df$run, trc_num, df$trc_id), , drop = FALSE]
  write.table(df, out_csv, sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
  log_msg(sprintf("Wrote %d TE-derived TRC row(s) -> %s", nrow(df), basename(out_csv)))
  invisible()
}
# ---- end verbatim copy ------------------------------------------------------

fail <- function(msg) { cat("FAIL:", msg, "\n"); quit(status = 1L) }
eq   <- function(got, want, what)
  if (!identical(as.character(got), as.character(want)))
    fail(sprintf("%s: got '%s' want '%s'", what, got, want))

gr <- function(sn, s, e) GRanges(sn, IRanges(s, e), strand = "*")

# ---- 1. .order_domains ------------------------------------------------------
eq(paste(.order_domains(c("RT","GAG","INT","RH")), collapse="|"), "GAG|INT|RT|RH", "domain order")
eq(paste(.order_domains(c("RT","ZZZ","GAG","ABC")), collapse="|"), "GAG|RT|ABC|ZZZ", "domain order unknown")

# ---- 2. read_trc_periods: tarean -> kite -> prevalent, drop empty -----------
tf <- tempfile(fileext = ".tsv")
writeLines(c("TRC_ID\tmonomer_tarean\tmonomer_kite\tprevalent_founder",
             "TRC_1\t13134\t79\t18591",   # TAREAN wins over the SSR kite peak
             "TRC_2\t\t2455\t",           # no tarean -> kite founder
             "TRC_3\t\t\t900",            # only prevalent_founder
             "TRC_4\t\t\t"), tf)          # nothing -> dropped
p <- read_trc_periods(tf)
eq(p[["TRC_1"]], 13134, "period: TAREAN wins")
eq(p[["TRC_2"]], 2455,  "period: kite founder fallback")
eq(p[["TRC_3"]], 900,   "period: prevalent_founder fallback")
eq("TRC_4" %in% names(p), FALSE, "period: all-empty row dropped")
eq(length(read_trc_periods("/no/such/file.tsv")), 0, "period: missing file -> empty")

# ---- 3. te_domain_rhythm: derived vs insertion vs multi-array ---------------
# derived: 1000 bp array, P=100 -> 10 windows, a domain in each -> occ 1.0
arr1 <- gr("c1", 1, 1000)
dense <- gr("c1", seq(20, 920, by = 100), seq(45, 945, by = 100))   # one per window
r <- te_domain_rhythm(arr1, dense, 100)
eq(sprintf("%.3f", r$occ), "1.000", "rhythm derived occ")
eq(r$frac, 1, "rhythm derived frac")
# insertion: same array, domains only in a single 2-window block -> occ 0.2
block <- gr("c1", c(20, 120), c(45, 145))
r2 <- te_domain_rhythm(arr1, block, 100)
eq(sprintf("%.3f", r2$occ), "0.200", "rhythm insertion occ")
eq(r2$frac, 0, "rhythm insertion frac (single array < 0.5)")
# no domains at all -> occ 0 (the TRC_421/469 embedded-satellite case)
r0 <- te_domain_rhythm(arr1, GRanges(), 100)
eq(r0$occ, 0, "rhythm no-domains occ 0")
# multi-array: one fully occupied, one empty -> frac 0.5
arr2 <- c(gr("c1", 1, 1000), gr("c2", 1, 1000))
r3 <- te_domain_rhythm(arr2, dense, 100)   # dense only overlaps c1
eq(r3$frac, 0.5, "rhythm multi-array frac 0.5")

# ---- 4/5. write_te_derived_trc_table: element count + period + occupancy ----
mk_te <- function(sn, s, e, name, cls, tool, struct = NA_character_) {
  g <- GRanges(sn, IRanges(s, e), strand = "+")
  g$Name <- name; g$TE_origin <- cls; g$TE_origin_structure <- struct; g$source_tool <- tool
  g
}
# TRC_5 default: satellite c1:1-1000 over 2 standalone Ale, domains in every 100bp window
te5 <- mk_te("c1", 1, 1000, "TRC_5", "Class_I/LTR/Ty1_copia/Ale", "TideCluster_default")
# TRC_9 default: satellite c2:1-1200 over an LTR_RT_TR container with 3 members
te9 <- mk_te("c2", 1, 1200, "TRC_9", "Class_I/LTR/Ty3_gypsy/chromovirus/CRM",
             "TideCluster_default", "tandem_LTR_RT")
level1 <- suppressWarnings(c(te5, te9))

mk_t1 <- function(sn, s, e, struct = NA_character_) {
  g <- GRanges(sn, IRanges(s, e), strand = "+"); g$structure <- struct; g }
t1 <- suppressWarnings(c(
  mk_t1("c1", 100, 400), mk_t1("c1", 600, 900),       # TRC_5: 2 standalone Ale
  mk_t1("c2", 101, 700, "LTR_RT_TR")                  # TRC_9: container (excluded from count)
))
t1_members <- GRanges("c2", IRanges(c(101, 301, 501), c(300, 500, 700)), strand = "+")  # 3 copies
# domains: one per 100bp window on c1 (TRC_5) and c2 (TRC_9) -> occupancy 1.0
d_c1 <- GRanges("c1", IRanges(seq(20, 920, by = 100), seq(45, 945, by = 100)), strand = "+")
d_c1$domain <- rep(c("GAG","PROT","INT","RT","RH"), length.out = length(d_c1))
d_c2 <- GRanges("c2", IRanges(seq(20, 1120, by = 100), seq(45, 1145, by = 100)), strand = "+")
d_c2$domain <- rep(c("GAG","PROT","RT","RH","INT"), length.out = length(d_c2))
t2 <- suppressWarnings(c(d_c1, d_c2))

period_default <- c(TRC_5 = 100L, TRC_9 = 200L)   # expected monomers: 10, 6
out <- tempfile(fileext = ".csv")
write_te_derived_trc_table(level1, t1, t1_members, t2, period_default, integer(0), out)
df <- read.table(out, header = TRUE, sep = ",", stringsAsFactors = FALSE,
                 quote = "", colClasses = "character")
r5 <- df[df$trc_id == "TRC_5", ]; r9 <- df[df$trc_id == "TRC_9", ]

# scattered case
eq(r5$n_complete_elements, 2, "TRC_5 n_complete (2 standalone)")
eq(r5$monomer_length_bp, 100, "TRC_5 monomer from period map")
eq(r5$n_expected_monomers, 10, "TRC_5 expected monomers")
eq(r5$domain_occupancy, "1.0000", "TRC_5 occupancy 1.0")
eq(r5$frac_arrays_in_rhythm, "1.0000", "TRC_5 frac in rhythm")

# container case: count the 3 members, NOT the container
eq(r9$n_complete_elements, 3, "TRC_9 n_complete (members, not container)")
eq(r9$te_origin_structure, "tandem_LTR_RT", "TRC_9 structure")
eq(r9$monomer_length_bp, 200, "TRC_9 monomer from period map")
eq(r9$domain_occupancy, "1.0000", "TRC_9 occupancy 1.0")

# header-only when no TE_origin
out2 <- tempfile(fileext = ".csv")
write_te_derived_trc_table(GRanges(), t1, t1_members, t2, period_default, integer(0), out2)
if (length(readLines(out2)) != 1L) fail("empty input should write header-only")

cat("test_te_derived_trc_table: PASSED\n")
