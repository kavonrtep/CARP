#!/usr/bin/env Rscript
# Regression test for the LTR density panel's track discovery
# (discover_lineage_bw_files in make_repeat_report.R).
#
# The panel filters discovered BigWigs against a fixed list of REXdb lineage
# names, taking the last dot-component of the filename as the label. Under
# `dante_ltr_mode: core` an element's classification routinely stops at an
# INTERNAL node -- Class_I/LTR/Ty3_gypsy, or .../chromovirus -- so its track's
# last component is not a lineage name and the fixed list DISCARDED it. On a
# Draparnaldia run that is 86% of complete elements, which would render a panel
# showing two lineages next to a composition table reporting 3.3% LTR content.
#
# The widening is gated on the mode on purpose: internal-node tracks also occur
# in lineage mode (from DANTE domain calls and RepeatMasker), but at
# 0.001-0.13% of the genome, and the canonical report must not move. This test
# pins both halves -- lineage mode unchanged, core mode inclusive.
#
# discover_lineage_bw_files() is copied VERBATIM from
# scripts/make_repeat_report.R (its main is unguarded, so it cannot be sourced).
# Keep in sync with the source.

REPO <- normalizePath(file.path(dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])), ".."))
source(file.path(REPO, "scripts", "classification.R"))   # provides get_canonical()

# ---- copied verbatim from scripts/make_repeat_report.R ----------------------
discover_lineage_bw_files <- function(rm_dir, bin_width, mode = "lineage") {
  suffix <- if (bin_width == 100000L) "_100k.bw" else "_10k.bw"
  known  <- c("Ale","Alesia","Angela","Bianca","Bryco","Gymco-I","Gymco-II","Gymco-III",
              "Gymco-IV","Ikeros","Ivana","Lyco","Osser","SIRE","TAR","Tork","Chlamyvir",
              "chromo-unclass","CRM","Galadriel","Reina","Tcn1","Tekay","Athila","Ogre",
              "Retand","TatI","TatII","TatIII","Phygy","Selgy")
  if (!dir.exists(rm_dir)) return(list())
  pat    <- paste0("LTR\\.Ty.*", gsub("\\.", "\\\\.", suffix), "$")
  fnames <- list.files(rm_dir, pattern = pat, full.names = FALSE)
  if (length(fnames) == 0) return(list())

  if (identical(mode, "core")) {
    # Accept any node under Class_I/LTR/Ty* in the vocabulary, in vocabulary
    # order, and mark the internal ones so a bucket does not read as a lineage.
    vocab <- tryCatch(get_canonical(), error = function(e) character(0))
    ltr   <- vocab[startsWith(vocab, "Class_I/LTR/Ty")]
    if (length(ltr) > 0) {
      paths <- gsub("\\.", "/", sub("_100k\\.bw$|_10k\\.bw$", "", fnames))
      keep  <- paths %in% ltr
      fnames <- fnames[keep]; paths <- paths[keep]
      if (length(fnames) == 0) return(list())
      idx   <- match(paths, ltr)
      internal <- vapply(ltr, function(p) any(startsWith(ltr, paste0(p, "/"))),
                         logical(1))
      leaf  <- vapply(strsplit(paths, "/", fixed = TRUE),
                      function(x) x[length(x)], character(1))
      lnames <- ifelse(internal[idx], paste0(leaf, " (unspecified)"), leaf)
      ord   <- order(idx)
      return(setNames(as.list(file.path(rm_dir, fnames[ord])), lnames[ord]))
    }
    # No usable vocabulary -> fall through to the fixed list rather than nothing.
  }

  lnames <- gsub("_100k\\.bw$|_10k\\.bw$", "", fnames)
  lnames <- gsub(".*\\.", "", lnames)
  keep   <- lnames %in% known
  fnames <- fnames[keep]; lnames <- lnames[keep]
  if (length(fnames) == 0) return(list())
  ord    <- order(match(lnames, known))
  setNames(as.list(file.path(rm_dir, fnames[ord])), lnames[ord])
}
# -----------------------------------------------------------------------------

ok <- 0L
check <- function(what, got, want) {
  if (!isTRUE(all.equal(got, want)))
    stop(sprintf("FAIL %s:\n  got  %s\n  want %s", what,
                 paste(got, collapse = ", "), paste(want, collapse = ", ")), call. = FALSE)
  ok <<- ok + 1L
}

# A directory holding both lineage-leaf and internal-node tracks, exactly as
# make_bigwig_density emits them from the per-class split.
d <- file.path(tempdir(), "bwpanel"); dir.create(d, showWarnings = FALSE)
on.exit(unlink(d, recursive = TRUE), add = TRUE)
files <- c("Class_I.LTR.Ty1_copia.Ale_100k.bw",
           "Class_I.LTR.Ty1_copia.SIRE_100k.bw",
           "Class_I.LTR.Ty1_copia_100k.bw",                       # internal
           "Class_I.LTR.Ty3_gypsy_100k.bw",                       # internal
           "Class_I.LTR.Ty3_gypsy.chromovirus_100k.bw",           # internal
           "Class_I.LTR.Ty3_gypsy.chromovirus.Tekay_100k.bw",
           "Class_I.LTR.Ty3_gypsy.non-chromovirus.OTA.Tat.Ogre_100k.bw",
           "Class_II.Subclass_1.TIR.hAT_100k.bw")                 # not LTR at all
invisible(file.create(file.path(d, files)))

lin  <- discover_lineage_bw_files(d, 100000L)                     # default mode
core <- discover_lineage_bw_files(d, 100000L, mode = "core")

# 1. Lineage mode is unchanged: leaves only, in the fixed order, no buckets.
check("lineage labels", names(lin), c("Ale", "SIRE", "Tekay", "Ogre"))
check("lineage has no (unspecified) rows", any(grepl("unspecified", names(lin))), FALSE)

# 2. Core mode adds the internal nodes, labelled so they cannot read as lineages.
check("core adds the three internal nodes",
      sort(setdiff(names(core), names(lin))),
      sort(c("Ty1_copia (unspecified)", "Ty3_gypsy (unspecified)",
             "chromovirus (unspecified)")))
check("core keeps every lineage leaf", all(names(lin) %in% names(core)), TRUE)

# 3. Neither mode picks up a non-LTR track.
check("no Class_II track in either mode",
      any(grepl("hAT", c(names(lin), names(core)))), FALSE)

# 4. Core ordering follows the vocabulary: a parent precedes its children.
p <- match("Ty3_gypsy (unspecified)", names(core))
c1 <- match("chromovirus (unspecified)", names(core))
c2 <- match("Tekay", names(core))
check("parent precedes child in core order", p < c1 && c1 < c2, TRUE)

# 5. Every returned value is a real path in the directory.
check("paths resolve", all(file.exists(unlist(core))), TRUE)

cat(sprintf("test_lineage_bw_panel.R: %d checks passed\n", ok))
