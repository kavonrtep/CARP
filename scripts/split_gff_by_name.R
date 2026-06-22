#!/usr/bin/env Rscript
# Split a GFF3 into one file per distinct `Name` attribute.
#
# Used to produce per-family inputs for density BigWig generation. With
# --name-prefix TRC_ it extracts only the tandem-repeat families (Name=TRC_<n>)
# from a larger annotation (e.g. Repeat_Annotation_Unified.gff3) and writes one
# GFF3 per family; feeding those to calculate_density_batch.R yields one density
# BigWig per family. See docs/bigwig_redesign_plan.md.
#
# Robust to empty input: a missing, 0-byte, or header-only GFF3 (a legitimate
# state — e.g. RM_on_TideCluster_Library.gff3 is truncated to 0 bytes when there
# is no dimer library) produces no output files and exits 0, so callers under
# `set -euo pipefail` do not fail.

suppressPackageStartupMessages({
  library(optparse)
  library(rtracklayer)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input GFF3"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL,
              help="Output directory (one <Name>.gff3 per distinct Name)"),
  make_option(c("-p", "--name-prefix"), type="character", default=NULL,
              dest="name_prefix",
              help="If set, only split features whose Name starts with this prefix (e.g. TRC_)")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input) || is.null(opt$output_dir)){
  stop("Please provide --input and --output_dir")
}

dir.create(opt$output_dir, showWarnings=FALSE, recursive=TRUE)

# Guard empty input *before* import(): rtracklayer's GFF sniffer dies on a
# 0-byte file ("argument is of length zero" in .sniffGFFVersion).
if (!file.exists(opt$input) || file.info(opt$input)$size == 0){
  message("Input ", opt$input, " is missing or empty — nothing to split")
  quit(save="no", status=0)
}

g <- import(opt$input, format="gff3")

if (length(g) == 0 || is.null(g$Name)){
  message("No features (or no Name attribute) in ", opt$input, " — nothing to split")
  quit(save="no", status=0)
}

if (!is.null(opt$name_prefix)){
  keep <- startsWith(as.character(g$Name), opt$name_prefix)
  g <- g[!is.na(g$Name) & keep]
  if (length(g) == 0){
    message("No features with Name prefix '", opt$name_prefix, "' in ", opt$input,
            " — nothing to split")
    quit(save="no", status=0)
  }
}

groups <- split(g, g$Name)
for (n in names(groups)){
  if (is.na(n) || n == "") next
  # sanitise the Name for use as a filename (TRC_<n> needs no change)
  n_clean <- gsub("[/ ]", "_", n)
  export(groups[[n]], file.path(opt$output_dir, paste0(n_clean, ".gff3")),
         format="gff3")
}
