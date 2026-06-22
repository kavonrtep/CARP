#!/usr/bin/env Rscript
# Split a GFF3 into one file per distinct `Name` attribute.
#
# Used to produce per-family inputs for density BigWig generation: the
# RepeatMasker tandem pass (Tandem_repeats_RepeatMasker.gff3) carries
# Name=TRC_<n>, and feeding the per-family GFF3s to calculate_density_batch.R
# yields one density BigWig per family — mirroring the per-cluster set the
# TideCluster pass already produces. See docs/carp_feature_requests.md FR-2b.

suppressPackageStartupMessages({
  library(optparse)
  library(rtracklayer)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input GFF3"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL,
              help="Output directory (one <Name>.gff3 per distinct Name)")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input) || is.null(opt$output_dir)){
  stop("Please provide --input and --output_dir")
}

dir.create(opt$output_dir, showWarnings=FALSE, recursive=TRUE)

g <- import(opt$input, format="gff3")

if (length(g) == 0 || is.null(g$Name)){
  message("No features (or no Name attribute) in ", opt$input, " — nothing to split")
  quit(save="no", status=0)
}

groups <- split(g, g$Name)
for (n in names(groups)){
  if (is.na(n) || n == "") next
  # sanitise the Name for use as a filename (TRC_<n> needs no change)
  n_clean <- gsub("[/ ]", "_", n)
  export(groups[[n]], file.path(opt$output_dir, paste0(n_clean, ".gff3")),
         format="gff3")
}
