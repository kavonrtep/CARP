#!/usr/bin/env Rscript
library(optparse)

# Resolve this script's own directory so we can source the shared
# density helper that sits beside it. Works whether the script is run
# from the repo checkout or from /opt/pipeline/scripts in the container,
# because --file= always carries the resolved path to this script.
.density_script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
source(file.path(.density_script_dir, "density_utils.R"))

# get input arguments
# bed file
# window size
# output bigwig file
option_list <- list(
  make_option(c("-b", "--bed"), type="character", default=NULL, help="BED or GFF file"),
  make_option(c("-w", "--window"), type="integer", default=1000000, help="Window size"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output BigWig file"),
  make_option(c("-f", "--format"), type="character", default="gff3", help="Input format (gff3 or bed)"),
  make_option(c("-m", "--merge"), type="logical", action="store_true", default=FALSE, help="[deprecated/no-op] overlaps are now ALWAYS merged into a union before density (kept for backward-compatible invocation)"),
  make_option(c("-g", "--genome"), type="character", default=NULL, help="Genome file in fasta format")


)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# check mandatory arguments
if (is.null(opt$bed) || is.null(opt$output)){
  stop("Please provide bed file and output file")
}

suppressPackageStartupMessages({
  library(rtracklayer)
})
g <- import(opt$bed, format=opt$format)

# check if gff file is not empty
if (length(g)==0){
  # exit normally - create empty bigwig file
  print("No regions found in the input file")
  write.table(data.frame(), file=opt$output, quote=FALSE, sep="\t", row.names=FALSE)
  quit()
}

print(opt)
chr_size_all <- readRDS(opt$genome)

window_size <- opt$window/10 # 10 bins per window

# density_track():
#  - union coverage (overlaps reduced away, so a score is a fraction in [0,1] —
#    the unified annotation deliberately tolerates overlap and feeding that
#    straight to coverage() counts stacking depth, which pushed this total track
#    to ~3.5x; the legacy --merge flag is now implied);
#  - keep_unoccupied = TRUE: this is the genome-wide TOTAL track, so scaffolds
#    with no annotation stay in the output as zero-score intervals;
#  - FR-1 run-length merge (adjacent equal-value tiles, incl. zero runs,
#    collapsed into one interval) — lossless, non-zero values unchanged at
#    every position.
d <- density_track(g, chr_size_all, window_size, N_for_mean = 10,
                   keep_unoccupied = TRUE)

export(d, opt$output, format="bigwig")
