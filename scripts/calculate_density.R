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

get_density <- function(x, chr_size=NULL, tw=1000000){
  cvg <- coverage(x)
  bins <- tileGenome(chr_size, tilewidth = tw)
  d <- binnedAverage(unlist(bins), cvg, "score")
  d
}
max_chr_length <- function(g){
  x <- split(g, ~seqnames)
  L <- sapply(x, function(x)max(end(x), na.rm=TRUE))
  L
}


get_density2 <- function(x, chr_size=NULL, N_for_mean = 10,step_size=100000){
  cvg <- coverage(x)
  bins <- tileGenome(chr_size, tilewidth = step_size)
  d <- binnedAverage(unlist(bins), cvg, "score")
  s_part <- split(d$score, seqnames(d))
  s_part_smooth <- lapply(s_part, function(x)smooth_score2(x, N_for_mean))
  d$score <- unlist(s_part_smooth)
  d
}


smooth_score <- function(x, N_for_mean = 10){
  # extend the score in each direction by N_for_mean-1 zeros
  sc <- c(rep(0, N_for_mean-1), x$score, rep(0, N_for_mean-1))
  sc_smooth <- filter(sc, rep(1/N_for_mean, N_for_mean), sides=2)
  # remove the first N_for_mean-1 and the last N_for_mean-1 elements
  x$score <- sc_smooth[(N_for_mean):(length(sc_smooth)-N_for_mean+1)]
  x
}

smooth_score2 <- function(x, N_for_mean = 10){
  # extend the score in each direction by N_for_mean-1 zeros
  sc <- c(rep(0, N_for_mean-1), x, rep(0, N_for_mean-1))
  sc_smooth <- filter(sc, rep(1/N_for_mean, N_for_mean), sides=2)
  # remove the first N_for_mean-1 and the last N_for_mean-1 elements
  out <- sc_smooth[(N_for_mean):(length(sc_smooth)-N_for_mean+1)]
  out
}


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

# Density tracks are UNION coverage: every base is either repeat-covered or
# not, so a score must never exceed 1.0. The unified annotation deliberately
# tolerates overlap — L1 Simple_repeat / Low_complexity sitting on top of a TE,
# and ALL L2 nested children (tandem-array members inside an LTR_RT_TR
# container, simple repeats nested in a satellite). Feeding those straight to
# coverage() counts stacking depth (>1) and pushed the total track to ~3.5x.
# reduce() to a non-overlapping union first. ignore.strand=TRUE is essential: a
# '*'-strand simple repeat must merge with the +/- element it overlaps. The
# legacy --merge flag is now implied (always on) and kept only so existing
# invocations that pass it keep working.
g <- reduce(g, ignore.strand = TRUE)
print(opt)
chr_size_all <- readRDS(opt$genome)

# add missing seqlevels to g
chr_size <- chr_size_all[seqlevels(g)]
not_used <- setdiff(names(chr_size_all), names(chr_size))
chr_size_not_used <- chr_size_all[not_used]
seqlevels(g) <- c(seqlevels(g), names(chr_size_not_used))
chr_size_in_order <- chr_size_all[seqlevels(g)]


window_size <- opt$window/10 # 10 bins per window
d <- get_density2(g, chr_size_in_order, N_for_mean = 10, step_size = window_size)

# FR-1: write the BigWig run-length-merged (adjacent equal-value tiles,
# incl. zero runs, collapsed into one interval) instead of one entry per
# window. Lossless; non-zero values unchanged at every position.
export(rle_merge_granges(d), opt$output, format="bigwig")
