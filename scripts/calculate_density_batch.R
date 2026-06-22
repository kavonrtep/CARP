#!/usr/bin/env Rscript
library(optparse)
suppressPackageStartupMessages({
  library(rtracklayer)
})

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


max_chr_length <- function(g){
  x <- split(g, ~seqnames)
  L <- sapply(x, function(x)max(end(x), na.rm=TRUE))
  L
}

option_list <- list(
  make_option(c("-d", "--dir"), type="character", default=NULL, help="Directory of GFF3 files"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, help="Output directory for BigWig files"),
  make_option(c("-g", "--genome"), type="character", default=NULL, help="Genome file in fasta format")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

files <- list.files(opt$dir, pattern="\\.gff3$", full.names=TRUE)


directory_for_10k <- paste0(opt$output_dir, "/10k")
directory_for_100k <- paste0(opt$output_dir, "/100k")

dir.create(directory_for_10k, showWarnings=FALSE, recursive=TRUE)
dir.create(directory_for_100k, showWarnings=FALSE, recursive=TRUE)
print("genome:")
print(opt$genome)
chr_size_all <- readRDS(opt$genome)
# exclude seqlevels not included in g



for (f in files) {
  base <- basename(f)


  base_noext <- sub("\\.gff3$", "", base)

  base_bw10k <- paste0(directory_for_10k, "/", base_noext, "_10k.bw")
  base_bw100k <- paste0(directory_for_100k, "/", base_noext, "_100k.bw")
  g <- import(f, format="gff3")

  if (length(g)==0){
    print(paste("No regions found in the input file:", f))
    next
  }

  # add missing seqlevels to g
  chr_size <- chr_size_all[seqlevels(g)]
  not_used <- setdiff(names(chr_size_all), names(chr_size))
  chr_size_not_used <- chr_size_all[not_used]
  seqlevels(g) <- c(seqlevels(g), names(chr_size_not_used))
  chr_size_in_order <- chr_size_all[seqlevels(g)]

  # FR-1: write run-length-merged BigWigs (adjacent equal-value tiles,
  # incl. zero runs, collapsed into one interval) instead of one entry
  # per window. Lossless; non-zero values unchanged at every position.
  #d10k <- get_density(g, chr_size, 10000)
  d10k_smooth <- get_density2(g, chr_size_in_order, N_for_mean = 10, step_size=1000)
  export(rle_merge_granges(d10k_smooth), base_bw10k, format="bigwig")
  #d100k <- get_density(g, chr_size, 100000)
  d100k_smooth <- get_density2(g, chr_size_in_order, N_for_mean = 10, step_size=10000)
  export(rle_merge_granges(d100k_smooth), base_bw100k, format="bigwig")
}

