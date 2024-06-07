#!/usr/bin/env Rscript
library(optparse)
suppressPackageStartupMessages({
  library(rtracklayer)
})

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
  # calculate sliding window mean for N_for_mean bins

  # split by chromosome
  d_part <- split(d, ~seqnames)
  d_part <- lapply(d_part, function(x)smooth_score(x, N_for_mean))
  d <- unlist(GRangesList(d_part))
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

dir.create(directory_for_10k, showWarnings=FALSE)
dir.create(directory_for_100k, showWarnings=FALSE)
print("genome:")
print(opt$genome)
chr_size <- readRDS(opt$genome)

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

  #d10k <- get_density(g, chr_size, 10000)
  d10k_smooth <- get_density2(g, chr_size, N_for_mean = 10, step_size=1000)
  export(d10k_smooth, base_bw10k, format="bigwig")
  #d100k <- get_density(g, chr_size, 100000)
  d100k_smooth <- get_density2(g, chr_size, N_for_mean = 10, step_size=10000)
  export(d100k_smooth, base_bw100k, format="bigwig")
}

