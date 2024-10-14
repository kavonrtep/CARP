#!/usr/bin/env Rscript
library(optparse)

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




# get input arguments
# bed file
# window size
# output bigwig file
option_list <- list(
  make_option(c("-b", "--bed"), type="character", default=NULL, help="BED or GFF file"),
  make_option(c("-w", "--window"), type="integer", default=1000000, help="Window size"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output BigWig file"),
  make_option(c("-f", "--format"), type="character", default="gff3", help="Input format (gff3 or bed)"),
  make_option(c("-m", "--merge"), type="logical", action="store_true", default=FALSE, help="Merge overlapping regions"),
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

if (opt$merge){
  g <- reduce(g)
}
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

export(d, opt$output, format="bigwig")
