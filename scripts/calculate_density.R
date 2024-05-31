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

# get input arguments
# bed file
# window size
# output bigwig file
option_list <- list(
  make_option(c("-b", "--bed"), type="character", default=NULL, help="BED or GFF file"),
  make_option(c("-w", "--window"), type="integer", default=1000000, help="Window size"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output BigWig file"),
  make_option(c("-f", "--format"), type="character", default="gff3", help="Input format (gff3 or bed)"),
  make_option(c("-m", "--merge"), type="logical", action="store_true", default=FALSE, help="Merge overlapping regions")


)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# check mandatory arguments
if (is.null(opt$bed) || is.null(opt$output)){
  stop("Please provide bed file and output file")
}

suppressPackageStartupMessages(library(rtracklayer))
g <- import(opt$bed, format=opt$format)
if (opt$merge){
  g <- reduce(g)
}
chr_size <- max_chr_length(g)
d <- get_density(g, chr_size, opt$window)

export(d, opt$output, format="bigwig")
