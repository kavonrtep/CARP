#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(rtracklayer)
  library(Biostrings)
  library(optparse)
})

# get input arguments
# inputs:
# dante gff file
# fasta file
# output - filtered cleaned fasta:

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Input GFF3 file from DANTE"),
  make_option(c("-f", "--fasta"), type="character", default=NULL, help="Input genome FASTA file"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output filtered FASTA file"),
)


# good LINE elements criteria:
# sequence name must contain "#LINE"
# presence of at least ENDO and RT domains or ENDO RT RH domains for LINE elements
# domains must be in corect order and orientation
# no other domain which would indicate other type of transposon


opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

dante <- import(opt$input)
lib <- readDNAStringSet(opt$fasta)
names(lib) <- gsub(" .+", "", names(lib))
# split dante by seqnames and analyze each sequence
dante_split <- split(dante, seqnames(dante))

number_of_unique_classes  <-  sapply(dante_split, function(x)length(unique(x$Final_Classification)))
dante_split <- dante_split[number_of_unique_classes==1]
# keep only those with LINE in the classification
is_line <- sapply(dante_split, function(x) x$Final_Classification[1] == "Class_I|LINE")

dante_split <- dante_split[is_line]

is_valid <- function(x){
  if (length(x) %in% 2){
    x <- x[order(start(x))]
    ok <- all(x$Name == c("ENDO", "RT")) & all(strand(x) == c("+", "+"))  || all(x$Name == c("RT", "ENDO")) & all(strand(x) == c("-", "-"))
    return(ok)
  }
  if (length(x) %in% 3){
    x <- x[order(start(x))]
    ok <- all(x$Name == c("ENDO", "RT", "RH")) & all(strand(x) == c("+", "+", "+")) || all(x$Name == c("RH", "RT", "ENDO")) & all(strand(x) == c("-", "-", "-"))
    return(ok)
  }
  return(FALSE)
}

valid_LINE <- sapply(dante_split, is_valid)

lib_valid_LINE <- lib[names(valid_LINE)[valid_LINE]]

# use compatible lables
names(lib_valid_LINE) <- paste0(gsub("#.+", "", names(lib_valid_LINE)), "#Class_I/LINE")

# export
writeXStringSet(lib_valid_LINE, filepath = opt$output)


nchar(lib_valid_LINE)