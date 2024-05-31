#!/usr/bin/env Rscript
library(optparse)

# get input arguments

option_list <- list(
  make_option(c("-r", "--rm"), type="character", default=NULL, help="repeats_NoSat GFF3"),
  make_option(c("-s", "--sat"), type="character", default=NULL, help="repeats_Sat GFF3"),
  make_option(c("-g", "--genome"), type="character", default=NULL, help="Genome FASTA"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output csv table"),
  make_option(c("-d", "--dir"), type="character", default="repeats", help="Output directory")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

suppressPackageStartupMessages({
  library(rtracklayer)
  library(Biostrings)
})

g <- readDNAStringSet(opt$genome)
rm <- import(opt$rm)
sat <- import(opt$sat)


genome_size <- sum(width(g))
# calculate size of intervals, split by Name attribute
SR <- grepl("Simple_repeat", rm$Name)
if (any(SR)){
  rm$Name[SR] <- "Simple_repeat"
}
group_size <- sapply(split(rm, ~Name), function(x)sum(width(x)))
group_size <- c(group_size, Satellites=sum(width(sat)))


group_perc <- group_size/genome_size*100

out <- data.frame("Repeat_type"=names(group_size),
                  "Total_size [bp]"=group_size,
                  "Proportion [%]"=group_perc,
                  check.names=FALSE)



# write groups to separate files
dir.create(opt$dir, showWarnings = FALSE)
groups <-  split(rm, ~Name)
for (n in names(groups)){
  # clean name
  n_clean <- gsub("/", ".", n)
  export(groups[[n]], paste0(opt$dir,"/",n_clean,".gff3"), format="gff3")
}

write.table(out, file=opt$output, row.names=FALSE, sep="\t")
