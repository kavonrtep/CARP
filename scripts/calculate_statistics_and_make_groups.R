#!/usr/bin/env Rscript
library(optparse)

# get input arguments

option_list <- list(
  make_option(c("-r", "--rm"), type="character", default=NULL, help="repeats_NoSat GFF3"),
  make_option(c("-s", "--sat"), type="character", default=NULL, help="repeats_Sat GFF3"),
  make_option(c("-g", "--genome"), type="character", default=NULL, help="Genome FASTA"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output csv table"),
  make_option(c("-d", "--dir"), type="character", default="repeats", help="Output directory"),
  make_option(c("-R", "--rDNA"), type="character", default=NULL, help="rDNA GFF3"),
  make_option(c("-S", "--simple_repeats"), type="character", default=NULL, help="Simple repeats GFF3"),
  make_option(c("-L", "--low_complexity"), type="character", default=NULL, help="Low complexity GFF3"),
  make_option(c("-M", "--mobile_elements"), type="character", default=NULL, help="Mobile elements GFF3")
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



# make additional groups:

mobile_elements <- rm[grepl("^Class", rm$Name)]
simple_repeats <- rm[grepl("^Simple_repeat", rm$Name)]
low_complexity <- rm[grepl("^Low_complexity", rm$Name)]
rDNA <- rm[grepl("^rDNA", rm$Name)]

# export but check if there are any elements, otherwise make empty gff3 file with
# header and comment
if (length(mobile_elements) > 0){
  export(mobile_elements, paste0(opt$dir,"/Mobile_elements.gff3"), format="gff3")
} else {
  writeLines(c("##gff-version 3", "##", "## No mobile elements found"), paste0(opt$dir,"/Mobile_elements.gff3"))
}

if (length(simple_repeats) > 0){
  export(simple_repeats, paste0(opt$dir,"/Simple_repeats.gff3"), format="gff3")
} else {
  writeLines(c("##gff-version 3", "##", "## No simple repeats found"), paste0(opt$dir,"/Simple_repeats.gff3"))
}

if (length(low_complexity) > 0){
  export(low_complexity, paste0(opt$dir,"/Low_complexity.gff3"), format="gff3")
} else {
  writeLines(c("##gff-version 3", "##", "## No low complexity found"), paste0(opt$dir,"/Low_complexity.gff3"))
}

if (length(rDNA) > 0){
  export(rDNA, paste0(opt$dir,"/rDNA.gff3"), format="gff3")
} else {
  writeLines(c("##gff-version 3", "##", "## No rDNA found"), paste0(opt$dir,"/rDNA.gff3"))
}

# table output is at the end - it server as checkpoint in snakemake
write.table(out, file=opt$output, row.names=FALSE, sep="\t")