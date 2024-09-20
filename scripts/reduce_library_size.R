#!/usr/bin/env Rscript
library(optparse)
# input - library as FASTA
# output - library as FASTA with reduced size
# require cap3 to be installed and blastn in the PATH

## FUNCTIONS
analyze_blast <- function (blfile){
  # some blast file could be empty
  if (file.size(blfile) == 0){
    return(0)
  }
  bl <- read.table(blfile, header = FALSE, sep = "\t", as.is = TRUE, comment.char = "")
  colnames(bl) <- c("qseqid", "sseqid", "pident", "length", "qstart", "qend", "evalue", "bitscore", "qlen", "slen", "qcovs")
  bl_parts <- split(bl, bl$qseqid)
  qcov <- sapply(bl_parts, calculate_total_coverage)
  qcov
}

calculate_total_coverage <- function (bl_table){
  # for each line it take start and end of the alignment to get covered region
  # then it merges the regions and calculate the total length of covered region
  # do not count overlaps!!!
  # the total length is divided by the length of the query sequence
  region <- rep(FALSE, bl_table$qlen[1])
  for (i in 1:nrow(bl_table)){
    # minimum length of region is 50 bp
    if (bl_table$length[i] < 50){
      next
    }
    region[bl_table$qstart[i]:bl_table$qend[i]] <- TRUE
  }
  if (FALSE){
    # just for testing
    region_segments <- rle(region)
    # get starts and end of the regions, where region is FALSE
    boundaries <- c(1, cumsum(region_segments$lengths) - 1)
    starts <- boundaries[1:(length(boundaries) - 1)]
    ends <- boundaries[2:length(boundaries)] + 1
    qcov_df <- data.frame(start = starts, end = ends, region = region_segments$values)
  }
  qcov <- sum(region) / bl_table$qlen[1]
  qcov
}



opt_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Input library as FASTA", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output library as FASTA with reduced size", metavar="character"),
  make_option(c("-t", "--threads"), type="numeric", default=4, help="number of threads", metavar="numeric"),
  make_option(c("-d", "--directory"), type="character", default=NULL, help="directory for temporary files", metavar="character")
)

opt <- parse_args(OptionParser(option_list = opt_list))

# input and output are mandatory
if (is.null(opt$input) | is.null(opt$output)){
  stop("Both input and output are mandatory")
}

if (is.null(opt$directory)){
  opt$directory <- tempdir()
}
dir.create(opt$directory)
# for testing
if (FALSE){
  opt <- list(
    input = "/mnt/raid/users/petr/workspace/assembly_repeat_annotation_pipeline/output_abblast/Libraries/combined_library.fasta",
    output = "/mnt/raid/users/petr/workspace/assembly_repeat_annotation_pipeline/output_abblast/Libraries/test5.fasta",
    threads = 4,
    directory = "/mnt/raid/users/petr/workspace/assembly_repeat_annotation_pipeline/output_abblast/Libraries/tmp_dir"
  )
}
cap3 <- "/mnt/raid/opt/tgicl_novy-cap/bin/cap3"

suppressPackageStartupMessages({
  library(Biostrings)
  library(parallel)
})


s <- readDNAStringSet(opt$input)
classification <- gsub(".+#", "", names(s))
unique_classification <- unique(classification)
# split by classification
s_split <- split(s, classification)

message("Library loaded:")
print(s)
message("Classification of sequences in the library:")
print(unique_classification)



dirs <- paste0(opt$directory, "/", seq_along(s_split))
tm <- sapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
# write FASTA files to the directories
for (i in seq_along(s_split)){
  writeXStringSet(s_split[[i]], file.path(dirs[i], "input.fasta"))
}
input_fasta <- paste0(dirs, "/input.fasta")
output_contigs <- paste0(dirs, "/input.fasta.cap.contigs")
output_singlets <- paste0(dirs, "/input.fasta.cap.singlets")
output_aln <- paste0(dirs, "/input.fasta.cap.aln")
blast_output <- paste0(dirs, "/input.fasta.cap.singlets.blastn")
# run cap3 in parallel


cmds <- paste0(cap3, " ", input_fasta, " -p 80 -o 50 > ", output_aln) # 0.72783 reduction
#cmds <- paste0(cap3, " ", input_fasta, " > ", output_aln) # 0.846

message("Running cap3")
tm <- mclapply(cmds, function(x) system(x, intern = TRUE), mc.cores = opt$threads)
message("Cap3 finished")

#  read the contigs and singletons to check the size
contigs <- lapply(output_contigs, readDNAStringSet)
singlets <- lapply(output_singlets, readDNAStringSet)
names(singlets) <- output_singlets
size_contigs <- sapply(contigs, function(x) sum(nchar(x)))
size_singlets <- sapply(singlets, function(x) sum(nchar(x)))
size_input <- sapply(s_split, function(x) sum(nchar(x)))
# ratio of contigs  + singletons to input
ratio <- (size_contigs + size_singlets) / size_input

total_ratio <- sum(size_contigs + size_singlets) / sum(size_input)

run_blast <- size_contigs > 0 &  size_singlets > 0
if (any(run_blast)){
  # run blastn - singletons against contigs
  cmd_create_db <- paste0("makeblastdb -in ", output_contigs[run_blast], " -dbtype nucl")
  tm <- mclapply(cmd_create_db, function(x) system(x, intern = TRUE), mc.cores = opt$threads)
  # run blastn
  cmd_blastn <- paste0("blastn -query ", output_singlets[run_blast],
                       " -db ", output_contigs[run_blast],
                       " -outfmt '6 qseqid sseqid pident length qstart qend evalue bitscore qlen slen qcovs'",
                       " -evalue 1e-20 -perc_identity 80 -word_size 9 -max_target_seqs 20",
                       " -gapextend 1 -gapopen 2 -reward 1 -penalty -1 ",
                       " -num_threads ", opt$threads,
                       " -out ", blast_output[run_blast])
  message("Running blastn")
  lapply(cmd_blastn, function(x) system(x, intern = TRUE))
  message("Blastn finished")
  qcov <- lapply(blast_output[run_blast], analyze_blast)
  names(qcov) <- output_singlets[run_blast]
  message("parsing blastn results")
  for (i in seq_along(qcov)){
    if (any(qcov[[i]]> 0.98)){
      seq_id <- names(qcov[[i]][qcov[[i]] > 0.95])
      s_name <- gsub(".blastn$", "", names(qcov)[i])
      singlets_filtered <- singlets[[s_name]][!names(singlets[[s_name]]) %in% seq_id]
      s_name_filtered <- paste0(s_name, "_filtered")
      writeXStringSet(singlets_filtered, s_name_filtered)
    }
  }
  message("Exporting library")
  for (i in seq_along(contigs)){
    if (length(contigs[[i]]) > 0){
      names(contigs[[i]]) <- paste0(names(contigs[[i]]), "#",classification[i])
    }
  }
  new_singlets_files <- output_singlets
  i <- 0
  for (n in names(singlets)){
    i <- i + 1
    if (file.exists(paste0(n, "_filtered"))){
      new_singlets_files[i] <- paste0(n, "_filtered")
    }
  }
  singlets2 <- sapply(new_singlets_files, readDNAStringSet)

  # unlist singlets2 and contig and save to output
  out <- c( do.call(c, contigs), do.call(c, unname(singlets2)))
  writeXStringSet(out, opt$output)



}else{
  message("Skipping blastn step")
  out <- c( do.call(c, contigs), do.call(c, unname(singlets)))
  message("Exporting library")
  writeXStringSet(out, opt$output)
}

output_size <- sum(nchar(out))
total_reduction <- round(output_size / sum(size_input) * 100,2)
message("-----------------------------------------------")
message("Input library size: ", sum(size_input), " bp")
message("Output library size: ", output_size, " bp")
message("Reduction: ", total_reduction, "%")
message("-----------------------------------------------")

save.image("tmp.RData")

