#!/usr/bin/env Rscript
library(optparse)
# input - library as FASTA
# output - library as FASTA with reduced size
# require cap3, mmseqs and blastn to be installed and in the PATH

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
  for (i in seq_len(nrow(bl_table))){
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
dir.create(opt$directory, recursive = TRUE, showWarnings = FALSE)
cap3='cap3'
# for testing
if (FALSE){
  opt <- list(
    input = "/mnt/raid/users/petr/workspace/assembly_repeat_annotation_pipeline/output_abblast/Libraries/combined_library.fasta",
    output = "/mnt/raid/users/petr/workspace/assembly_repeat_annotation_pipeline/output_abblast/Libraries/test5.fasta",
    threads = 4,
    directory = "/mnt/raid/users/petr/workspace/assembly_repeat_annotation_pipeline/output_abblast/Libraries/tmp_dir"
  )
  cap3 <- "/mnt/raid/opt/tgicl_novy-cap/bin/cap3"
}

suppressPackageStartupMessages({
  library(Biostrings)
  library(parallel)
})

read_dna_or_empty <- function(path){
  if (!file.exists(path) || file.size(path) == 0){
    return(DNAStringSet())
  }
  readDNAStringSet(path)
}

should_use_cap3 <- function(classification_name){
  startsWith(classification_name, "Class_I/LTR")
}

run_cap3 <- function(input_fasta, output_aln, cap3_bin){
  status <- system2(
    cap3_bin,
    args = c(input_fasta, "-p", "80", "-o", "50"),
    stdout = output_aln,
    stderr = output_aln
  )
  if (is.null(status)){
    status <- 0L
  }
  as.integer(status)
}

run_mmseqs <- function(input_fasta, output_prefix, tmp_dir, log_file){
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  status <- system2(
    "mmseqs",
    args = c("easy-cluster", input_fasta, output_prefix, tmp_dir, "--threads", "1"),
    stdout = log_file,
    stderr = log_file
  )
  if (is.null(status)){
    status <- 0L
  }
  status <- as.integer(status)
  rep_seq <- paste0(output_prefix, "_rep_seq.fasta")
  list(
    status = status,
    rep_seq = rep_seq,
    success = status == 0L && file.exists(rep_seq)
  )
}

process_classification_group <- function(index, classification_name, input_fasta, output_contig,
                                         output_singlet, output_aln, cap3_bin){
  mmseqs_prefix <- file.path(dirname(input_fasta), "mmseqs_cluster")
  mmseqs_tmp <- file.path(dirname(input_fasta), "mmseqs_tmp")
  mmseqs_log <- file.path(dirname(input_fasta), "mmseqs.log")
  stale_mmseqs_files <- list.files(dirname(input_fasta), pattern = "^mmseqs_cluster",
                                   full.names = TRUE)

  unlink(c(output_contig, output_singlet, output_aln), force = TRUE)
  if (length(stale_mmseqs_files) > 0){
    unlink(stale_mmseqs_files, recursive = TRUE, force = TRUE)
  }
  unlink(c(mmseqs_tmp, mmseqs_log), recursive = TRUE, force = TRUE)

  if (should_use_cap3(classification_name)){
    message("[", classification_name, "] Using CAP3")
    cap3_status <- run_cap3(
      input_fasta = input_fasta,
      output_aln = output_aln,
      cap3_bin = cap3_bin
    )
    if (cap3_status != 0L){
      stop("[", classification_name, "] CAP3 failed with status ", cap3_status)
    }
    return(list(
      index = index,
      classification = classification_name,
      mode = "cap3",
      reason = "classification matches Class_I/LTR",
      contig_file = output_contig,
      singlet_file = output_singlet,
      mmseqs_rep_file = NA_character_
    ))
  }

  message("[", classification_name, "] Using MMseqs2 easy-cluster")
  mmseqs_result <- run_mmseqs(
    input_fasta = input_fasta,
    output_prefix = mmseqs_prefix,
    tmp_dir = mmseqs_tmp,
    log_file = mmseqs_log
  )
  if (!mmseqs_result$success){
    stop("[", classification_name, "] MMseqs2 failed with status ",
         mmseqs_result$status)
  }
  list(
    index = index,
    classification = classification_name,
    mode = "mmseqs",
    reason = "classification does not match Class_I/LTR",
    contig_file = output_contig,
    singlet_file = output_singlet,
    mmseqs_rep_file = mmseqs_result$rep_seq
  )
}


s <- readDNAStringSet(opt$input)
classification <- gsub(".+#", "", names(s))
# split by classification
s_split <- split(s, classification)
unique_classification <- names(s_split)

message("Library loaded-:")
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
unlink(c(blast_output, paste0(output_singlets, "_filtered")), force = TRUE)
# run CAP3 only for LTR classes; all other classes use MMseqs2
message("Running per-class reduction")
group_results <- mclapply(
  seq_along(s_split),
  function(i){
    process_classification_group(
      index = i,
      classification_name = unique_classification[i],
      input_fasta = input_fasta[i],
      output_contig = output_contigs[i],
      output_singlet = output_singlets[i],
      output_aln = output_aln[i],
      cap3_bin = cap3
    )
  },
  mc.cores = opt$threads
)
message("Reduction step finished")
message("Reduction mode summary:")
for (result in group_results){
  message("  ", result$classification, ": ", result$mode, " (", result$reason, ")")
}

# read CAP3 outputs for LTR classes and MMseqs representatives for all other classes
group_modes <- vapply(group_results, function(x) x$mode, character(1))
contigs <- lapply(seq_along(output_contigs), function(i){
  if (group_modes[i] != "cap3"){
    return(DNAStringSet())
  }
  read_dna_or_empty(output_contigs[i])
})
singlets <- lapply(seq_along(output_singlets), function(i){
  if (group_modes[i] != "cap3"){
    return(DNAStringSet())
  }
  read_dna_or_empty(output_singlets[i])
})
names(singlets) <- output_singlets
mmseqs_representatives <- lapply(seq_along(group_results), function(i){
  if (group_results[[i]]$mode != "mmseqs"){
    return(DNAStringSet())
  }
  read_dna_or_empty(group_results[[i]]$mmseqs_rep_file)
})
size_contigs <- sapply(contigs, function(x) sum(nchar(x)))
size_singlets <- sapply(singlets, function(x) sum(nchar(x)))
size_mmseqs <- sapply(mmseqs_representatives, function(x) sum(nchar(x)))
size_input <- sapply(s_split, function(x) sum(nchar(x)))
size_output_by_group <- ifelse(group_modes == "mmseqs", size_mmseqs, size_contigs + size_singlets)
ratio <- size_output_by_group / size_input
total_ratio <- sum(size_output_by_group) / sum(size_input)

run_blast <- size_contigs > 0 & size_singlets > 0 & group_modes == "cap3"
if (any(run_blast)){
  # run blastn - singletons against contigs
  cmd_create_db <- paste0("makeblastdb -in ", output_contigs[run_blast], " -dbtype nucl")
  tm <- mclapply(cmd_create_db, function(x) system(x, intern = TRUE), mc.cores = opt$threads)
  # run blastn
  cmd_blastn <- paste0("blastn -query ", output_singlets[run_blast],
                       " -db ", output_contigs[run_blast],
                       " -outfmt '6 qseqid sseqid pident length qstart qend evalue bitscore qlen slen qcovs'",
                       " -evalue 1e-20 -perc_identity 95 -word_size 9 -max_target_seqs 20",
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
    if (length(qcov[[i]]) > 0 && any(qcov[[i]] > 0.98)){
      seq_id <- names(qcov[[i]][qcov[[i]] > 0.98])
      s_name <- gsub(".blastn$", "", names(qcov)[i])
      singlets_filtered <- singlets[[s_name]][!names(singlets[[s_name]]) %in% seq_id]
      s_name_filtered <- paste0(s_name, "_filtered")
      writeXStringSet(singlets_filtered, s_name_filtered)
    }
  }
}else{
  message("Skipping blastn step")
}

message("Exporting library")
for (i in seq_along(contigs)){
  contig_suffix <- paste0("_", i)
  if (length(contigs[[i]]) > 0){
    names(contigs[[i]]) <- paste0(names(contigs[[i]]), contig_suffix, "#", unique_classification[i])
  }
}
new_singlets_files <- output_singlets
for (i in seq_along(output_singlets)){
  if (file.exists(paste0(output_singlets[i], "_filtered"))){
    new_singlets_files[i] <- paste0(output_singlets[i], "_filtered")
  }
}
singlets2 <- lapply(seq_along(new_singlets_files), function(i){
  if (group_results[[i]]$mode != "cap3"){
    return(DNAStringSet())
  }
  read_dna_or_empty(new_singlets_files[i])
})

out <- c(do.call(c, contigs), do.call(c, unname(singlets2)),
         do.call(c, unname(mmseqs_representatives)))
writeXStringSet(out, opt$output)

output_size <- sum(nchar(out))
total_reduction <- round(output_size / sum(size_input) * 100,2)
message("-----------------------------------------------")
message("Input library size: ", sum(size_input), " bp")
message("Output library size: ", output_size, " bp")
message("Reduction: ", total_reduction, "%")
message("-----------------------------------------------")
