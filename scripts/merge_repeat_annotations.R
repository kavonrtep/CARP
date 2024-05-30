#!/usr/bin/env Rscript
resolve_name <- function(x){
  if (length(x)==1){
    # no conflict
    return(x)
  } else{
    y <- sapply(x, strsplit, split="/", fixed = TRUE)
    ny <- table(unlist(sapply(y, function(x)paste(seq_along(x), x))))
    if (max(ny)<length(x)){
      return("Unknown")
    }else{
      k <- which(ny==length(x))
      r <- max(as.numeric((gsub(" .+", "", names(k)))))
      out <- paste(y[[1]][1:r], collapse="/")
      return(out)
    }
  }
}

convert_names <- function(n, old_sep = "|" , new_sep = "\""){
  # remove all characters which are new_sep with -
  n_new <- gsub(old_sep, new_sep,
                gsub(new_sep,"_", n, fixed = TRUE),
                fixed = TRUE)
  return(n_new)
}

gff_cleanup <- function(gff){
  ## remove overlapin annotation track - assign new annot
  gff_disjoin <- disjoin(gff, with.revmap=TRUE)
  ## append annotation:
  # get number of cores
  num_cores <- detectCores()
  gff_names <- mclapply(as.list(gff_disjoin$revmap), FUN = function(x)gff$Name[x], mc.cores = round(num_cores *0.8))
  gff_strands <- mclapply(as.list(gff_disjoin$revmap), FUN = function(x)strand(gff[x]), mc.cores = round(num_cores *0.8))
  new_annot <- sapply(sapply(gff_names, unique), paste, collapse="|")
  new_annot_uniq <- unique(new_annot)
  lca_annot <- sapply(strsplit(new_annot_uniq, "|", fixed = TRUE), resolve_name)
  names(lca_annot) <- new_annot_uniq
  new_annot_lca <- lca_annot[new_annot]
  #new_annot_lca = sapply(sapply(gff_names, unique), resolve_name)
  strand_attribute <- sapply(sapply(gff_strands, unique), paste, collapse="|")
  gff_disjoin$source <- "RM"
  gff_disjoin$type <- "repeat"
  gff_disjoin$score <- NA
  gff_disjoin$phase <- NA
  gff_disjoin$Name <- new_annot_lca
  gff_disjoin$Original_names <- new_annot
  gff_disjoin$strands <- strand_attribute
  gff_disjoin$revmap <- NULL
  return(gff_disjoin)
}

N <- length(commandArgs(TRUE))
if (N < 3) {
  message("expecting at least 2 input files and 1 output file")
  stop("Usage: merge_repeat_annotations.R <infile1> <infile2> ... <outfile>")
}
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(parallel))



# output is last one
out_gff <- commandArgs(TRUE)[N]
# input is all but last
in_gff <- commandArgs(TRUE)[-N]

# concatenate all input files by cat and read in as gff
tmp_gff <- tempfile()
# purpose of concatenation is get everithing into one file to avoid mising seqlevels in individual files
system(paste("cat", paste(in_gff, collapse = " "), ">", tmp_gff))
gff <- import.gff(tmp_gff)
unlink(tmp_gff)
# set strand for all records to *
strand(gff) <- "*"

if (any(grepl("|", gff$Name, fixed = TRUE))){
  to_clean <- grepl("|", gff$Name, fixed = TRUE)
  gff1 <- gff[!to_clean,]
  gff2 <- gff[to_clean,]
  message('replacing c lassification separator character "|" with "/"')
  gff2$Name <- convert_names(gff2$Name, old_sep = "|", new_sep = "/")
  gff <- c(gff1, gff2)
}

result <- unlist(reduce(split(gff, gff$Name)))
result$Name <- names(result)
result_clean <-  gff_cleanup(result)

gff_out <-  sortSeqlevels(result_clean)
gff_out <- sort(gff_out)
gff_out$type <- "repeat_region"
gff_out$source <- "RepeatMasker_parsed"
gff_out$ID <- paste0(gff_out$Name, "_", seq_along(gff_out$Name))
export(gff_out,  format = "gff3", con=out_gff)

