#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(parallel))

gff_cleanup <- function(gff){
  ## remove overlapin annotation track - assign new annot
  gff_disjoin <- disjoin(gff, with.revmap=TRUE)
  ## append annotation:
  gff_names <- mclapply(as.list(gff_disjoin$revmap), FUN = function(x)gff$Name[x], mc.cores = 8)
  gff_strands <- mclapply(as.list(gff_disjoin$revmap), FUN = function(x)strand(gff[x]), mc.cores = 8)
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
                gsub(new_sep,"-", n, fixed = TRUE),
                fixed = TRUE)
  return(n_new)
}


infile <- commandArgs(T)[1]
outfile <- commandArgs(T)[2]


## infile = "./test_data/raw_rm.out"

rm_out <- read.table(infile, as.is=TRUE, sep="", skip = 2, fill=TRUE, header=FALSE, col.names=paste0("V", 1:16))

gff <- GRanges(seqnames = rm_out$V5, ranges = IRanges(start = rm_out$V6, end=rm_out$V7))

# repeat class after # symbol - syntax 1
# detect separator
# if "|" is present replace "|" -> "/" and "/" -> "-"
if (any(grepl("|", rm_out$V11, fixed = TRUE))){
  gff$Name <- convert_names(rm_out$V11, old_sep = "|", new_sep = "/")
  message('replacing classification separator character "|" with "/"')
  print(gff)
}else{
  gff$Name <- rm_out$V11
}
# if rm_out$V11 is Simple_repeat, add column V10 to Name
new_name <- ifelse(rm_out$V11 == "Simple_repeat", paste0(rm_out$V11,"", rm_out$V10), rm_out$V11)
gff$Name <- new_name


## is repeat type is specifies by double underscore:
## then rm_out$V11 is unspecified
if (any(rm_out$V11 == "Unspecified")){
  ## set Name from prefix
  ## TODO
  inc <- rm_out$V11 == "Unspecified"
  Name <- gsub("__.+", "", rm_out$V10)
  # chanche Usnpsecified to new name
  gff$Name[inc] <- Name
}


## join neighbors with the same annotation, disregard strand!
result <- unlist(reduce(split(gff, gff$Name)))

result$Name <- names(result)

result_clean <- gff_cleanup(result)

## TODO
## identify conflicting annotation, replace by LCA but keep origin list of classifications

gff_out <- sortSeqlevels(result_clean)
gff_out <- sort(gff_out)
gff_out$type <- "repeat_region"
gff_out$source <- "RepeatMasker_parsed"
gff_out$ID <- paste0(gff_out$Name, "_", seq_along(gff_out$Name))
export(gff_out,  format = "gff3", con=outfile)


