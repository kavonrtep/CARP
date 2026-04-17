#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(parallel))

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "classification.R"))

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

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]

rm_out <- read.table(infile, as.is=TRUE, sep="", skip = 2, fill=TRUE, header=FALSE, col.names=paste0("V", 1:16))

gff <- GRanges(seqnames = rm_out$V5, ranges = IRanges(start = rm_out$V6, end=rm_out$V7))

# Canonicalise V11 separators (pipe-encoded libraries become slash-encoded);
# validation is off because this is raw RepeatMasker output and may contain
# entries like "Simple_repeat" that we decorate with V10 on the next line.
pipe_source <- any(grepl("|", rm_out$V11, fixed = TRUE))
if (pipe_source) message('normalising classification separator "|" -> "/"')
v11_norm <- canonicalise(rm_out$V11,
                         source = if (pipe_source) "DANTE" else "RepeatMasker",
                         validate = FALSE)
# Simple_repeat entries carry the monomer sequence in V10; RepeatMasker
# consumers expect that motif in Name (e.g. "Simple_repeat(TAAACCC)n").
gff$Name <- ifelse(v11_norm == "Simple_repeat",
                   paste0(v11_norm, rm_out$V10),
                   v11_norm)


## is repeat type is specifies by double underscore:
## then rm_out$V11 is unspecified
if (any(rm_out$V11 == "Unspecified")){
  ## set Name from prefix
  ## TODO
  inc <- rm_out$V11 == "Unspecified"
  Name <- gsub("__.+", "", rm_out$V10[inc])
  # change Unspecified to new name
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


