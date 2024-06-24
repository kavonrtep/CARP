#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Biostrings)
})
sl <- seqlengths(readDNAStringSet(commandArgs(TRUE)[1]))
names(sl) <- gsub(" .*", "", names(sl))
saveRDS(sl, commandArgs(TRUE)[2])
