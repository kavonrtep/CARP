#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Biostrings)
})
sl <- seqlengths(readDNAStringSet(commandArgs(TRUE)[1]))
saveRDS(sl, commandArgs(TRUE)[2])
