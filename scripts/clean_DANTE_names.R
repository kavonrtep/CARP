#!/usr/bin/env Rscript
suppressPackageStartupMessages(
  library(rtracklayer)
)

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "classification.R"))

args <- commandArgs(trailingOnly = TRUE)
g <- import(args[1])
g$Name <- canonicalise(g$Final_Classification, source = "DANTE")
export(g, format = "gff3", args[2])
