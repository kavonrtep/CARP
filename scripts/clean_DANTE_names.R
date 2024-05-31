#!/usr/bin/env Rscript
suppressPackageStartupMessages(
  library(rtracklayer)
)
g <- import(commandArgs(T)[1])
new_name <- gsub("|", "/",
                 gsub("/", "_",
                      g$Final_Classification, fixed = TRUE
                 ),
                 fixed = TRUE
)
g$Name <- new_name
export(g, format = 'gff3', commandArgs(T)[2])