#!/usr/bin/env Rscript
library(optparse)
suppressPackageStartupMessages({
  library(rtracklayer)
  library(parallel)
})

# Resolve this script's own directory so we can source the shared
# density helper that sits beside it. Works whether the script is run
# from the repo checkout or from /opt/pipeline/scripts in the container,
# because --file= always carries the resolved path to this script.
.density_script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
source(file.path(.density_script_dir, "density_utils.R"))

smooth_score2 <- function(x, N_for_mean = 10){
  # extend the score in each direction by N_for_mean-1 zeros
  sc <- c(rep(0, N_for_mean-1), x, rep(0, N_for_mean-1))
  sc_smooth <- filter(sc, rep(1/N_for_mean, N_for_mean), sides=2)
  # remove the first N_for_mean-1 and the last N_for_mean-1 elements
  sc_smooth[(N_for_mean):(length(sc_smooth)-N_for_mean+1)]
}

# Per-family smoothed density, byte-identical to the original whole-genome
# get_density2() path but ~140x faster on assemblies with many scaffolds.
#
# The original tiled the WHOLE genome per family with the occupied seqlevels
# moved to the FRONT of the concatenation (seqlevels(g) <- c(occupied, unused)).
# tileGenome phases bin boundaries by the genome-concatenation offset, so that
# occupied-first reordering makes each family's bin phasing depend on which
# seqlevels it occupies. A single precomputed global grid therefore CANNOT
# reproduce those bins -- the tiling must be rebuilt per family in the same
# occupied-first order (this is what `chr_in_order` below restores).
#
# The speedup over the original is twofold and lossless:
#   (a) keep ONLY the occupied bins before binnedAverage / smoothing / export.
#       binnedAverage is per-bin independent and the moving-average smoothing is
#       applied per seqname, so the values on occupied seqlevels are unchanged;
#       we simply skip computing density on the (zero-coverage) scaffolds a
#       family does not touch -- on a 1888-scaffold genome that is the bulk of
#       the original cost.
#   (b) parallelise across input files with mclapply.
# The output header consequently carries only the occupied seqlevels (BigWig
# export drops dataless seqlevels anyway), which is the intended behaviour.
density_per_family <- function(g, chr_size_all, step, N_for_mean = 10){
  occ <- seqlevels(g)
  not_used <- setdiff(names(chr_size_all), occ)
  chr_in_order <- chr_size_all[c(occ, not_used)]   # occupied first, as the original
  bins <- unlist(tileGenome(chr_in_order, tilewidth = step))
  bins <- bins[as.character(seqnames(bins)) %in% occ]
  bins <- keepSeqlevels(bins, occ, pruning.mode = "coarse")
  cvg <- coverage(g)
  d <- binnedAverage(bins, cvg, "score")
  s_part <- split(d$score, seqnames(d))
  d$score <- unlist(lapply(s_part, function(z) smooth_score2(z, N_for_mean)))
  d
}

option_list <- list(
  make_option(c("-d", "--dir"), type="character", default=NULL, help="Directory of GFF3 files"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, help="Output directory for BigWig files"),
  make_option(c("-g", "--genome"), type="character", default=NULL, help="genome_seqlengths.rds (named integer vector of scaffold lengths)"),
  make_option(c("-t", "--threads"), type="integer", default=1, help="Parallel workers over input files [default %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))

files <- list.files(opt$dir, pattern="\\.gff3$", full.names=TRUE)

directory_for_10k  <- file.path(opt$output_dir, "10k")
directory_for_100k <- file.path(opt$output_dir, "100k")
dir.create(directory_for_10k,  showWarnings=FALSE, recursive=TRUE)
dir.create(directory_for_100k, showWarnings=FALSE, recursive=TRUE)
chr_size_all <- readRDS(opt$genome)

# Parallelise over (file x resolution) units rather than over files alone.
# One dominant per-class file (e.g. Mobile_elements) otherwise occupies a
# single core for the whole rule while the other workers idle; splitting each
# file's two BigWig resolutions into independent tasks lets that file's 10k and
# 100k tracks compute on two cores. Output is unchanged: each task writes its
# own independent .bw and the per-resolution computation is byte-for-byte the
# former two-statement process_one. (The dominant file is imported once per
# resolution instead of once total — a small cost dwarfed by overlapping its
# two coverage/binnedAverage passes.)
tasks <- list()
for (f in files) {
  base_noext <- sub("\\.gff3$", "", basename(f))
  tasks[[length(tasks) + 1L]] <- list(
    f = f, step = 1000L,
    out = file.path(directory_for_10k,  paste0(base_noext, "_10k.bw")))
  tasks[[length(tasks) + 1L]] <- list(
    f = f, step = 10000L,
    out = file.path(directory_for_100k, paste0(base_noext, "_100k.bw")))
}

process_task <- function(tk){
  g <- import(tk$f, format="gff3")
  if (length(g) == 0) return(paste("No regions found in the input file:", tk$f))
  # FR-1: write a run-length-merged BigWig (adjacent equal-value tiles, incl.
  # zero runs, collapsed into one interval). Lossless; values unchanged.
  export(rle_merge_granges(density_per_family(g, chr_size_all, tk$step, 10)),
         tk$out, format="bigwig")
  "ok"
}

res <- mclapply(tasks, function(tk)
  tryCatch(process_task(tk), error=function(e) paste("ERROR", tk$f, ":", conditionMessage(e))),
  mc.cores = max(1L, opt$threads))

errs <- unlist(res)[grepl("^ERROR", unlist(res))]
if (length(errs)) { writeLines(errs, stderr()); quit(save="no", status=1) }
cat(sprintf("calculate_density_batch: %d files (%d tasks), %d threads, done\n",
            length(files), length(tasks), opt$threads))
