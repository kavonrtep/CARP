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
source(file.path(.density_script_dir, "mem_utils.R"))

.script_start <- proc.time()[3]
log_msg <- function(...) message(sprintf("[+%6.1fs] ", proc.time()[3] - .script_start), ...)
log_mem <- function(label) {
  s <- mem_str()
  if (!is.na(s)) log_msg(sprintf("[mem] %-40s %s", label, s))
}

option_list <- list(
  make_option(c("-d", "--dir"), type="character", default=NULL, help="Directory of GFF3 files"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, help="Output directory for BigWig files"),
  make_option(c("-g", "--genome"), type="character", default=NULL, help="genome_seqlengths.rds (named integer vector of scaffold lengths)"),
  make_option(c("-t", "--threads"), type="integer", default=1, help="Parallel workers over input files [default %default]"),
  make_option(c("--max_workers"), type="integer", default=0,
              help="Hard ceiling on concurrent workers; 0 = no ceiling (the memory gate still applies) [default %default]"),
  make_option(c("--mem_budget_gb"), type="double", default=0,
              help="Memory budget in GB (the job's allocation) for sizing the worker pool; 0 = auto-detect: AGENT_MEMORY, then the scheduler environment (PBS_RESC_MEM, SLURM_MEM_PER_NODE, ...), then the tightest of the cgroup limit and MemAvailable [default %default]")
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

# Longest-processing-time first: heaviest input, and within an input the finer
# (step=1000) resolution, are dispatched first so a dominant class file cannot
# start late and strand the pool. Dispatch order only — every task writes its own
# independent file, so the outputs cannot depend on it.
if (length(tasks) > 0) {
  sizes <- vapply(tasks, function(tk) as.numeric(file.size(tk$f)), numeric(1))
  steps <- vapply(tasks, function(tk) as.numeric(tk$step), numeric(1))
  tasks <- tasks[order(-sizes, steps)]
}

process_task <- function(tk){
  g <- import(tk$f, format="gff3")
  if (length(g) == 0) return(paste("No regions found in the input file:", tk$f))
  # density_track() streams the tileGenome grid one sequence at a time and
  # returns the run-length-merged track (FR-1: adjacent equal-value tiles, incl.
  # zero runs, collapsed into one interval). Lossless; values unchanged.
  export(density_track(g, chr_size_all, tk$step, 10), tk$out, format="bigwig")
  "ok"
}

# Each task reports its own peak RSS so the pool can be sized from a measurement
# rather than a guess, and so a run that gets close to the limit says so.
run_task <- function(tk) {
  msg <- tryCatch(process_task(tk),
                  error = function(e) paste("ERROR", tk$f, ":", conditionMessage(e)))
  list(msg = msg, peak_mb = mem_hwm_mb())
}

# ── worker pool sizing ───────────────────────────────────────────────────────
# History: this script ran mc.cores = threads (= workflow.cores = 96 on the HPC
# profile) with no memory bound, and every task built the WHOLE-genome tile grid.
# On run-000156 (94.26 Gbp, 5,096 sequences) that peaked at 568.5 GB — 74% of a
# 768 GB host — for 1,603 tiny per-family tracks. density_track() removed the
# grid blow-up; this bounds what is left, which is dominated by import() of the
# GFF3 and by the size of the merged output track. Both vary by orders of
# magnitude between a 1-feature family and a 52 M-feature class, so instead of
# guessing a per-task cost we MEASURE it: run the heaviest task first, on its
# own, and divide the budget by its peak.
n_workers <- max(1L, min(opt$threads, length(tasks)))
log_msg(sprintf("%d file(s), %d task(s), up to %d worker(s) requested",
                length(files), length(tasks), n_workers))
log_mem("before probe task")

probe <- NULL
if (length(tasks) > 0 && n_workers > 1L) {
  # mcparallel/mccollect rather than mclapply: mclapply(mc.cores = 1) does not
  # fork, and running the probe in the parent would leave its heap inflated for
  # every subsequent fork.
  probe <- mccollect(mcparallel(run_task(tasks[[1]])))[[1]]
}
probe_peak <- if (is.list(probe) && !is.null(probe$peak_mb)) probe$peak_mb else NA_real_
if (!is.na(probe_peak))
  log_msg(sprintf("[mem] probe task %s (step=%d) peaked at %.2fG",
                  basename(tasks[[1]]$f), tasks[[1]]$step, probe_peak/1024))

if (opt$max_workers > 0 && n_workers > opt$max_workers) {
  log_msg(sprintf("[mem] capping workers %d -> %d (--max_workers)",
                  n_workers, opt$max_workers))
  n_workers <- opt$max_workers
}
{
  bud <- mem_budget_mb(opt$mem_budget_gb)
  warn_if_host_budget(bud$src)   # host-wide budget under PBS/Slurm -> say so now
  # NOTE: mem_budget_mb() already applies MEMORY_HEADROOM to an ALLOCATION
  # source (max_memory_gb / AGENT_MEMORY / scheduler vars); the 0.8 below is
  # this gate's own headroom and is kept for the measured sources (cgroup,
  # MemAvailable), which arrive unfactored. An explicit budget is therefore
  # spent conservatively — the gate caps concurrency only, never the result.
  if (!is.na(bud$mb) && !is.na(probe_peak) && probe_peak > 0) {
    # 20% headroom: the probe is the heaviest task, but the parent also holds the
    # seqlengths vector and the OS needs room for page cache during export.
    cap <- max(1L, as.integer(floor((bud$mb * 0.8) / probe_peak)))
    log_msg(sprintf("[mem] budget %.1fG (%s), probe peak %.2fG -> room for %d worker(s)",
                    bud$mb/1024, bud$src, probe_peak/1024, cap))
    if (cap < n_workers) {
      log_msg(sprintf("[mem] capping workers %d -> %d (memory gate)", n_workers, cap))
      n_workers <- cap
    }
  } else {
    log_msg(sprintf("[mem] no memory gate applied (budget=%s, probe peak=%s)",
                    ifelse(is.na(bud$mb), "unknown", sprintf("%.1fG", bud$mb/1024)),
                    ifelse(is.na(probe_peak), "unknown", sprintf("%.2fG", probe_peak/1024))))
  }
}
log_msg(sprintf("── Worker pool: %d concurrent worker(s) for %d remaining task(s) ─────",
                n_workers, length(tasks) - as.integer(!is.null(probe))))

# mc.preschedule = FALSE so tasks are handed out dynamically: with prescheduling
# mclapply splits the list into contiguous blocks, which would hand every heavy
# task to the first worker and defeat the LPT ordering above.
rest <- if (is.null(probe)) tasks else tasks[-1]
res <- if (length(rest) > 0) {
  mclapply(rest, run_task, mc.cores = n_workers, mc.preschedule = FALSE)
} else list()
all_res <- c(if (!is.null(probe)) list(probe), res)
done    <- if (is.null(probe)) tasks else c(tasks[1], rest)

log_mem("after worker pool")
peaks <- vapply(all_res,
                function(r) if (is.list(r) && !is.null(r$peak_mb)) r$peak_mb else NA_real_,
                numeric(1))
if (any(!is.na(peaks)))
  log_msg(sprintf("[mem] worker peak RSS: max %.2fG, median %.2fG over %d task(s)",
                  max(peaks, na.rm = TRUE)/1024, median(peaks, na.rm = TRUE)/1024,
                  sum(!is.na(peaks))))

# Validate positively. A signal-killed mclapply child returns NULL rather than a
# "try-error", so scanning the returned messages for "^ERROR" cannot see it — and
# because these rules are checkpointed by a .done marker rather than by the .bw
# files themselves, a silently missing track would have gone unnoticed.
errs <- character(0)
for (i in seq_along(done)) {
  r  <- all_res[[i]]
  tk <- done[[i]]
  if (!is.list(r) || is.null(r$msg) || !is.character(r$msg)) {
    errs <- c(errs, sprintf("ERROR %s (step=%d): worker produced no result (killed?)",
                            tk$f, tk$step))
  } else if (grepl("^ERROR", r$msg)) {
    errs <- c(errs, r$msg)
  } else if (identical(r$msg, "ok") && !file.exists(tk$out)) {
    errs <- c(errs, sprintf("ERROR %s (step=%d): reported ok but %s is missing",
                            tk$f, tk$step, tk$out))
  }
}
if (length(errs)) { writeLines(errs, stderr()); quit(save="no", status=1) }
cat(sprintf("calculate_density_batch: %d files (%d tasks), %d threads, done\n",
            length(files), length(tasks), opt$threads))
