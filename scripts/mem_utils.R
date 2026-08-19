# mem_utils.R — memory instrumentation and budget detection, shared by the R
# scripts that fork worker pools (make_unified_annotation.R,
# calculate_density_batch.R via density_utils.R).
#
# Added after run-000156 (94 Gbp) died in make_unified_annotation with 15 of 55
# mclapply workers killed and no memory figures to size a fix from: a rule's
# benchmark file is only written on success, so a failing run leaves nothing
# behind. These read /proc directly so every run records what it actually used.
#
#   VmRSS  current resident set
#   VmHWM  peak resident set for THIS process. A forked mclapply child inherits
#          the parent's counters, so a child's peak is its whole footprint
#          (shared parent pages + its own) — which is the number to divide the
#          host budget by when sizing the pool.
#   MemAvailable  host headroom, so a run heading for the OOM killer is visible
#          in the log before it dies rather than only in dmesg.
#
# All return NA off Linux; every call site degrades to printing nothing.
#
# NOTE: this file must stay dependency-free (base R only) — it is sourced by
# scripts with very different library stacks.

.mem_field_mb <- function(path, field) {
  if (!file.exists(path)) return(NA_real_)
  ln <- grep(paste0("^", field, ":"), readLines(path, warn = FALSE), value = TRUE)
  if (length(ln) == 0) return(NA_real_)
  suppressWarnings(as.numeric(sub("^[^0-9]*([0-9]+).*$", "\\1", ln[1])) / 1024)
}
mem_rss_mb   <- function() .mem_field_mb("/proc/self/status", "VmRSS")
mem_hwm_mb   <- function() .mem_field_mb("/proc/self/status", "VmHWM")
mem_avail_mb <- function() .mem_field_mb("/proc/meminfo",     "MemAvailable")

mem_str <- function() {
  r <- mem_rss_mb(); h <- mem_hwm_mb(); a <- mem_avail_mb()
  if (is.na(r)) return(NA_character_)
  sprintf("rss=%.1fG peak=%.1fG host_avail=%.1fG", r/1024, h/1024, a/1024)
}

# ── memory budget (scheduler- and cgroup-aware) ──────────────────────────────
# Python mirror: scripts/mem_utils.py — same chain, same source labels, so the
# run_pipeline.py banner, these R gates and TideCluster 1.20's own gates all
# report the same number from the same source. See that file's header for the
# full rationale; the short version is that /proc/meminfo reports the HOST (the
# kernel does not namespace it), so inside a .sif under PBS/Slurm it is the
# node's memory and not the job's limit — TideCluster issue #6, a 128 GB job
# that believed it had 1.6 TB and was OOM-killed hours in.
#
# First hit wins:
#   1. explicit --mem_budget_gb (config max_memory_gb)     x MEMORY_HEADROOM
#   2. AGENT_MEMORY (GB)                                   x MEMORY_HEADROOM
#   3. scheduler env: PBS_RESC_MEM (bytes), SLURM_MEM_PER_NODE (MB),
#      LSB_MAX_MEM_RUSAGE (KB), SLURM_MEM_PER_CPU x SLURM_CPUS_ON_NODE (MB)
#                                                          x MEMORY_HEADROOM
#   4/5. the TIGHTEST of the cgroup headroom (leaf -> root, v1+v2,
#        limit - current usage) and MemAvailable — both are already
#        availability readings, so no headroom factor is applied, and taking
#        the min is right because exceeding either one kills the job
#   6. nothing readable -> NA; callers then skip gating rather than guess
#
# Sources 1-3 name an ALLOCATION (what the job was granted), so MEMORY_HEADROOM
# leaves room for the parent process, the OS cache and tool children.

MEMORY_HEADROOM <- 0.8

.MEM_UNIT_MB <- c(b = 1/1024/1024, k = 1/1024, kb = 1/1024, m = 1, mb = 1,
                  g = 1024, gb = 1024, t = 1024*1024, tb = 1024*1024)

# (variable, unit assumed when the value carries no suffix)
.SCHEDULER_MEM_VARS <- list(c("PBS_RESC_MEM",       "b"),
                            c("SLURM_MEM_PER_NODE", "m"),
                            c("LSB_MAX_MEM_RUSAGE", "k"))
.SCHEDULER_JOB_VARS <- c("PBS_JOBID", "SLURM_JOB_ID", "SLURM_JOBID", "LSB_JOBID")

# Sources that name an allocation (and so carry MEMORY_HEADROOM).
.ALLOCATION_SOURCES <- c("max_memory_gb", "AGENT_MEMORY", "SLURM_MEM_PER_CPU",
                         vapply(.SCHEDULER_MEM_VARS, `[`, character(1), 1))

# Schedulers emit both bare numbers (PBS_RESC_MEM=137438953472) and suffixed
# ones (128gb), so accept either; NA when unparseable.
.parse_mem_to_mb <- function(raw, default_unit) {
  if (is.null(raw) || length(raw) == 0 || is.na(raw) || !nzchar(raw))
    return(NA_real_)
  x <- tolower(trimws(raw))
  m <- regmatches(x, regexec("^([0-9]+(\\.[0-9]+)?)[[:space:]]*([kmgt]?b?)$", x))[[1]]
  if (length(m) == 0) return(NA_real_)
  unit <- if (nzchar(m[4])) m[4] else default_unit
  factor <- unname(.MEM_UNIT_MB[unit])
  if (is.na(factor)) return(NA_real_)
  v <- suppressWarnings(as.numeric(m[2])) * factor
  if (is.na(v) || v <= 0) NA_real_ else v
}

.getenv <- function(name, environ) {
  if (is.null(environ)) Sys.getenv(name, "") else {
    v <- environ[[name]]
    if (is.null(v)) "" else as.character(v)
  }
}

.read_cgroup_num <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  v <- tryCatch(trimws(readLines(path, warn = FALSE)[1]),
                error = function(e) NA_character_, warning = function(w) NA_character_)
  if (is.na(v) || !nzchar(v)) return(NA_real_)
  if (identical(v, "max")) return(Inf)                       # cgroup v2 unlimited
  n <- suppressWarnings(as.numeric(v))
  if (is.na(n)) return(NA_real_)
  # cgroup v1 spells "unlimited" as a huge sentinel (PAGE_COUNTER_MAX).
  if (n > 1e15) return(Inf)
  n
}

# Headroom (MB) inside a single cgroup directory: limit - current usage.
.cgroup_dir_avail_mb <- function(dir, lim_file, cur_file) {
  lim <- .read_cgroup_num(file.path(dir, lim_file))
  if (is.na(lim) || !is.finite(lim)) return(NA_real_)         # absent or unlimited
  cur <- .read_cgroup_num(file.path(dir, cur_file))
  if (is.na(cur) || !is.finite(cur)) cur <- 0
  max(0, (lim - cur) / 1048576)
}

# Every ancestor path of `rel` under `base`, leaf first, then `base` itself.
.cgroup_chain <- function(base, rel) {
  parts <- Filter(nzchar, strsplit(rel, "/", fixed = TRUE)[[1]])
  paths <- character(0)
  for (i in rev(seq_along(parts)))
    paths <- c(paths, file.path(base, paste(parts[seq_len(i)], collapse = "/")))
  c(paths, base)
}

# Under PBS/Slurm the limit is typically set on an ANCESTOR scope (the job
# scope), not on the leaf the process sits in, so walk leaf -> root and take the
# tightest. The bare mount root is probed too: inside a cgroup-namespaced
# container the job's own limits appear there.
cgroup_avail_mb <- function(sysfs_root = "/sys/fs/cgroup",
                            proc_cgroup = "/proc/self/cgroup") {
  lines <- if (!file.exists(proc_cgroup)) character(0) else
    tryCatch(readLines(proc_cgroup, warn = FALSE),
             error = function(e) character(0), warning = function(w) character(0))
  rel_v2 <- ""; rel_v1 <- ""
  v2 <- grep("^0::", lines, value = TRUE)
  if (length(v2) > 0) rel_v2 <- sub("^0::", "", v2[1])
  v1 <- grep(":memory:", lines, value = TRUE)
  if (length(v1) > 0) rel_v1 <- sub("^.*:memory:", "", v1[1])

  cands <- numeric(0)
  for (d in .cgroup_chain(sysfs_root, rel_v2))
    cands <- c(cands, .cgroup_dir_avail_mb(d, "memory.max", "memory.current"))
  for (d in .cgroup_chain(file.path(sysfs_root, "memory"), rel_v1))
    cands <- c(cands, .cgroup_dir_avail_mb(d, "memory.limit_in_bytes",
                                           "memory.usage_in_bytes"))
  cands <- cands[!is.na(cands)]
  if (length(cands) == 0) NA_real_ else min(cands)
}

# Budget in MB. Returns list(mb = <MB or NA>, src = <source label>).
mem_budget_mb <- function(explicit_gb = 0, environ = NULL,
                          sysfs_root = "/sys/fs/cgroup",
                          proc_cgroup = "/proc/self/cgroup",
                          meminfo = "/proc/meminfo") {
  if (!is.na(explicit_gb) && explicit_gb > 0)
    return(list(mb = explicit_gb * 1024 * MEMORY_HEADROOM, src = "max_memory_gb"))

  agent <- .parse_mem_to_mb(.getenv("AGENT_MEMORY", environ), "g")
  if (!is.na(agent)) return(list(mb = agent * MEMORY_HEADROOM, src = "AGENT_MEMORY"))

  for (v in .SCHEDULER_MEM_VARS) {
    mb <- .parse_mem_to_mb(.getenv(v[1], environ), v[2])
    if (!is.na(mb)) return(list(mb = mb * MEMORY_HEADROOM, src = v[1]))
  }
  per_cpu <- .parse_mem_to_mb(.getenv("SLURM_MEM_PER_CPU", environ), "m")
  ncpu <- suppressWarnings(as.integer(.getenv("SLURM_CPUS_ON_NODE", environ)))
  if (!is.na(per_cpu) && !is.na(ncpu) && ncpu > 0)
    return(list(mb = per_cpu * ncpu * MEMORY_HEADROOM, src = "SLURM_MEM_PER_CPU"))

  cands <- c(cgroup = cgroup_avail_mb(sysfs_root, proc_cgroup),
             MemAvailable = .mem_field_mb(meminfo, "MemAvailable"))
  cands <- cands[!is.na(cands)]
  if (length(cands) == 0) return(list(mb = NA_real_, src = "none"))
  list(mb = min(cands), src = names(cands)[which.min(cands)])
}

# One-line summary, identical in wording to the Python mirror's describe_budget.
describe_budget <- function(mb, src) {
  if (is.na(mb))
    return("[mem] no memory budget could be detected (source: none) - worker pools will not be memory-gated")
  detail <- ""
  if (src %in% .ALLOCATION_SOURCES)
    detail <- sprintf(" (%.1f GB allocation x %g)", mb / MEMORY_HEADROOM / 1024,
                      MEMORY_HEADROOM)
  sprintf("[mem] budget %.0f MB from %s%s", mb, src, detail)
}

# Warn when a host-wide budget is in use under a batch scheduler: MemAvailable
# (or nothing at all) while $PBS_JOBID / $SLURM_JOB_ID is set means the gates are
# sizing pools against the NODE's free memory instead of the job's limit — the
# issue-#6 failure, otherwise invisible until the OOM killer arrives.
warn_if_host_budget <- function(src, environ = NULL) {
  if (!src %in% c("MemAvailable", "none")) return(invisible(FALSE))
  for (var in .SCHEDULER_JOB_VARS) {
    value <- .getenv(var, environ)
    if (nzchar(value)) {
      message(sprintf(
        paste0("WARNING: memory budget came from %s, but %s=%s is set. ",
               "MemAvailable is not cgroup-aware, so this is probably the whole ",
               "node's memory, not this job's limit — worker pools may ",
               "over-allocate and be OOM-killed. Pass the job's allocation with ",
               "'run_pipeline.py -m <GB>' (or set max_memory_gb in the config)."),
        if (identical(src, "MemAvailable")) "/proc/meminfo MemAvailable" else "nowhere",
        var, value))
      return(invisible(TRUE))
    }
  }
  invisible(FALSE)
}
