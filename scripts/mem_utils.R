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

# ── memory budget (cgroup-aware) ─────────────────────────────────────────────
# /proc/meminfo reports the HOST. Under PBS/Slurm (or any container) the limit
# that will actually kill this job is the cgroup's, and it is typically set on an
# ANCESTOR scope — the job scope — not on the leaf the process sits in. So walk
# from the leaf up to the root and take the tightest effective limit. Falls back
# to MemAvailable when no cgroup limit applies, and returns NA when neither can
# be read (callers then skip gating rather than guess).
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

cgroup_avail_mb <- function() {
  lines <- tryCatch(readLines("/proc/self/cgroup", warn = FALSE),
                    error = function(e) character(0))
  cands <- numeric(0)

  # cgroup v2: "0::/user.slice/.../job.scope"
  v2 <- grep("^0::", lines, value = TRUE)
  if (length(v2) > 0) {
    parts <- Filter(nzchar, strsplit(sub("^0::", "", v2[1]), "/", fixed = TRUE)[[1]])
    # leaf -> root, so an ancestor job-scope limit is picked up
    for (i in rev(seq_along(parts))) {
      d <- file.path("/sys/fs/cgroup", paste(parts[seq_len(i)], collapse = "/"))
      cands <- c(cands, .cgroup_dir_avail_mb(d, "memory.max", "memory.current"))
    }
  }
  # Namespaced container: the job's own limits appear at the mount root.
  cands <- c(cands, .cgroup_dir_avail_mb("/sys/fs/cgroup", "memory.max", "memory.current"))

  # cgroup v1
  v1 <- grep(":memory:", lines, value = TRUE)
  if (length(v1) > 0) {
    rel <- sub("^.*:memory:", "", v1[1])
    parts <- Filter(nzchar, strsplit(rel, "/", fixed = TRUE)[[1]])
    for (i in rev(seq_along(parts))) {
      d <- file.path("/sys/fs/cgroup/memory", paste(parts[seq_len(i)], collapse = "/"))
      cands <- c(cands, .cgroup_dir_avail_mb(d, "memory.limit_in_bytes",
                                             "memory.usage_in_bytes"))
    }
  }
  cands <- c(cands, .cgroup_dir_avail_mb("/sys/fs/cgroup/memory",
                                         "memory.limit_in_bytes", "memory.usage_in_bytes"))

  cands <- cands[!is.na(cands)]
  if (length(cands) == 0) NA_real_ else min(cands)
}

# Budget in MB. Explicit override wins; otherwise the tightest signal, because
# exceeding ANY of them kills the job.
mem_budget_mb <- function(explicit_gb = 0) {
  if (!is.na(explicit_gb) && explicit_gb > 0)
    return(list(mb = explicit_gb * 1024, src = "--mem_budget_gb"))
  cg   <- cgroup_avail_mb()
  host <- mem_avail_mb()
  cands <- c(cgroup = cg, host = host)
  cands <- cands[!is.na(cands)]
  if (length(cands) == 0) return(list(mb = NA_real_, src = "unavailable"))
  list(mb = min(cands), src = names(cands)[which.min(cands)])
}
