#!/usr/bin/env Rscript
# Unit test for the R memory-budget resolver (scripts/mem_utils.R).
#
# The case list mirrors tests/test_mem_utils.py one-for-one, because the two
# resolvers MUST agree: the same job is resolved by mem_utils.py for the
# run_pipeline.py startup banner and by mem_utils.R for the make_unified /
# BigWig worker-pool gates. If they diverge, the banner reports a budget the
# gates do not use — the diagnostic then lies about the failure it exists to
# prevent (TideCluster issue #6: a 128 GB PBS job that believed it had 1.6 TB,
# because /proc/meminfo is not namespaced inside a .sif).
#
# Base R only, like mem_utils.R itself — it is sourced by scripts with very
# different library stacks and must not drag in dependencies.
#
# Run: Rscript tests/test_mem_utils.R

root <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=",
             commandArgs(trailingOnly = FALSE), value = TRUE)[1]))))
source(file.path(root, "scripts", "mem_utils.R"))

GIB_MB     <- 1024
BYTES_128G <- 137438953472   # PBS emits bytes
KB_128G    <- 134217728      # LSF emits KB
MB_128G    <- 131072         # Slurm emits MB

ok <- function(cond, what) if (!isTRUE(cond)) stop("FAIL: ", what, call. = FALSE)
near <- function(a, b) !is.na(a) && !is.na(b) && abs(a - b) < 0.5

# Build a synthetic sysfs/proc trio; any file may be omitted to model
# "not readable in here". `files` maps a path under the fake /sys/fs/cgroup to
# its content, so a limit can sit on a leaf, on an ancestor, or nowhere.
tree <- function(cgroup_lines = NULL, files = list(), mem_available_kb = NULL) {
  td <- tempfile("memtree"); dir.create(td)
  sysfs <- file.path(td, "cgroup"); dir.create(sysfs)
  for (rel in names(files)) {
    path <- file.path(sysfs, rel)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(files[[rel]], path)
  }
  proc_cgroup <- file.path(td, "self_cgroup")
  if (!is.null(cgroup_lines)) writeLines(cgroup_lines, proc_cgroup)
  meminfo <- file.path(td, "meminfo")
  if (!is.null(mem_available_kb))
    writeLines(c("MemTotal:       131789948 kB",
                 sprintf("MemAvailable:   %.0f kB", mem_available_kb)), meminfo)
  list(sysfs_root = sysfs, proc_cgroup = proc_cgroup, meminfo = meminfo)
}

budget <- function(explicit = 0, env = list(), t = tree())
  mem_budget_mb(explicit, environ = env, sysfs_root = t$sysfs_root,
                proc_cgroup = t$proc_cgroup, meminfo = t$meminfo)

test_explicit_wins_and_carries_headroom <- function() {
  env <- list(AGENT_MEMORY = "500", PBS_RESC_MEM = as.character(BYTES_128G),
              SLURM_MEM_PER_NODE = as.character(MB_128G))
  b <- budget(128, env, tree(mem_available_kb = 1.6e9))
  ok(identical(b$src, "max_memory_gb"), "explicit source label")
  ok(near(b$mb, 128 * GIB_MB * MEMORY_HEADROOM), "explicit budget carries headroom")
}

test_agent_memory_then_scheduler_vars <- function() {
  t <- tree(); expected <- MB_128G * MEMORY_HEADROOM

  b <- budget(0, list(AGENT_MEMORY = "128"), t)
  ok(identical(b$src, "AGENT_MEMORY") && near(b$mb, expected), "AGENT_MEMORY rung")

  vars <- list(PBS_RESC_MEM = as.character(BYTES_128G),
               SLURM_MEM_PER_NODE = as.character(MB_128G),
               LSB_MAX_MEM_RUSAGE = as.character(KB_128G))
  for (v in names(vars)) {
    env <- list(); env[[v]] <- vars[[v]]
    b <- budget(0, env, t)
    ok(identical(b$src, v) && near(b$mb, expected), paste("scheduler rung", v))
  }

  # Suffixed spellings must parse to the same budget as bare numbers.
  for (value in c("128gb", "128GB", "131072mb", " 128g ")) {
    b <- budget(0, list(PBS_RESC_MEM = value), t)
    ok(identical(b$src, "PBS_RESC_MEM") && near(b$mb, expected),
       paste("suffixed spelling", value))
  }

  # Slurm's per-CPU form is a product, not a value.
  b <- budget(0, list(SLURM_MEM_PER_CPU = "4096", SLURM_CPUS_ON_NODE = "32"), t)
  ok(identical(b$src, "SLURM_MEM_PER_CPU") && near(b$mb, expected),
     "SLURM_MEM_PER_CPU x SLURM_CPUS_ON_NODE")

  # Unparseable / meaningless values must not be trusted as a budget.
  for (value in c("unlimited", "", "0", "lots")) {
    b <- budget(0, list(PBS_RESC_MEM = value), t)
    ok(identical(b$src, "none"), paste("unparseable value ignored:", value))
  }
}

test_cgroup_leaf_and_ancestor_scope <- function() {
  # Leaf carries the limit: 8 GiB granted, 1 GiB already used -> 7 GiB headroom.
  t <- tree("0::/user.slice/pbs_jobs.service/job-42.scope", list(
    "user.slice/pbs_jobs.service/job-42.scope/memory.max"     = format(8 * 1024^3, scientific = FALSE),
    "user.slice/pbs_jobs.service/job-42.scope/memory.current" = format(1024^3, scientific = FALSE)))
  b <- budget(0, list(), t)
  ok(identical(b$src, "cgroup") && near(b$mb, 7 * 1024), "cgroup v2 leaf limit")

  # The PBS case: nothing on the leaf, the job scope two levels up holds it.
  t <- tree("0::/user.slice/job-42.scope/leaf", list(
    "user.slice/job-42.scope/memory.max"     = format(8 * 1024^3, scientific = FALSE),
    "user.slice/job-42.scope/memory.current" = "0"))
  b <- budget(0, list(), t)
  ok(identical(b$src, "cgroup") && near(b$mb, 8 * 1024), "cgroup v2 ancestor scope")
}

test_cgroup_v1_and_unlimited_spellings <- function() {
  t <- tree("7:memory:/pbs_jobs/42", list(
    "memory/pbs_jobs/42/memory.limit_in_bytes" = format(8 * 1024^3, scientific = FALSE),
    "memory/pbs_jobs/42/memory.usage_in_bytes" = "0"))
  b <- budget(0, list(), t)
  ok(identical(b$src, "cgroup") && near(b$mb, 8 * 1024), "cgroup v1 layout")

  # v2 "max" and the v1 PAGE_COUNTER_MAX sentinel both mean unlimited, so the
  # resolver must fall through rather than report a nonsense budget.
  t <- tree(c("0::/scope", "7:memory:/scope"), list(
    "scope/memory.max"                    = "max",
    "memory/scope/memory.limit_in_bytes"  = "9223372036854710784"),
    mem_available_kb = 10 * 1024 * 1024)
  b <- budget(0, list(), t)
  ok(identical(b$src, "MemAvailable") && near(b$mb, 10 * 1024),
     "unlimited cgroup falls through to MemAvailable")
}

test_tightest_measured_reading_wins <- function() {
  # Both are availability readings and exceeding either kills the job, so a
  # generous cgroup limit on a busy node must not license a big pool.
  t <- tree("0::/scope", list("scope/memory.max"     = format(100 * 1024^3, scientific = FALSE),
                              "scope/memory.current" = "0"),
            mem_available_kb = 10 * 1024 * 1024)
  b <- budget(0, list(), t)
  ok(identical(b$src, "MemAvailable") && near(b$mb, 10 * 1024), "host tighter than cgroup")

  t <- tree("0::/scope", list("scope/memory.max"     = format(10 * 1024^3, scientific = FALSE),
                              "scope/memory.current" = "0"),
            mem_available_kb = 100 * 1024 * 1024)
  b <- budget(0, list(), t)
  ok(identical(b$src, "cgroup") && near(b$mb, 10 * 1024), "cgroup tighter than host")
}

test_nothing_readable <- function() {
  b <- budget(0, list(), tree())
  ok(identical(b$src, "none") && is.na(b$mb), "nothing readable -> none, never a guess")
}

test_warning_only_fires_under_a_scheduler_on_a_host_budget <- function() {
  fired <- function(src, env) {
    msgs <- character(0)
    res <- withCallingHandlers(
      warn_if_host_budget(src, env),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
      })
    ok(isTRUE(res) == (length(msgs) > 0), "warning return value tracks the message")
    list(fired = isTRUE(res), text = paste(msgs, collapse = ""))
  }

  r <- fired("MemAvailable", list(PBS_JOBID = "42.pbs-m1"))
  ok(r$fired && grepl("42.pbs-m1", r$text, fixed = TRUE) &&
       grepl("-m <GB>", r$text, fixed = TRUE), "warns and names the job + the fix")
  ok(fired("none", list(SLURM_JOB_ID = "9"))$fired, "warns when nothing was readable")
  # No scheduler: a host budget is simply the truth on a workstation.
  ok(!fired("MemAvailable", list())$fired, "silent without a scheduler job id")
  # A real limit was found: nothing to warn about, scheduler or not.
  ok(!fired("cgroup", list(PBS_JOBID = "42"))$fired, "silent on a cgroup budget")
  ok(!fired("max_memory_gb", list(PBS_JOBID = "42"))$fired, "silent on an explicit budget")
}

test_describe_budget_distinguishes_allocation_from_reading <- function() {
  b <- budget(128, list(), tree())
  ok(identical(describe_budget(b$mb, b$src),
               "[mem] budget 104858 MB from max_memory_gb (128.0 GB allocation x 0.8)"),
     "allocation sources report the allocation and the factor")
  ok(identical(describe_budget(10240, "cgroup"), "[mem] budget 10240 MB from cgroup"),
     "measured readings report no factor")
  ok(grepl("not be memory-gated", describe_budget(NA_real_, "none")), "none is explicit")
}

test_explicit_wins_and_carries_headroom()
test_agent_memory_then_scheduler_vars()
test_cgroup_leaf_and_ancestor_scope()
test_cgroup_v1_and_unlimited_spellings()
test_tightest_measured_reading_wins()
test_nothing_readable()
test_warning_only_fires_under_a_scheduler_on_a_host_budget()
test_describe_budget_distinguishes_allocation_from_reading()
cat("OK  test_mem_utils.R: R resolver matches the Python mirror case-for-case\n")
