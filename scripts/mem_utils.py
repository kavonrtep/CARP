#!/usr/bin/env python3
"""Memory-budget resolution, shared by the Python side of the pipeline.

The Python mirror of ``scripts/mem_utils.R``. Both resolve the same chain and
return the same source labels, so every layer of a run — the ``run_pipeline.py``
banner, the R worker-pool gates, and TideCluster's own gates — reports the same
number from the same source.

WHY THIS EXISTS
---------------
CARP normally runs from a ``.sif`` under PBS/Slurm. ``/proc/meminfo``
``MemAvailable`` reports the **host**, not the job — the kernel does not
namespace it — and the cgroup limit that will actually kill the job is usually
set on an ANCESTOR scope (the job scope), which may not be reachable from inside
the container at all: controllers not delegated, ``/sys/fs/cgroup`` not mounted,
or a scheduler that enforces by polling rather than by cgroup. A pool sized from
host memory then over-allocates by an order of magnitude and the run is
OOM-killed hours in, with no earlier warning (TideCluster issue #6: a 128 GB PBS
job believed it had 1.6 TB).

The chain below is deliberately identical to TideCluster 1.20.0's
``tc_utils.memory_budget_mb()``. First hit wins:

  1. explicit budget      -- config ``max_memory_gb`` / ``--mem_budget_gb``
  2. ``AGENT_MEMORY``     -- GB, sandbox override
  3. scheduler env        -- ``PBS_RESC_MEM`` (bytes), ``SLURM_MEM_PER_NODE`` (MB),
                             ``LSB_MAX_MEM_RUSAGE`` (KB),
                             ``SLURM_MEM_PER_CPU`` x ``SLURM_CPUS_ON_NODE`` (MB).
                             These survive the container boundary, which is what
                             makes them the reliable signal under a scheduler.
  4/5. the tightest measured reading -- cgroup headroom (v2 ``memory.max`` /
                             v1 ``memory.limit_in_bytes``, minus current usage)
                             and ``/proc/meminfo`` ``MemAvailable``
  6. nothing readable     -- ``(None, "none")``; callers skip gating rather than guess

Sources 1-3 are ALLOCATIONS (what the job was granted), so ``MEMORY_HEADROOM``
is applied: a pool sized to exactly the job limit leaves nothing for the parent
process, the OS cache or tool children. Sources 4-5 are already
measured-availability numbers (cgroup ``limit - current``, ``MemAvailable``) and
pass through unchanged.

Three deliberate differences from upstream, all retained from mem_utils.R: the
cgroup rung reports ``limit - current`` rather than the raw limit; it walks
leaf -> root so a limit set on an ancestor job scope is picked up; and the two
measured rungs are combined as ``min()`` rather than first-hit, because
exceeding either one kills the job.

Dependency-free (stdlib only): ``run_pipeline.py`` imports it before any conda
environment exists. Every filesystem root is injectable so the tests can build
fake ``/sys/fs/cgroup`` and ``/proc`` trees.
"""
from __future__ import annotations

import os
import re
import sys

# Fraction of an ALLOCATION a worker pool may plan to use. Matches
# TideCluster 1.20.0's MEMORY_HEADROOM, so both agree on the same job.
MEMORY_HEADROOM = 0.8

# cgroup v1 spells "unlimited" as a huge sentinel (PAGE_COUNTER_MAX).
_CGROUP_UNLIMITED = 1e15

_MEM_UNIT_MB = {
    "b": 1.0 / 1024 / 1024,
    "k": 1.0 / 1024, "kb": 1.0 / 1024,
    "m": 1.0, "mb": 1.0,
    "g": 1024.0, "gb": 1024.0,
    "t": 1024.0 * 1024, "tb": 1024.0 * 1024,
}

# (variable, unit when the value carries no suffix)
_SCHEDULER_MEM_VARS = (
    ("PBS_RESC_MEM", "b"),        # PBS/Torque: bytes
    ("SLURM_MEM_PER_NODE", "m"),  # Slurm: MB
    ("LSB_MAX_MEM_RUSAGE", "k"),  # LSF: KB
)

_SCHEDULER_JOB_VARS = ("PBS_JOBID", "SLURM_JOB_ID", "SLURM_JOBID", "LSB_JOBID")

_MEM_RE = re.compile(r"(\d+(?:\.\d+)?)\s*([kmgt]?b?)")


def _parse_mem_to_mb(raw, default_unit: str):
    """Parse a scheduler memory string to MB. ``None`` if unparseable.

    Accepts a bare number in ``default_unit`` (``PBS_RESC_MEM=137438953472``)
    or a suffixed one (``128gb``, ``2048mb``) — schedulers emit both.
    """
    if raw is None:
        return None
    m = _MEM_RE.fullmatch(str(raw).strip().lower())
    if not m:
        return None
    factor = _MEM_UNIT_MB.get(m.group(2) or default_unit)
    if factor is None:
        return None
    value = float(m.group(1)) * factor
    return value if value > 0 else None


def _read_cgroup_num(path: str):
    """Numeric content of a cgroup control file.

    Returns ``None`` when absent/unreadable and ``float("inf")`` for the two
    spellings of "unlimited" (v2 ``max``, v1 sentinel).
    """
    try:
        with open(path) as fh:
            raw = fh.read().strip()
    except (OSError, ValueError):
        return None
    if not raw:
        return None
    if raw == "max":
        return float("inf")
    try:
        value = float(raw)
    except ValueError:
        return None
    return float("inf") if value > _CGROUP_UNLIMITED else value


def _cgroup_dir_avail_mb(directory: str, lim_file: str, cur_file: str):
    """Headroom (MB) inside one cgroup directory: ``limit - current usage``."""
    limit = _read_cgroup_num(os.path.join(directory, lim_file))
    if limit is None or limit == float("inf"):
        return None                                  # absent or unlimited
    current = _read_cgroup_num(os.path.join(directory, cur_file))
    if current is None or current == float("inf"):
        current = 0.0
    return max(0.0, (limit - current) / 1048576)


def _chain(base: str, rel: str) -> list:
    """``base`` plus every ancestor path of ``rel`` under it, leaf first."""
    paths, path = [], rel
    while path and path != "/":
        paths.append(base + path)
        path = os.path.dirname(path)
    paths.append(base)
    return paths


def cgroup_avail_mb(sysfs_root: str = "/sys/fs/cgroup",
                    proc_cgroup: str = "/proc/self/cgroup"):
    """Tightest effective cgroup headroom in MB, or ``None`` if there is none.

    Walks leaf -> root because the limit that kills a batch job is typically set
    on an ancestor scope, not on the leaf the process sits in. The bare mount
    root is probed too: inside a cgroup-namespaced container the job's own
    limits appear there.
    """
    rel_v2, rel_v1 = "", ""
    try:
        with open(proc_cgroup) as fh:
            for line in fh:
                parts = line.strip().split(":", 2)
                if len(parts) != 3:
                    continue
                if parts[0] == "0" and not rel_v2:                     # v2 unified
                    rel_v2 = parts[2]
                elif "memory" in parts[1].split(",") and not rel_v1:   # v1
                    rel_v1 = parts[2]
    except OSError:
        pass

    candidates = []
    for path in _chain(sysfs_root, rel_v2):
        mb = _cgroup_dir_avail_mb(path, "memory.max", "memory.current")
        if mb is not None:
            candidates.append(mb)
    for path in _chain(os.path.join(sysfs_root, "memory"), rel_v1):
        mb = _cgroup_dir_avail_mb(path, "memory.limit_in_bytes",
                                  "memory.usage_in_bytes")
        if mb is not None:
            candidates.append(mb)
    return min(candidates) if candidates else None


def mem_available_mb(meminfo: str = "/proc/meminfo"):
    """``MemAvailable`` in MB — the HOST's, inside a container. ``None`` off Linux."""
    try:
        with open(meminfo) as fh:
            for line in fh:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) / 1024
    except (OSError, ValueError, IndexError):
        pass
    return None


def memory_budget_mb(explicit_gb=None, environ=None,
                     sysfs_root: str = "/sys/fs/cgroup",
                     proc_cgroup: str = "/proc/self/cgroup",
                     meminfo: str = "/proc/meminfo"):
    """Resolve the memory budget. Returns ``(budget_mb, source)``.

    ``source`` is one of ``max_memory_gb``, ``AGENT_MEMORY``, a scheduler
    variable name, ``cgroup``, ``MemAvailable`` or ``none`` — see the module
    docstring for the full chain. ``budget_mb`` is ``None`` only for ``none``.
    """
    env = os.environ if environ is None else environ

    if explicit_gb:
        return float(explicit_gb) * 1024 * MEMORY_HEADROOM, "max_memory_gb"

    agent = _parse_mem_to_mb(env.get("AGENT_MEMORY"), "g")
    if agent:
        return agent * MEMORY_HEADROOM, "AGENT_MEMORY"

    for var, unit in _SCHEDULER_MEM_VARS:
        mb = _parse_mem_to_mb(env.get(var), unit)
        if mb:
            return mb * MEMORY_HEADROOM, var

    per_cpu = _parse_mem_to_mb(env.get("SLURM_MEM_PER_CPU"), "m")
    ncpu = (env.get("SLURM_CPUS_ON_NODE") or "").strip()
    if per_cpu and ncpu.isdigit() and int(ncpu) > 0:
        return per_cpu * int(ncpu) * MEMORY_HEADROOM, "SLURM_MEM_PER_CPU"

    # Rungs 4-5 are both measured-availability readings, so take the TIGHTEST
    # rather than the first hit: exceeding either one kills the job. The label
    # names whichever won (cgroup on a tie).
    measured = [(mb, src) for mb, src in (
        (cgroup_avail_mb(sysfs_root=sysfs_root, proc_cgroup=proc_cgroup), "cgroup"),
        (mem_available_mb(meminfo=meminfo), "MemAvailable")) if mb]
    if measured:
        return min(measured, key=lambda pair: pair[0])

    return None, "none"


def scheduler_job_id(environ=None):
    """``(variable, value)`` of the first batch-scheduler job id set, else ``None``."""
    env = os.environ if environ is None else environ
    for var in _SCHEDULER_JOB_VARS:
        value = (env.get(var) or "").strip()
        if value:
            return var, value
    return None


# Sources that name an ALLOCATION (and therefore carry MEMORY_HEADROOM), as
# opposed to a measured-availability reading (cgroup headroom, MemAvailable).
_ALLOCATION_SOURCES = frozenset(
    ("max_memory_gb", "AGENT_MEMORY", "SLURM_MEM_PER_CPU")
    + tuple(var for var, _unit in _SCHEDULER_MEM_VARS))


def describe_budget(budget_mb, source) -> str:
    """The one-line ``[mem] …`` summary printed at startup and by the gates."""
    if budget_mb is None:
        return ("[mem] no memory budget could be detected (source: none) - "
                "worker pools will not be memory-gated")
    detail = ""
    if source in _ALLOCATION_SOURCES:
        detail = (F" ({budget_mb / MEMORY_HEADROOM / 1024:.1f} GB allocation"
                  F" x {MEMORY_HEADROOM:g})")
    return F"[mem] budget {budget_mb:.0f} MB from {source}{detail}"


def warn_if_host_budget(source, environ=None, stream=None) -> bool:
    """Warn when a host-wide budget is in use under a batch scheduler.

    ``MemAvailable`` (or nothing at all) while ``$PBS_JOBID`` / ``$SLURM_JOB_ID``
    is set means the memory gates are sizing pools against the NODE's free
    memory instead of the job's limit — the issue-#6 failure, which otherwise
    only becomes visible when the OOM killer arrives. Returns whether it warned.
    """
    if source not in ("MemAvailable", "none"):
        return False
    job = scheduler_job_id(environ)
    if job is None:
        return False
    var, value = job
    print(F"WARNING: memory budget came from "
          F"{'/proc/meminfo MemAvailable' if source == 'MemAvailable' else 'nowhere'}, "
          F"but {var}={value} is set. MemAvailable is not cgroup-aware, so this is "
          F"probably the whole node's memory, not this job's limit — worker pools "
          F"may over-allocate and be OOM-killed. Pass the job's allocation with "
          F"'run_pipeline.py -m <GB>' (or set max_memory_gb in the config).",
          file=sys.stderr if stream is None else stream)
    return True


def main() -> int:
    """Report what this machine/job resolves to — a debugging aid.

    Useful inside a container to see what the pipeline will believe:
    ``singularity exec image.sif /opt/pipeline/scripts/mem_utils.py``
    """
    import argparse
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--max_memory", "--max-memory", type=float, default=None,
                        metavar="GB", help="explicit allocation to resolve against")
    args = parser.parse_args()
    budget_mb, source = memory_budget_mb(args.max_memory)
    print(describe_budget(budget_mb, source))
    warn_if_host_budget(source)
    job = scheduler_job_id()
    print(F"[mem] scheduler job: {job[0]}={job[1]}" if job
          else "[mem] no batch-scheduler job id in the environment")
    return 0


if __name__ == "__main__":
    sys.exit(main())
