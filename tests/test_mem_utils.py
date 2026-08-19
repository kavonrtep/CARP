#!/usr/bin/env python3
"""Unit test for the memory-budget resolver (scripts/mem_utils.py).

Why this matters: CARP runs from a .sif under PBS/Slurm, where /proc/meminfo
reports the HOST rather than the job (the kernel does not namespace it) and the
cgroup limit usually sits on an ANCESTOR job scope that may not be reachable
from inside the container. Sizing a worker pool from the wrong number is the
TideCluster issue-#6 failure: a 128 GB job that believed it had 1.6 TB and was
OOM-killed hours in. The resolver is the one place that decides what the whole
pipeline believes, and its inputs (cgroup trees, scheduler variables) cannot be
reproduced on a CI runner — so every rung is exercised against a synthetic
/sys/fs/cgroup + /proc tree here.

The R mirror (scripts/mem_utils.R) is tested case-for-case by
tests/test_mem_utils.R; the two MUST agree, or the startup banner and the R
worker-pool gates would disagree about the same job.

Run: python3 tests/test_mem_utils.py
"""
import io
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from mem_utils import (  # noqa: E402
    MEMORY_HEADROOM, describe_budget, memory_budget_mb, scheduler_job_id,
    warn_if_host_budget,
)

GIB_MB = 1024.0
# Deliberately awkward but real spellings: PBS emits bytes, Slurm MB, LSF KB.
BYTES_128G = 137438953472
KB_128G = 134217728
MB_128G = 131072


def _tree(td, cgroup_lines="", cgroup_files=None, mem_available_kb=None):
    """Build a synthetic sysfs/proc trio and return it as resolver kwargs.

    ``cgroup_files`` maps a path relative to the fake /sys/fs/cgroup to its
    content, so a test can put a limit on a leaf, on an ancestor, or nowhere.
    Any of the three files may be omitted to model "not readable in here".
    """
    sysfs = os.path.join(td, "cgroup")
    os.makedirs(sysfs, exist_ok=True)
    for rel, content in (cgroup_files or {}).items():
        path = os.path.join(sysfs, rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            fh.write(content + "\n")

    proc_cgroup = os.path.join(td, "self_cgroup")
    if cgroup_lines:
        with open(proc_cgroup, "w") as fh:
            fh.write(cgroup_lines if cgroup_lines.endswith("\n")
                     else cgroup_lines + "\n")

    meminfo = os.path.join(td, "meminfo")
    if mem_available_kb is not None:
        with open(meminfo, "w") as fh:
            fh.write("MemTotal:       131789948 kB\n"
                     F"MemAvailable:   {mem_available_kb} kB\n")

    return {"sysfs_root": sysfs, "proc_cgroup": proc_cgroup, "meminfo": meminfo}


def _budget(explicit=None, env=None, **tree):
    return memory_budget_mb(explicit, environ=(env or {}), **tree)


def test_explicit_wins_and_carries_headroom():
    """An explicit allocation outranks every other signal, and is not spent whole."""
    with tempfile.TemporaryDirectory() as td:
        env = {"AGENT_MEMORY": "500", "PBS_RESC_MEM": str(BYTES_128G),
               "SLURM_MEM_PER_NODE": str(MB_128G)}
        mb, src = _budget(128, env, **_tree(td, mem_available_kb=1_600_000_000))
        assert src == "max_memory_gb", src
        assert mb == 128 * GIB_MB * MEMORY_HEADROOM, mb


def test_agent_memory_then_scheduler_vars():
    """Rungs 2 and 3, each in its own unit, all landing on the same 128 GB job."""
    with tempfile.TemporaryDirectory() as td:
        tree = _tree(td)
        expected = MB_128G * MEMORY_HEADROOM

        mb, src = _budget(None, {"AGENT_MEMORY": "128"}, **tree)
        assert (src, mb) == ("AGENT_MEMORY", expected), (src, mb)

        for var, value in (("PBS_RESC_MEM", str(BYTES_128G)),
                           ("SLURM_MEM_PER_NODE", str(MB_128G)),
                           ("LSB_MAX_MEM_RUSAGE", str(KB_128G))):
            mb, src = _budget(None, {var: value}, **tree)
            assert (src, mb) == (var, expected), (var, src, mb)

        # Suffixed spellings must parse to the same budget as bare numbers.
        for value in ("128gb", "128GB", "131072mb", " 128g "):
            mb, src = _budget(None, {"PBS_RESC_MEM": value}, **tree)
            assert (src, mb) == ("PBS_RESC_MEM", expected), (value, src, mb)

        # Slurm's per-CPU form is a product, not a value.
        mb, src = _budget(None, {"SLURM_MEM_PER_CPU": "4096",
                                 "SLURM_CPUS_ON_NODE": "32"}, **tree)
        assert (src, mb) == ("SLURM_MEM_PER_CPU", expected), (src, mb)

        # Unparseable / meaningless values must not be trusted as a budget.
        for value in ("unlimited", "", "0", "lots"):
            mb, src = _budget(None, {"PBS_RESC_MEM": value}, **tree)
            assert src == "none", (value, src)


def test_cgroup_leaf_and_ancestor_scope():
    """The PBS case: the limit lives on an ancestor scope, not on the leaf."""
    with tempfile.TemporaryDirectory() as td:
        lines = "0::/user.slice/pbs_jobs.service/job-42.scope"
        # Leaf carries the limit: 8 GiB granted, 1 GiB already used -> 7 GiB.
        leaf = _tree(td, cgroup_lines=lines, cgroup_files={
            "user.slice/pbs_jobs.service/job-42.scope/memory.max": str(8 * 1024**3),
            "user.slice/pbs_jobs.service/job-42.scope/memory.current": str(1024**3),
        })
        mb, src = _budget(None, {}, **leaf)
        assert (src, round(mb)) == ("cgroup", 7 * 1024), (src, mb)

    with tempfile.TemporaryDirectory() as td:
        # Nothing on the leaf; the job scope two levels up holds the limit.
        anc = _tree(td, cgroup_lines="0::/user.slice/job-42.scope/leaf",
                    cgroup_files={
                        "user.slice/job-42.scope/memory.max": str(8 * 1024**3),
                        "user.slice/job-42.scope/memory.current": "0",
                    })
        mb, src = _budget(None, {}, **anc)
        assert (src, round(mb)) == ("cgroup", 8 * 1024), (src, mb)


def test_cgroup_v1_and_unlimited_spellings():
    """v1 layout, plus the two ways a cgroup says 'no limit'."""
    with tempfile.TemporaryDirectory() as td:
        v1 = _tree(td, cgroup_lines="7:memory:/pbs_jobs/42",
                   cgroup_files={
                       "memory/pbs_jobs/42/memory.limit_in_bytes": str(8 * 1024**3),
                       "memory/pbs_jobs/42/memory.usage_in_bytes": "0",
                   })
        mb, src = _budget(None, {}, **v1)
        assert (src, round(mb)) == ("cgroup", 8 * 1024), (src, mb)

    with tempfile.TemporaryDirectory() as td:
        # v2 "max" and the v1 PAGE_COUNTER_MAX sentinel both mean unlimited, so
        # the resolver must fall through rather than report a nonsense budget.
        unl = _tree(td, cgroup_lines="0::/scope\n7:memory:/scope",
                    cgroup_files={
                        "scope/memory.max": "max",
                        "memory/scope/memory.limit_in_bytes": "9223372036854710784",
                    },
                    mem_available_kb=10 * 1024 * 1024)
        mb, src = _budget(None, {}, **unl)
        assert (src, round(mb)) == ("MemAvailable", 10 * 1024), (src, mb)


def test_tightest_measured_reading_wins():
    """cgroup and MemAvailable are both availability readings: take the min.

    Exceeding either kills the job, so a generous cgroup limit on a node whose
    free memory has been eaten by neighbours must not license a big pool.
    """
    with tempfile.TemporaryDirectory() as td:
        tree = _tree(td, cgroup_lines="0::/scope", cgroup_files={
            "scope/memory.max": str(100 * 1024**3),
            "scope/memory.current": "0",
        }, mem_available_kb=10 * 1024 * 1024)
        mb, src = _budget(None, {}, **tree)
        assert (src, round(mb)) == ("MemAvailable", 10 * 1024), (src, mb)

        # ... and the other way round, with no headroom factor on either.
        tree = _tree(td, cgroup_lines="0::/scope", cgroup_files={
            "scope/memory.max": str(10 * 1024**3),
            "scope/memory.current": "0",
        }, mem_available_kb=100 * 1024 * 1024)
        mb, src = _budget(None, {}, **tree)
        assert (src, round(mb)) == ("cgroup", 10 * 1024), (src, mb)


def test_nothing_readable():
    """Off Linux / in a stripped container: report 'none', never a guess."""
    with tempfile.TemporaryDirectory() as td:
        mb, src = _budget(None, {}, **_tree(td))
        assert (mb, src) == (None, "none"), (mb, src)


def test_warning_only_fires_under_a_scheduler_on_a_host_budget():
    """The diagnostic that turns a five-hour OOM into a first-second warning."""
    def warned(src, env):
        stream = io.StringIO()
        fired = warn_if_host_budget(src, environ=env, stream=stream)
        text = stream.getvalue()
        assert fired == bool(text), (src, env, fired, text)
        return fired, text

    fired, text = warned("MemAvailable", {"PBS_JOBID": "42.pbs-m1"})
    assert fired and "42.pbs-m1" in text and "-m <GB>" in text, text
    assert warned("none", {"SLURM_JOB_ID": "9"})[0]
    # No scheduler: a host budget is simply the truth on a workstation.
    assert not warned("MemAvailable", {})[0]
    # A real limit was found: nothing to warn about, scheduler or not.
    assert not warned("cgroup", {"PBS_JOBID": "42"})[0]
    assert not warned("max_memory_gb", {"PBS_JOBID": "42"})[0]

    assert scheduler_job_id({"LSB_JOBID": "7"}) == ("LSB_JOBID", "7")
    assert scheduler_job_id({"PBS_JOBID": "  "}) is None


def test_describe_budget_distinguishes_allocation_from_reading():
    """The banner must say where the number came from — that is its whole job."""
    line = describe_budget(*memory_budget_mb(128, environ={}))
    assert line == "[mem] budget 104858 MB from max_memory_gb (128.0 GB allocation x 0.8)", line
    assert describe_budget(10240.0, "cgroup") == "[mem] budget 10240 MB from cgroup"
    assert "not be memory-gated" in describe_budget(None, "none")


def main():
    test_explicit_wins_and_carries_headroom()
    test_agent_memory_then_scheduler_vars()
    test_cgroup_leaf_and_ancestor_scope()
    test_cgroup_v1_and_unlimited_spellings()
    test_tightest_measured_reading_wins()
    test_nothing_readable()
    test_warning_only_fires_under_a_scheduler_on_a_host_budget()
    test_describe_budget_distinguishes_allocation_from_reading()
    print("OK  test_mem_utils: explicit > AGENT_MEMORY > scheduler > cgroup "
          "(leaf+ancestor, v1+v2) > MemAvailable > none; tightest measured "
          "reading wins; scheduler warning fires only when it should")


if __name__ == "__main__":
    main()
