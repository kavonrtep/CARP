# DANTE_LTR: the chunk pool's memory gate reads the host, not the job

> Feature request for `github.com/kavonrtep/dante_ltr` — observed on **0.5.3.0**.
>
> **Resolved in dante_ltr 0.5.4.0** (upstream issue #13), pinned by CARP. Every
> proposal below was implemented as written: `--max_memory <GB>` (alias
> `--max-memory`), the five-rung chain, the source named in the pool-sizing line,
> and the scheduler warning — including the detail that the helper returns the
> *raw* budget so the existing `0.8` at the call site is not applied twice.
>
> This is the DANTE_LTR counterpart of
> [TideCluster issue #6](https://github.com/kavonrtep/TideCluster/issues/6),
> resolved there in **TideCluster 1.20.0** with a `--max_memory <GB>` flag and a
> five-rung resolution chain. The same class of bug is present in DANTE_LTR, and
> DANTE_LTR is now the only tool in the CARP stack with **no way at all** to be
> told the real budget — no flag, no environment variable.

## Summary

0.5.3.0 parallelised the per-chunk detection loop and, correctly, gated the pool
by memory as well as by cores:

```python
pool_size = min(args.cpu, 0.8 * MemAvailable / per-chunk peak RSS)
```

The per-chunk term is measured properly — the largest chunk runs first as a probe
and `ru_maxrss` gives its exact footprint. The **budget** term is the problem:
`MemAvailable` comes from `/proc/meminfo`, which the kernel does **not**
namespace. Inside a container, and under a batch scheduler, it reports the whole
node's free memory rather than what this job is allowed to use.

So the one safeguard that stands between a large `-c` and an OOM is fed a number
that can be an order of magnitude too large, and it silently authorises a pool
that cannot fit. The failure is invisible until the OOM killer arrives, hours in.

## Environment

```
tool        dante_ltr 0.5.3.0  (share/dante_ltr/dante_ltr, version.py '0.5.3.0')
caller      CARP (Assembly Repeat Annotation Pipeline), run from a Singularity
            .sif under PBS, one `dante_ltr` invocation per assembly
invocation  dante_ltr -o <prefix> -s genome_cleaned.fasta -g DANTE.gff3 \
                      -c 96 -M 1 -S 50000000
genome      GCA_963277665.1 (drVisAlbu1.1), 94,261,041,113 bp — 1,889 chunks
host        96 cores, 768 GB RAM  (job allocations are a fraction of this)
```

## Evidence

### 1. The budget is host-wide

`dante_ltr:296`, the only memory-reading code in the tool:

```python
def _mem_available_kb():
    """Return MemAvailable from /proc/meminfo in kB, or None if it cannot be read."""
    try:
        with open('/proc/meminfo') as f:
            for line in f:
                if line.startswith('MemAvailable:'):
                    return int(line.split()[1])
    except (OSError, ValueError):
        pass
    return None
```

and its single consumer, `dante_ltr:1333`:

```python
        avail_kb = _mem_available_kb()
        pool_size = min(args.cpu, len(remaining)) if remaining else 0
        if remaining and avail_kb and peak_kb > 0:
            mem_cap = max(1, int(avail_kb * 0.8) // peak_kb)
            if mem_cap < pool_size:
                ...
                pool_size = mem_cap
```

`grep -n "os.environ\|getenv" dante_ltr` returns nothing, and there is no CLI
option for a budget, so a caller that *knows* the allocation has no way to say so.

### 2. Why `MemAvailable` is the wrong number in the environment DANTE_LTR runs in

* `/proc/meminfo` is not namespaced. Inside a `.sif` it is the host's.
* The cgroup limit that will actually kill the job is normally set on an
  **ancestor** scope (the job scope), not on the leaf the process sits in — so
  even a naive `/sys/fs/cgroup/memory.max` read at the mount root usually misses
  it.
* Some sites enforce `mem` by polling rather than by cgroup, in which case no
  cgroup limit exists at all and only the scheduler's own variables
  (`PBS_RESC_MEM`, `SLURM_MEM_PER_NODE`, …) carry the truth. Those *do* survive
  the container boundary.
* Even with no limit anywhere, `MemAvailable` on a shared node counts memory that
  neighbouring jobs are about to take.

`pool_size` is **linear** in the budget, so the error passes straight through:
the measured TideCluster #6 case was a 128 GB PBS job on a node reporting
~1.6 TB — a **12.5× oversubscription** of every memory-gated pool.

### 3. Why the gate matters here in particular

On the 94 Gbp genome above, the pre-0.5.3.0 **serial** run — exactly one chunk in
flight — peaked at **108.9 GB** RSS (snakemake benchmark, `-S 50000000`). Chunks
of that class are not cheap, which is precisely why 0.5.3.0 added the gate. With a
budget inflated 12.5× the gate computes a `mem_cap` roughly 12× too permissive
and effectively stops binding, leaving `-c` to set the pool — the situation the
gate was introduced to prevent.

## Proposed change

Mirror what TideCluster 1.20.0 did, so the two tools resolve the same budget from
the same signals and report it the same way.

**1. A `--max_memory` flag** (GB, float, alias `--max-memory`), next to `-c`:

```python
    parser.add_argument(
        '--max_memory', '--max-memory', default=None, required=False, type=float,
        metavar='GB',
        help="Memory available to this run, in GB — the allocation the scheduler "
             "granted. Used to size the chunk pool. Set this on a cluster or in a "
             "container: without it the budget falls back to the scheduler "
             "environment, then the cgroup limit, then /proc/meminfo MemAvailable, "
             "which reports the whole node rather than this job's limit."
        )
```

**2. A resolution chain in place of the bare `MemAvailable` read**, first hit
wins — identical to `tc_utils.memory_budget_mb()`:

| # | source | unit |
|---|--------|------|
| 1 | `--max_memory` | GB |
| 2 | `AGENT_MEMORY` | GB |
| 3 | `PBS_RESC_MEM` (bytes), `SLURM_MEM_PER_NODE` (MB), `LSB_MAX_MEM_RUSAGE` (KB), `SLURM_MEM_PER_CPU` × `SLURM_CPUS_ON_NODE` (MB) | per var |
| 4 | cgroup v2 `memory.max` / v1 `memory.limit_in_bytes`, **walking `/proc/self/cgroup` leaf → root** so an ancestor job scope is found | bytes |
| 5 | `/proc/meminfo MemAvailable` | kB |

Returning `(budget_kb, source)` keeps the change local:

```python
        avail_kb, mem_source = _memory_budget_kb(args.max_memory)
        ...
            print('memory-gating chunk pool to {} concurrent chunks '
                  '(~{} MB/chunk, {} MB available from {})'.format(
                      mem_cap, peak_kb // 1024, avail_kb // 1024, mem_source))
```

**Keep the existing `* 0.8` where it is** and have the helper return the raw
budget — TideCluster's helper applies its headroom internally, so a straight copy
would apply the factor twice here.

**3. Warn when the budget is host-wide under a scheduler.** If the source is
`MemAvailable` and any of `PBS_JOBID` / `SLURM_JOB_ID` / `SLURM_JOBID` /
`LSB_JOBID` is set, print one warning saying so and pointing at `--max_memory`.
This is the part that converts a silent OOM into a first-second diagnostic, and
it is useful even for users who never pass the flag.

## Compatibility and determinism

* On a plain host with no scheduler, container or `AGENT_MEMORY`, every new rung
  misses and the budget is `MemAvailable` exactly as today — behaviour unchanged.
* Where a new source does hit, the only effect is a **smaller pool**. 0.5.3.0
  already keys per-chunk results by index and concatenates in index order, so
  completion order cannot reach the output: this is a scheduling change, and the
  GFF3 stays byte-identical.
* Nothing about the fd-budget guard or the largest-first probe ordering changes.

## Testing

The inputs cannot be reproduced on an ordinary CI runner, so make the helper's
filesystem roots and environment injectable (`environ`, `sysfs_root`,
`proc_cgroup`, `meminfo`) and unit-test it against a synthetic tree: explicit
flag wins; each scheduler variable in its own unit; cgroup v2 on the leaf and on
an **ancestor** scope; cgroup v1 including the `PAGE_COUNTER_MAX` "unlimited"
sentinel and v2 `"max"`; `MemAvailable` fallback; nothing readable; and the
warning firing only under a scheduler on a host-wide budget. That is the shape
CARP uses for its own mirror of this chain
(`scripts/mem_utils.py`, `tests/test_mem_utils.py`), which is available if a
reference implementation is useful.

## Why this exact shape

CARP resolves one budget per run and hands it to every tool that will take it —
`--cap3_max_memory` for DANTE_TIR 0.3.0, `--max_memory` for TideCluster 1.20.0 —
from a single user-facing knob (`max_memory_gb` / `run_pipeline.py -m <GB>`).
DANTE_LTR is the one stage where that value currently has nowhere to go, so it is
also the one stage where a 96-core run on a memory-limited job still sizes its
pool from the node. Matching TideCluster's flag name, units and chain means the
whole stack reports the same number from the same source, and a user who sets the
allocation once gets it honoured everywhere.

## Secondary note (not part of this request)

`utils/mmseq_clustering.R` and `utils/refine_boundaries.py` run `mmseqs
easy-cluster` with `--threads` but no `--split-memory-limit`. mmseqs sizes its own
splits from host RAM and has the same blind spot, but it was not the failure
observed here and is better handled separately if it ever bites.
