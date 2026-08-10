# TideCluster: two parallelism collapses on a 94 Gbp genome

> Feature request for `github.com/kavonrtep/TideCluster` — observed on **1.18.0**.
>
> **Resolved in TideCluster 1.19.0** (2026-08-10), pinned by CARP. All three
> issues below were addressed: `tc_reannotate` now parses and remaps hits in the
> chunk workers with a k-way merge in the parent; TAREAN dispatches longest-first;
> and TAREAN gives large TRCs multiple threads (memory-gated) instead of always
> `-n 1`. Upstream reports the changes byte-identical to 1.18.0 on the
> deterministic outputs. The secondary observation at the end of this document
> (the repeated `there must be an error` message) is **not** covered by 1.19.0 and
> remains open.

## Summary

On a 94.26 Gbp assembly run with `-c 96` on a 96-core / 768 GB host, two
TideCluster stages spend almost all of their wall time single-threaded while 90+
cores idle:

| stage | wall | cores used (mean) | efficiency |
|---|---|---|---|
| `TideCluster.py run_all --long` | **75.3 h** | 9.8 / 96 | 10.2 % |
| `tc_reannotate.py` | **43.3 h** | 2.3 / 96 | 2.4 % |

Neither is a bad pool size — **both pools are correctly sized to 96 at runtime**
(the 1.18.0 memory gate is not the constraint). The loss is elsewhere:

1. **TAREAN is gated by one straggler job.** 441 of 445 TAREAN jobs finish in
   ~2 h; a single TRC then runs **~55 h alone**, because the per-job thread count
   is hardcoded to `1`. That one job is 74 % of the stage's wall time.
2. **`tc_reannotate`'s chunked RepeatMasker is bracketed by serial phases.** The
   96-worker pool accounts for ~1 h of CPU; the other ~42 h is a serial split and
   a serial single-threaded merge over millions of hits.

Together these cost ~10,500 idle core-hours on one genome.

## Environment

```
genome     GCA_963277665.1 (drVisAlbu1.1)  94,261,041,113 bp in 5,096 sequences
host       96 cores, 768 GB RAM, local NVMe scratch
TideCluster 1.18.0   (TideCluster_pipeline_stats.json: "tidecluster_version": "1.18.0")
invocation TideCluster.py run_all -pr TideCluster -c 96 -f genome_cleaned.fasta --long
           tc_reannotate.py -s genome_cleaned.fasta -f ..._dimer_library_reduced.fasta \
                            -o RM_on_TideCluster_Library.gff3 -c 96 --sensitivity rush \
                            --chunk_size 50000000
result     n_trcs_total 1603, n_trcs_above_threshold 453, total_tr_length 305,723,448
```

---

# Issue 1 — TAREAN: one straggler job serialises the stage

## Evidence

Phase timeline reconstructed from output-artifact mtimes (`TideCluster/default/`):

| phase | duration |
|---|---|
| split + TideHunter (3 rounds) | 7.5 h |
| clustering + consensus | 2.3 h |
| kite + stats | 0.5 h |
| **TAREAN** | **56.0 h** |
| dotplots | 2.5 h |
| rDNA + legacy report | 6.4 h |
| **total** | **75.3 h** |

**TideHunter is not the problem** — it parallelised exactly as designed:

```
TideHunter: 2073 parts, pool_size 96 (per-part peak ~1331 MB, budget ~595753 MB, cap 96)
TideHunter: 2073 parts, pool_size 96 (per-part peak ~3187 MB, budget ~595125 MB, cap 96)
TideHunter: 2073 parts, pool_size 96 (per-part peak ~2674 MB, budget ~596147 MB, cap 96)
```

6,219 single-threaded instances (3 rounds × 2,073 parts), ~93 s each, 7.5 h total.
The 1.18.0 memory gate behaved correctly and never bound.

**TAREAN job completion times** (`TideCluster_tarean/*_tarean` directory mtimes,
445 jobs):

| completed | jobs |
|---|---|
| 2026-08-01 02:xx | 210 |
| 2026-08-01 03:xx | 231 |
| 2026-08-01 14:19 | 1 &nbsp;(`TRC_8`) |
| **2026-08-03 10:32** | **1 &nbsp;(`TRC_1`)** |

441 jobs finish inside ~2 h. `TRC_1` then runs **alone for ~55 hours** with 95 of
96 cores idle. Its input is a extreme outlier:

```
TRC_1_dimers.fasta   21,677,078 B      <- largest
median over 3,212 consensus fastas     7,194 B      (TRC_1 is ~3,013x the median)
```

## Root cause (`TideCluster.py`, `tarean()`)

The per-job thread count is a hardcoded literal — the function's `cpu` argument
never reaches `tarean.R`:

```python
# TideCluster.py, tarean(prefix, gff, fasta=None, cpu=4, ...)   line ~118
tarean_out = F"{tarean_dir}/{v_basename}_tarean"
cmd = F"{script_path}/tarean/tarean.R -i {v} -s 0 -n {1} -o {tarean_out}"
#                                                    ^^^^^  always 1
cmd_list.append(cmd)
...
# line ~335
with Pool(cpu) as p:                       # pool IS 96 - that part is fine
    for _ in p.imap(tc.run_cmd, cmd_list):
```

`Pool(96)` is correct while there is work to spread. Once the queue drains to the
last job, the stage degrades to one core, and `-n 1` means that job cannot use the
95 free ones. Because TAREAN cost scales steeply with array size and TRC sizes span
~3,000×, the largest TRC dominates the whole stage on any large genome.

`imap` also dispatches in `fasta_dict` order rather than largest-first, so the
straggler may start late — but ordering alone cannot fix a job that is longer than
all others combined.

## Proposed fix (layered — recommend 1 + 2)

**1. Pass the real thread count through (the actual fix).** Give `tarean.R` a
share of the cores, so the straggler is not stuck on one:

```python
cmd = F"{script_path}/tarean/tarean.R -i {v} -s 0 -n {threads_per_job} -o {tarean_out}"
```

The simplest correct policy is a two-phase run: pool the small jobs at `-n 1`, then
run the few largest **sequentially with `-n cpu`**. Even a static
`threads_per_job = max(1, cpu // pool_size)` helps once the queue is short. This
requires `tarean.R` to honour `-n`; if it does not scale internally, see 3.

**2. Dispatch longest-job-first (cheap, no behaviour change).** Sort `cmd_list` by
the TRC's total array length descending before `imap`. Output is per-job and
independent, so ordering affects only scheduling and the progress counter. This is
the standard LPT fix and it removes the "straggler starts last" case at near-zero
cost. (CARP does exactly this in `repeatmasker_wrapper.py` and it holds 96/96
utilisation.)

**3. If `tarean.R` cannot be parallelised internally**, expose a size cap /
subsample for pathological TRCs (e.g. `--tarean-max-input-bp`, analysed on a
representative subsample with the cap recorded in the report), so one megasatellite
cannot set the runtime of the whole stage. Silent truncation would be worse than
the current behaviour, so this must be logged.

Expected effect on the reported run: TAREAN **56 h → a few hours**, i.e.
`run_all --long` from 75.3 h to under 20 h.

---

# Issue 2 — `tc_reannotate`: serial split and serial merge around a correct pool

## Evidence

```
Chunked RepeatMasker: 1855 chunk file(s), 6876 piece(s), 96 worker(s)
Chunked RepeatMasker: 1855 chunk(s) completed -- 1855 with hits, 0 hit-less, 0 failed
```

The pool is correctly sized to 96 and every chunk succeeded. Yet:

```
wall            43.3 h
cpu_time        365,743 s  (101.6 core-hours)
mean cores      2.3 / 96          (2.4 %)
io              2.10 TB written, 499 GB read
```

All 101.6 core-hours of CPU are consistent with the pool phase alone, which at 96
workers is **~1.1 h of wall**. The remaining **~42 h consumed almost no CPU** — it
is serial, I/O- and single-thread-bound work before and after the pool.

## Root cause (`tc_utils.py`, `run_repeatmasker_genome_chunked`)

The function is `serial split → serial warmup → parallel pool → serial merge`:

```python
# line ~1219  SERIAL: splits 94 Gbp into 1,855 chunk files (writes ~94 GB single-threaded)
chunk_files, matching_table = split_fasta_to_chunk_files(
    ref_seq, work_dir, chunk_size=chunk_size, overlap=overlap)

_repeatmasker_warmup(rm_library, sensitivity_flag, work_dir)   # serial, and necessary

pool_size = max(1, min(cpu, len(tasks)))                       # = 96, correct
with Pool(pool_size) as pool:
    results = list(pool.map(_repeatmasker_chunk_worker, tasks))

# line ~1266  SERIAL: re-reads every chunk .out in the parent, one record at a time
records = []
for status, _chunk, out_file in results:
    with open(out_file) as f:
        ...
        for line in f:
            items = line.split()
            row = token_row.get(items[4])
            records.append((ori_header, start, end, strand, name))

records.sort(...)
with open(out_gff3, "w") as gff3_out:
    for seqid, start, end, strand, name in records:
        gff_line = (F"{seqid}\tRepeatMasker\ttandem repeat\t{start}...")
        gff3_out.write(Gff3Feature(gff_line).print_line())   # object per record
```

Two serial bottlenecks scale linearly with genome size:

- **The split** writes ~94 GB through one thread before any worker starts.
- **The merge** re-reads all 1,855 chunk `.out` files in the parent, accumulates
  every hit in one list, sorts it, and writes the GFF3 constructing a
  `Gff3Feature` object per record. On a genome with 305 Mb of tandem arrays this
  is tens of millions of records on one core.

## Proposed fix

**1. Parse in the workers, not the parent.** `_repeatmasker_chunk_worker` already
owns its chunk's `.out`; have it do the `token_row` lookup and coordinate offset
and return the parsed tuples. The parent then only concatenates, sorts and writes.
This moves the dominant serial pass into the existing 96-way pool with no change
to output — the final `records.sort()` already makes the result order-independent.

**2. Avoid the per-record object in the writer.** `Gff3Feature(gff_line).print_line()`
constructs and re-parses an object per hit; formatting the line directly is
equivalent and removes tens of millions of allocations.

**3. Parallelise or stream the split.** `split_fasta_to_chunk_files` can write
chunks from a worker pool, or the pool can start on chunk *k* as soon as it is
written rather than waiting for all 1,855.

Expected effect: `tc_reannotate` **43.3 h → ~2–4 h** on this genome.

### Note on the `-num_threads 1` interaction

The caller (CARP, `tidecluster_reannotate_culling_limit: 2`) installs an `rmblastn`
shim that appends `-culling_limit 2` and forces `-num_threads 1`, because BLAST
culling is non-deterministic across threads. That is deliberate and is *not* the
problem here — the chunk pool is the intended replacement for intra-RepeatMasker
threading and it works. It does mean the pool is the **only** source of parallelism
in this stage, so the serial phases around it are fully exposed.

---

# Secondary observation (not performance)

`tidecluster_long.err` contains, **115 times**:

```
there must be an error: 11634 deleted from 22991 that now is empty, but not assigned to a cluster
```

The run completed and produced plausible output (1,603 TRCs, 305.7 Mb of tandem
repeat), so this is not fatal, but a message that self-describes as an error and
repeats 115 times on a large genome is worth either fixing or downgrading to an
explained warning. Happy to open this separately with more context if useful.
