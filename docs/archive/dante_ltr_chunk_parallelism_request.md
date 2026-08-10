# DANTE_LTR: the chunk loop is serial — 1 core of 96 on a 94 Gbp genome

> Feature request for `github.com/kavonrtep/dante_ltr` — observed on **0.5.2.0**.
>
> **Resolved in dante_ltr 0.5.3.0**, pinned by CARP. All four proposals below were
> implemented: the chunk loop runs in a `Pool` with `imap_unordered` and `-c 1`
> children; chunks are ordered largest-first; the largest chunk is run serially
> first as a memory probe and the pool is then capped at
> `min(cpu, 0.8 * MemAvailable / per-chunk peak RSS)`; and the `-S` help text now
> states that it controls pool granularity as well as per-chunk memory. The change
> is confined to that loop — the fd-budget guard from the too-many-open-files fix
> is retained, and per-chunk outputs are still concatenated in index order, so
> completion order cannot reach the output.

## Summary

On a 94.26 Gbp assembly, `dante_ltr -c 96` ran for **24.2 hours using 1.0 core**.
The genome is split into 1,889 chunks and the chunks are then consumed by a
**strictly serial `for` loop** that blocks on `subprocess.check_call` per chunk.
`-c` is forwarded to each `detect_putative_ltr.R` invocation, but at this chunk
size that does not translate into concurrency — the measured average is 0.98
cores across the whole run.

This is the single cheapest large-genome win in the pipeline: the work is already
partitioned into 1,889 independent units, and nothing about them is order- or
state-dependent (they are concatenated afterwards). Running the existing loop in a
process pool should take the stage from **24.2 h to well under 1 h**.

There is an irony worth flagging: the chunk count comes from `-S`, which callers
raised to work around the `Too many open files` bug
(`dante_ltr_too_many_open_files_request.md`). **That workaround is what multiplies
the serial iterations** — the smaller the chunk, the more serial passes.

## Evidence (DANTE_LTR 0.5.2.0, 94.26 Gbp, 96 cores / 768 GB)

Invocation (from the CARP `dante_ltr` rule):

```
dante_ltr -o DANTE_LTR/DANTE_LTR -s genome_cleaned.fasta -g DANTE_filtered.gff3 \
          -c 96 -M 1 -S 50000000
```

Snakemake benchmark for the rule:

```
wall        87,118.97 s   = 24.2 h
cpu_time    85,134.53 s   = 23.6 core-hours
mean cores  85,134 / 87,119 = 0.98  of 96 available   (1.0 % efficiency)
max_rss     108.9 GB
io_out      434 GB
```

`cpu_time ≈ wall` is the signature of a single-threaded run: over 24 hours the
process never averaged more than one core.

The log is a flat serial sequence of 1,889 chunks, each a complete
read-gff → read-fasta → cluster → detect cycle:

```
Running retrotransposon detection on file 0
...
Running retrotransposon detection on file 1888
DANTE LTR-RTs finished
```

Every per-chunk stage line appears exactly 1,889 times, confirming one full pass
per chunk with no interleaving:

```
   1889  reading gff...done
   1889  reading fasta...done
   1889  Total number of domain clusters
   1889  Running retrotransposon detection on file
   1889  Number of putative TE with identified LTR
   1889  Number of complete TE found
```

24.2 h / 1,889 chunks ≈ **46 s per chunk, one at a time**.

## Root cause (`dante_ltr`, `main()`, ~line 1244)

The chunks are produced in parallel-ready form and then consumed serially:

```python
# run retrotransposon detection on each temporary file
output_files = make_temp_files(number_of_temp_files)
print('number of chunks:": ' + str(number_of_temp_files))
for i in range(number_of_temp_files):                       # <-- strictly serial
    print('Running retrotransposon detection on file ' + str(i))
    cmd = [f'{tool_path}/utils/detect_putative_ltr.R', '-s', temp_files_fasta[i],
         '-g', temp_files_gff[i], '-o', output_files[i], '-c', str(args.cpu), '-M',
         str(args.max_missing_domains), '-L', str(args.min_relative_length)]
    if constrains_path:
        cmd += ['--te_constrains', constrains_path]
    if args.no_ambiguous_domains:
        cmd += ['--no_ambiguous_domains']
    subprocess.check_call(cmd)                              # <-- blocks per chunk
```

`args.cpu` *is* passed to the R detector as `-c`, so the intent is clearly that
parallelism happens inside `detect_putative_ltr.R`. The measurement says it does
not materialise: whatever `-c` parallelises inside the R script, it is not enough
to move the average above one core on 50 Mb chunks. Two candidate explanations —
the maintainer is better placed to say which:

- the R detector's parallel section is a small fraction of per-chunk work
  (Amdahl), so 96 threads buy nothing at this chunk size; or
- `-c` only reaches a sub-step that is not on the hot path here.

Either way, the outer loop is where the 1,889-fold parallelism actually exists,
and it is unused.

The chunks are fully independent: each gets its own FASTA/GFF pair and its own
output file, and results are only combined afterwards by plain concatenation of
the per-chunk suffixed outputs (`_D.fasta`, `_DL.fasta`, `.gff3`,
`_statistics.csv`, …).

## Proposed fix

**1. Run the chunk loop in a process pool with single-threaded children.**
This is the whole fix — the loop body is already a self-contained subprocess call:

```python
from multiprocessing import Pool

def _chunk_worker(task):
    subprocess.check_call(task)

cmds = []
for i in range(number_of_temp_files):
    cmd = [f'{tool_path}/utils/detect_putative_ltr.R', '-s', temp_files_fasta[i],
           '-g', temp_files_gff[i], '-o', output_files[i], '-c', '1', '-M',
           str(args.max_missing_domains), '-L', str(args.min_relative_length)]
    ...
    cmds.append(cmd)

pool_size = max(1, min(args.cpu, len(cmds)))
with Pool(pool_size) as p:
    for _ in p.imap_unordered(_chunk_worker, cmds):
        pass
```

Output is unchanged: the concatenation step downstream already walks
`range(number_of_temp_files)` in index order, so completion order is irrelevant
and the result stays byte-identical to the serial run.

**2. Dispatch longest-first.** Chunks are built by round-robin over
size-sorted sequence ids, so they are roughly equal — but on assemblies with a few
huge scaffolds they will not be. Sorting `cmds` by input FASTA size descending
before `imap_unordered` costs nothing and removes the straggler tail.

**3. Gate the pool on memory, not just cores.** This run peaked at 108.9 GB RSS
*serially*; N concurrent chunks cost ~N × per-chunk peak. A pool sized purely from
`-c` could OOM on a smaller host. The sibling project already has the pattern worth
copying — TideCluster 1.18.0's `_tidehunter_memory_budget_mb()` measures the first
part's peak RSS, reads `MemAvailable`, and sizes the pool as
`min(cpu, budget // per_part_peak)`. Something equivalent here would make `-c 96`
safe by default.

**4. Consider decoupling `-S` from the parallel width.** `-S` is documented as the
per-chunk *memory* control for the R detector, but it currently also sets the
serial iteration count, so tuning memory down makes the run slower rather than just
leaner. Once the loop is pooled this mostly resolves itself, but it is worth stating
in the help text that `-S` now controls chunk memory *and* pool granularity.

## Expected effect

At 23.6 core-hours of real work and 96 cores, a pooled loop bounded by memory
rather than by the loop gives roughly:

| | wall |
|---|---|
| today (serial) | 24.2 h |
| pooled at 96, perfect scaling | ~0.25 h |
| pooled, memory-gated to ~8–16 concurrent chunks (108.9 GB peak, 768 GB host) | **~1.5–3 h** |

Even the conservative memory-gated figure is an **8–16× reduction** in wall time
for this stage, with no change to output.

## Pipeline-side stopgap

None available. The serialisation is inside `dante_ltr`; a caller can only trade
chunk size against memory via `-S`, and lowering `-S` makes it worse (more serial
iterations). Raising `-S` reduces iterations but raises per-chunk RSS above the
108.9 GB already observed and re-approaches the open-files limit this flag was
lowered to avoid. This one needs the upstream change.
