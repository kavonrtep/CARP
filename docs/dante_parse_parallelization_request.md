# DANTE: single-threaded post-search merge dominates wall-time on large genomes

> Feature request for `github.com/kavonrtep/dante` — observed on **0.2.5**.

## Summary

On a large genome (~90 Gbp, `-c 32`), DANTE spends a first phase with all cores
busy (the per-chunk `lastal` search + parse, run in a `Pool`) followed by a long
**single-threaded phase** that recombines the per-chunk results. The serial
phase — per-chunk coordinate recalculation, a global `sort`, and a dedup pass —
runs on **one core** and approaches the search time, so overall CPU utilisation
collapses in the back half of the run. Two design choices amplify it: the
per-chunk `lastal` search is run **twice** (TAB + MAF), and the chunk size is
hard-coded to **500 kb**, producing ~180,000 chunks on a 90 Gbp genome and thus
a 180,000-iteration serial merge.

## Where the time goes (`dante` wrapper + `dante.py`, 0.2.5)

**Parallel phase (good).** For a genome exceeding the chunk size the wrapper
splits into `number_of_temp_files = int(total_size / chunk_size) + 1` chunks and
runs the full per-chunk pipeline in a process pool:

```python
# dante  (wrapper), line 519
chunk_size = 500000                       # 500 kb  ← ~180,000 chunks on 90 Gbp
...
# line 571
with Pool(args.cpu) as p:
    p.map(dante.main, args_list)          # each worker: lastal search + streaming parse
```

Each `dante.main` runs `lastal` and parses its stdout **in the same worker**
(`dante.py:domain_search`, streaming parse at `dante.py:505-533`), so the parse
*is* parallel across chunks. That is not the bottleneck.

**Serial phase (the problem).** After the pool returns, the wrapper does all of
this on one core (`dante`, lines 576-607):

```python
# 1) per-chunk coordinate recalculation — serial loop over ALL ~180k chunks
for args_part in args_list:
    tmp_gff = recalculate_gff3_back_to_original_coordinates(args_part.domain_gff, ...)
    # ... concatenate into one big unfiltered gff
# 2) single-threaded global sort of the whole concatenated domain GFF
p = subprocess.Popen(['sort', '-k1,1', '-k4,4n'], stdin=f_tmp, stdout=subprocess.PIPE)
# 3) serial dedup pass over the sorted stream
for line in p.stdout:
    ... drop line if (seqid,start,end) == previous ...
```

On a repeat-dense 90 Gbp genome the concatenated domain GFF holds tens of
millions of records; recalculating coordinates for 180k chunk files, then a
single-threaded `sort`, then a serial Python dedup, is the single-core tail the
user observes.

**Amplifier — `lastal` is run twice per chunk** (`dante.py:491-503`): one
invocation with `-f TAB` (scores/positions) and one with `-f MAF` (alignment
sequence), each a full search over the same chunk. `lastal` is given no thread
flag (the `-p` is the scoring matrix), so parallelism comes only from the chunk
pool — and the search cost is paid twice.

## Proposed fixes (in impact order)

> **Status:** fixes 1 and 2 are **implemented and verified** (branch
> `fix/parallelize-merge-tail`): the recalc now runs across `args.cpu` workers
> and the merge `sort` uses `--parallel`. A/B on a 4 Mb chunked run (8 chunks)
> is **byte-identical** to stock with `PYTHONHASHSEED` fixed. Fixes 3–4
> (chunk-size scaling, single `lastal`) are left as suggestions — bigger
> behavioural changes needing their own validation. **Aside:** DANTE's
> `Region_Hits_Classifications` attribute is built from a Python `set`, so its
> internal order varies run-to-run with `PYTHONHASHSEED` (same domains/coords,
> just reordered list items) — pre-existing, harmless, but it makes DANTE output
> non-reproducible byte-for-byte across invocations; sorting that list before
> emit would make runs reproducible.

1. **Do coordinate recalculation inside the pool workers.** Have each
   `dante.main` (or a thin wrapper around it) emit its GFF already in original
   coordinates, so the serial `for args_part in args_list:
   recalculate_...` loop (line 579) disappears entirely. This is the largest,
   cleanest win — it moves the 180k-iteration loop onto all cores.

2. **Parallelise / bound the merge sort.** Pass `--parallel={cpu}` and a large
   `-S` buffer to `sort` (line 594), and point `TMPDIR` at fast scratch. Better:
   duplicates can only occur in the small `overlap` (2 kb) region at chunk
   boundaries, so a **global** sort+dedup over the whole genome is unnecessary —
   dedup can be restricted to boundary intervals, removing the whole-genome sort.

3. **Scale chunk size with genome size instead of the fixed 500 kb** (line 519).
   Target a bounded chunk count (e.g. a few × `cpu`, or a floor like 10–50 Mb)
   so a 90 Gbp genome makes thousands, not ~180k, chunks. Fewer chunks ⇒ far
   less process-spawn / DB-mmap overhead, fewer temp files (see note), and a much
   shorter serial merge. (Chunks are independent, so larger chunks don't change
   results; only per-chunk peak memory grows, which is modest for `lastal`.)

4. **Run `lastal` once, not twice** (`dante.py:491-503`). If both TAB and MAF are
   needed, derive one from the other or use a single format that carries both
   the coordinates and the aligned sequence, halving search time.

## Notes

- **Not an fd bug.** DANTE already avoids the `dante_ltr`-style "too many open
  files" failure — line 544 explicitly switched the chunk split to open-append
  per sequence (`with open(path, 'a')`, line 552). Good. But
  `make_temp_files(180000)` still creates ~180k temp files in `TMPDIR`; fix 3
  (bigger chunks) also relieves that inode/space pressure.
- **Pooled workers keep `args.cpu = 32`** (deepcopy at `dante:561` doesn't reduce
  it). Harmless today because `lastal` isn't given `-P`, but if any per-chunk
  step later uses `args.cpu` it would oversubscribe `cpu²` threads — worth
  setting `args_part.cpu = 1` for pooled calls.

## Pipeline-side stopgap

Limited options without an upstream change. The pipeline can raise the wrapper's
chunk size via the tool's window/db mechanism only indirectly; the practical
lever now is ensuring `TMPDIR` is on a large, fast filesystem so the 180k temp
files and the global `sort` spill efficiently. The durable fixes (1–4) are
upstream.
