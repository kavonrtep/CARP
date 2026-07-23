# Pipeline Performance / Parallelization Report

**Scope:** Snakemake rules that call **repo-owned scripts** (`scripts/*.py`, `scripts/*.R`)
and can be sped up **without changing outputs**. External tools wrapped by the
pipeline (DANTE, DANTE_LTR, DANTE_TIR, TideCluster, RepeatMasker itself) are only
discussed where *our* wrapper/scheduling controls their efficiency.

**Hard constraint (from the request):** every change must be output-preserving
(byte-identical GFF3 / FASTA / masked-bp). All recommendations below are scheduling,
parallelism, memory, or redundant-work changes — **no algorithmic changes** that
alter results. Validation guidance is in the last section.

**This is analysis only — no code was modified.**

---

## 1. Evidence base — two real runs

Benchmarks were read directly from the pipeline's own `benchmarks/*.tsv`.

| Run | Genome | Scaffolds | `--cores` | Source |
|-----|--------|-----------|-----------|--------|
| CARP | 781 Mb | 152 | 16 | `tmp/carp_output/benchmarks/` |
| **Pisum Pearl** | **3.93 Gb** | **3317** | **50** | `/mnt/ceph/.../PBA-Pearl_v1.2.1_2026_05_27/.../output/benchmarks/` |

The Pisum run is the most informative (largest, most cores). Wall-time + the
`mean_load` column (snakemake's avg CPU%, where 100 = 1 core busy) on **50 cores**:

| Rule | Wall | mean_load (cores) | max_rss | Backed by repo script? |
|------|-----:|:--:|--------:|:--|
| tidecluster_reannotate | 16 h 47 m | ~1.1 / 50 | 28 GB | partly (`reduce_dimer_library.py`) + external `tc_reannotate` |
| **repeatmasker** | **15 h 15 m** | **10.9 / 50** | 195 GB | **`repeatmasker_wrapper.py`** ✅ |
| tidecluster_long | 14 h 01 m | ~1.5 | 59 GB | external (TideCluster) |
| dante_ltr | 7 h 08 m | ~0.7 | 137 GB | external |
| **dante_tir_fallback** | **4 h 00 m** | **15 / 50** | 19 GB | **`dante_tir_fallback.py`** ✅ |
| dante_tir | 3 h 20 m | ~1.9 | 281 GB | external |
| dante | 1 h 55 m | ~0.2 | 0.6 GB | external |
| **reduce_library** | **44 m 45 s** | **1.2 / 50** | 8 GB | **`reduce_library_size.py`** ✅ |
| **make_bigwig_density** | **38 m 08 s** | **~1.0** | 4 GB | **`calculate_density_batch.R`** ✅ |
| **make_repeat_report** | **24 m 13 s** | **~1.0** | 2 GB | **`make_repeat_report.R`** ✅ |
| **merge_rm_and_dante** | **18 m 15 s** | **4.4 / 50** | **296 GB** ⚠️ | **`merge_repeat_annotations.R`** ✅ |
| make_summary_plots | 10 m 13 s | ~1.0 | 2 GB | `make_summary_plots.R` ✅ |
| make_unified_annotation | 6 m 39 s | ~0.9 | 4 GB | `make_unified_annotation.R` ✅ |
| calculate_bigwig_density | 2 m 30 s | 0.4 | 2.5 GB | `calculate_density.R` ✅ |
| make_summary_statistics… | 2 m 22 s | ~1.0 | 8.8 GB | `calculate_statistics_and_make_groups.R` ✅ |

**The recurring story is the `mean_load` column.** On a 50-core box, our heaviest
repo script (`repeatmasker`) averaged **~11 cores**, and almost every other repo
script averaged **~1 core**. The cores are there; the code isn't using them. That is
the single biggest, safest source of speed-up in this pipeline.

---

## 2. Top findings (ranked by actionable impact)

### 🥇 FINDING 1 — `repeatmasker_wrapper.py`: tiny scaffolds are not packed → thousands of redundant RepeatMasker launches
**Rule:** `repeatmasker` (15 h 15 m, **10.9 / 50 cores**, the dominant actionable cost)
**File:** `scripts/repeatmasker_wrapper.py`

What happens now (`split_fasta_to_chunks`, lines 46–94):
- Sequences **>10 Mb** (`min_chunk_size = chunk_size*2`, chunk_size=5 Mb at line 278) are split into ~5 Mb chunks. ✔ correct.
- Sequences **≤10 Mb** each become **their own one-sequence chunk** (lines 79–81). The
  docstring promises *"If sequences are shorter, chunk with multiple sequences are
  created"* (line 48) — **but that packing is not implemented.**

Consequence on Pisum (verified from the log — `Number of fasta parts: 4068`):
- 7 chromosomes → ~750 chunks of real work.
- **3310 small scaffolds → 3310 separate chunks** (3179 are <100 kb; the smallest are 1 kb).
- `split_fasta_to_files` (line 238) then writes **one FASTA file per chunk → one
  `RepeatMasker` process per chunk** (`multiprocessing.Pool`, line 295).
- **Every one of those 4068 processes re-runs RepeatMasker's full library setup**
  (re-formats / `makeblastdb` the 3569-seq library) before masking. ~3300 of those
  preps are pure waste on kilobase scaffolds. This is the bulk of why 50 cores
  average 10.9: most processes spend their life in single-threaded library prep, and
  a long tail of repeat-dense 5 Mb chunks runs on a handful of cores for hours.

**Output-preserving fixes (RepeatMasker treats each sequence in a multi-FASTA
independently, so batching ≠ result change):**
1. **Pack small sequences into ~5 Mb multi-sequence chunks** (deliver the documented
   behavior). ~3310 tiny chunks collapse to ~100–200 packed chunks → thousands of
   redundant library preps eliminated and far better core packing. **Highest impact.**
   `recalculate_rm_out_coordinates` already maps results per original header, so
   coordinate handling is unaffected — but a packed (multi-seq) chunk needs the
   matching-table/offset logic to treat each packed sequence with offset 0 (each small
   scaffold is whole), which is straightforward.
2. **Adaptive chunk size** = `max(5 Mb, ceil(total_bp / (k·threads)))` (e.g. k≈3) so
   the number of *big-scaffold* chunks tracks the core count: enough for load balance,
   not 10× more than needed. On few-core machines this avoids over-chunking the big
   chromosomes; on Pisum/50-core it stays ≈ current.
3. **Format the library once and reuse** across chunks instead of per-process. Higher
   effort (RepeatMasker's `-lib` setup is per-run); investigate `-lib` with a
   pre-built `*.nsq`/consensus next to the lib, or a single warm-up run, then verify
   byte-identical `.out`. Mark as *investigate*, not *assume-safe*.
4. **Scalability hygiene:** the per-header filter `[x for x in matching_table if
   x[0]==header]` (line 90) is O(headers × chunks) and the genome is read fully into a
   dict **twice** (lines 86 and 245, the latter inside `split_fasta_to_files`). Build a
   `header → rows` dict once and stream the split. Negligible today, but on a
   100k-scaffold assembly the O(N²) scan and double in-RAM genome bite hard.

**Safety:** (1),(2),(4) are pure batching/scheduling → byte-identical `.out`. (3) needs
a masked-bp parity check before trusting.
**Expected impact:** Large. Removing ~3300 redundant library preps and packing the
tail should both cut wall-time and lift mean_load well above 11; the floor becomes the
genuine search cost of the big chromosome chunks.

---

### 🥈 FINDING 2 — `merge_repeat_annotations.R`: `detectCores()` override forks ~all physical cores → 296 GB RSS
**Rule:** `merge_rm_and_dante` (18 m, **296 GB max_rss**, OOM risk)
**File:** `scripts/merge_repeat_annotations.R`, function `gff_cleanup`

```r
33  num_cores <- as.integer(Sys.getenv("CPU_COUNT"))   # rule exports CPU_COUNT={threads}
34  if (is.na(num_cores)) num_cores <- detectCores()
35      num_cores <- detectCores()                     #  ← unconditional override
37  num_cores <- detectCores()                         #  ← AGAIN: kills CPU_COUNT entirely
38  gff_names   <- mclapply(..., mc.cores = round(num_cores*0.8))
39  gff_strands <- mclapply(..., mc.cores = round(num_cores*0.8))
```

- Line 37 throws away the `CPU_COUNT` the Snakefile carefully exports
  (`Snakefile:1326 export CPU_COUNT={threads}`). The env plumbing is **dead code.**
- `detectCores()` returns the **physical machine** core count, not the cores the job
  was given. On the HPC node this forks `0.8 × (all node cores)` workers — far more
  than the 50 the job owns → oversubscription + context-thrash.
- `mclapply` **forks** the parent; each worker COW-shares the large RM `GRanges`. With
  ~100 forks over a 3.9 Gb-genome annotation, reported RSS explodes to **296 GB**
  (the exact copy-on-write over-count that motivated rewriting `reduce_library` from R
  to Python — see `CLAUDE.md`). Real risk: OOM-kill on a shared node or a bigger genome.
- The work is done **twice** (lines 38 *and* 39: one `mclapply` for names, one for
  strands, over the *same* `gff_disjoin$revmap`).

**Output-preserving fixes:**
- Honor the budget: delete the line-35/37 overrides; use
  `num_cores <- CPU_COUNT (env) else detectCores()`. `mc.cores` never affects results.
- Collapse the two `mclapply` passes into **one** that returns `list(name=…, strand=…)`
  per disjoint range → half the forks, half the COW pressure.
- Consider `mc.cores` capped to the real budget (e.g. `min(CPU_COUNT, …)`); for this
  largely vectorizable disjoin/revmap step even a modest cap cuts the memory blow-up.

**Safety:** `mc.cores` / loop fusion are scheduling only → identical GFF3.
**Expected impact:** High on memory (296 GB → a few GB, removes OOM risk); moderate on
wall-time (less thrash, one pass instead of two).

---

### 🥉 FINDING 3 — `clean_rm_output.R`: hard-coded `mc.cores = 8` + double `mclapply`, hidden inside the 15 h RepeatMasker rule
**Rule:** `repeatmasker` (runs **after** the wrapper, inside the same rule)
**File:** `scripts/clean_rm_output.R`, lines 12–13

```r
12 gff_names   <- mclapply(as.list(gff_disjoin$revmap), function(x) gff$Name[x],   mc.cores = 8)
13 gff_strands <- mclapply(as.list(gff_disjoin$revmap), function(x) strand(gff[x]), mc.cores = 8)
```

- `mc.cores = 8` is hard-coded — ignores the rule's `{threads}` (16 on CARP, 50 on
  Pisum). Same COW-fork pattern as Finding 2, same redundant **double** pass.
- Also parses the (huge) RepeatMasker `.out` with base-R `read.table(..., fill=TRUE,
  col.names=paste0("V",1:16))` (line 54) — slow and memory-heavy for millions of rows.
- This time is **invisible** because it is folded into the `repeatmasker` benchmark
  (15 h), but it sits on the critical path right after masking.

**Output-preserving fixes:** take `mc.cores` from an arg/env (`CPU_COUNT`) instead of
`8`; fuse the two `mclapply` passes into one; optionally parse the `.out` with
`data.table::fread` (if available in the env) for a large parse speed-up.
**Safety:** scheduling / parser swap → identical GFF3 (verify fread parity).
**Expected impact:** Medium; shaves the post-mask tail of the most expensive rule.

---

### FINDING 4 — `reduce_library_size.py`: Phase 1 is strictly one-class-at-a-time
**Rule:** `reduce_library` (44 m 45 s, **1.2 / 50 cores**)
**File:** `scripts/reduce_library_size.py`

```python
630  big_done   = _run_phase(big_jobs,   threads=1,            …)  # Phase 1: ONE class at a time
634  small_done = _run_phase(small_jobs, threads=args.threads, …)  # Phase 2: small CAP3 in parallel
649  run_blast_filter(jobs, threads=args.threads)                  # Phase 3: serial loop
```

Why mean_load ≈ 1.2 (the math matches the run): Phase 1 (`ProcessPoolExecutor(max_workers=1)`,
line 578+630) runs the memory-heavy classes **sequentially** —
- big **CAP3** classes are **single-threaded by nature** (CAP3 has no thread flag,
  lines 68/238) → ~1 core each, and they run one after another;
- **mmseqs** classes each get the full thread budget (lines 622–624) but still run
  one at a time.
So during most of the 44 min, exactly one class is active. Phase 2 (parallel) and the
small classes finish fast; **Phase 3's BLAST filter is a serial `for` loop** (line 648–649
→ `run_blast_filter`), one `makeblastdb`+`blastn` at a time.

The two-pass design exists for a real reason — mmseqs has a ~8 GB nucleotide memory
floor and 4× parallel mmseqs OOM-killed 16 GB CI runners (documented at lines 549–557).
So **do not** parallelize mmseqs. But:

**Output-preserving fixes (parity test `tests/test_reduce_library_parity.sh` must pass):**
- **Run the big *CAP3* classes in parallel** within a memory budget (CAP3 peaks at a
  few hundred MB, not GB) while keeping **mmseqs serial**. Split Phase 1 into
  "big-CAP3 (parallel, N workers)" + "mmseqs (sequential)". This is where most of the
  44 min goes when an assembly has several large LTR classes.
- **Parallelize Phase 3's BLAST filter** across classes with a small pool
  (`max_workers≈4`, `threads = budget/4`). BLAST is deterministic → identical drops.
- Auto-tune / raise `--max-parallel-bp` on high-RAM machines so fewer CAP3 classes are
  forced into the sequential bucket (it's a pure scheduling knob).

**Safety:** all are scheduling — deterministic tools, separate workdirs, reassembled in
original class order (lines 638–644). Parity test guards it.
**Expected impact:** Medium-high (Phase 1 is the bulk of the 44 m and is ~1 core).

---

### FINDING 5 — `dante_tir_fallback.py`: per-subtype work is serialized
**Rule:** `dante_tir_fallback` (4 h, **15 / 50 cores**)
**Files:** `scripts/dante_tir_fallback.py`, `scripts/global_local_aln.py`

- The main loop processes each TIR subtype **sequentially** (`for subtype in
  sorted(...)`, ~line 1065), passing the full thread budget to each `process_subtype`.
  Subtypes are independent (separate output dirs, merged afterward) → they can run
  concurrently with a **divided** thread budget (`threads // n_subtypes`). The 15/50
  average is consistent with one subtype's inner alignment being partially parallel
  while the rest of the pipeline idles.
- `global_local_aln.py` parallelizes the all-vs-all alignment with a
  **`ThreadPoolExecutor`** (line 829). The inner work is `parasail` SIMD alignment.
  **Verify whether this parasail build releases the GIL** during alignment: if yes,
  threads are fine; if not, CPU-bound alignment is GIL-serialized and a
  `ProcessPoolExecutor` would unlock real parallelism. (Don't assume — measure.)
- O(N²) hot spot: per-anchor flank computation re-filters all features by seqname and
  scans the list **four times** per anchor (`create_prime_bed_files`, ~lines 243–288).
  Pre-build a `seqname → [features]` dict once and merge the passes. Output identical.
- Minor: the 5′ and 3′ `seqkit subseq` calls run back-to-back (~lines 358–370); they
  are independent and could overlap.

**Safety:** subtype parallelism + dict pre-grouping + thread/process choice are all
output-neutral (results are merged/serialized after).
**Expected impact:** Medium-high (4 h rule; subtype parallelism + a correct
thread/process choice could roughly halve it on multi-subtype genomes).

---

### FINDING 6 — `make_bigwig_density` / `calculate_density_batch.R`: one dominant class file starves the file-level parallelism
**Rule:** `make_bigwig_density` (38 m, **~1 core**)
**File:** `scripts/calculate_density_batch.R`, line 92

`mclapply(files, …, mc.cores = threads)` parallelizes **across input files** — but the
per-class split dir has only ~10–15 files and **one of them (`Mobile_elements`, the
bulk of the genome's repeats) dominates**. The 14 small files finish in seconds; the
one giant file runs ~38 min on a single core → mean_load ≈ 1.

**Output-preserving fix:** parallelize **within** the dominant file. `density_per_family`
(lines 48–60) is per-seqname independent after `coverage()`; split the per-seqname
`binnedAverage`+smoothing across cores, or shard the big GFF by chromosome, compute,
and concatenate. The merge/RLE step (`density_utils.R::rle_merge_granges`) is per-seqname
and lossless, so re-concatenation is byte-identical.
**Expected impact:** Medium-high (38 m at 1/50 cores → potentially a few minutes).

---

### FINDING 7 — `make_unified_annotation.R`: batches balanced by base-pairs, not feature count
**Rule:** `make_unified_annotation` (6 m 39 s, **~1 core** despite `mclapply`)
**File:** `scripts/make_unified_annotation.R`, lines 849–870

It *does* batch and `mclapply` (good intent), but batches are sized by **genome bp**
(`ceiling(genome_bp / threads)`, line 851). Repeat **features cluster on the largest
chromosomes** — the code even says so (line 860: *"features cluster on largest"*). So
the batch holding the big chromosomes carries most of the features and dominates the
`mclapply`, while the other workers finish early → ~1 effective core.

**Output-preserving fix:** assign **whole scaffolds** to batches by a quick
**feature-count** estimate (count features/seqname from the input GFFs first, then
bin-pack), instead of by bp. Whole-scaffold assignment means tier resolution is
unchanged (no feature is split across a batch boundary) → identical output.
**Expected impact:** Medium (better balance on the 6.6 m; larger relative win on
feature-dense genomes).

---

### FINDING 8 — Single-threaded report/plot scripts (`make_repeat_report.R`, `make_summary_plots.R`)
**Rules:** `make_repeat_report` (24 m), `make_summary_plots` (10 m), both **~1 core**
**Files:** `scripts/make_repeat_report.R`, `scripts/make_summary_plots.R`

Neither uses any parallelism (no `mclapply`/`parallel`). The report's cost is many
**independent** file reads/bins: per-family/per-class BigWig `import()` and `bin_bw`
loops (e.g. `make_repeat_report.R` lines 487–545, `lapply(paths, bin_bw …)` line 493),
per-cluster GFF imports (lines 541–545). These are embarrassingly parallel and feed an
HTML assembled at the end.

**Output-preserving fix:** wrap the independent file-read/bin loops in `mclapply`
(threaded via the rule's `{threads}`); assemble the report from the collected pieces
exactly as now. Same bytes, fewer minutes.
**Expected impact:** Medium (24 m + 10 m of single-core work that overlaps poorly today).

---

### FINDING 9 — `calculate_bigwig_density` runs four `calculate_density.R` processes back-to-back
**Rule:** `calculate_bigwig_density` (2 m 30 s, **0.4 cores**)
**File:** `Snakefile` lines 1603–1606 → `scripts/calculate_density.R`

The rule launches `calculate_density.R` **4× sequentially** (unified@10k, unified@100k,
tc@10k, tc@100k). Each call re-spawns R, re-imports the GFF, and recomputes
`coverage()`. The two unified calls re-import the **same** GFF and recompute the same
coverage twice; likewise tc.

**Output-preserving fix:** run the 4 independent calls concurrently (background +
`wait`, or one R invocation looping over the 4 outputs and reusing one `coverage()` per
input). Identical BigWigs.
**Expected impact:** Low-medium (small rule, but trivially parallel).

---

### FINDING 10 — `reduce_dimer_library.py` processes TRC groups sequentially
**Rule:** `tidecluster_reannotate` (the dimer-reduction portion)
**File:** `scripts/reduce_dimer_library.py`, line 218 (`for trc_id in sorted(groups)`)

Each TRC group runs its own `mmseqs easy-search` (line 124+), passing `--threads`, but
groups run **one at a time**. Most groups are tiny (a handful of dimers) so mmseqs
barely uses the threads — per-group startup dominates, serialized over all groups.

> ⚠️ Caveat: `tidecluster_reannotate` is 16.8 h, but that wall-time is dominated by the
> **external** `tc_reannotate.py` RepeatMasker step (mean_load ~1.1 = the same
> "RepeatMasker not parallelized" problem as Finding 1, but inside TideCluster's code,
> not ours). `reduce_dimer_library.py` is a small fraction. Parallelizing across TRC
> groups (with a divided thread budget) is a clean, output-neutral win for the dimer
> step, but **won't move the 16.8 h** — that needs the chunked-RM approach applied
> upstream in `tc_reannotate`, which is out of repo scope.

**Output-preserving fix:** parallelize the per-group loop with a process pool, dividing
threads per worker. Per-group reps are independent → identical output.
**Expected impact:** Low for the rule total (overshadowed by external RM); nice for the
reduction step itself and for satellite-rich genomes.

---

### FINDING 11 — `calculate_statistics_and_make_groups.R` loads the whole genome just to sum lengths
**Rule:** `make_summary_statistics_and_split_by_class` (2 m 22 s, 8.8 GB RSS, ~1 core)
**File:** `scripts/calculate_statistics_and_make_groups.R`

`readDNAStringSet(opt$genome)` (line 28) pulls the entire 3.9 Gb FASTA into memory only
to compute `genome_size <- sum(width(g))` (line 32) — that's the 8.8 GB RSS. The
`.fai` already has the lengths. The per-group `export()` loop (line 72) is serial but
the groups are independent.

**Output-preserving fix:** read lengths from the `.fai` (sum column 2) instead of the
FASTA; optionally `mclapply` the per-group export. Identical CSV/GFF3.
**Expected impact:** Low on time (2 m), good on memory (8.8 GB → tens of MB).

---

## 3. Cross-cutting issue — `threads: workflow.cores` on single-threaded rules blocks the scheduler

`run_pipeline.py` runs `snakemake --cores {threads}` (line 164). **19 rules declare
`threads: workflow.cores`** — i.e. each grabs the *entire* core pool. For the heavy
external tools (dante, dante_ltr, repeatmasker, tidecluster_*) that is correct; they
are internally parallel. But several rules that claim **all** cores are effectively
**single-threaded**, so they hold 50 cores while using ~1 and **prevent any other rule
from running concurrently**:

| Rule claiming all cores | Actual mean_load | Effect |
|---|---|---|
| `make_bigwig_density` | ~1 / 50 | blocks the other density/per-family bigwig rules for 38 m |
| `make_unified_annotation` | ~0.9 / 50 | blocks for 6.6 m |
| `make_unified_tandem_per_family_bigwig` | low | independent of the above but can't overlap |
| `make_tidecluster_tandem_per_family_bigwig` | low | same |
| `merge_rm_and_dante` | 4.4 / 50 | (also Finding 2) |
| `make_tir_combined_library`, `build_fallback_tir_library` | ~0 | fast, but still serialize |

The four density/bigwig rules and the two report rules are mutually independent and sit
near the end of the DAG. With **honest `threads:` declarations** (small/realistic
numbers for the single-threaded ones), snakemake could overlap them, collapsing a long
serial tail toward the duration of the longest single rule.

> Note: rules with **no** `threads:` directive (e.g. `make_repeat_report`,
> `make_summary_plots`, `calculate_bigwig_density`, `make_summary_statistics…`) already
> default to `threads=1`, so they *can* be co-scheduled today — the blocker is the
> all-cores rules running between them.

**Safety:** `threads:` affects scheduling and the `{threads}` value passed to scripts
only — never the result. The only care needed: don't starve a rule that genuinely uses
its threads internally (lower a rule's threads *and* improve its internal parallelism
together, per Findings 6–9).

---

## 4. Prioritized roadmap

**Tier A — biggest, safest wins**
1. **Finding 1** — pack small scaffolds in `repeatmasker_wrapper.py` (+ adaptive chunk
   size). Dominant rule, ~3300 redundant library preps removed, core utilization up.
2. **Finding 2** — fix the `detectCores()` override + double-`mclapply` in
   `merge_repeat_annotations.R`. Kills the 296 GB RSS / OOM risk.
3. **Finding 3** — same fork pattern in `clean_rm_output.R` (and parser swap).

**Tier B — strong wins on already-parallel-but-imbalanced code**
4. **Finding 4** — parallel big-CAP3 + parallel Phase-3 BLAST in `reduce_library_size.py`.
5. **Finding 6** — within-file parallelism in `calculate_density_batch.R`.
6. **Finding 5** — per-subtype parallelism + parasail GIL check in `dante_tir_fallback.py`.
7. **Cross-cutting §3** — honest `threads:` on single-threaded rules.

**Tier C — smaller / scalability**
8. Finding 7 (feature-count batching, `make_unified_annotation.R`)
9. Finding 8 (parallel report/plot file reads)
10. Findings 9, 10, 11 (4× density calls; per-TRC dimer parallelism; `.fai` instead of
    genome load)

---

## 5. How to prove "same output" for each change

Because outputs must not change, pair every change with a parity check on a small but
non-trivial genome (the CARP 781 Mb run, or `tests/fixtures/output_small`):

- **Library reductions** (`reduce_library_size.py`, `reduce_dimer_library.py`,
  `containment_reduce_library.py`): `tests/test_reduce_library_parity.sh` and a
  `diff`/`md5sum` of the output FASTA vs the pre-change run.
- **RepeatMasker wrapper** (Finding 1): compare the final `RM_on_combined_library.out`
  / `.gff3` and **total masked bp** before/after (the repo already uses a ±0.15%
  masked-bp bar for losslessness elsewhere — reuse it). Packing/chunking should be
  exactly identical, not just within tolerance.
- **R GFF3 producers** (`merge_repeat_annotations.R`, `clean_rm_output.R`,
  `make_unified_annotation.R`, `calculate_statistics_and_make_groups.R`): sort + `diff`
  the GFF3 / CSV. `mc.cores`, loop fusion, and batch rebalancing must give identical
  records.
- **Density / report** (`calculate_density*.R`, `make_repeat_report.R`,
  `make_summary_plots.R`): `md5sum` the BigWigs; the HTML/PDF are presentation outputs
  (compare the underlying numbers).
- **Classification correctness** is already gated by the `validate_classifications`
  rule and `tests/test_classification.{py,R}` — keep them green.
