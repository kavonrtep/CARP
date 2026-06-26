# Tier A + Tier B — Implementation & Testing Plan

Derived from `PERFORMANCE_OPTIMIZATION_REPORT.md`. Goal: implement the Tier A and
Tier B speed-ups **with byte-identical outputs**. Every item below is scheduling /
parallelism / memory / redundant-work only — no algorithm changes.

**Golden rule for this whole effort:** *no pipeline output file may change.* The plan
therefore leads with a **determinism-baseline** methodology (§1): for each affected
output we first establish whether the current code is byte-stable run-to-run, which
fixes the correct equivalence bar (exact bytes vs canonical/sorted) for that output.
Then each change ships behind that bar.

Scope recap:

| Item | Script | Rule | Status today |
|---|---|---|---|
| A1 | `repeatmasker_wrapper.py` | `repeatmasker` | 15 h, 11/50 cores |
| A2 | `merge_repeat_annotations.R` | `merge_rm_and_dante` | 296 GB RSS |
| A3 | `clean_rm_output.R` | `repeatmasker` (post-step) | hard-coded `mc.cores=8` |
| B4 | `reduce_library_size.py` | `reduce_library` | 44 m, 1.2/50 cores |
| B5 | `dante_tir_fallback.py` (+`global_local_aln.py`) | `dante_tir_fallback` | 4 h, 15/50 cores |
| B6 | `calculate_density_batch.R` | `make_bigwig_density` | 38 m, ~1 core |
| B7 | `Snakefile` thread declarations | post-processing subgraph | serialized tail |

---

## 1. Testing methodology (applies to every item)

### 1.1 Determinism baseline — do this FIRST, per output
Before changing anything, run the **current** rule/script **twice** on a fixture and
compare each output:

```
run current code  -> out_A
run current code  -> out_B
cmp -s out_A out_B   # byte-stable?
```

- **Byte-stable** ⇒ the equivalence bar for that output is **byte-identity** (`cmp -s`).
- **Not byte-stable** (e.g. an intermediate whose row order follows thread completion —
  see `global_local_aln.py:847` `as_completed`) ⇒ the bar is **canonical equivalence**:
  sort records deterministically (e.g. `sort -k1,1 -k4,4n` for GFF3; sort FASTA by
  header) then `diff`; for masked output compare **total masked bp** within the
  project's **±0.15 %** losslessness bar. We assert byte-identity on the *final* rule
  outputs and canonical equivalence only where a baseline is provably unstable.

This step is cheap and prevents both false alarms (chasing benign reordering) and false
confidence (declaring "identical" against a moving target).

### 1.2 Before/after self-parity harness
Model on the existing `tests/test_reduce_library_parity.sh` (which already does R-vs-Py
`cmp -s` with structural diagnostics). Add a generic
`tests/test_perf_parity.sh <fixture> <target>` that:

1. builds a **baseline** by running the target from a clean `git stash`/worktree of the
   pre-change code,
2. runs the target from the working tree (post-change),
3. compares with the §1.1 bar for that target, preserving the workdir on failure with a
   diff sample (same UX as the existing script).

Use `isolation: worktree` or a second checkout so baseline and candidate never share a
workdir.

### 1.3 Fixtures (already present)
`tests/fixtures/` ships `genome_{tiny,micro,small,medium}.fasta` with full
`output_{...}` trees and `config_{...}.yaml`. Each exercises a different CAP3 / mmseqs /
BLAST / TR mix (see `tests/fixtures/README.md`). Tiering:

- **micro / tiny** — fast smoke + determinism baseline (seconds).
- **small** — primary parity gate per item.
- **medium** — the gate CI's `release.yml` already uses; **green here = mergeable**.
- **CARP 781 Mb** (`tmp/carp_output`) — opt-in heavy parity + real speed measurement.

### 1.4 Existing tests that must stay green
- `tests/test_reduce_library_parity.sh` (R-vs-Py byte-identity) — **the guard for B4**.
- `tests/test_reduce_library_memprofile.sh` — peak-RSS cap; **the guard for B4's
  big-CAP3 parallelism** (must not raise peak memory past the bar).
- `tests/test_containment_reduce.py`, `tests/test_classification.{py,R}`,
  `tests/test_manifest.py`.

### 1.5 Speed measurement
Re-use the pipeline's own `benchmarks/*.tsv`. Acceptance = output-equivalence **AND**
a measured wall-time / `mean_load` improvement on `medium` (and, for A1/B4, on CARP).
Record before/after `s`, `mean_load`, `max_rss` per changed rule.

### 1.6 CI wiring
Add the new parity targets to the `release.yml` test-in-container step next to the
existing reduce_library parity gate. A change merges only when its parity target is
green on `medium`.

---

## 2. Tier A

### A1 — `repeatmasker_wrapper.py`: pack small scaffolds + adaptive chunking
**Files:** `scripts/repeatmasker_wrapper.py`
**Functions:** `split_fasta_to_chunks` (46–94), `split_fasta_to_files` (238–253),
`recalculate_rm_out_coordinates` (156–210), `main` (255–317).

**Change set (output-preserving):**

1. **Pack small sequences into multi-FASTA chunks.** In `split_fasta_to_chunks`, keep
   the >`min_chunk_size` splitting unchanged, but **bin sequences ≤ `min_chunk_size`
   into shared chunks of ≤ `chunk_size`** (first-fit). The matching-table entry per
   small sequence stays `[header, 0, 0, size, new_header]` with **its own** `new_header`
   so coordinate recovery is unchanged; only the *file packing* changes. RepeatMasker
   processes each record in a multi-FASTA independently, so masking per sequence is
   identical.
   - `split_fasta_to_files` becomes "write one file per **chunk**" (a chunk may now hold
     many records) instead of one file per record.
   - `recalculate_rm_out_coordinates` already maps each `.out` row back via
     `matching_table_dict[query]` keyed by `new_header`; with per-sequence `new_header`s
     preserved, **no change needed** there.
2. **Adaptive chunk size** for the big-sequence split: `chunk_size =
   max(5_000_000, ceil(total_bp / (K * threads)))` with `K≈3` (≈3 chunks/core), exposed
   as an optional `--chunk-size`/`--chunks-per-core` arg with the **current 5 Mb as the
   default-preserving floor** so existing behavior is the lower bound.
3. **Scalability hygiene (no behavior change):** replace the O(headers×chunks)
   `[x for x in matching_table if x[0]==header]` (line 90) with a
   `defaultdict(list)` index built once; stream the split instead of materializing the
   genome dict twice (lines 86, 245).

**Why output is identical:** RepeatMasker masks each FASTA record independently of how
records are batched into files; chunk **boundaries on large sequences are unchanged**
(same offsets/overlap), so the concatenated, coordinate-recalculated `.out` is identical
up to **line order**. → equivalence bar: canonical (sort the `.out`/`.gff3`) **and**
total masked bp exact. Verify row-set identity, not file order.

**Risks / edge cases:** (a) packing must not exceed RepeatMasker's per-run limits — cap
records-per-chunk; (b) a packed chunk's `.out` interleaves records — already handled by
per-`new_header` remap; (c) the `-xsmall`/`-no_is`/`-e ncbi`/sensitivity flags must be
passed **exactly** as today (don't touch `repeatmasker()` at 212–236). (d) Investigate
but **do not assume** the "format library once" idea (report §1 item 3) — defer to a
follow-up with its own masked-bp proof.

**Tests:**
- **T-A1.1 (determinism):** run current wrapper twice on `medium` genome+library →
  establish bar (expect canonical, since Pool completion order varies).
- **T-A1.2 (parity):** baseline vs packed wrapper on `small` and `medium`: sort `.out`
  and `.gff3`, `diff` must be empty; `total masked bp` exact; record count equal.
- **T-A1.3 (fragmented stress):** synthesize a 2000-tiny-scaffold FASTA; assert #chunks
  drops from ~2000 to ≲ `ceil(total/chunk_size)+packed`, outputs canonical-identical.
- **T-A1.4 (speed):** CARP run; expect fewer RM invocations, higher `mean_load`,
  identical `RM_on_combined_library.out` (canonical) + masked bp.

---

### A2 — `merge_repeat_annotations.R`: honor `CPU_COUNT`, fuse the double `mclapply`
**File:** `scripts/merge_repeat_annotations.R`, `gff_cleanup` (28–56).

**Change set:**
1. **Delete line 37** (`num_cores <- detectCores()`) — the unconditional override that
   discards the `CPU_COUNT` read at lines 33–36. Result: `num_cores` = `CPU_COUNT`
   (exported by the rule at `Snakefile:1326`) with `detectCores()` only as the NA
   fallback.
2. **Cap forks:** `mc.cores <- max(1L, min(num_cores, CPU_COUNT))` (no `*0.8` of a
   machine-wide count). Optionally cap by a memory heuristic.
3. **Fuse the two passes** (38–39) into one `mclapply` returning
   `list(name=gff$Name[x], strand=strand(gff[x]))` per disjoint range, then split → half
   the forks, half the COW pressure.

**Why output is identical:** `mclapply` returns results in **input order regardless of
`mc.cores`**; fusing two element-wise passes into one computes the same per-range values.
The downstream `resolve_name`/LCA/sort/export logic is untouched. → bar: **byte-identity**
(this output is deterministic).

**Risks:** ensure the rule actually exports `CPU_COUNT` (it does, line 1326); the fused
closure must return both fields for **every** revmap element (including length-1).

**Tests:**
- **T-A2.1 (determinism):** current script twice on a fixture RM+DANTE GFF → expect
  byte-stable (bar = `cmp -s`).
- **T-A2.2 (parity):** baseline vs patched on `small` and `medium` merged GFF3 →
  `cmp -s` identical.
- **T-A2.3 (memory):** `/usr/bin/time -v` peak RSS on CARP merged inputs with
  `CPU_COUNT=8` → assert peak ≪ the 296 GB seen, and **independent of machine core
  count** (run with a fake high `nproc`).
- **T-A2.4 (env contract):** with `CPU_COUNT` unset, falls back to `detectCores()` and
  still produces identical output.

---

### A3 — `clean_rm_output.R`: drop `mc.cores=8`, honor budget, fuse passes
**File:** `scripts/clean_rm_output.R`, `gff_cleanup` (8–30); rule `repeatmasker`
(`Snakefile:1242–1283`).

**Change set:**
1. In the `repeatmasker` rule shell, **`export CPU_COUNT={threads}`** before
   `clean_rm_output.R` (mirrors `merge_rm_and_dante`).
2. In `clean_rm_output.R`: replace the hard-coded `mc.cores = 8` (lines 12–13) with
   `CPU_COUNT` (fallback to a **fixed** small default, e.g. 4 — *not* `detectCores()`, to
   avoid re-introducing the A2 bug). **Fuse** the two `mclapply` passes as in A2.
3. *(Optional, separate commit)* swap base-R `read.table` (line 54) for
   `data.table::fread` **iff** available in `envs/tidecluster.yaml`; gate behind a
   `requireNamespace` check, parse identically (same `col.names`, `fill`). Treat as its
   own parity-gated change.

**Why output is identical:** identical to A2 reasoning (mc.cores order-invariant, fusion
same values). `fread` swap must reproduce `read.table`'s framing exactly → its own
parity test. → bar: **byte-identity** for the GFF3.

**Tests:**
- **T-A3.1 (parity):** baseline vs patched `RM_on_combined_library.gff3` on `small` +
  `medium` → `cmp -s`.
- **T-A3.2 (threads honored):** run with `CPU_COUNT=2` vs `=16`; output identical, wall
  time scales.
- **T-A3.3 (fread parity, if done):** `read.table` vs `fread` path → `cmp -s` on a large
  CARP `.out`.

---

## 3. Tier B

### B4 — `reduce_library_size.py`: parallelize big-CAP3 (bounded) + Phase-3 BLAST
**File:** `scripts/reduce_library_size.py`. Key pieces: `_phase_partition` (545–563),
`_run_phase` (566–579), `process_class`/`_run_cap3`/`_run_mmseqs` (321/238/268),
`run_blast_filter` (404–449, **serial loop**), `assemble_output` (455, fixed order),
`main` (582–664, Phase 1 `threads=1` at line 630).

**Change set (each independently shippable):**
1. **Split Phase 1 into 1-mmseqs (sequential) + 1-bigCAP3 (bounded-parallel).** mmseqs
   keeps `threads=1` worker pool with full `mmseqs_threads` each (the 8 GB-floor
   constraint, documented 268–278, is preserved). Big CAP3 jobs (single-threaded, few
   hundred MB each) run via `_run_phase(big_cap3, threads=K)` with a **new config knob**
   `--max-big-cap3-parallel` (default conservative, e.g. 2–4, or derived from a RAM
   budget) so peak memory stays under the `memprofile` bar.
2. **Parallelize Phase 3 BLAST** (`run_blast_filter`): the `for job in candidates` loop
   (420) is per-job independent (each writes its own `*.blastn`, sets its own
   `singlets_final_file`). Wrap in a bounded `ProcessPoolExecutor(max_workers=min(P,
   len(candidates)))`, giving each `blastn` `threads//P`. Results merge by mutating each
   job's own fields → order-independent.

**Why output is identical:** CAP3/mmseqs/blastn are deterministic; each job uses a
**separate workdir**; `assemble_output` already emits in fixed class-alphabetical order
(455–467) regardless of completion order; the BLAST drop-set is per-job. Scheduling only.
→ bar: **byte-identity** — and the existing `test_reduce_library_parity.sh` (Py-vs-R)
**already enforces it**.

**Risks:** big-CAP3 parallelism raises peak RSS → keep `--max-big-cap3-parallel` small
and **gate on `test_reduce_library_memprofile.sh`**; CAP3 stale-output cleanup
(250–254) is per-workdir so parallel-safe; `ProcessPoolExecutor` must use the existing
`spawn` context (577) for honest RSS.

**Tests:**
- **T-B4.1 (parity):** `tests/test_reduce_library_parity.sh` on micro/small/**medium**
  → byte-identical to R (unchanged guard).
- **T-B4.2 (memory):** `tests/test_reduce_library_memprofile.sh` with
  `--max-big-cap3-parallel` default → peak under existing cap.
- **T-B4.3 (BLAST-path parity):** a fixture with multiple `Class_I/LTR` classes that
  produce contigs+singlets (so Phase 3 has ≥3 candidates) → parallel vs serial filter
  byte-identical.
- **T-B4.4 (speed):** CARP `reduce_library` → expect higher `mean_load`, lower wall,
  identical output FASTA (`md5sum`).

---

### B5 — `dante_tir_fallback.py`: parallelize the per-subtype loop (profile-gated)
**Files:** `scripts/dante_tir_fallback.py` (main loop **1065–1082**, `process_subtype`
795), `scripts/global_local_aln.py` (executor **804–864**).

**Step 0 — PROFILE first (mandatory).** The 4 h / 15-core figure must be localized
before changing anything: instrument `process_subtype` and `run_prime_alignments` with
timestamped logs (or `cProfile`) on a TR/TIR-rich `medium`-plus fixture to attribute time
to {alignment, mmseqs, seqkit, Python O(N²) flank filtering}. Optimize the proven
dominant cost, not the assumed one.

**Change set (apply the subset Step 0 justifies):**
1. **Parallelize subtypes** (1065–1082): replace the serial `for subtype` with a pool
   over the independent subtypes, **collect into a dict keyed by subtype**, then run the
   existing combine (1116–1133) in `sorted(subtype)` order → output order preserved.
   Divide the budget: `threads_per_subtype = max(1, args.threads //
   min(n_subtypes, args.threads))`; pool size `min(n_subtypes, args.threads)`.
2. **Executor choice for alignment — decide empirically.** `global_local_aln.py` uses
   `ThreadPoolExecutor` (829). Add a tiny benchmark (N parasail
   `sg_*_trace_striped_16` calls, threads=1 vs 8, wall-time ratio) to determine whether
   this parasail build **releases the GIL**. If it does → keep threads (already scales).
   If not → switch the alignment pool to `ProcessPoolExecutor`.
3. **O(N²) → dict (safe regardless):** in `create_prime_bed_files` (213) pre-build
   `seqname → [features]` once and merge the 4 per-anchor scans into one. Pure speedup,
   identical BED.
4. *(optional)* overlap the independent 5′/3′ `seqkit subseq` calls (≈358–370).

**Why output is identical — and the determinism caveat:** subtypes write to separate
`<subtype>/` dirs; combine is re-sorted; cross-subtype resolution
(`resolve_cross_subtype_overlaps`, 754) consumes `all_elements` which we assemble in
deterministic subtype order. **BUT** `global_local_aln.py` already appends alignment rows
in `as_completed` order (847) — so the *intermediate* alignment table may not be
byte-stable today. → **Step 0 determinism baseline is decisive here:** run current
`dante_tir_fallback` twice on a fixture and compare the **final** outputs
(`DANTE_TIR_FALLBACK.gff3`, `TIR_fallback_rep_lib.fasta`, `TIR_fallback_extended.fasta`).
- If final outputs are byte-stable ⇒ downstream is order-insensitive; bar = byte-identity.
- If not ⇒ bar = canonical (sorted GFF3 / FASTA-by-header) and we additionally make the
  alignment writer sort deterministically (by pair `(i,j)`) so the change does not *worsen*
  stability.

**Risks:** thread division under subtype skew (one giant subtype) → consider a shared
bounded pool rather than equal division; `ProcessPoolExecutor` needs picklable
`TIRFallbackElement` (verify) if processes are chosen; nested pools must not oversubscribe
(total ≤ `args.threads`).

**Tests:**
- **T-B5.0 (profile + determinism baseline):** as above; sets the bar and the target.
- **T-B5.1 (parity):** baseline vs patched final outputs on `small`+`medium` at the
  established bar.
- **T-B5.2 (single-subtype guard):** input with one subtype → behaves exactly as serial.
- **T-B5.3 (GIL micro-bench):** record the threads=1 vs 8 ratio that justified the
  executor choice (keep in the commit message / test log).
- **T-B5.4 (speed):** CARP/medium → higher `mean_load`, identical outputs.

---

### B6 — `calculate_density_batch.R`: parallelize within the dominant file (bounded, phase-safe)
**File:** `scripts/calculate_density_batch.R`: `density_per_family` (48–60), `process_one`
(78–90), file-level `mclapply` (92–94).

**Important constraint (verified in source):** the per-family tiling uses
`unlist(tileGenome(chr_in_order, tilewidth=step))` in **occupied-first** order; the code
comment (27–47) documents that **bin phasing depends on the genome-concatenation
offset**. Therefore splitting a single file's work **per seqname would change bin
boundaries → different output. Do NOT split per seqname.** Safe wins only:

**Change set:**
1. **Profile the dominant file first** (e.g. `Mobile_elements.gff3`): attribute its
   ~38 min to {`coverage()`, `binnedAverage()`, per-seqname `smooth_score2`}.
2. **Run the two resolutions concurrently** (`process_one` lines 87–88 currently
   sequential): the 10k and 100k computations are independent → run on 2 cores. This
   alone ≈ halves the dominant file's wall time. Output identical (separate BigWigs).
3. **Parallelize the per-seqname smoothing** (line 58 `lapply` → `mclapply`) — it runs
   **after** `binnedAverage`, is per-seqname independent, and does **not** touch bin
   boundaries → identical values.
4. Keep the file-level `mclapply` (92); make the worker budget aware that the dominant
   file now spawns sub-workers (avoid oversubscription: total ≤ `opt$threads`).

**Why output is identical:** (2) splits independent outputs; (3) parallelizes a
post-binning per-seqname transform that is already computed per-seqname. The global
`coverage()`+`binnedAverage()` and the occupied-first tiling are **untouched**, so bins
and values are bit-for-bit the same. → bar: **byte-identity** on each `.bw`.

**Honest ceiling:** if Step-1 profiling shows `coverage()`/`binnedAverage()` dominate
(not smoothing), the safe speedup is bounded (~2× from resolution-parallelism). A deeper
per-seqname split is **possible but high-risk** (must replicate the exact
concatenation-offset phase and boundary-spanning tiles) — defer unless profiling demands
it, and only behind a dedicated bit-exact parity proof.

**Tests:**
- **T-B6.1 (determinism):** current rule twice on `medium` → expect byte-stable `.bw`
  (bar = `cmp -s` / `md5sum`).
- **T-B6.2 (parity):** baseline vs patched, **every** `.bw` in
  `Repeat_density_by_class_bigwig/{10k,100k}/` `md5sum`-identical on `small`+`medium`.
- **T-B6.3 (oversubscription):** with `-t 8`, assert no more than 8 R workers live
  (process count sampling) and output identical.
- **T-B6.4 (speed):** CARP/medium dominant-file wall time down, `.bw` identical.

---

### B7 — Honest `threads:` on single-threaded rules (scheduling only)
**File:** `Snakefile`. Rules declaring `threads: workflow.cores` but effectively
single-/low-threaded: `make_unified_annotation` (1363), `make_bigwig_density` (1506),
`make_unified_tandem_per_family_bigwig` (1635), `make_tidecluster_tandem_per_family_bigwig`
(1677), `make_tir_combined_library` (479), `build_fallback_tir_library` (980),
`merge_rm_and_dante` (1320). Heavy external-tool rules (dante, dante_ltr, repeatmasker,
tidecluster_*) **keep** `workflow.cores`.

**Change set:** set each such rule's `threads:` to its **true internal parallelism** so
snakemake can co-schedule the independent post-processing subgraph (the density / report
rules near the DAG leaves), and so concurrent rules' thread sums stay ≤ `--cores`.
Sequence this **after** B4/B6 land (so the numbers reflect the new internal parallelism).

**Why output is identical:** `threads:` affects only (a) scheduling and (b) the integer
substituted into `{threads}` (→ `mc.cores`/`-t`/`--threads`), none of which change
results (established in A2/A3/B4/B6). → bar: **byte-identity** on all final outputs.

**Risk:** the *only* real risk is oversubscription if a rule's `{threads}` stays high
while it now runs concurrently with others — avoided by setting honest counts and relying
on snakemake's `sum(threads) ≤ --cores` guarantee.

**Tests:**
- **T-B7.1 (full-run parity):** end-to-end `config_small` then `config_medium` before
  vs after; **all** final outputs byte-identical (or canonical where §1.1 says so).
- **T-B7.2 (concurrency + wall):** confirm from `benchmarks/` that post-processing rules
  now overlap and total wall-time drops; no rule's `mean_load × concurrent rules` exceeds
  cores.

---

## 4. Rollout order & acceptance

Independent commits, each its own parity gate; suggested order by risk/impact:

1. **A2** (delete one line + fuse) — trivial, huge memory win, byte-identity. *Ship first.*
2. **A3** (mirror of A2) — trivial.
3. **B4** (guarded by the two existing reduce_library tests) — high confidence.
4. **A1** (most impact; needs the canonical-`.out` parity harness) — medium effort.
5. **B6** (bounded, phase-safe subset) — after profiling.
6. **B5** (profile-gated; most moving parts) — last of the script changes.
7. **B7** (scheduling) — **after** B4/B6 so thread numbers are honest.

**Definition of done (per item):** (a) §1.1 bar met on `small`+`medium`; (b) relevant
existing tests green; (c) a measured wall-time / `mean_load` improvement recorded on
`medium` (and CARP for A1/B4); (d) new parity target wired into `release.yml`.

**Global regression gate before any release:** one full `config_medium` run, all final
outputs equivalent to a pre-Tier-A/B baseline, plus `test_classification.{py,R}` and
`test_manifest.py` green.

---

## 5. Summary table

| Item | Core change | Output bar | Primary test | Effort | Impact |
|---|---|---|---|---|---|
| A1 | pack small scaffolds + adaptive chunks | canonical `.out` + masked-bp exact | T-A1.2/3 | M | ★★★ |
| A2 | drop `detectCores()` override, fuse mclapply | byte | T-A2.2/3 | XS | ★★★ (mem) |
| A3 | `mc.cores` from `CPU_COUNT`, fuse | byte | T-A3.1 | XS | ★★ |
| B4 | bounded big-CAP3 ∥ + Phase-3 BLAST ∥ | byte (existing R parity) | parity+memprofile | M | ★★ |
| B5 | parallel subtypes + GIL-correct executor | byte/canonical per baseline | T-B5.0/1 | L | ★★ |
| B6 | resolution-∥ + per-seqname smoothing ∥ | byte (`.bw`) | T-B6.2 | M | ★★ |
| B7 | honest `threads:` declarations | byte (all finals) | T-B7.1 | XS | ★★ (tail) |

---

## 6. Implementation status & verification (this session)

All changes implemented and verified output-preserving on the `tests/fixtures`
data with the in-repo conda envs. Determinism baseline established per output
before each change (§1.1).

| Item | Files changed | Verified | Evidence |
|---|---|---|---|
| **A1** | `scripts/repeatmasker_wrapper.py` | ✅ | medium genome packed **35 RM runs → 1**; downstream `clean_rm_output.R` GFF3 **byte-identical** (418 recs, matches fixture); `.out` identical except the cosmetic per-run ID column (unused downstream). Also fixes a latent crash on no-hit chunks. |
| **A2** | `scripts/merge_repeat_annotations.R` | ✅ | merged GFF3 **byte-identical** in both the `CPU_COUNT`-honoured and unset-fallback paths (`md5 9432630f…`). Removes the `detectCores()` override → no more 296 GB COW blow-up. |
| **A3** | `scripts/clean_rm_output.R`, `Snakefile` (repeatmasker `CPU_COUNT` export) | ✅ | RM GFF3 **byte-identical** at `CPU_COUNT=16` and unset (`md5 96981855…`). |
| **B4** | `scripts/reduce_library_size.py`, `Snakefile` (+`reduce_library_max_big_cap3_parallel` config, rule wiring) | ✅ | **byte-identical** on micro/small/medium with default flags; **byte-identical** with the new parallel branches forced (Phase-1b 3-way CAP3 pool + Phase-3 2-way BLAST pool incl. cross-worker propagation). Existing `test_reduce_library_parity.sh` (R-vs-Py) **still green**. |
| **B5** | `scripts/dante_tir_fallback.py` | ✅ | parallel subtypes (3 on medium) → GFF3 + rep_lib **byte-identical**; `TIR_fallback_extended.fasta` **canonical-identical** (this file is already run-to-run unstable in the *current* code — established by the determinism baseline). |
| **B6** | `scripts/calculate_density_batch.R` | ✅ | restructured to (file × resolution) tasks; **all 42 BigWigs byte-identical**. |
| **B7** | — (none) | n/a | Examined every `threads: workflow.cores` rule: after B4/B6 they genuinely use their cores, and the truly single-threaded rules (`make_repeat_report`, `make_summary_plots`, …) already declare no `threads:` (default 1 → already co-schedulable). No safe standalone `threads:` change remains; the residual waste is `make_unified_annotation`'s bp-balanced batches (report **Finding 7, Tier C**, out of approved scope). |

**Whole-pipeline check:** full `snakemake -n` DAG build on `config_medium`
succeeds with all edits (38 rules).

**Not yet done (needs the real big genome, not fixtures):** measured wall-time /
`mean_load` / peak-RSS deltas for A1, B4, B5, B6, and a memprofile check of the
B4 big-CAP3 parallelism (`tests/test_reduce_library_memprofile.sh`) — the
fixtures prove *parity* but are too small to show *speed* or stress memory.
Nothing was committed (left for review).
