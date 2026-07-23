# CARP large-genome (~90 Gbp) readiness audit

Scope: CARP's **own** code (`scripts/*`, `Snakefile`) — not the upstream DANTE/
TideCluster tools (covered in the `docs/*_request.md` set). Target: a ~90 Gbp
assembly, >10,000 contigs, N50 ~2.13 Gb, a unified annotation of tens of
millions of features, and a de-novo library of up to ~10⁵ consensi.

Audited 2026-07-15 (dante_line/dante_tir_fallback read directly; the R
annotation, density/report, library, and Snakefile subsystems via parallel
sub-audits). **fd-exhaustion was NOT found in any CARP script** — the only fd
risks are in the Snakefile rules that shell out to the upstream tools (already
patched with `ulimit` stopgaps).

---

## 0. TL;DR — two things will stop a 90 Gbp run cold

1. **The scheduler will OOM-kill almost every rule** — the Snakefile declares
   **no `resources: mem_mb`/`runtime` on any rule**, so on SLURM/PBS the heavy
   tools get the 1–4 GB default. This is the single highest-leverage fix.
   (Deferred by maintainer for now.)
2. **Several stages load the whole genome or whole annotation into one process**
   — a hard OOM even with correct scheduling (see §2).

**Not a blocker (resolved by the maintainer):** the 32-bit IRanges/BigWig 2³¹
limit (§4) — assemblies with chromosomes ≥ 2.147 Gb are split into sub-2³¹
chunks upstream before CARP, so no sequence reaches the R/BigWig layer over the
limit. Kept in §4 as a documented invariant to preserve, not an action item.

---

## 1. Fixed in this session

| Fix | File | What |
|-----|------|------|
| **Flank-clipping O(N²)** | `scripts/dante_line.py`, `scripts/dante_tir_fallback.py` | Per-pattern/anchor full rescans of `all_features` + `mask_features` (millions of TideHunter records) replaced by a shared `FeatureIndex` (bisect); O(patterns × features) → O(patterns × log features). **Byte-identical**, proven by `tests/test_flank_index.py` (primitive + both scripts, both strands, ± mask). |
| **Whole-genome OOM** | `scripts/calculate_statistics_and_make_groups.R` | `readDNAStringSet(genome)` (~90 GB, used only to sum widths) → `fasta.seqlengths()`. |
| **fd stopgaps** | `Snakefile` (`tidecluster_long`, `tidecluster_short`, `repeatmasker`) | `ulimit -n "$(ulimit -Hn)"` added (mirrors `dante_ltr`/`tidecluster_reannotate`) — chunked runs + >10k contigs open >1024 handles. |
| **Whole-genome dict ×2 OOM** | `scripts/repeatmasker_wrapper.py` | `split_fasta_to_chunks` (:98) and `split_fasta_to_files` (:272) slurped the whole genome into a dict (~90 GB each) → new `iter_fasta_records()` streams one record at a time (peak = largest single sequence). Also fixes the two leaked handles. **Byte-identical**, proven by `tests/test_repeatmasker_wrapper_streaming.py`. |
| **Unified-annotation hot paths** | `scripts/make_unified_annotation.R` | (A) `subset_seqs` full-Rle→character coercion per batch → `seqnames(gr) %in% seqs`; (B) provenance header no longer `readLines()` the whole multi-GB GFF3 (streaming prepend); (C) tier-1 resolver O(N²) → O(k²) (only overlapping features enter the greedy loop; no per-iteration whole-set append). All **output-identical** — A/B verified in R; C by `tests/test_resolve_tier1_overlaps.R` (60 overlap-dense trials). |
| **RM-layer 296 GB OOM** | `scripts/clean_rm_output.R`, `scripts/merge_repeat_annotations.R` | `gff_cleanup` did `as.list(gff_disjoin$revmap)` (tens of millions of tiny R vectors) inside an `mclapply` that COW-forked the whole GRanges (296 GB max_rss at 3.9 Gb) → single vectorized `extractList` + `unstrsplit` in the compressed S4Vectors domain (no fork). Plus `clean_rm_output.R` reads only the 5 used `.out` columns via `colClasses`. **Output-identical**, proven by `tests/test_gff_cleanup.R` (+ colClasses verified on ragged rows). |
| **Density hang** | `scripts/density_utils.R` | `rle_merge_granges` scanned all tiles once per seqlevel (`d[as.logical(seqnames(d)==sn)]`), O(tiles × seqnames) ≈ 9×10¹¹ → hang. Rewritten to a single pass over the sorted seqnames-Rle runs. **Output-identical**, `tests/test_density_utils.R`. |

---

## 2. OOM blockers still open (must fix before a 90 Gbp run)

Ranked by how surely they crash. All are "load whole genome/annotation into one
process" unless noted.

| # | File:line | Mechanism | Fix |
|---|-----------|-----------|-----|
| B1 | `Snakefile` (all heavy rules) | **No `mem_mb`/`runtime`** → scheduler default OOM-kills `dante`, `dante_ltr`, `dante_tir*`, `tidecluster_*`, `repeatmasker`, `merge_rm_and_dante`, `make_unified_annotation`, the density/stats/report R rules | Add `resources: mem_mb=lambda wc,input,attempt: …, runtime=…` with an `attempt` multiplier; 64–256 GB floors for whole-genome tools. **Needs your scheduler policy** (SLURM/PBS/local) to pick numbers. |
| ~~B2~~ | ~~`repeatmasker_wrapper.py:98, :272`~~ | **FIXED this session** — now streams via `iter_fasta_records()` (byte-identical; test added). | — |
| B3 | `make_unified_annotation.R` (imports ~14 layers; `:1285` combine) | Loads every annotation layer + the multi-GB RM GFF into one R process, then a tens-of-millions-feature combined GRanges | Biggest `mem_mb` floor (B1); longer term, process per-seqname and stream the export. |
| ~~B4~~ | ~~`merge_repeat_annotations.R`; `clean_rm_output.R`~~ | **FIXED this session** — `gff_cleanup` vectorized (`extractList`+`unstrsplit`, no `as.list`/`mclapply`); `tests/test_gff_cleanup.R`. | — |
| B5 (partial) | `clean_rm_output.R:70` | **read.table narrowed** to the 5 used columns via `colClasses` (dependency-free) — a `data.table::fread` would be faster still but is not in the env. `merge_repeat_annotations.R:94` `cat`+`import` of the whole RM+DANTE layer remains (bounded by feature count; needs `mem_mb`, not an algorithm change). | Optional: add `data.table` to the env for `fread`. |
| B6 | `calculate_density_batch.R:58`; `calculate_density.R:117` | `tileGenome` over the **whole 90 Gbp** per task (90M 1 kb bins) + `as.character` filter, run concurrently by `mclapply` and re-done per TRC family | Tile only occupied seqlevels (`tileGenome(chr_size[occ], …)`) — byte-identical, drops the 90M allocation. |
| ~~B7~~ | ~~`make_unified_annotation.R:1035`~~ | **FIXED this session** — streaming prepend (verified identical). | — |

Off-path (fix or delete before anyone wires them in): `calculate_coverage_from_blast.py` (per-base dict-of-dicts, catastrophic — currently orphaned) and `reduce_library_size.R` (dead; `.py` is the live rule).

---

## 3. Quadratic / wall-time hazards (won't OOM, but can hang or add hours)

| # | File:line | Mechanism | Fix |
|---|-----------|-----------|-----|
| ~~Q1~~ | ~~`density_utils.R:28-29`~~ | **FIXED this session** — single pass over the sorted seqnames-Rle runs; `tests/test_density_utils.R`. | — |
| B6 (memory, open) | `calculate_density_batch.R:58`, `calculate_density.R:30` | `tileGenome` still materialises a whole-genome (~90M-tile at 10 kb-window) grid per task. **Occupied-only tiling was investigated and REJECTED**: `tileGenome` derives the tile *width* from `sum(seqlengths)`, so tiling occupied-only changes every bin (verified). With Q1 fixed this is a **memory** cost (a few GB transient/task, not a hang), not an algorithm bug — needs `mem_mb` (B1). | `mem_mb` floor on the density rules; optionally cap that rule's `--threads` so concurrent 90M-tile grids don't sum. |
| ~~Q2~~ | ~~`make_unified_annotation.R:630-642`~~ | **FIXED this session** — greedy loop restricted to the overlapping subset (O(N²)→O(k²)); `tests/test_resolve_tier1_overlaps.R`. | — |
| ~~Q3~~ | ~~`make_unified_annotation.R:130`~~ | **FIXED this session** — `seqnames(gr) %in% seqs` (no full character coercion). Deeper `split()`-once refactor still possible if profiling shows the per-batch subset itself dominates. | — |
| Q4 | `make_summary_plots.R:81-98` | `xx[xx$seqnames==sn,]` per track per contig, no min-length/max-tracks cap → O(rows × 10⁴ contigs); draws 10⁴ separators | `split(xx, seqnames)` once; cap to largest N seqs (mirror `make_repeat_report`'s 50). |
| Q5 | `make_unified_annotation.R:649, :1076, :1110` | `findOverlaps(gr, gr)` self-joins materialize huge Hits only to test emptiness | `isDisjoint()` / `countOverlaps()`; run the L1 overlap report per-seqname. |
> **Scope note (maintainer):** the repeat *library* scales with library size, not
> genome size (~80 MB at 14 Gbp → maybe a few hundred MB at 90 Gbp), so the
> library reducers are **wall-time**, not OOM/hang, concerns. Q6/Q7 below are
> throughput items, not blockers; the CAP3 O(N²) is bounded by per-class library
> size and would only be a *slow* stage on the biggest LTR lineage, not a hard
> failure — so the risky, byte-identity-breaking CAP3 replacement (Q6) is
> deliberately **not** pursued.
>
| Q6 (deferred) | `reduce_library_size.py:256, :573` | CAP3 on an LTR leaf class is ~O(N²), single-threaded. **Byte-locked to the R reference** (`test_reduce_library_parity.sh`), so a scalable fix (pre-cluster then CAP3 within sub-clusters) would change the library — needs a masked-bp acceptance bar, not pursued (bounded by library size). | (if ever needed) mmseqs/linclust pre-cluster + masked-bp validation; or just raise `--max-big-cap3-parallel`. |
| Q7 (open) | `containment_reduce_library.py:64, :151, :85` | Whole-library self all-vs-all blastn with all hits + records in RAM — bounded by library size (hundreds of MB), so time not OOM. Class-sharding is *mostly* output-identical but can find MORE containers where `-max_target_seqs` truncates → needs masked-bp validation. | Shard by `#classification`; stream hits per query. |
| ~~Q8~~ | ~~`reduce_dimer_library.py:218`~~ | **FIXED this session** — `ProcessPoolExecutor` over TRC groups (each mmseqs single-threaded), per-group `rmtree`. Output-identical (mmseqs thread-invariant, verified `-t1==-t4` on a real 18-group fixture; `tests/test_reduce_dimer_parallel.sh`). | — |

---

## 4. The 32-bit IRanges / BigWig ceiling (invariant to preserve, NOT a blocker)

`make_repeat_report.R`, `calculate_density*.R`, `make_summary_plots.R`,
`make_unified_annotation.R` — every GRanges/`coverage`/`tileGenome`/`Views`/
BigWig operation uses 32-bit signed positions (max **2,147,483,647 ≈ 2.147 Gb**).
A sequence ≥ that length cannot be represented and those stages would overflow.

**Resolved upstream (maintainer):** assemblies whose chromosomes exceed 2³¹ bp
are split into sub-2³¹ chunks *before* CARP runs, so no sequence reaching the
R/BigWig layer is over the limit. No code change needed. Recorded here as an
invariant to keep in mind: do not add a stage that reconstitutes full oversized
chromosomes in a GRanges.

---

## 5. Verdict — safe vs. needs-work

- **Safe at 90 Gbp as-is:** `build_fallback_tir_library.py` (default off), the
  rmblast culling shim, `resolve_ltr_tandems.py`, `classification.py/.R`,
  `manifest.py`, `validate_classifications`, the streaming `clean_genome_fasta`/
  `make_track_for_Ns` (slow but low-mem), and — after this session's fixes —
  `dante_line.py` / `dante_tir_fallback.py` (flank step) and
  `calculate_statistics_and_make_groups.R`.
- **Near-safe with one small fix:** `validate_unified_gff3.py` (`readlines()` →
  streaming; relax `ID_RE` to `\d{8,}` — IDs exceed 8 digits past 10⁸ features).
- **Needs real work before 90 Gbp:** everything in §2 (OOM) and §3 (quadratic).
- **Hard ceiling to decide on:** §4 (2^31) — determines whether the track/report
  layer can run at all on this genome.

## 6. Suggested order of attack

1. **`Snakefile` `resources:`** (B1) — unblocks the scheduler; nothing else
   matters if jobs are OOM-killed. Needs your scheduler details.
2. **`repeatmasker_wrapper.py` streaming** (B2) — mandatory path, direct
   genome-size axis.
3. **`density_utils.R` Q1 + `calculate_density*` B6/Q3** — convert hang/OOM into
   finite runs for the track layer.
4. **`make_unified_annotation.R` Q2/Q3/B7 + `merge/clean` B4/B5** — the unified
   annotation is the product; these make it finish.
5. **Library `reduce_*` Q6/Q7/Q8** — throughput.
6. **Decide on §4 (2^31)** — possibly gate the BigWig/report stages by max
   sequence length.
