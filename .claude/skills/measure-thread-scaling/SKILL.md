# Measure all-vs-all alignment thread scaling + peak RSS (PR-C / item 5)

`measure_thread_scaling.py` profiles the shared all-vs-all alignment engine
(`run_all_vs_all_alignment` in `scripts/global_local_aln.py`, used by both the
`dante_line` and `dante_tir_fallback` rules) at several thread counts on a real
input, and reports peak RSS alongside. It answers item 5 of
`docs/alignment_engine_hardening_plan.md` — *does adding threads actually add
throughput here?* — and doubles as the memory (RSS) confirmation for the PR-B
streaming/heap changes (items 1/2/4).

## When to run
On a genome with a **genuinely high-copy** LINE (or TIR) family — the O(N²) step
only matters at scale. The small CI fixtures are too tiny to be informative
(they just confirm the harness works).

## Input
Point `-f` at an actual all-vs-all input FASTA from a real run:
- dante_line: `<output>/DANTE_LINE/ENDO_RT_5prime.fasta`, `ENDO_RT_3prime.fasta`,
  `ENDO_RT_RH_5prime.fasta`, `ENDO_RT_RH_3prime.fasta`
- fallback:   `<output>/DANTE_TIR_FALLBACK/<subtype>/TPase_5prime.fasta`,
  `TPase_3prime.fasta`
Pick the one with the most sequences (that is the expensive family).

## Run (needs the dante_line conda env — parasail; mmseqs only with the flags below)
```bash
# find the env (parasail + mmseqs + seqkit), then:
conda run -p <dante_line_env> python \
  .claude/skills/measure-thread-scaling/measure_thread_scaling.py \
  -f <output>/DANTE_LINE/ENDO_RT_5prime.fasta \
  -t 1,4,8,16,32 --end 5 --reps 2
```
- `--end 5|3` — must match how the file is aligned (5prime files use `--end 5`).
- `--max-group-size N` — also exercise the grouped/streaming path (this is the
  configuration whose **peakRSS** validates PR-B on the high-copy family).
- `--prefilter` — include the mmseqs prefilter (default off, to isolate the
  parasail alignment scaling).
- Each (thread-count, rep) runs in its own subprocess so peak RSS is per-config.

## Reading the output
Columns: `wall_s`, `cpu_self` (in-process CPU), `eff_cores = cpu_self/wall`,
`speedup` + `effic%` vs the smallest thread count, `peakRSS_MB`, `records`.
- **eff_cores ≈ threads and effic% high** → the alignment loop scales; keep/raise
  `-t`. (Expected: parasail's Python binding is ctypes and releases the GIL, so
  unlike the blastn case that motivated this, threads should help.)
- **eff_cores plateaus near 1 (effic% collapses as threads rise)** → the loop is
  GIL/Python-overhead-bound (common when flanks are very short so the tiny C
  alignment is dwarfed by per-pair dict/filter work). The fix is then
  process-level concurrency (several groups at once / a ProcessPoolExecutor),
  NOT more threads per group.
- **peakRSS_MB** should stay roughly flat as the family grows (that is the PR-B
  streaming/heap win). A figure that climbs with family size means something
  still accumulates.
- **determinism** line must say OK — output is byte-identical across all thread
  counts (guards item 3). A WARNING there is a real regression.

## Notes
- Lives under `.claude/` so it is not copied into the release SIF and is not a
  pipeline rule — run it by hand.
- `cpu_self`/`eff_cores` exclude mmseqs subprocess CPU (grouping/prefilter), so
  they measure the alignment loop itself even with `--max-group-size`/`--prefilter`.
