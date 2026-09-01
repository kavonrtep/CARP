# Flank-extension boundary benchmark

Scores any boundary estimator for `dante_line.py` / `dante_tir_fallback.py`
against **exact** boundaries, with **no pipeline run** (~1 min per small genome).

Ground truth comes from tools that already know the answer:
- **LTR** — DANTE_LTR elements whose `Rank` contains `L` (both LTRs found).
- **TIR** — DANTE_TIR primary elements carrying `tir_seq5` + `tir_seq3` + `tsd`.

Both have a domain core, so the true core→boundary distance is known. The
harness hands the algorithm only the core edge and asks it to rediscover it.

```bash
PY=<a python with parasail>          # e.g. any .snakemake/conda env with parasail
$PY eval_diagonal.py --mode tir --tir-group class --group-by rexdb \
    --gff <out>/DANTE_TIR/DANTE_TIR_final.gff3 \
    --fasta <out>/genome_cleaned.fasta --out result.tsv --repo <repo>
python3 purity.py "label=result.tsv"      # reliable-core purity + contamination
python3 analyze.py "label=result.tsv"     # per-policy error / recovery tables
python3 sf_bounds.py <DANTE_TIR_final.gff3>...   # per-superfamily true extents
python3 fallback_contam.py "label=<carp_output>" # real-population contamination
```

## Score purity, not recall

The purpose of these layers is a **reliable core**, not a complete element.
Under-extension costs nothing; only overshoot does. `purity.py` is therefore the
right scorer — `analyze.py`'s "recovered" column measures complete-element
recall and will make a correctly-conservative estimator look broken.

## `--group-by rexdb`

Group by the core domain's REXdb protein hit. Core-only, so it is non-circular
with respect to flanks — unlike DANTE_TIR's `Cluster_ID`, which clusters whole
elements and partly assumes the answer. Grouping dominates the result: on
Boechera LTR it moved "recovered" from 12.6 % to 56.3 % with no code change.

## Results that came out of it (2026-09-01)

- Anchored-**diagonal/drift** boundary detection: evaluated and **rejected**.
  Tested on both element types, 3 genomes, 3 groupings — it lost every time.
  Sensitivity also collapses where it is needed (84 % sparse → 48–63 % dense).
- TIR purity 99.93–99.96 %, **zero** sides overshooting >50 bp.
- The TPase→TIR interval is genuinely less conserved than an LTR internal
  region (73.2 % vs 85.3 % within-family identity), so TIR extension repays an
  indel at two-thirds the rate — which is why simple extension underperforms there.
- Per-superfamily p95 extents from `sf_bounds.py` are the source of
  `max_extension_per_side` in `classification_vocabulary.yaml`.

**Caveat:** ground truth is DANTE_LTR/DANTE_TIR *successes*, i.e. exactly the
population where extension works. The fallback runs where DANTE_TIR FAILED, so
these numbers are an upper bound.

Data files (`*.tsv`, `wheat_chr.gff3`) are regenerable and deliberately not committed.
