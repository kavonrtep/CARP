# End-to-end validation of the new LINE boundary bounds

## Commands

```bash
cp -r <this dir> /nfsroot/projects/darwin/runs/tmp3/bounds_e2e
cd /nfsroot/projects/darwin/runs/tmp3/bounds_e2e
./RUN.sh --check     # seconds; verifies image, patch, genomes, baselines. RUNS NOTHING.
./RUN.sh             # both genomes, then the comparison
```

Also `./RUN.sh --compare` (re-print without running) and `./RUN.sh <tag>` (one
genome). Resumable: finished genomes are skipped via `results/DONE_<tag>`; a
failed one prints its log tail and clears its marker so a re-run retries it.

Knobs: `THREADS` (default 32), `SIF`, `SCRATCH`.

Expect roughly **3 hours** if the two run together — the baselines took 2.8 h
each at 0.43 and 0.52 Gb.

## What is being tested

1.7.0rc1 caps a LINE's flank extension at **2000 bp (5')** and **800 bp (3')**.
Measured on 218 LINE loci across 7 genomes where BOTH ends were observed (poly-A
tail plus a TSD clearing a measured chance floor), those caps truncate real
element. The 3' extension turned out not to be a stable distance at all — it
varies **13.8x** (p90/p10), because it is only what is left over after the domain
annotation stops. Pearson r between core length and 3' extension is **−0.936**:
the two are nearly a constant sum. What *is* stable is the span from the START of
the domain core to the element's 3' end, which varies **1.16x**.

So `main` replaces the flat 3' cap with a **4500 bp bound on that span** and
raises the 5' cap to **3400 bp**.

**This is the only behaviour change on main that has never been through a full
pipeline run.**

## The two genomes

| tag | Gb | baseline | why this one |
|---|---:|---|---|
| GCA_973357735.1 | 0.52 | run-000232 | **sensitive**: 494 of 550 LINE elements hit a cap under rc1, so if the change does anything it shows here |
| GCA_937616625.2 | 0.43 | run-000230 | **control**: only 8 of 633 capped, so it should barely move |

Same design as the earlier control batch: one case where the change should bite,
one where it should not. A control is only a control if it can fail.

## How the container is patched

The shipped image is 1.7.0rc1 and does not carry the change. Rather than rebuild
it, `patch/scripts/` and `patch/classification_vocabulary.yaml` are bind-mounted
over `/opt/pipeline/`. **The Snakefile is byte-identical between the tag and
main**, so the pipeline structure is untouched — only the vocabulary and three
scripts differ.

The run goes through the container runscript (`run_pipeline.py`), not bare
snakemake, because provenance-gated code is skipped otherwise — that gap shipped
the 1.1.0 `make_unified` crash.

Settings mirror the baseline runs exactly (`rush`, both culling limits 2,
`reduce_library`, `cleanup_intermediates: maximal`), so the bounds are the only
difference.

## What to expect, and what would be alarming

The change lets LINE elements grow, so `Class_I/LINE` should rise **slightly**.
A local A/B on one small genome moved **4 of 343 elements (+0.25% extension)`**,
so single-digit-percent moves are the expectation.

| observation | reading |
|---|---|
| small LINE rise, total annotated roughly flat | as designed |
| **large LINE rise, or LTR/Tandem falling to pay for it** | **investigate** — that is the displacement pattern the whole boundary investigation started from (wheat read 15.93% LINE against a real ~1.4–2.5%) |
| `n_capped` falls sharply | expected; the caps bind less often now |
| `n_at_flank_ceiling` rises above ~0 | **red flag** — extensions running to the 10 kb search limit is the runaway signature |

## Known limitation of this test

The bounds were calibrated on confirmed loci — LINEs still carrying a poly-A tail
and an intact TSD, i.e. young insertions. On one genome the extension layer
extended 19.0% of all elements but **none** of the 27 confirmed ones, so that
population may not represent the genome. The bounds are *caps*, and calibrating a
cap on the longest elements is the safe direction (it binds rarely), but this run
measures the genome-wide effect precisely because the calibration set might not.

## Contents

| path | what |
|---|---|
| `RUN.sh` | the runner |
| `patch/scripts/` | full patched CARP scripts (dante_line.py pulls in modules via importlib, so a subset does not work) |
| `patch/classification_vocabulary.yaml` | the new bounds — the change under test |
| `genomes.tsv` | tag, assembly, baseline run, rationale |
| `results/` | per-genome output, logs, DONE markers |

## Note

`--check` must run on the exec host. The agent sandbox only mounts
`/nfsroot/projects/darwin/runs`, not `genomes/`, so the two assembly paths could
not be verified from there — they were taken verbatim from the baseline runs'
`job.sh`.
