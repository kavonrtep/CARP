# unified_multibatch — make_unified_annotation.R batch-invariance fixture

A tiny, synthetic set of `make_unified_annotation.R` inputs used by
`tests/test_unified_multibatch_determinism.sh` to prove the unified GFF3 is
**byte-identical across thread counts** (single-batch vs multi-batch).

## Why it exists

`make_unified_annotation.R` splits the genome into per-thread batches whose
*composition* depends on the thread count. Several per-batch operations were
batch-composition-dependent — most severely the TE-derived-satellite (`te_sat`)
pre-pass, which trimmed still-self-overlapping `t1` against the satellite regions
and thereby **disjoined `t1`'s internal overlaps**, but only in batches that
contained a `te_sat` array. On a large genome the same sequence then resolved
differently at `--threads 1` vs `--threads N` (816 / 2.82 M features on the
run116 test genome). Full write-up:
`docs/archive/tier1_resolution_determinism_audit.md`.

The full-pipeline determinism gate (`.github/workflows/pipeline.yml`) cannot catch
this: its `small`/`medium` fixtures are far below the batch-target floor, so
`make_unified` always takes the single-batch path regardless of thread count.

## Design

`make_unified_annotation.R` reads only **sequence lengths** from the `.fai`
(never the FASTA), so the fixture declares three 2 Mb sequences while the GFF3s
stay tiny. At `--batch_size 1000000 --threads 3`, `make_batches` puts each 2 Mb
sequence in its own batch (3 batches); at `--threads 1` it is one batch.

- **seq_te** — a TideCluster clustering array `TRC_1` (100001–140000) tiled by 4
  complete CRM `DANTE_LTR` elements + one DANTE domain per 10 kb monomer window
  (`trc_table` period 10000). This clears the coverage, LCA-depth and
  domain-rhythm tests, so `TRC_1` is tagged `TE_origin` → `te_sat` is non-empty.
- **seq_ovl** — two *overlapping* complete CRM `DANTE_LTR` elements
  (100000–107030 and 105557–115187), far from any `te_sat`. In the single batch
  they share the batch with seq_te's `te_sat` (pre-fix: pre-disjoined → 3 pieces);
  in the multi-batch run they are alone (no `te_sat` → correct greedy 2 pieces).
  This is the exact shape that differed on the real genome.
- **seq_pad** — a lone RepeatMasker hit, padding to a third batch.

The test asserts the fixture actually triggers `te_sat` and a real >1-batch split
(so a future edit that stops exercising the path fails loudly rather than
silently passing), then diffs the two unified GFF3s. It PASSES on the fixed code
and FAILS (9 differing lines) against the pre-fix code (commit 67acf4c).

## Regenerating

`python3 generate.py <outdir>` rewrites every file here. Edit `generate.py`, not
the GFF3s, and keep this README in sync.

## Out of scope

The fixture deliberately avoids isolated same-tier features that overlap *each
other* but no higher tier (e.g. two overlapping DANTE domains far from any
tier-1). Those trip a separate, rarer, pre-existing strand (`*` vs source)
non-determinism in the `trim_to_nonoverlap` / `resolve_within_tier` interaction
that is noted in the audit doc and not addressed here.
