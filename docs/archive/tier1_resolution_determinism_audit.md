# Pre-existing determinism gap: tier-1 structural-TE overlap resolution (large multi-sequence genomes)

**Status:** open — separate scoped effort. Surfaced during the TE-derived-TRC
domain-rhythm work (which is itself deterministic); flagged here so it is not lost.

## Symptom

On a large, many-sequence genome the unified GFF3 is **not** byte-identical
run-to-run across thread counts. Measured on run116 (GCA of a 2.82 Gbp,
2246-sequence assembly), `make_unified_annotation.R` at `--threads 4` vs
`--threads 1`:

- `Repeat_Annotation_Unified.te_derived_trc.csv` — **identical** (the TE_origin
  layer is deterministic after the fixes below).
- `Repeat_Annotation_Unified.gff3` — **816 of 2.82 M features differ** (0.03 %),
  all `DANTE_LTR` / `DANTE_TIR` `transposable_element` (tier-1 structural)
  overlap-resolution artifacts, concentrated on the sequences that carry the
  TE_origin TRCs (i.e. around the te_sat-trimmed regions).

The differences are split-vs-kept-whole outcomes of the greedy longest-first
tier-1 trim (e.g. one run emits `5557–15187` whole + `375–5556` trimmed; the
other emits `375–5556`, `5557–7405`, `7406–15187`).

## Root cause

`make_unified_annotation.R` splits the genome into **batches whose composition
depends on the thread count** — `--threads 1` (or a single-sequence genome)
takes the `batches <- list(all_seqs)` path; `--threads > 1` calls
`make_batches()` with a per-thread target. Any per-batch operation whose result
depends on which features are grouped together, or on their arrival order, is
therefore thread-dependent. The tier-1 structural-overlap resolution still has
order-sensitivity beyond the one spot already fixed.

## Already fixed (in the TE-derived-TRC change)

- **`resolve_tier1_overlaps` tie-break** — the greedy longest-first order broke
  width ties by input position; now ordered by `(seqname, start, end, strand)`.
  Unit-tested for order-invariance (`tests/test_resolve_tier1_overlaps.R`).
- **TE-derived TRC decision made global** — `identify_te_derived_trcs` ran per
  batch; now runs once over each TRC's arrays genome-wide (+ stray-tolerant LCA).
  This removed the TE_origin-layer non-determinism (CSV now identical).

Together these reduced the differing-feature count but did **not** eliminate it,
so at least one more order-sensitive operation remains in the tier-1 path.

## Candidate spots to audit

- `resolve_within_tier()` — `disjoin(with.revmap)` then, per disjoint bin,
  `source_tool <- gr$source_tool[idx[1]]` and `source_tier <- min(...)`. The
  `idx[1]` pick depends on feature order within the bin → order-sensitive when a
  bin has ≥2 equal-tier features from different tools. (Used on t2/t4/t5, not t1
  directly, but worth checking for tier-1 interactions.)
- `trim_to_nonoverlap()` chains and the Step-1..7 ordering in `process_batch` —
  verify each is a pure function of (this batch's features, the reduced higher
  set) and independent of how features were grouped into the batch.
- Any remaining `order()` / greedy without a full coordinate tie-break.
- Identical-interval duplicates in tier-1 (same seqname/start/end/strand from
  DANTE_LTR complete+partial vs TIR/LINE) — these tie even under the coordinate
  tie-break and fall back to input order.

## Why CI doesn't catch it

The determinism gate (`.github/workflows/pipeline.yml`) runs the **small** and
**medium** fixtures, both of which are small enough to take the single-batch
path regardless of thread count — so per-batch composition never varies and the
bug cannot manifest. It only appears on large multi-sequence genomes, i.e. the
manual full-size check scenario CLAUDE.md already flags.

## Reproduction

```bash
IN=<run output dir>
for TH in 4 1; do
  make_unified_annotation.R ... \
    --tc_trc_table_default $IN/TideCluster/default/TideCluster_report/data/trc_table.tsv \
    --tc_trc_table_short   $IN/TideCluster/short_monomer/TideCluster_report/data/trc_table.tsv \
    --output out_t$TH.gff3 --threads $TH
done
diff <(grep -v '^#' out_t4.gff3 | cut -f1-8 | sort) \
     <(grep -v '^#' out_t1.gff3 | cut -f1-8 | sort) | grep -c '^[<>]'
```

## Recommended approach

The robust fix is to make the tier-resolution operations **batch-composition-
invariant** — either (a) canonicalise the per-batch input order and give every
greedy/disjoin resolver a full deterministic tie-break, or (b) hoist the tier-1
structural-overlap resolution out of the per-batch loop and run it once globally
(as was done for the TE_origin decision). Add a large-multi-sequence fixture (or
a threads-varying full-size check) to the determinism gate so a regression is
caught automatically.
