# RESOLVED: tier-1/2 overlap resolution non-determinism (large multi-sequence genomes)

**Status:** resolved (2026-07). Root cause was a `trim_to_nonoverlap` side-effect,
not a tie-break gap. Surfaced during the TE-derived-TRC domain-rhythm work (which is
itself deterministic); fixed in `scripts/make_unified_annotation.R`.

## Symptom

On a large, many-sequence genome the unified GFF3 was **not** byte-identical
run-to-run across thread counts. Measured on run116 (a 2.82 Gbp, 2246-sequence
assembly), `make_unified_annotation.R` at `--threads 1` (one batch) vs a fragmented
multi-batch run:

- `Repeat_Annotation_Unified.te_derived_trc.csv` — identical (the TE_origin layer
  was already deterministic).
- `Repeat_Annotation_Unified.gff3` — **816 of 2.82 M features differed** (0.03 %),
  all `DANTE_LTR` / `DANTE_TIR` / `DANTE_LINE` (tier-1) and a few `DANTE` (tier-2)
  `transposable_element` overlap-resolution artifacts, concentrated on sequences
  that carry TE-derived-TRC arrays.

The differences were split-vs-kept-whole outcomes: one run emitted the *greedy
longest-first* result (`5557–15187` whole `+`, `375–5556` trimmed `*`); the other
emitted a naive *disjoin* of the same two overlapping elements (`375–5556` `+`,
`5557–7405` `+`, `7406–15187` `+`, all untrimmed).

## Root cause

`process_batch`'s TE-derived-satellite pre-pass trimmed the tier tracks against the
satellite regions **before** tier-1 overlaps were resolved:

```r
t1 <- trim_to_nonoverlap(t1, te_sat_r, min_len)   # t1 still self-overlapping here
...
t1 <- resolve_tier1_overlaps(t1, min_len)          # greedy — but t1 already disjoined
```

`trim_to_nonoverlap` computes `disjoin(c(lower, higher))`. When `lower` (t1) still
carries its **internal** overlaps, that disjoin decomposes them into disjoint pieces
as a side effect — both overlapping elements kept, split at the overlap, the contested
span assigned to the lower **input index** (`min(revmap)`). The subsequent
`resolve_tier1_overlaps` then sees an already-non-overlapping `t1`, so `findOverlaps`
returns nothing and the greedy is a silent no-op.

The trigger is batch-scoped: `te_sat` is non-empty only in a batch that **contains a
TE-derived-TRC array**, and `trim_to_nonoverlap`'s `any(lower_overlaps_higher)` guard
is evaluated over the **whole batch**. So:

- `--threads 1` (one batch): `te_sat` always present → the entire genome's `t1` is
  pre-disjoined → greedy defeated everywhere (even on sequences far from any te_sat).
- multi-batch: only batches that happen to contain a te_sat array pre-disjoin their
  `t1`; a sequence whose batch has no te_sat is resolved correctly.

Same input, different thread count → different batching → different `t1` → the 816
differing features. It was **also an annotation-quality bug**, not merely cosmetic:
whenever a batch contained a te_sat array, overlapping structural TEs in that batch
were over-fragmented (naive disjoin, arbitrary metadata) instead of resolved greedy
longest-first.

The tier-2 twin: `t2 <- trim_to_nonoverlap(t2, te_sat_r)` pre-disjoined `t2` the same
way, which made the Step-2 `resolve_within_tier` a no-op in te_sat batches (strand
preserved) but not in others (disjoined → strand `*`) — the residual 26 strand-only
`DANTE` differences.

## Fix (`scripts/make_unified_annotation.R`)

1. **Resolve tier-1 overlaps BEFORE the te_sat trim.** `resolve_tier1_overlaps(t1)`
   now runs first, so `t1` is internally non-overlapping when trimmed against
   `te_sat_r` — the trim only carves, never disjoins, and the greedy is applied
   consistently regardless of batching.
2. **Fold the te_sat carve of `t2` into Step 2.** The separate pre-pass
   `trim_to_nonoverlap(t2, te_sat_r)` is gone; `t2` is now trimmed against
   `reduce(level1)` (= te_sat ∪ tier-1) inside Step 2 and disjoined once by
   `resolve_within_tier`, batch-invariantly. Identical to the old behaviour when
   there are no TE-derived TRCs (then `level1 == t1` at Step 2).
3. **Full content tie-break in the greedy.** `resolve_tier1_overlaps`'s longest-first
   order now breaks residual ties on `(seqname,start,end,strand,classification,
   source_tool)`, so genuinely identical intervals from different tools/classes
   resolve independently of input position.
4. **Canonical mcols column order before export.** rtracklayer emits column-9
   attributes in mcols order, which was inherited from whichever batch combined
   first (thread-dependent — `ID` first under multi-batch, later under one batch).
   `finalise_output` now reorders mcols to a fixed order before `export`, so the
   GFF3 is byte-identical across thread counts, not just data-identical.

## Validation

Same-machine A/B on the full run116 data (single batch vs 34 batches), current code:

- geometry (`cut -f1-8`, sorted) diff: **816 → 0**
- full record incl. col-9 attributes: **0** (after fix 4)
- `te_derived_trc.csv`: identical; feature counts equal (2 822 010 either way)

Regression test: `tests/test_resolve_tier1_overlaps.R` now asserts (a) the greedy is
order-invariant with the extended tie-break and (b) the te_sat trim must run **after**
resolve — a batch that merely *contains* a te_sat array must not change a far-away
sequence's tier-1 result (fixed = 2 greedy pieces; the old trim-before-resolve
over-fragments to 3).

## Still worth doing (separate, optional)

The determinism CI gate (`.github/workflows/pipeline.yml`) runs the **small** and
**medium** fixtures, both single-batch regardless of thread count, so this class of
bug could never manifest there. Adding a small **multi-sequence** fixture large enough
to split into >1 batch (or a threads-varying full-size check) would catch a regression
of any of the four fixes above automatically. This bug was found only by a manual
full-size A/B.
