# TPRT-signature boundary detection for DANTE_LINE — exploratory probes

Testing whether a LINE boundary can be found from the **TPRT signature**
(poly-A tail + target-site duplication) instead of from flank alignment.
Motivation: every failure mode diagnosed in the 1.7.0 boundary work came from
"the flanks agree", which a shared neighbouring repeat produces just as readily
as a real element. A TSD pair cannot be manufactured by an adjacent LTR-RT.

## `polya_probe.py` — WRONG, kept as the record of two mistakes

Searched 600 bp for "a 20 bp window >=70% A". Both wrong:
- 600 bp is too short. The RT domain end is NOT the element end — the ORF2
  C-terminal domain plus the 3'UTR lie between them (measured median distance
  to the tail: ~1.6 kb).
- "70% A in 20 bp" is satisfied by coding AT-richness. The 3' side scored
  *below* the 5' negative control (0.77x, 0.00x) — the test found ORF, not tail.

## `polya_probe2.py` — corrected

Genuine homopolymer run (>=N A, <=1 mismatch) over a realistic range, with the
element's own 5' side as the NEGATIVE CONTROL: a real tail must be enriched 3'
and not 5'.

```
genome     window  minrun     3'   5'ctl  enrich   median dist
boechera     1000      15   4.7%    0.6%   8.00x           134
boechera     2500      15  22.2%    7.0%   3.17x          1597
lycopus      2500      15   5.7%    0.0%     inf          1520
boechera     2500      10  77.8%   66.8%   1.17x             —
```

**Conclusions carried into the implementation plan:**
- The anchor is real: 3–8x enriched on the correct side at minrun 15.
- minrun must be >=15. At 10 the enrichment collapses to ~1.2x — background
  A-runs swamp it.
- Search out to ~2.5 kb, not "a few hundred bp".
- **Yield is the binding constraint**: 22% (boechera) / 5.7% (lycopus) of loci.
  TSD confirmation is a subset of that, so a fully structural boundary will be
  available for a minority of loci, not the majority.

Always report against a negative control — an absolute A-fraction is meaningless
in an AT-rich genome, which is exactly how v1 got it backwards.

---

# Phase 0 + 0b results (2026-09-01)

## `tsd_finder.py` — direct-repeat finder with a decoy-window null

Null = the same query searched in windows drawn from distant loci on the SAME
sequence. Preserves query and realistic genomic composition/repeat structure; a
shuffle would destroy exactly the local structure that makes chance hits likely
and would understate the false-positive rate.

**0b validation, wheat DANTE_TIR (true TSDs 8-12 bp, the LINE-range case):**
```
window=6000  k=8-20  n=28
  GEOMETRY self-check   28/28 (100%)   <- coordinates verified independently
  RECALL                25/28 (89.3%)
  null: median 6  p95 9  p99 10
  loci beating EVERY decoy: 19/28 (67.9%)
```
The finder works. The difficulty is **discrimination**, not detection.

## `null_curve.py` — the chance floor (the real Phase 0 deliverable)

A perfect TSD of length L scores exactly L (3*matches - 2k = k when perfect), so
a TSD is separable from chance only if L exceeds the chance floor of the window
it is searched in. Measured on Lycopus, no ground truth needed:

```
  window    p99   min TSD length separable at p99
     500      7                                  8
    1000      9                                 10
    2000     10                                 11
    4000     10                                 11
    8000     11                                 12
```

**This is why the domain anchor matters.** Bounding the 5' search to ~2 kb (the
LINE per-side vocabulary cap) sets the bar at 11 bp; an unbounded multi-kb search
pushes it to 12+. Published LINE TSDs run ~9-18 nt with the consensus clearest at
>=12 nt, so a substantial but not overwhelming fraction clears the bar.

**Corollary — this method cannot be validated on canonical LTR TSDs.** A perfect
5 bp TSD scores 5, below every floor above. Measured: wheat DANTE_LTR at k>=5
gives 12.7% recall, and correct hits score BELOW wrong ones (median 6 vs 7),
because a chance 12-mer with mismatches outscores a perfect 5-mer. Use only the
>=8 bp ground truth (wheat DANTE_TIR: 739 elements at 8-12 bp).

## Two errors this phase caught — both now guarded

1. **`TSD=not_found` is a 9-character SENTINEL, not a 9 bp TSD.** A naive
   `len()` counts it as a LINE-range length, and it is the single most abundant
   "TSD" in the file (wheat: 108,581). The first validation run was silently
   asking the finder to locate repeats that do not exist -> 0% recall.
   `load_truth()` now requires the value to be pure ACGT.
2. **The search window is measured back from the ANCHOR (element 3' end).** If it
   is shorter than the element, the true 5' TSD copy is simply not inside it, and
   the finder scores 0% for a reason unrelated to the finder. Elements longer
   than the window are now skipped and counted.

The **geometry self-check** (do the two TSD copies actually match in the genome?)
is permanent, and is what separates "the coordinates are wrong" from "the finder
is wrong". It should read ~100%.
