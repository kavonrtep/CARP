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
