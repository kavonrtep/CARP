# How LINE boundaries are decided, and what was tried and rejected

CARP finds a LINE from its protein domains outward. The domains are observed;
everything past them is inferred. This page explains how that inference works,
which alternatives were tested, and why they were not adopted — with the
measurements behind each decision, so the reasoning can be re-checked rather
than taken on trust.

For the short user-facing summary see **LINE identification** in the README.

---

## What the pipeline does

1. **Core.** DANTE locates the protein domains a LINE encodes — endonuclease
   (ENDO) and reverse transcriptase (RT), sometimes with RNase H (RH). Two or
   three of these close together on the same strand define an element, and the
   span they cover is its **core**. This is the only part read directly from the
   sequence.

2. **Extension.** Up to 10 kb of flank is taken on each side, stopping early at
   any tandem array, LTR retrotransposon or DNA transposon already annotated, so
   a LINE cannot run through its neighbour. The flanks of many copies of the same
   family are compared with each other, and the element is extended to the
   distance where at least half of them still match. An element with fewer than
   five comparable partners gets no extension.

3. **Limits.** If the agreement between copies fades on its own, the boundary is
   evidence and is kept — up to 3,400 bp on the 5′ side, and 4,500 bp measured
   from the *start* of the core to the 3′ end. If it instead runs to the edge of
   the search window, no boundary was found and the conservative limits from
   1.7.0rc1 apply (2,000 / 800 bp). Values live in
   `classification_vocabulary.yaml`.

4. **What is trusted where.** The full element goes into `DANTE_LINE.gff3` and
   the unified annotation. **The RepeatMasker library is built from the cores
   only.** A boundary that is slightly wrong costs one locus; a wrong consensus
   in the search library mislabels every sequence it matches, genome-wide.

---

## Why the 3′ limit is measured from the core's start

The distance from the end of the domains to the end of the element is **not** a
stable quantity. It is only what is left over after the domain annotation stops,
and it varies **13.8-fold** (p90/p10) across elements. Whether RNase H was
annotated changes it completely.

What is stable is the distance from the **start** of the core to the element's 3′
end: **1.16-fold** variation, per-genome medians of 3,779–4,247 bp across seven
genomes. Correlation between core length and leftover 3′ extension is
**r = −0.936** — they are two halves of one constant.

| element type | n | median core | median 3′ leftover | core start → 3′ end |
|---|---:|---:|---:|---:|
| ENDO–RT (2 domains) | 157 | 2,148 bp | 1,728 bp | 3,899 bp |
| ENDO–RT–RH (3 domains) | 61 | 3,980 bp | 162 bp | 4,158 bp |

Same biology, measured from different starting points. A flat 3′ cap treats the
two arrangements completely differently for no biological reason.

---

## Rejected: finding boundaries from the insertion scar

When a LINE inserts by target-primed reverse transcription it leaves two marks: a
**poly-A tail** at the 3′ end and a **target-site duplication (TSD)** — a short
run of host DNA copied on both sides. Unlike sequence similarity, neither can be
produced by a neighbouring repeat, so this looked like a way to replace inference
with observation.

It was implemented, measured, and **rejected on yield**.

| genome | LINE loci | poly-A found | TSD confirmed | end-to-end |
|---|---:|---:|---:|---:|
| *Boechera stricta* | 343 | 22.2 % | 22.4 % of those | **5.0 %** |
| *Lycopus europaeus* | 53 | 5.7 % | 33.3 % | **1.9 %** |
| *Triticum aestivum* | 22,050 (300 sampled) | 3.7 % | 18.2 % | **0.7 %** |

The marks survive only in **young** insertions, so the method is weakest exactly
where the problem is worst: wheat, whose LINE annotation motivated this whole
investigation, has the most decayed elements and the lowest yield.

It is not a tuning problem. Doubling the tail search window (2,500 → 5,000 bp)
added candidate tails but **zero** additional confirmed loci — the extra tails are
background A-runs with no TSD behind them.

### But the two tests together are extremely specific

Low yield is not the same as low reliability, and the combination is worth
recording because it is unusually clean.

Measured by running the identical chain on the element's **5′ side**, where a
LINE has no poly-A tail so anything found is spurious by construction:

| | 3′ (real) | 5′ (decoy) |
|---|---:|---:|
| poly-A run found | 359 / 2,743 (13.1 %) | 77 / 2,743 (2.8 %) |
| **+ TSD confirmed** | **72 (2.6 %)** | **2 (0.07 %)** |

- spurious poly-A rate: **2.8 %**
- chance TSD given a spurious anchor: **2.6 %** (2 of 77)
- product, if the two are independent: **0.073 %**
- **observed joint false-positive rate: 0.073 %** (2 of 2,743)

The product matches the observation exactly, which is the interesting part: the
two tests are **statistically independent**. Finding a poly-A run tells you
nothing about whether a matching direct repeat will be found upstream, so their
error rates multiply rather than overlap.

The practical consequence: a locus confirmed by *both* is right about **97 %** of
the time, and confirmed calls are enriched **36-fold** over the decoy side. A TSD
alone would not carry that — a 9–18 bp repeat searched over kilobases of AT-rich
DNA is found by chance often enough that in a 2 kb window a repeat must be at
least **11 bp** to stand out at all, which many decayed TSDs are not.

So the scar is **too rare to find elements with, but reliable enough to measure
them with**. That is precisely how it is used: the limits in step 3 above come
from 218 elements across 7 genomes where a poly-A tail and an intact TSD made the
true end directly visible. The pipeline does not search for these features at run
time; it inherits the numbers they produced.

---

## Rejected: reading the shape of the flank alignment

A second idea was to use the *geometry* of the alignment rather than its length —
a real boundary should show a clean break, while running into a neighbouring
repeat should show a characteristic disturbance (a gap burst, or the alignment
drifting off the diagonal).

The signal is real but **weakest where it is needed**:

| genome | over-extensions | caught by drift |
|---|---:|---:|
| *Boechera* (sparse) | 25 | **84 %** |
| *Lycopus* (dense) | 107 | **48 %** |
| wheat (dense) | 158 | **63 %** |

The predicted failure mode is the observed one: if a LINE inserts at a consistent
position inside another repeat, the neighbour is colinear across copies and there
is no disturbance to detect.

Applying it anyway cost **30–40 bp of legitimate element for every 1 bp of
over-extension avoided**. On TIR elements, where over-extension was already
0.0 %, it was pure loss.

---

## Rejected: loosening the flank-pair prefilter

Most LINE elements get no extension because the mmseqs prefilter admits very few
flank pairs. Two levers were tested under a decision rule fixed in advance:

- **`--prefilter-identity`** is nearly inert. Lowering it 0.8 → 0.6 moved admitted
  pairs 342 → 344, because mmseqs requires a shared exact k-mer *before* it ever
  evaluates identity.
- **`-k 7`** works mechanically (342 → 2,014 pairs, and the extra pairs are
  *enriched* for genuine family members), but its benefit is unpredictable:
  coverage gain was **+3.0, +20.0 and +7.2** percentage points on three genomes.
  One genome alone would have given either an obvious yes or an obvious no.

Rejected for variance, not for harm — purity stayed 98.4–100 %, so the extra
extension was real element. Full record in
`tmp/prefilter_eval_2026-09-02/` (`PREREGISTERED.md`, `VERDICT.md`).

---

## Reproducing any of this

The harnesses are kept because each answers a question that will recur:

| directory | what it measures |
|---|---|
| `tmp/diagonal_eval_2026-09-01/` | scores any boundary estimator against exact boundaries (DANTE_LTR / DANTE_TIR elements whose ends are known), ~1 min per genome, no pipeline run |
| `tmp/line_tprt_2026-09-01/` | poly-A and TSD detection, the chance floor, the joint false-positive rate |
| `tmp/prefilter_eval_2026-09-02/` | the pre-registered prefilter experiment |
| `tmp/bounds_ab_2026-09-02/` | bounds A/B at the `dante_line` step |
| `tmp/unified_ab_2026-09-02/` | the same, at the unified-annotation level |

**Score purity, not recall.** These layers exist to produce a reliable core, not a
complete element, so under-extension is not a defect. A scorer that rewards
recovering the full element will make a correctly conservative estimator look
broken.
