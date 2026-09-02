# Unified-annotation effect of the LINE bounds (rc1 vs gated vs ungated)

## What was actually run, and why it is valid

Not the full pipeline — this sandbox has no container runtime and the input
assemblies are not mounted. Instead `make_unified_annotation` was run three times
per genome on the **baseline run's own inputs**, changing only which
`DANTE_LINE.gff3` it consumed.

That substitution is sound because **the LINE library is unchanged by the
bounds**: `dante_line_library_source: core` means the library is built from the
domain cores, which the bounds do not touch. Verified per genome — same headers,
same sequences, same total bp in both arms (one genome's file differs only in
record order, which `canonical_sort_fasta.py` normalises before RepeatMasker).
So `DANTE_LINE.gff3` is the only input to this step that the change affects, and
every other layer is byte-identical by construction.

**What this therefore does NOT prove:** that a real end-to-end run reproduces it.
It should, since the only changed input is held fixed correctly, but the
RepeatMasker and TideCluster layers were reused rather than recomputed.

## GCA_973357735.1 — the sensitive genome (63% of inferences unconverged)

| | rc1 | **gated** | ungated |
|---|---:|---:|---:|
| Class_I/LINE bp | 2,908,290 | **3,451,616** | 3,936,016 |
| LINE as % of genome | 0.56% | **0.66%** | 0.76% |
| total annotated bp | 297,536,067 | **297,958,412** | 298,194,570 |
| total annotated change | — | **+0.14%** | +0.22% |

Where the LINE gain comes from:

| | LINE gain | from unannotated | **taken from other classes** |
|---|---:|---:|---:|
| gated | +543,326 | 422,345 | **120,981 (22%)** |
| ungated | +1,027,726 | 658,503 | **369,223 (36%)** |

Classes that lose ground:

| class | rc1 | gated | ungated |
|---|---:|---:|---:|
| rDNA/45S_rDNA/18S | 1,328,359 | 1,317,244 | 1,195,240 |
| TRC_33 (satellite) | 5,730,307 | 5,662,338 | 5,615,395 |
| TRC_4 | 2,210,769 | 2,192,597 | 2,170,197 |
| Low_complexity | 1,556,260 | 1,555,232 | 1,546,684 |

## GCA_937616625.2 — the control

Total annotated +3,390 bp (+0.00%). **Class_I/LINE is the only class that
changes at all**, and the Level-1 feature count is identical (336,919 both
ways). A clean no-op.

## Reading

The gate roughly halves everything: the LINE gain (+543 kb vs +1,028 kb), the
share of that gain taken from other classes (22% vs 36%), and the rDNA loss
(11 kb vs 133 kb, a 92% reduction).

**The displacement is small and edge-shaped, not wholesale.** Each affected
class loses ~1% of itself (TRC_33 1.2%, rDNA 18S 0.8%), which is what boundary
adjustment at contact points looks like, rather than a LINE consuming an array.
Total annotation moves +0.14% and LINE moves +0.10 points of genome.

## The rDNA question — checked, and it is benign

The 11 kb of displaced 45S rDNA was the one implausible-looking item. It is not
what it looked like.

**Every base pair of it is at an array edge.** Of 11,085 bp of overlap between
gated LINE features and rc1 rDNA arrays, 100% lies within 100 bp of an array end
(distance 0 for every large case) and **0 bp is more than 1 kb inside an array**.

And the affected "arrays" are not arrays:

| rc1 rDNA array size | arrays | total bp | bp lost to LINE | % lost |
|---|---:|---:|---:|---:|
| under 200 bp | 352 | 21,164 | **9,511** | 44.9% |
| 200 bp – 1 kb | 462 | 309,990 | 408 | 0.13% |
| 1 – 10 kb | 799 | 3,898,575 | 1,166 | 0.03% |
| **over 10 kb (the real arrays)** | 5 | 55,277 | **0** | **0.00%** |

86% of the loss comes from rDNA fragments **under 200 bp** — isolated 44–585 bp
calls, of which the largest single overlap is 585 bp. All five genuine rDNA
arrays (10.5–11.9 kb) lost **exactly zero**.

So no LINE is extending into ribosomal DNA. LINE features are absorbing short,
isolated rDNA-labelled fragments that sit at their own boundaries. A 46 bp "rDNA"
call is a RepeatMasker hit to the rDNA library, not a ribosomal gene, and folding
it into an adjacent LINE is at worst neutral. (Stated as observation: the
fragments' small size and isolation are measured, their spuriousness is
inferred.)

## Verdict

The gated bounds are safe to ship on this evidence. The control genome is a
clean no-op; the sensitive genome gains 0.10 points of genome as LINE, 78% of it
from previously unannotated sequence, and the displacement is confined to array
edges and short fragments with the real satellite and rDNA arrays untouched.
