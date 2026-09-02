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

**The one thing worth a look is the rDNA.** A LINE extending into 45S rDNA is
not biologically expected, and 11 kb of it survives the gate. It may be genuine
contact at array edges — LINEs do insert near rDNA — but it is the least
plausible of the displaced classes and the reason to look before shipping,
rather than a reason not to ship.
