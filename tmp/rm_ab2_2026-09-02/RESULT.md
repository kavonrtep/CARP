# Second genome: the benefit is real but strongly genome-dependent

Same A/B as *Boechera*, on a genome with less than half the confirmation rate.

**Target:** GCA_986270105 chromosome OZ488822.1 (474 Mb). 250 confirmed elements
from a full extraction over all 11,900 LINE loci — a **2.1 % confirmation rate**
against *Boechera*'s 5.0 %.

| library | hits | masked bp | % of chromosome |
|---|---:|---:|---:|
| A: cores only (2,701 seqs) | 157,334 | 86,503,186 | 18.25 % |
| B: + 250 confirmed | 158,904 | 87,126,283 | 18.38 % |
| **difference** | **+1,570** | **+623,097** | **+0.13 pts** |

**Relative gain: +0.7 %** — against **+28.2 %** on *Boechera*.

## Still clean, just not needed here

| | contamination |
|---|---:|
| A: cores only | 11.40 % |
| B: + confirmed | 11.32 % |
| **the new masking alone** | **0.71 %** |

Same picture as *Boechera* (0.91 %): the added elements are far cleaner than the
cores they join, so they lower the library's overall rate. The change is never
harmful — it is simply not always useful.

## Why the two genomes differ

| | *Boechera* | GCA_986270105 |
|---|---:|---:|
| final `Class_I/LINE` | 1.33 % | 8.75 % |
| LINE loci | 343 | 11,900 |
| cores in the library | 189 | 2,701 |
| cores alone mask | 2.17 % | 18.25 % |
| **gain from confirmed** | **+28.2 %** | **+0.7 %** |

The benefit tracks how much LINE the **core library already misses**, not the
confirmation rate. On *Boechera* 189 cores mask 2.17 % of the genome and leave
most of each element's UTRs and ORF1 uncovered, so 25 full-length elements add a
great deal. On GCA_986270105, 2,701 cores drawn from 11,900 loci already cover
the population densely — the families are represented many times over, and
full-length copies have little new ground to reach.

This is the more common situation for a LINE-rich genome, and *Boechera* is
probably the more common situation for a LINE-poor one.

## Reading

The proposal is **safe everywhere and valuable somewhere**. It cannot be sold as
a general improvement on the strength of *Boechera* alone: on a second genome the
same procedure bought 0.7 %.

What would make it worth shipping is that the cost is small (250 sequences on
2,701) and the failure mode is benign (0.71 % contamination, better than the
library's own 11.4 %). A knob that helps 28 % on some genomes, 0.7 % on others
and never hurts is defensible — but it should be presented that way, not as a
uniform gain.

**Unresolved:** both tests used the LINE library alone. Inside the full combined
library the confirmed elements also compete with LTR and TIR consensi, and pass
through `reduce_library` and containment reduction, either of which could remove
them as redundant with their own cores.
