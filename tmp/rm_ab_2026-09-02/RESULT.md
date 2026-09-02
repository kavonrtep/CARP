# RepeatMasker A/B: do confirmed full-length elements improve LINE masking?

**Genome:** *Boechera stricta* (183 Mb, 343 LINE loci — the one genome where every
locus was examined, so its 25 confirmed elements are the complete set, not a sample).

**Arm A:** the LINE library the pipeline ships — 189 domain cores, 492,729 bp.
**Arm B:** the same 189 cores **plus** the 25 confirmed full-length elements
(4,759–7,372 bp, median 6,106 bp), 644,937 bp.

A confirmed element is one where BOTH ends were observed — a poly-A tail and a
TSD clearing the measured chance floor. Its boundaries are not inferred, so it
cannot be chimeric in the way that motivated the whole boundary investigation.

RepeatMasker 4.1.2-p1, `-qq -nolow -no_is`, identical settings both arms.

## Result

| library | hits | masked bp | % of genome |
|---|---:|---:|---:|
| A: cores only | 5,120 | 3,968,368 | 2.17 % |
| B: cores + 25 confirmed | 7,482 | 5,086,016 | 2.78 % |
| **difference** | **+2,362** | **+1,117,648** | **+0.61 pts** |

**+28.2 % more LINE sequence masked**, from 25 added sequences — a 13 % increase
in library bp producing a 28 % increase in masking.

## And the new masking is cleaner than what is already shipped

Contamination measured against independent structural annotation from the same
run (DANTE_LTR elements, DANTE_TIR primary elements, TideCluster arrays —
21.6 Mb of the genome):

| | masked bp | on other elements | rate |
|---|---:|---:|---:|
| A: cores only | 3,968,368 | 349,263 | **8.80 %** |
| B: cores + confirmed | 5,086,016 | 359,473 | **7.07 %** |
| **the new masking alone** | 1,118,663 | 10,210 | **0.91 %** |

The 1.12 Mb of newly masked sequence is **99.09 % clean**. Adding these elements
does not merely avoid making things worse — it lowers the library's overall
contamination rate, because the added masking is an order of magnitude cleaner
than the cores' own.

That is the opposite of what inferred extensions do: measured earlier, inferred
flank contamination ran 1.3–4.2 % and *rose* with extension length.

## Why it works

The cores are ~2.1 kb of ENDO…RT. A real LINE is 4–7 kb. Masking with cores alone
finds the domain region of each copy and stops; the UTRs and ORF1 go unmasked
even though they are real element. A confirmed full-length element carries those
regions with **observed** boundaries, so it masks them without the risk that made
inferred extensions unusable.

## What this does NOT establish

- **One genome**, and the most favourable one: *Boechera*'s 5.0 % confirmation
  rate is the highest of those tested (others: 1.2–4.1 %).
- **Family coverage limits the ceiling.** These 25 elements belong to 19 of 64
  LINE families, covering 62.7 % of loci. Families with no young copy gain
  nothing — this adds to the core library, it cannot replace it.
- **"Clean" means "not overlapping DANTE_LTR / DANTE_TIR / TideCluster."**
  Same-class or unannotated contamination is invisible to this test.
- **LINE library only.** The interaction with the full combined library —
  `reduce_library`, containment reduction, and competition with LTR/TIR consensi
  during masking — is untested. The genome-wide `Class_I/LINE` figure would not
  move by the full +0.61 points, because the unified annotation resolves LINE
  against higher-tier calls.

## Reading

The idea works on this genome, and for a reason that generalises: an element whose
ends were measured is a better library entry than one whose ends were guessed.
The honest next step is a second genome with a lower confirmation rate, and then
the same test inside the full library rather than the LINE library alone.
