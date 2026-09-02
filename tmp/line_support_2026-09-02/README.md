# What actually stops LINE extension: the mmseqs prefilter, not the support rule

## The question

71.4% of LINE elements (245 of 343 on Boechera stricta) get ZERO extension on
both sides. The suspect was the support rule -- `min_group_alignments=5` plus
`support_fraction=0.5`, which sets the kept length to the k-th largest of a
group's alignment lengths.

## It is not the support rule

Measured on the per-flank alignment tables:

| flank set          | groups | median partners | dropped: <5 partners | k-th largest = 0 |
|--------------------|-------:|----------------:|---------------------:|-----------------:|
| ENDO_RT_5prime     |    142 |               2 |          125 (88%)   |            **0** |
| ENDO_RT_3prime     |     93 |               1 |           89 (96%)   |            **0** |
| ENDO_RT_RH_5prime  |     77 |              12 |            5 (6%)    |            **0** |
| ENDO_RT_RH_3prime  |     70 |               8 |           21 (30%)   |            **0** |

**No group that clears the size gate ends up with a k-th largest of zero.** The
loss is entirely upstream: most groups never reach 5 partners.

## What gates the partners

A pair is only ever aligned if `mmseqs easy-search` matches the **30 nt adjacent
to the domain core** at **>=80% identity**
(`global_local_aln.run_mmseqs_prefilter`, `prefilter_identity=0.8`). Admission:

| flank set          | seqs | possible pairs | aligned pairs | admitted |
|--------------------|-----:|---------------:|--------------:|---------:|
| ENDO_RT_5prime     |  264 |         34,716 |           171 |   0.49%  |
| ENDO_RT_3prime     |  264 |         34,716 |            85 |   0.24%  |
| ENDO_RT_RH_5prime  |   79 |          3,081 |           500 |  16.23%  |
| ENDO_RT_RH_3prime  |   79 |          3,081 |           292 |   9.48%  |

**That 0.8 is hardcoded in `run_all_vs_all_alignment` and exposed nowhere** --
no CLI flag in `dante_line.py` or `dante_tir_fallback.py`, no config key in the
Snakefile -- yet it gates the entire extension layer for both tools.

## Is 0.8 defensible?

`prefilter_identity_probe.py` measures the identity of exactly the 30 nt the
prefilter compares, with families defined from the CORE (the REXdb protein each
domain hit) so the test is independent of the flanks being judged.

```
ENDO_RT_5prime      SAME family n=811  median 0.63    DIFF family n=3575  median 0.37
ENDO_RT_RH_5prime   SAME family n= 40  median 0.63    DIFF family n=3461  median 0.43

threshold   admits SAME    admits DIFF
   0.6         60.2%           8.1%
   0.7         32.8%           3.7%
   0.8         10.0%           0.8%   <- shipped
   0.9          2.8%           0.2%
```

**The shipped 0.8 admits only 10% of genuine same-family pairs.** The median
same-family identity over that window is 0.63, so 0.8 sits well above the
middle of the distribution it is meant to accept. The signal is real -- same
family 0.63 vs different family 0.37, against ~0.25 for random DNA -- so a lower
threshold has room to work: 0.6 admits 6x more same-family pairs while still
rejecting 92% of different-family ones.

## Caveats

- **"Same family" is a proxy.** Two elements hitting the same REXdb protein share
  a lineage, which may be broader than a repeat family. If REXdb LINE proteins are
  coarse, some "same family" pairs are genuinely unrelated and their low identity
  is correct -- which would mean 0.8 is less wrong than it looks. The same/diff
  separation (0.63 vs 0.37) says the grouping captures something real, but not
  how sharply.
- One genome, two pattern sets.
- **More partners does not automatically mean better boundaries.** Lowering the
  threshold admits more pairs into the k-th largest, which will produce more and
  longer extensions. Whether those are correct is a separate question and must be
  measured against ground truth (tmp/diagonal_eval_2026-09-01/) before any change.

## A correction made during this investigation

I first reported that the prefilter samples the DISTAL 30 nt of the flank rather
than the core-adjacent end -- reading `extract_end_region`'s `end="5"` branch
(`seq[:30]`) without checking which `end` value the 5prime file is given. It is
given `end='3'` (dante_line.py:1445, "5' seqs -> fix 3' end"), so it takes
`seq[-30:]`, the core-adjacent window. Verified empirically: the core-adjacent
30 nt is at the END of the 5' flank in 134/134 elements, and 171/171 aligned
segments sit at that same end. **There is no orientation bug.**
