# Candidate genomes for the rc4 end-to-end tests

What rc4 adds over rc2 is the **`Status` work**: LINE elements with both ends
observed (poly-A tail + TSD) are marked `Status=complete` and given their
measured span; everything else is `inferred` or `core`.

The size of that annotation change is driven by the **confirmation rate**, so
that is what these candidates are ranked by. Rates were measured with the
shipped `scripts/line_complete_elements.py` against each genome's existing
`DANTE_LINE.gff3` — seconds per genome, no pipeline run.

`candidates.tsv` has the full table. Survey was still running when this was
written; more genomes may be added.

## The range found

| accession | Gb | LINE loci | complete | rate | run |
|---|---:|---:|---:|---:|---|
| GCA_054642775.1 | 0.76 | 1,598 | 247 | **15.5 %** | run-000150 |
| GCA_978657495.1 | 0.14 | 155 | 21 | 13.5 % | run-000126 |
| GCA_937616625.2 | 0.45 | 633 | 32 | 5.1 % | run-000164 |
| GCA_018361405.1 | 0.19 | 343 | 17 | 5.0 % | run-000165 |
| GCA_965151555.2 | 0.32 | 480 | 17 | 3.5 % | run-000160 |
| GCA_050576095.1 | 1.39 | 587 | 18 | 3.1 % | run-000104 |
| GCA_964197575.1 | 1.75 | 822 | 13 | 1.6 % | run-000128 |
| GCA_042453785.1 | 0.65 | 1,788 | 16 | 0.9 % | run-000140 |
| GCA_964200675.1 | 0.48 | 431 | 1 | 0.2 % | run-000134 |
| GCA_973357735.1 | 0.55 | 550 | 0 | **0.0 %** | run-000154 |
| GCA_055776105.1 | 0.40 | 264 | 0 | 0.0 % | run-000163 |

Rates span **15.5 % to 0 %**, which is wider than expected and covers every case
worth testing.

## Suggested set

| role | accession | why |
|------|-----------|-----|
| **maximum effect** | GCA_054642775.1 | 247 complete elements, by far the most — if the change breaks anything it breaks here. Also the genome the original chimera fix moved most (LINE 12.43 % → 2.39 %), so it has history. |
| **high rate, tiny genome** | GCA_978657495.1 | 13.5 % at 0.14 Gb — fastest possible run at a high rate, good for a quick signal |
| **middle** | GCA_937616625.2 | 5.1 %; already has rc1 and rc2 baselines (run-000230, run-000237), so it gives a three-way comparison |
| **do-no-harm control** | GCA_973357735.1 | 0 % complete — the annotation must be identical to rc2. It is also the genome where the bounds bind hardest (494 of 550 elements capped), so it doubles as the bounds stress test |

A control is only a control if it can fail: GCA_973357735.1 must come out
**byte-identical** to its rc2 result for `Class_I/LINE`. Any movement there means
something changed that should not have.

## What to check in each run

- `Status` distribution in `DANTE_LINE.gff3` — `complete` / `inferred` / `core`
- No `complete` element whose `TSD=` is periodic (the telomere failure mode; the
  shipped guards should make this zero)
- `Class_I/LINE` versus the genome's rc2 baseline — expect a small rise
  proportional to the confirmation rate, and **nothing else moving much**
- `DANTE_LINE/LINE_complete_elements.fasta` present and its record count matching
  the `Status=complete` count

## Caveat

These rates come from each genome's EXISTING `DANTE_LINE.gff3`, produced by
whichever version that run used. rc4 also carries the bounds and convergence-gate
changes, which alter element spans slightly before the confirmation step sees
them, so the realised rate may differ a little from the table.
