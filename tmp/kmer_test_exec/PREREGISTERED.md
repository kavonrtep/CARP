# Pre-registration — mmseqs prefilter identity for DANTE_LINE

**Written and committed BEFORE any threshold was run.** Six genomes is a small
sample; a threshold tuned to look good on them can easily fail elsewhere. This
file fixes the split, the metric and the decision rule in advance so the result
cannot be chosen after the fact.

## What is being tested

`prefilter_identity` (hardcoded 0.8 until now, newly exposed as
`--prefilter-identity`) gates which flank PAIRS are aligned at all. It compares
the 30 nt adjacent to the domain core. At 0.8 only 0.2-16% of possible pairs are
admitted, most groups then fall below `--min-group-alignments 5`, and 71% of
LINE elements get no extension.

Thresholds to try: **0.8 (shipped, the null hypothesis), 0.7, 0.6.**

## Ground truth

LINE loci where BOTH boundaries were observed — a poly-A tail found and a TSD
clearing the measured chance floor — from `tmp/line_tprt_2026-09-01/w6_*.tsv`,
filtered to score < 20 (saturated matches are not TSDs) and element <= 7500 bp
(longer is not a plant LINE).

## Split — fixed now, by accession, NOT by result

| role | genomes | usable confirmed loci |
|------|---------|----------------------:|
| **CALIBRATION** | boechera (GCA_018361405), GCA_977018505, GCA_986270105 | 25 + 25 + 27 = 77 |
| **VALIDATION** | GCA_986264755, GCA_965641495, GCA_984574005 | 105 + 20 + 13 = 138 |

Validation deliberately holds the largest single set (105 loci). Wheat is
excluded: 3 usable loci is too few and it is 14 Gb.

**The validation genomes are not looked at until a threshold has been chosen on
calibration alone. No retuning afterwards.**

## Metrics

For each confirmed locus, compare the extension DANTE_LINE emits against the
true extension:

- **PURITY** = emitted bp lying inside the true element / emitted bp.
  The primary metric. This is the quantity the project cares about: an
  over-extended consensus mislabels whatever its tail matches genome-wide,
  a truncated one only loses a little masking.
- **COVERAGE** = emitted bp inside the true element / true extension bp.
  Secondary. How much real element is recovered.

Reported **per genome**, and judged on the **worst** genome, never the pooled
mean — a pooled mean lets one good genome hide a bad one.

## Decision rule

Lower the threshold only if **all** of the following hold.

On **every** calibration genome:
1. coverage increases by **>= 10 percentage points** versus 0.8, and
2. purity stays **>= 95%**.

Then, on the **validation** genomes at that single chosen threshold, with no
further tuning:
3. purity stays **>= 95%** on every one, and
4. coverage is still higher than at 0.8 on every one.

**If any condition fails, keep 0.8.** A threshold that helps on average but hurts
one genome is overfitting, and is rejected.

The 95% purity floor is set from earlier measurements on the same kind of data:
TIR extension purity was 99.9%, LTR 94-99%. Below 95% we would be shipping worse
boundaries than the layers we already trust.

## What this cannot settle

Confirmed loci are the youngest, least-decayed elements — the easiest cases and
probably the longest. Purity measured here is an upper bound for the whole
population. A threshold that passes is not thereby proven safe for degraded
copies; it is only not-rejected.

---

# Amendment 1 — the k-mer test (written BEFORE running it)

The `--prefilter-identity` arm above was **rejected**: on the first calibration
genome all three thresholds were byte-identical, so the >=10 point coverage bar
failed immediately. Cause established afterwards: mmseqs requires a shared exact
k-mer *before* it evaluates identity, and on the 30 nt window compared here its
default k discards most related pairs unseen.

So the real parameter is the k-mer length, newly exposed as `--prefilter-kmer`.

**Everything in the original design carries over unchanged** — same calibration
and validation genomes, same purity/coverage definitions, same worst-genome
judging, same accept/reject thresholds, same 95% purity floor. Only the
parameter under test changes. Recording it as an amendment rather than a fresh
document so that the identity arm's failure stays on the record: this is the
second knob tried, not the first, and that matters when reading a positive
result.

Values to try: **0 (mmseqs default, the null hypothesis) and 7.**

7 is not tuned here. It is the value `scripts/reduce_dimer_library.py` already
passes for the same reason (short nucleotide queries), so it is a transplant
from existing practice in this codebase, not a number picked to fit these
genomes. **No search over k is performed.** If 7 fails the rule, the answer is
that the default stays -- not that some other k should be tried, which would be
exactly the overfitting this document exists to prevent.

## Additional accept condition for this arm

Because lowering the k-mer admits ~5.5x more pairs, and more pairs is the lever
that produced the chimeric consensi, one condition is added to the four already
listed:

5. On every calibration AND validation genome, purity must not fall by more than
   **1 percentage point** versus the default. Absolute purity >= 95% is not
   enough on its own here: the baseline measures 100.0% on boechera, so a drop
   to 95% would be a real degradation hiding inside a passing threshold.
