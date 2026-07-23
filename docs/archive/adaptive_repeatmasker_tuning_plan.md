# Adaptive RepeatMasker tuning (culling / rush) — design notes

**Status:** EXPLORATORY — *not implemented*. This document captures the design
discussion for auto-selecting `repeatmasker_culling_limit` (and, separately,
`repeatmasker_sensitivity: rush`) from a measurable property of the library, so
the reasoning is on record when we build it.
**Decision (2026-07-01):** the preferred path is the **realized probe**
(§2d) because it measures *both* the speedup and its masking cost on the actual
genome+library. Implementation is **deferred** until several genomes have been
annotated with `culling ∈ {0, 2}` (and rush on/off), so we can calibrate the
switch on real speedup-vs-cost data before wiring anything.
**Prereqs already shipped (1.0.0):** the two manual knobs this would automate —
`repeatmasker_culling_limit` (default 0) and
`tidecluster_reannotate_culling_limit` (default 0). See `docs/configuration.md`
and the `## culling` discussion in `CLAUDE.md`.

---

## 0. Goal

Turn culling (and possibly rush) on **only when it pays off**, automatically, so
a user never has to reason about library redundancy. The intuition to encode:
on a small genome with a non-redundant library that isn't the runtime
bottleneck, culling should stay **off**; on a redundant de-novo library it
should be **2**. That requires a *redundancy criterion* we can compute and
threshold.

## 1. What "redundancy" must mean here (tie it to the mechanism)

Culling does **not** exploit redundancy in the abstract — it collapses **HSP
multiplicity per genomic locus**. A locus explodes into hundreds of HSPs when
many *different* library consensi share the same **seed words** with that locus,
so `rmblastn` seeds and extends each one. `-culling_limit 2` keeps ~2 HSPs per
locus instead of *F*, so the wall-time win is roughly **~F/2** (matches the
measured ~3× at *F*≈6; cf. `CLAUDE.md`).

So the quantity to estimate is:

> **F = expected number of distinct library sequences a random genomic seed
> matches** (≈ expected HSP multiplicity per seed).

Two distinctions follow, which most off-the-shelf "library redundancy" numbers
get wrong:

- **Removable vs. residual redundancy.** Clustering / CD-HIT / the pipeline's
  own `reduce_library` measure redundancy you can *delete* (near-identical
  full-length duplicates). Culling handles the redundancy that *survives*
  reduction — diverged sibling consensi (~55–80 % local identity) a reducer
  correctly keeps but that still co-seed a locus. A library can cluster tightly
  yet still have high *F*. **The metric must count similarity at the seed/local
  level, not the mergeable-full-length level.**
- **Intrinsic vs. realized.** A library can be hugely redundant in a family that
  is *absent* from this genome → zero HSP explosion → culling buys nothing.
  Intrinsic (library-only) metrics are cheap and genome-agnostic but can
  mis-predict; realized (library-vs-this-genome) metrics are accurate but cost a
  probe run.

## 2. Candidate metrics (cheapest → most faithful)

**(a) Seed fan-out, k-mer based — intrinsic, cheap, most mechanistically honest.**
Build `kmer → {distinct library seq IDs}` at a word size comparable to
`rmblastn`'s seed (~11–14), then
`F = mean over kmers (weighted by occurrence) of |sequences containing it|`.
Direct estimator of HSP multiplicity; predicted speedup ≈ *F*/limit falls
straight out. Make it cheap with **minimizer sampling** (every *m*-th k-mer) —
the average is preserved, memory drops. Lead candidate for an intrinsic switch.

**(b) Reduction compression ratio — free, already computed.**
`reduce_library` (CAP3/mmseqs per class), `reduce_library_containment` (blastn
longest-first) and `reduce_dimer_library` already self-compare the library;
`bp_in/bp_out` and `N_in/N_out` are a redundancy signal we already pay for.
Caveat: it measures *removable* redundancy → a lower-bound proxy for *F*.
Correlated but not the same axis. Useful as a free secondary feature and doubles
as a **QC metric for the reduction step itself**.

**(c) Cheap sketches — gzip ratio, distinct-kmers/bp, seqs/Mbp.** Trivial, but
conflate within-sequence low-complexity with inter-sequence fan-out. Only useful
as extra regression features, not the decision variable.

**(d) Realized HSP-per-bp probe — genome-aware, most robust. ← chosen direction.**
Because the pipeline already **chunks** the genome, run 1–2 chunks with
`culling=0`, count total HSPs and masked bp: `HSPs / masked_bp` *is* realized
*F*, measured on the actual genome+library. Optionally probe one chunk at both
`culling ∈ {0, 2}` to also get the **masking-delta** for free. Decide culling for
the remaining chunks from it. Immune to the absent-family bias of (a)/(b). A bit
more plumbing (probe → decide → run rest), but the chunking makes it cheap.

## 3. Decision policy — "no-regret", not a plain on/off

Key point: **culling's only downside is a small masking perturbation** (−0.7 %
main RM at limit 2; the real hazard is megasatellite under-masking on the
`tidecluster_reannotate` path, and limit 1 overshooting at −3 %). When *F*≈1 it
also has ~no *benefit*. So the principled rule is not "off on small genomes" but:

> **Enable culling iff the predicted speedup is large enough to *pay for* the
> masking perturbation.** Below threshold you accept a masking change you get
> nothing for; above it you buy real time.

This makes *F* (or realized HSP/bp) the single switch — e.g. `culling=2` when
*F* ≳ 3–4 else 0 (**threshold to be calibrated, not asserted**). It also
satisfies the "small non-redundant genome" intuition automatically: those have
*F*≈1 → auto-off, with no separate genome-size rule.

**Keep `rush` OUT of this criterion.** Rush (`-qq`) only enlarges `word_size`;
its benefit tracks the *fraction of spurious seeds*, which is *lower* in
repeat-dense genomes (seeds are mostly real) — so its benefit is roughly
orthogonal (even anti-correlated) to library redundancy, and it costs real
sensitivity. Better model: an **escalation ladder** by predicted absolute cost
(∝ `genome_bp × library_bp × F`): culling first (near-lossless), then rush only
under an explicit wall-time budget / user opt-in. The redundancy metric must not
flip rush.

## 4. Calibration experiment (do this before wiring anything)

When reduced libraries from distinct sources are available, hold the genome (or
a fixed multi-chunk sample) constant and, per library, record:

- **features:** *F* (k-mer fan-out), reduction ratio, gzip ratio, seqs/Mbp,
  class composition (LTR-dispersed fraction);
- **outcomes:** actual RM wall-time at `culling ∈ {0, 2}` (and optionally rush),
  and masked-bp delta.

Then scatter **speedup vs F** and **masked-bp-delta vs F**. Target: the *F* where
speedup crosses ~1.5× while masked-bp-delta stays inside the ±0.15 % bar. If one
feature (probably *F*) separates cleanly → ship a threshold; else a 2-feature
logistic (`F × LTR-fraction`) is still fully explainable. **Diversity of sources
matters more than count** — include one megasatellite-heavy genome so the
under-masking hazard shows up on the masked-bp axis.

Since the chosen direction is (d), the same probe that would drive the switch at
runtime is also the calibration instrument: annotate the candidate genomes with
both options, log `HSPs`, `masked_bp`, and wall-time per run, and fit the
speedup/cost curves from those.

## 5. Open questions / caveats to resolve at implementation time

- **Where to count HSPs** for the probe — from the `rmblastn` `.out`/tabular per
  chunk, before ProcessRepeats, so the count reflects pre-filter multiplicity
  (the thing culling attacks), not post-filter survivors.
- **`tidecluster_reannotate` is a separate switch.** Its culling target
  (phase-redundant tandem HSPs) is even stronger, but it carries the
  megasatellite under-masking risk. Do **not** couple the two auto-decisions;
  calibrate `tidecluster_reannotate_culling_limit` on a satellite-heavy genome
  independently.
- **Probe cost accounting.** The probe run must be counted against the
  optimization's own budget — if RM is already cheap (small genome), the probe
  itself may dominate; a cheap intrinsic pre-filter (§2a/§2b) could gate whether
  the probe is even worth running.
- **Determinism / reproducibility.** Any auto-decision must be logged in
  `record_provenance.py` output (the *chosen* culling value, the metric value,
  and the threshold) so a run is reproducible and the decision is auditable.
