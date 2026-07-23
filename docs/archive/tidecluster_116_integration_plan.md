# TideCluster 1.16.0 integration + unified-annotation overlap fix — implementation & test plan

Scope: the three items in [`code_updates_and_bugs.md`](code_updates_and_bugs.md).
Unlike the Tier A/B performance work, **these are intended functional changes** —
outputs *will* change (new rDNA labels, non-overlapping satellite TRCs, a new
`TE_origin` attribute, faster `tc_reannotate` with ≤0.15 % masked-bp drift at
chunk cuts). Byte-identity is therefore **not** the bar here; correctness on the
curated test genome is.

## Test genome is curated ground truth

`tmp/genome_cleaned.fasta.gz` (5 seqs, `OZ408683.1`–`OZ408687.1`) is not arbitrary:

- `OZ408684.1` (27.8 Mb) is TideCluster 1.16.0's **own rDNA + overlap-resolution
  calibration fixture** — a 5.4 Mb 45S array that fragments into TRC_1/3/4/5, plus
  single 5S genes. Expected: those four TRCs labelled `rDNA_type=45S` at coverage
  1.0; cross-TRC boundary overlaps (239 kb, 3.9 %) resolved to 0.
- `OZ408687.1:18,715,299..18,883,725` is the **TR-from-TE fixture**: a tandemly
  organised array of `Class_I/LTR/Ty1_copia/Ale` elements that TideCluster calls
  `TRC_13`. Today the unified annotation reports **both** the satellite and the
  individual Ale LTR-RTs over the same bases.

So every change below has a concrete expected result on this genome.

---

## Locked design decisions (from maintainer)

| # | Decision |
|---|----------|
| D1 | **TR-from-TE:** satellite wins, judged **only against DANTE structural annotation (Tier 1)** — never RepeatMasker. Trigger is a length-qualified clustering TRC covering **multiple same-type structural TE elements**. Record a **note** of the TE origin on the satellite; **drop** the underlying structural TEs from the unified file (they remain in `DANTE_*.gff3`). |
| D2 | **rDNA:** label array-level `rDNA_45S` / `rDNA_5S` (no internal 18S/ITS/… substructure). |
| D3 | **Overlap verification:** write a residual-overlap **TSV report, non-fatal** (warn, don't fail). |

"Minimal length already applied" (D1) is structural, not a new parameter: Tier-3
TRCs come from `TideCluster_clustering.gff3`, where TideCluster has already applied
`-m/--min_length` (per array) and `-M/--min_total_length` (**default 50 kb** total
per TRC; pipeline sets `-m 5000` on the short run). Sub-threshold tandems are in the
TideHunter residual tracks (Tier 6) and are excluded from the TE-origin logic by
construction.

---

## Item 1 — Upgrade to TideCluster 1.16.0

1.16.0 is on the `petrnovak` conda channel (verified). The upgrade alone delivers
two of the three doc items with **no Snakefile logic change**, because both new
behaviours are default-on in `run_all` and write into the existing
`TideCluster_clustering.gff3`:

- **`tc_reannotate` parallelism** (the requested optimization; your FR-1+FR-3):
  chunked, pooled single-threaded RepeatMasker + sub-quadratic `filter_intervals`.
  The `tidecluster_reannotate` rule already passes `-c {threads}`; new
  `--chunk_size`/`--overlap` default to 50 Mb/100 kb. No rule change needed for the
  speedup.
- **rDNA identification** (default-on): adds `rDNA_type=45S|5S;rDNA_coverage=<frac>`
  to rDNA TRC rows + a `<prefix>_rdna.tsv`. Uses TideCluster's bundled
  `data/rdna_library.fasta`.
- **Cross-TRC overlap resolution** (default-on): clustering GFF3 made
  non-overlapping *across satellite TRCs* (dominant-TRC-wins).

### Changes

1. `envs/tidecluster_run.yaml`: pin `tidecluster=1.15.2` → `1.16.0`; refresh the
   header comment (it currently asserts 1.15.2 specifics).
2. **(optional, recommended)** Expose the new knobs as pipeline config with
   safe defaults so behaviour is explicit and tunable:
   - `tidecluster_chunk_size` (default 50_000_000) → `tc_reannotate.py --chunk_size`
     in the `tidecluster_reannotate` rule.
   - `tidecluster_detect_rdna` (default `True`) → pass `--no_rdna` when `False`.
   - `tidecluster_rdna_library` (default empty → TideCluster bundled) → `--rdna_library`.
     Note: the pipeline's own richer `data/rdna_library.fasta` is a candidate override.
   - `tidecluster_keep_trc_overlaps` (default `False`) → pass `--keep_overlaps` when `True`.
   These thread into `tidecluster_long` / `tidecluster_short` (`run_all`) and
   `tidecluster_reannotate`.
3. The env-hash change triggers a one-time conda env rebuild on next run.

No `Singularity %files` impact (env yaml already covered).

---

## Item 2 — Surface rDNA in the unified annotation (D2)

Consume the new `rDNA_type` attribute in `scripts/make_unified_annotation.R`,
`load_tier3_tidecluster()` (and, for symmetry, `load_tier4_tc_rm()`):

```r
normalise_tc <- function(gr, tool) {
  if (length(gr) == 0) return(gr)
  gr$type <- "repeat_region"
  rdna <- if (!is.null(gr$rDNA_type)) as.character(gr$rDNA_type) else rep(NA, length(gr))
  cls  <- ifelse(!is.na(rdna) & rdna == "45S", "rDNA_45S",
          ifelse(!is.na(rdna) & rdna == "5S",  "rDNA_5S",
                 paste0("Satellite/TideCluster/", as.character(gr$Name))))
  set_meta(gr, cls, cls, 3L, tool)   # Name == classification, like RM rDNA rows
}
```

- `rDNA_45S` / `rDNA_5S` are already in `classification_vocabulary.yaml` (leaves
  135 & 142) and already match the `^rDNA` branch in `finalise_output()`'s type
  regex → emitted as `repeat_region`.
- `calculate_statistics_and_make_groups.R` routes `^rDNA` Names into `rDNA.gff3`
  (line 85), so a relabelled TRC is counted as rDNA, not generic Satellite — the
  doc's "clearly labelled as rDNA".
- Keep the originating `TRC` id in a side attribute (`TRC=TRC_1`) for traceability.
- **Sub-decision to confirm:** should rDNA arrays *also* count in the
  `Tandem_repeats` aggregation row (doc: "rDNA is both TR but also rDNA")? Default
  proposed: **no** — rDNA is its own category (cleaner totals); revisit
  `aggregation_buckets` if double-counting is wanted.
- rtracklayer note: `rDNA_type`/`TRC` survive only if propagated through
  `subset_seqs()` (which rebuilds mcols from `.META_COLS`). rDNA needs nothing new
  (it rides `classification`); `TRC` would need adding to `.META_COLS` if we keep it.

Reconciliation with RepeatMasker rDNA: RM still matches the subunit-level
`data/rdna_library.fasta` (Tier 5). Inside an rDNA TRC the RM subunit hits are
lower tier and trimmed away → the array shows as one `rDNA_45S`/`rDNA_5S` feature.
RM rDNA outside any TRC is unaffected.

---

## Item 3 — Fix TR-from-structural-TE double annotation (D1) — the core change

### Current behaviour (root cause, traced)

In `process_batch()`:
- Step 1: `level1 <- t1` (structural TE: DANTE_LTR/TIR/LINE).
- Step 3: TideCluster default clusters are added to `level1` **without trimming
  against Tier 1**; only clusters *fully within a single LTR* are demoted to
  Level-2 nested. A tandem **array of Ale** is the opposite topology — many Ale
  elements lie *within the TRC* — so `findOverlaps(type="within")` matches nothing,
  the whole TRC is appended to Level 1, and it coexists with the Ale elements.
- `sanity_check()` excludes `^Satellite` from its TE-overlap test, so the double
  annotation is **silently tolerated** today.

### New behaviour

Add a pre-pass in `process_batch()` (before Step 1 builds `level1`) that classifies
each Tier-3 clustering TRC as **TE-derived** or not, operating on `t3_def` + `t3_sho`
grouped by TRC `Name`, against `t1` (structural only):

A TRC is **TE-derived** when **all** hold:
1. its arrays overlap **≥ 2** distinct Tier-1 structural TE elements;
2. those structural TEs share a classification lineage of depth ≥ 3
   (`lca_classification()` of their classes is at least e.g. `Class_I/LTR/Ty1_copia`);
3. structural TEs cover **≥ `TE_ORIGIN_MIN_COVERAGE`** (proposed 0.5) of the TRC's
   total array bp in the batch.

Thresholds = named constants (also surfaceable as CLI opts):
`TE_ORIGIN_MIN_ELEMENTS=2`, `TE_ORIGIN_MIN_LCA_DEPTH=3`, `TE_ORIGIN_MIN_COVERAGE=0.5`.

For **TE-derived** TRCs:
- the satellite arrays go to Level 1 and **win** the region;
- set `gr$TE_origin <- <LCA class of the covered structural TEs>` (e.g.
  `Class_I/LTR/Ty1_copia/Ale`) and `gr$TE_origin_note <- "satellite over tandem TE array"`;
- **remove** the covered Tier-1 structural TE elements from `level1` (they stay in
  `DANTE_LTR/TIR/LINE.gff3`). "Covered" = structural TE whose reciprocal overlap
  with the TE-derived TRC arrays ≥ 0.5 (avoids deleting a TE that merely clips the
  array edge).

For **non-TE-derived** TRCs: unchanged from today (Step 3 logic preserved).

`TE_origin` is a new metadata column → add to `.META_COLS` and initialise `NA` on
every tier so `c()`/`rbind` across batches stays schema-consistent; rtracklayer
omits per-feature `NA` on export, so only TE-derived satellites carry it.

### Why this respects every constraint

- "must be structural annotation, not RepeatMasker" → trigger reads `t1` only.
- "covering multiple same TE-like elements" → criteria 1+2.
- "minimal length already applied" → Tier-3 input only (50 kb `-M` etc.).
- "add note … but do not keep [the TE] in annotation" → `TE_origin` note +
  structural TEs dropped from unified (still in `DANTE_*`).

---

## Item 3b — Overlap verification report (D3)

Extend `sanity_check()` (or a small post-step) to emit
`Repeat_Annotation_Unified.overlaps.tsv`: every Level-1 vs Level-1 overlap whose
**neither** partner is `Simple_repeat`/`Low_complexity`, with columns
`seqid,start,end,name_a,class_a,tier_a,name_b,class_b,tier_b,overlap_bp`.
Print a one-line WARNING with the count; **do not fail** the run. After the Item-3
fix this file should contain no satellite-vs-structural-TE rows (the bug class);
any residual rows (e.g. same-tool partial TEs) are surfaced for review. Declare it
as an output of `make_unified_annotation`.

---

## Files touched (summary)

| File | Change |
|---|---|
| `envs/tidecluster_run.yaml` | pin 1.16.0; refresh comment |
| `Snakefile` | (opt.) new config + flags on the 3 TideCluster rules; add `overlaps.tsv` output to `make_unified_annotation` |
| `config.yaml` + validation block | (opt.) defaults for the new `tidecluster_*` keys |
| `scripts/make_unified_annotation.R` | rDNA relabel (Item 2); TE-derived-TRC pre-pass + `TE_origin` + structural-TE drop (Item 3); overlaps TSV (Item 3b); `.META_COLS` += `TE_origin`(,`TRC`) |
| `CLAUDE.md` | document new behaviour + config keys |
| `docs/…` | this plan |

No new runtime files → no `Singularity %files` edit (hook stays green).

---

## Testing plan

Baseline-then-change on the curated genome (no prior run exists — fresh baseline).

1. **Baseline (1.15.2 + current code):** run the full pipeline on
   `tmp/genome_cleaned.fasta.gz`; archive `Repeat_Annotation_Unified.gff3`,
   `summary_statistics.csv`, the split GFFs, and a snapshot of the `OZ408687.1`
   TRC_13 window and the `OZ408684.1` rDNA window. Confirm the bug is present
   (satellite + Ale both reported over 18.715–18.884 Mb).
2. **After upgrade + code:**
   - **rDNA (OZ408684.1):** TRC_1/3/4/5 carry `rDNA_type=45S`; unified shows
     `rDNA_45S` features there; `rDNA.gff3` non-empty; 5S genes handled. Compare to
     TideCluster's `_rdna.tsv`.
   - **Overlap resolution:** TideCluster `clustering.gff3` multiply-covered bp = 0
     across TRCs (matches the 239 kb→0 fixture); union bp preserved.
   - **TR-from-TE (OZ408687.1):** the 18.715–18.884 Mb window is a single
     `Satellite/TideCluster/TRC_13` feature carrying
     `TE_origin=Class_I/LTR/Ty1_copia/Ale`; **no** Ale LTR features remain in
     unified there; the Ale elements **are** still in `DANTE_LTR.gff3`.
   - **Verification TSV:** zero satellite-vs-structural-TE rows genome-wide.
   - **`tc_reannotate` speedup:** benchmark `mean_load` ≫ 1 (was ~1.1/50 on Pearl);
     masked-bp within ±0.15 % of a `--keep`-style whole-genome reference on one seq.
   - **No regressions:** pipeline reaches all leaf outputs (report HTMLs, BigWigs,
     summary) exit 0; classification validation passes; `summary_statistics.csv`
     biologically sane (Ale array now under Tandem/rDNA, not double-counted).
3. **Targeted unit checks** for `make_unified_annotation.R`: a synthetic batch with
   (a) a TRC over 3 same-class LTRs → TE-derived, TEs dropped, `TE_origin` set;
   (b) a TRC clipping one LTR edge → **not** TE-derived, both kept;
   (c) a TRC over mixed-class TEs (LCA too shallow) → not TE-derived;
   (d) an rDNA TRC → `rDNA_45S`, never TE-derived.

## Rollout / commit order

1. env pin + comment (Item 1) — isolated, rebuilds env.
2. rDNA relabel (Item 2) — small, independently testable.
3. TE-derived pre-pass + `TE_origin` + TE drop (Item 3).
4. overlaps TSV (Item 3b).
5. optional config knobs + `config.yaml` + validation + `CLAUDE.md`.
6. full validation run on the curated genome; archive before/after.

Per-item commits on `main` (sole-developer policy), after the validation run passes.

## Open sub-decisions to confirm before coding

- **S1** `TE_origin` thresholds (`MIN_ELEMENTS=2`, `MIN_LCA_DEPTH=3`,
  `MIN_COVERAGE=0.5`) — accept or tune?
- **S2** Attribute naming: `TE_origin=<class>` (+ `TE_origin_note=…`). OK, or prefer
  e.g. `derived_from` / `Note=`?
- **S3** rDNA double-counting in the `Tandem_repeats` aggregation row — keep
  separate (proposed) or count in both?
- **S4** Expose the new TideCluster knobs as config (recommended) or hard-wire
  defaults?
