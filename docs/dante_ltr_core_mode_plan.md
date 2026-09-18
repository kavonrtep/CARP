# Implementation plan — optional DANTE_LTR core-domain mode

Adds `dante_ltr_mode: core` as an **opt-in** alternative to DANTE_LTR's default
lineage-based element detection, for genomes that REXdb covers poorly. Upstream
work is specified separately in
[`dante_ltr_core_library_policy_request.md`](dante_ltr_core_library_policy_request.md).

**Scope note (read first).** This plan spans two repositories. CARP items are
`C*`; the `dante_ltr` items they depend on are `U*` and live in the request
document above. C3 cannot land before U1 ships.

**Invariants to preserve.**

* (a) **The canonical run is byte-identical.** `dante_ltr_mode: lineage` is the
  default; every behaviour below is reachable only from `core`. The one
  unavoidable change to the default path is the `dante_ltr` version bump (C1),
  which must be verified as output-neutral rather than assumed.
* (b) No manifest or `OUTPUT_SCHEMA_VERSION` change — see "Manifest" below.
* (c) No `classification_vocabulary.yaml` change — all 25 classifications core
  mode emits already validate (verified with `classification.py`).
* (d) New order-producing code carries a `PYTHONHASHSEED`-invariance test, per
  the determinism rules in `CLAUDE.md`.

---

## Why

DANTE_LTR's default mode requires the protein domains of an element to agree on
a REXdb *lineage*. On a genome far from REXdb they disagree, and the element is
discarded. Core mode instead seeds on the ordered RT/RH/INT triplet — the order
alone separates the superfamilies (`RT RH INT` = Ty3/gypsy, `INT RT RH` =
Ty1/copia) — and assigns the classification afterwards as the LCA of the
contained domains' calls, clipped at that superfamily. Elements therefore keep a
call even when the lineage is unresolvable, at a coarser depth.

Measured on a *Draparnaldia* assembly (397 Mb), DANTE_LTR 0.6.0.0:

| | lineage mode (today) | core mode |
|---|---|---|
| complete elements | **0** | **1,408** (13.17 Mb = 3.3% of assembly) |
| LTR library CARP builds | 13 seqs / 13.3 kb | 397 seqs / 436 kb |
| labels below lineage level | — | 1,320 of 1,530 TE records (86%) |

On a REXdb-covered genome (the DANTE_LTR sample assembly) core mode demotes
**nothing** — all 14 classification values stay at lineage level — so the coarse
labels appear only where REXdb is genuinely thin. Upstream measures core mode at
96% recall of lineage-mode elements on a covered genome, so it is *not* a
superset and must not become the default.

---

## Manifest: no action

`manifest.py::OUTPUTS` contains no per-class file entries — only fixed top-level
files plus four directory entries. The new per-class outputs land inside
`Repeat_density_by_class_bigwig/` (a directory entry, so
`assert_run_determinism.py` `rglob`s them automatically) and
`Repeat_Annotation_NoSat_split_by_class_gff3/`, which is not in `OUTPUTS` at all.
`OUTPUT_SCHEMA_VERSION` bumps only on a renamed or removed output; core mode does
neither, so it stays `"3"`.

These file names already occur today — `Class_I.LTR.Ty1_copia.gff3`,
`Class_I.LTR.Ty3_gypsy.gff3`, `Class_I.LTR.Ty3_gypsy.chromovirus.gff3` and
`...non-chromovirus.gff3` are present in existing runs, produced by DANTE tier-2
domain calls and RepeatMasker, which have always been able to stop at an internal
node. They are merely tiny: 0.0007%–0.13% of the genome across `tmp/gca_ath`,
`tmp/ath_1162`, `tmp/gca_dun`. **Both report defects in C5 and C7 are therefore
pre-existing latent bugs**; core mode takes them from negligible to dominant.

---

## C1 — Bump DANTE_LTR

`envs/tidecluster.yaml`: `dante_ltr=0.5.4.0` → **0.6.2.0** (verified on the
`petrnovak` channel 2026-09-17; the env dry-run-solves clean on the same r-base
4.1 stack). 0.6.0.0 provides `--mode core`, 0.6.1.0 the library policy (U1),
0.6.2.0 a clustering-determinism fix.

**Verify, don't assume:** run a fixture before and after and diff, since this
touches the default path. Upstream states lineage-mode output is unchanged, and
the detection executable is byte-identical from 0.6.0.0 to 0.6.2.0 (`utils/`
differs only by the new `library_policy.R`) — but 0.6.2.0 passes
`--spaced-kmer-mode 0` to `mmseqs easy-cluster`, which **changes the library
content** (~1.5% more representatives upstream; −0.4% on a CARP library here).
So the LTR library, and the RepeatMasker annotation downstream of it, WILL move
on this bump even in lineage mode. Expect a diff; confirm it is this and not
something else.

**Acceptance test for C3 must assert deltas, not absolutes.** The 397 → 459
figures in the request were measured under 0.6.0.0, before that flag. Run both
`strict` and `nested` under 0.6.2.0 and check the shape: ~+62 clusters, 5
conflicts still dropped, 69 promotions splitting 50 `chromovirus` / 19
`Chlamyvir`.

## C2 — Config knob

* `config.yaml`: `dante_ltr_mode: lineage`.
* `Snakefile`: validate against `{lineage, core}` at load, alongside the existing
  validation block near `Snakefile:266`.
* `dante_ltr` rule: pass `--mode {params.mode}` **unconditionally**. Upstream's
  default is `lineage`, so this is output-neutral and makes the log
  self-documenting.

## C3 — Library policy pass-through *(depends on U1)*

`make_library_of_ltrs`: pass `--annotation_conflict nested` to
`dante_ltr_to_library` only when `dante_ltr_mode == "core"`. One `params` entry
and one shell conditional.

**Test:** run the library step in lineage mode and assert the output matches the
`strict` result, proving this path is unreachable by default.

## C4 — Provenance *(blocks C5 and C6)*

`scripts/record_provenance.py::_filter_config` is an explicit 17-key allowlist;
`dante_ltr_mode` must be added or the report cannot see the mode.
`run_provenance.json` is excluded from the determinism gate (`_EXCLUDE_KEYS`), so
adding a field is free.

## C5 — Report: density panel, core mode only

`discover_lineage_bw_files()` (`scripts/make_repeat_report.R:582`) filters
discovered BigWigs against a hardcoded 31-name lineage list, labelling each by the
last dot-component of the filename. Simulated against the filenames core mode
produces:

```
Class_I.LTR.Ty1_copia_100k.bw                  -> Ty1_copia    DROPPED
Class_I.LTR.Ty1_copia.Ale_100k.bw              -> Ale          shown
Class_I.LTR.Ty3_gypsy_100k.bw                  -> Ty3_gypsy    DROPPED
Class_I.LTR.Ty3_gypsy.chromovirus_100k.bw      -> chromovirus  DROPPED
Class_I.LTR.Ty3_gypsy.chromovirus.Chlamyvir…bw -> Chlamyvir    shown
```

On the *Draparnaldia* run those three dropped tracks carry 1,320 of 1,530 elements
(86%), so the panel would render two lineages while the composition table above it
reports 3.3% LTR content — a contradiction, which is worse than an omission.

**Change.** Keep the hardcoded list as the lineage-mode behaviour. In core mode,
widen the accepted set to any node under `Class_I/LTR` in
`classification_vocabulary.yaml`, ordered by vocabulary order, and render internal
nodes as `chromovirus (unspecified)` so they do not read as lineages. Keep the
bare `Ty3_gypsy` row visually distinct from the `All_Ty3_Gypsy` roll-up, which
means something different.

**Mode source.** `load_provenance(outdir)$config$dante_ltr_mode`. **When
provenance is missing or the key is absent, fall back to lineage behaviour** — a
bare `snakemake` run writes no provenance, and that exact gap caused the 1.1.0
crash (memory `feedback_local_validation_provenance_gap`). Any test exercising
core-mode reporting must therefore go through `run_pipeline.py`, never bare
snakemake.

## C6 — Report: state the mode

Show `DANTE_LTR mode: core` in the report header, from the same provenance read.
Render nothing in lineage mode, so the canonical report is unchanged.

## C7 — Composition table: element count on internal nodes

`build_comp_tree()` attaches `dante_count` only to `row_type == "leaf"`; the
`total` and `unspecified` rows both pass `NA`. On the *Draparnaldia* run that
hides 1,320 complete elements, though their bp is correct. It is also
data-dependent: if no child of `chromovirus` happens to be annotated, the node
becomes a leaf and the count *does* appear — the same biology renders two ways.

**Change.** Carry the node's own `dante_count` onto the `unspecified` row,
parallel to how `own_bp` is already handled; optionally a subtree sum on `total`.
`html_comp_table()` already draws the cell.

**Not mode-gated** (unlike C5): it only adds a number to a row that already
exists, is correct in both modes, and moves the canonical report by <0.13% of the
genome. Flip to a mode conditional if that is preferred.

## C8 — Docs

`docs/configuration.md` + README table + `config*.yaml` via the `config-docs`
agent; `CHANGELOG.md` via the `changelog` agent. `tests/test_config_docs.py`
blocks the release otherwise.

---

## Order

C7 (standalone, testable against runs already in `tmp/`) → C1 + C2 → **U1** →
C3 → C4 + C5 + C6 → C8.

## Follow-up to raise upstream (not blocking)

`nested` is **not** a superset of `strict`. `utils/library_policy.R` keeps a
cluster only when `prop > threshold || chain`, so a cluster `strict` keeps via
`LCA == majority` that is not a single chain — `{Ty3/gypsy, chromovirus/Tekay,
chromovirus/Reina}` with `Ty3/gypsy` the majority — is dropped under `nested`.
Upstream documents this deliberately. It does not occur on the calibration
genome, so no published number moves, but on another genome switching policy can
silently remove library sequences as well as add them. The fix looks like one
condition: keep on `(majority || LCA == majority || chain)`, promote only on
`chain` — no loss, identical promotion safety.

## Cross-cutting guardrails

1. `dante_ltr_mode: lineage` byte-identical to today apart from C1 — verified with
   the `check-determinism` skill or a stored-run diff.
2. C3's test proves the `nested` path is unreachable by default.
3. **No core-mode CI fixture.** Core mode requires all three core domains, so the
   small fixtures would likely find zero elements and the test would assert
   nothing. Unit-test at the function level upstream (U3) instead.
4. Report changes cannot trip the determinism gate — `report_*` keys are excluded
   via `_EXCLUDE_KEY_PREFIXES`.
