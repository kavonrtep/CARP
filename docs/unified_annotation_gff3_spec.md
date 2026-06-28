# `Repeat_Annotation_Unified.gff3` — output specification

The authoritative, non-overlapping (tier-prioritised) repeat annotation produced
by `make_unified_annotation.R`. This document is the **contract** for the file's
shape and attribute values; `scripts/validate_unified_gff3.py` is its executable
mirror and `tests/test_unified_gff3_spec.py` guards against drift. Change all
three together.

Downstream consumers (this pipeline's own `split_gff_by_name.R`,
`calculate_statistics_and_make_groups.R`, the density/report rules, **and
external apps**) depend on these values — in particular the **satellite `Name`
invariant** below, whose breakage silently emptied the per-family BigWig outputs.

## File-level

- Standard GFF3. First line `##gff-version 3`; `rtracklayer` adds
  `##source-version` / `##date`; provenance lines `##pipeline-version`,
  `##git-sha`, `##run-started` may follow when `run_provenance.json` is present.
- One line per feature; tab-separated 9 columns. No sequence-region directives.
- Features are **Level 1** (top-level, non-overlapping except Simple/Low) or
  **Level 2** (nested children carrying a `Parent`).

## Columns 1–8

| Col | Field | Allowed values |
|-----|-------|----------------|
| 1 | seqid | genome sequence name (non-empty) |
| 2 | source | one of the nine `source_tool` values (table below); **equals the `source_tool` attribute** |
| 3 | type | `transposable_element` or `repeat_region` (rule below) |
| 4 | start | integer ≥ 1 |
| 5 | end | integer ≥ start |
| 6 | score | `.` |
| 7 | strand | `+`, `-`, or `.` |
| 8 | phase | `.` |

## Column 9 — attributes

| Attribute | Presence | Allowed values / format |
|-----------|----------|-------------------------|
| `ID` | **required** | `UA_L1_<8 digits>` (Level 1) or `UA_L2_<8 digits>` (Level 2); unique |
| `Name` | **required** | see *Name conventions* — for TideCluster satellites it is **always** `TRC_<n>` |
| `classification` | **required** | slash-separated path; semantic label (drives `type`) |
| `source_tier` | **required** | integer `1`–`6`, fixed per `source_tool` (table) |
| `source_tool` | **required** | one of the nine tools (table); equals column 2 |
| `element_type` | DANTE_LTR individual elements (not tandem containers) | `complete` or `partial` |
| `structure` | tandem LTR-RT **container** only (DANTE_LTR) | `LTR_RT_TR` |
| `copy_number` | with `structure=LTR_RT_TR` only | integer (member copies in the array) |
| `TE_origin` | optional; only `TideCluster_default`/`_short` | a slash classification path, e.g. `Class_I/LTR/Ty1_copia/Ale` |
| `TE_origin_structure` | optional; with `TE_origin` only | `tandem_LTR_RT` — the TE-derived satellite is a tandem of *full* LTR-RTs |
| `Parent` | iff Level 2 (`UA_L2_…`) | an existing `UA_L1_<8 digits>` ID |

A Level-1 feature (`UA_L1_…`) **must not** carry a `Parent`; a Level-2 feature
(`UA_L2_…`) **must**, and it must resolve to a Level-1 `ID` in the same file.

## `source_tool` ↔ `source_tier` (1:1) and provenance

| source_tool | source_tier | meaning |
|-------------|:-----------:|---------|
| `DANTE_LTR` | 1 | intact LTR-RT (structural) — carries `element_type` |
| `DANTE_TIR` | 1 | TIR DNA transposon (structural) |
| `DANTE_LINE` | 1 | LINE element (structural) |
| `DANTE` | 2 | DANTE protein domain (structural, trimmed under tier 1) |
| `TideCluster_default` | 3 | tandem-repeat cluster (default run) |
| `TideCluster_short` | 3 | tandem-repeat cluster (short-monomer run) — may be absent |
| `TideCluster_RM` | 4 | RepeatMasker-on-TideCluster-library re-annotation |
| `RepeatMasker` | 5 | similarity annotation (TEs, Simple/Low, rDNA subunits, …) |
| `TideHunter` | 6 | sub-threshold tandem residuals (`Satellite/Unknown`) |

## `type` ↔ `classification` rule

`type = repeat_region` **iff** `classification` begins with `Satellite`,
`Simple_repeat`, `Low_complexity`, `rDNA`, or `Unknown`; otherwise
`type = transposable_element`. (Applied uniformly in `finalise_output`.)

## Name conventions

| Feature kind | `Name` | `classification` |
|--------------|--------|------------------|
| DANTE_LTR / TIR / LINE / DANTE (tier 1–2) | = classification path | canonical `Class_I…` / `Class_II…` |
| RepeatMasker (tier 5) | = classification (RM name) | e.g. `Class_…`, `Simple_repeat(AT)n`, `Low_complexity`, `rDNA_45S/18S`, `Satellite/…`, `Unknown` |
| **TideCluster_default / _short / _RM (tier 3–4)** | **`TRC_<n>` (bare)** | `Satellite/TideCluster/TRC_<n>`, or `rDNA_45S` / `rDNA_5S` for rDNA arrays |
| TideHunter (tier 6) | `Satellite/Unknown` | `Satellite/Unknown` |

**Invariant (do not break):** every `TideCluster_*` satellite — *including* rDNA
arrays and TE-derived ones — keeps `Name = TRC_<n>`. The rDNA distinction lives
in `classification` (`rDNA_45S|5S`), not in `Name`; the TE-derived distinction
lives in the `TE_origin` attribute. `split_gff_by_name.R --name-prefix TRC_` and
external consumers select families by this `Name`. `calculate_statistics_…R`
routes rDNA by `classification`, so `Name` can stay `TRC_<n>` on disk.

## Special semantics

- **rDNA arrays** (`TideCluster_*`, `classification=rDNA_45S|5S`) are both a
  tandem family (they appear in the per-family split under their `TRC_<n>`
  Name) **and** rDNA (counted in the rDNA class via classification). No internal
  18S/ITS/5.8S/IGS/25S substructure is written here; that detail is only on the
  RepeatMasker `rDNA_45S/<subunit>` features. The rDNA call is taken from
  TideCluster's **authoritative** per-TRC table `<prefix>_rdna.tsv` (columns
  `TRC`/`rDNA_type`/`coverage`, written by `identify_rdna`), keyed by TRC id and
  applied to **all** TideCluster tiers — tier-3 clustering **and tier-4
  RM-on-TideCluster** (`TideCluster_RM`), which reannotates with the *default*
  dimer library and so shares the default run's TRC namespace. This is why a
  tier-4 array that the structural clustering did not cover is still labelled
  `rDNA_45S|5S` rather than `Satellite/TideCluster/<TRC>`. The clustering GFF3's
  per-feature `rDNA_type` attribute is a fallback used only when the TSV is
  unavailable (`--no_rdna`, detection failure, or older TideCluster); the
  RM-on-TideCluster GFF3 never carries `rDNA_type`.
- **TE-derived satellites** carry `TE_origin=<LCA class of the covered
  structural TEs>`; the underlying tier-1 TEs (and any tandem members) are absent
  from this file (they remain in `DANTE_*.gff3`). Two flavours occur:
  - **Full LTR-RTs in tandem** — the satellite coincides with a tandem LTR-RT
    array of *complete* elements (it overlaps an `LTR_RT_TR` container). Tagged
    `TE_origin_structure=tandem_LTR_RT`. This is the satellite-vs-tandem
    *conflict*: the satellite wins (one annotation for the region), the container
    and its member copies are dropped, and `TE_origin`/`TE_origin_structure`
    record the structural origin.
  - **Tandemised TE-derived sequence (incomplete)** — the satellite is built from
    degraded/partial TE monomers with no complete elements. Detected only when it
    still overlaps ≥2 complete same-lineage structural TEs (`TE_origin` set, no
    `TE_origin_structure`); a purely degraded array with no complete elements is
    *not* flagged TE-derived here and stays a plain `Satellite` (RepeatMasker
    similarity is deliberately not used as the trigger — structural only).
- **Tandem LTR-RT (`LTR_RT_TR`).** Head-to-tail, same-lineage LTR-RT arrays that
  share boundary LTRs (`scripts/resolve_ltr_tandems.py`, upstream of the unified
  annotation) are one **Level-1 container** (`structure=LTR_RT_TR`,
  `copy_number=N`, `classification` = the shared lineage, no `element_type`) with
  the member copies as **Level-2 children** (`Parent=<container>`, normal
  `element_type`). The container is counted once; the members (which legitimately
  overlap each other at their shared LTRs) are nested detail. All other
  overlapping tier-1 structural pairs (cross-lineage, cross-tool, partial) are
  resolved by trimming the shorter and keeping the longest complete, so **no two
  tier-1 features overlap** in this file.

## Validation

```bash
scripts/validate_unified_gff3.py <output_dir>/Repeat_Annotation_Unified.gff3   # exit 1 on any violation
python tests/test_unified_gff3_spec.py                                         # synthetic good + drift cases
```

The unit test runs in the cheap CI tier (`version`/`unit` jobs) and in the
release gate, so a schema change that isn't reflected here fails before a build.
