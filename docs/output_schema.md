# CARP output-contract schema versions

`carp_manifest.json` (written at the output root by `run_pipeline.py`) declares
the **output contract** a run conforms to: a `schema_version` plus an `outputs`
map of stable logical names → paths relative to the output root. Consumers
(e.g. the JBrowse portal) resolve files by logical name and key their extractor
on `schema_version` instead of sniffing the directory layout.

`schema_version` is bumped **only on a breaking change** to the layout or
contents of consumed outputs — **not** on every CARP release. A
renamed/removed/moved consumed output is a logical-name change **plus** a
`schema_version` bump, recorded here.

The manifest is implemented in `scripts/manifest.py`; the authoritative
`outputs` map lives there (`OUTPUTS`) and `OUTPUT_SCHEMA_VERSION` is the current
value.

## v3 — rDNA nested under a single parent (current)

**Breaking change to consumed outputs: the rDNA classification was restructured.**
The flat top-level classes `rDNA_45S` / `rDNA_5S` are now nested under a single
`rDNA` parent:

| Old (≤ v2) | New (v3) |
|---|---|
| `rDNA_45S` | `rDNA/45S_rDNA` |
| `rDNA_45S/18S`, `/25S`, `/5.8S`, `/IGS`, `/ITS1`, `/ITS2` | `rDNA/45S_rDNA/18S`, … |
| `rDNA_5S`, `rDNA_5S/5S` | `rDNA/5S_rDNA`, `rDNA/5S_rDNA/5S` |

This renames the **per-class rDNA outputs** (slash → dot in filenames):

- split GFF3s: `Repeat_Annotation_NoSat_split_by_class_gff3/rDNA_45S.gff3` →
  `rDNA.45S_rDNA.gff3`, `rDNA_45S.18S.gff3` → `rDNA.45S_rDNA.18S.gff3`, … and
  `rDNA_5S.5S.gff3` → `rDNA.5S_rDNA.5S.gff3`.
- per-class BigWigs: `Repeat_density_by_class_bigwig/{10k,100k}/rDNA_45S_10k.bw`
  → `rDNA.45S_rDNA_10k.bw`, etc.

**Unchanged:** the `rDNA` rollup outputs — `rDNA_RepeatMasker.gff3` (top-level
feature track), the split `rDNA.gff3`, and the `rDNA_{10k,100k}.bw` rollup
BigWig — all key on the `^rDNA` prefix, which still matches. The three REQUIRED
outputs (`genome_cleaned.fasta`, `Repeat_Annotation_Unified.gff3`,
`summary_statistics.csv`) and the marker dir `Repeat_density_by_class_bigwig/`
are **not renamed** (their *contents* change: the `classification` values and the
`summary_statistics.csv` row labels carry the new paths).

**No logical name in `OUTPUTS` points at an rDNA leaf file**, so manifest-aware
consumers are unaffected; consumers that glob `rDNA_45S*` / `rDNA_5S*` filenames
must update. **No backward-compat symlinks are emitted** — the schema bump is the
migration signal.

## v2 — CARP 0.9.0rc6 and rc7

The output **layout** contract is unchanged across rc6 and rc7; rc7 only *adds*
`carp_manifest.json` itself (additive, no layout change — see the note at the
bottom), so it stays `schema_version: "2"`.

Density BigWigs reworked and renamed, all sourced from
`Repeat_Annotation_Unified.gff3` and written sparse (run-length-merged):

- whole-genome total: `Repeat_density/Repeat_density_total_{10k,100k}.bw`
- per class / superfamily: `Repeat_density_by_class_bigwig/{10k,100k}/`
  (includes the `Class_II.Subclass_1.TIR` rollup — FR-2a)
- tandem per-family, two flavours (FR-2b):
  `Tandem_repeats_unified_split_by_family_bigwig/{10k,100k}/TRC_<n>_*.bw` and
  `Tandem_repeats_TideCluster_split_by_family_bigwig/{10k,100k}/TRC_<n>_*.bw`
- structural TideCluster totals: `Tandem_repeats_TideCluster_{10k,100k}.bw`

Replaces the pre-rc6 layout (see v1). GFF3 paths and the
`Repeat_Annotation_NoSat_split_by_class_gff3/` directory name are unchanged for
backward compatibility.

**Density value semantics (contract).** Every `*.bw` value is the *fraction of
the window covered* by that track's features, in `[0, 1]` (mean per-base
coverage per bin, 10-bin moving-average smoothed). Overlapping features are
merged into a strand-agnostic union before coverage, so a value never exceeds 1
even though the Unified annotation tolerates overlap (L1 `Simple_repeat`/
`Low_complexity` over a TE; nested L2 children). *(Fixed post-1.0.0rc2 — values
previously stacked above 1 in dense tandem/centromeric regions; additive, no
layout change, still `schema_version: "2"`.)*

**By-class set = partition + roll-ups (consumers: do not double-count).** The
`Repeat_density_by_class_bigwig/` exact-classification tracks are a **disjoint
partition** (each feature counted once at its leaf classification) and sum to
the total track. Alongside them are **cumulative roll-up** tracks that overlap
the partition: `All_Ty1_Copia`, `All_Ty3_Gypsy`, `Class_II.Subclass_1.TIR`, and
`Mobile_elements`. Use the partition for a summed/stacked view; use the roll-ups
as ready-made aggregates — never sum a roll-up together with its leaf tracks.

## v1 — pre-rc6 (≤ 0.9.0rc5), historical

- RM-only whole-genome total `RepeatMasker/Repeat_Annotation_NoSat_{10k,100k}.bw`
  (excluded tandems)
- per class: `Repeat_Annotation_NoSat_split_by_class_bigwig/`
- TideCluster per-cluster:
  `TideCluster/default/TideCluster_clustering_split_files_bigwig/`
- TideCluster totals: `TideCluster/default/TideCluster_clustering_{10k,100k}.bw`

No `carp_manifest.json` was emitted in v1.

---

**Note:** the manifest itself was introduced in CARP **0.9.0rc7**. Its presence
is additive (a new file at the output root) and did **not** bump the schema —
the rc6 layout it describes is unchanged. The first manifest therefore reports
`schema_version: "2"`.

**Note:** the tandem LTR-RT feature (post-1.0.0rc2) adds one new declared output,
`DANTE_LTR_tandems.gff3` (logical name `dante_ltr_tandems_gff3`) — a top-level
file of LTR_RT_TR container features, one per head-to-tail shared-LTR array, each
listing its member element IDs. **Always present; header-only when no tandems are
detected**, so consumers can open it unconditionally. Additive (a new logical
name; no existing output renamed/moved), so it does **not** bump the schema —
still `schema_version: "2"`.
