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

## v2 — CARP 0.9.0rc6 and rc7 (current)

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
