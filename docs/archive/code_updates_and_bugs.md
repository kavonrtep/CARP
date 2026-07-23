# Status

All items below are **RESOLVED** in the 1.0.0 release line. This file is kept
as the original request record; see the README changelog (1.0.0rc1–rc3) and
`CLAUDE.md` for the implemented behaviour.

# incorporate TideCluster 1.16.0  — ✅ DONE (1.0.0rc1)
- this version includes optimization - requested here
  → chunked/pooled `tc_reannotate` RepeatMasker (`tidecluster_chunk_size`, default 50 Mb)
- rDNA detection - this should be incorporated here to our pipeline to improve unified annotation - rDNA is both TR but also rDNA and again this should be clerly labeled as rDNA
  → TideCluster `rDNA_type=45S|5S` → array-level `classification=rDNA_45S`/`rDNA_5S` in the unified
    annotation (counted as rDNA via classification), while the satellite `Name` stays the bare
    `TRC_<n>` so the array remains a tandem family too. See `scripts/make_unified_annotation.R`
    (`normalise_tc_satellite`) and the "Unified annotation: rDNA labelling" section of `CLAUDE.md`.
- correction in reporting overlaping annotation
  → cross-TRC overlaps resolved by TideCluster 1.16.0; the unified annotation writes a non-fatal
    `Repeat_Annotation_Unified.overlaps.tsv` listing any residual Level-1 overlaps (empty on the
    calibration genomes).

# test genome:
`/home/petr/PycharmProjects/assembly_repeat_annotation_pipeline/tmp/genome_cleaned.fasta.gz`
(OZ408 TideCluster 1.16.0 ground-truth fixture; OZ408687.1 TRC_13 = the TE_origin calibration case.)

# bug to fix in unified annotation — ✅ DONE (TE_origin; 1.0.0rc1, refined 1.0.0rc3)

- there is a conflict in example genome in sttructural annotation - if tandem repeat originate from LTR retrotranspsons then both annotation are reporte
- such region can be found for example in OZ408687.1:18,715,299..18,883,725 which contain arrays of
Class_I|LTR|Ty1/copia|Ale element nicelly tandemly organized - here correct annotation should be Satellite/TideCluster/TRC_13 but with not that is is derived from TE  - something like TE origin
  → the satellite wins the region and carries `TE_origin=<LCA class of the covered TEs>`; the
    underlying structural elements are trimmed out of the unified file (they remain in `DANTE_*.gff3`).
    Two sub-types: a tandem of *complete* LTR-RTs gets `TE_origin_structure=tandem_LTR_RT`; a
    scattered/incomplete TE-derived array gets `TE_origin` only. See `identify_te_derived_trcs()`.
- run verificaion that there are not other overlaps in unified annotation (exceps simple repeat/low complexity)
  → `Repeat_Annotation_Unified.overlaps.tsv` does exactly this (non-fatal). Empty on the OZ408
    calibration genome and on all three full hardening assemblies (Boechera, Dunaliella, A. thaliana T2T).

# related follow-ups handled in the same line
- Tandem LTR-RT (`LTR_RT_TR`) head-to-tail shared-LTR arrays collapsed to container+members
  (`scripts/resolve_ltr_tandems.py`, new `DANTE_LTR_tandems.gff3`). (1.0.0rc3)
- Density BigWig tracks capped at `[0, 1]` via a strand-agnostic union before coverage. (1.0.0rc3)
- Unified-GFF3 output contract written down + drift guard (`scripts/validate_unified_gff3.py`,
  `docs/unified_annotation_gff3_spec.md`). (1.0.0rc2)
