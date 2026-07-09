# Changelog

## Unreleased
- **Pin CARP's own conda env dependencies to the validated 1.0.7 SIF versions
  (reproducibility).** The env files listed several deps unpinned
  (`envs/dante_line.yaml`: `seqkit`, `parasail-python`, `mmseqs2`, `blast`,
  `pyyaml`; `envs/tidecluster.yaml` / `envs/tidecluster_run.yaml`: the
  `r-jsonlite`/`r-yaml`/`r-optparse`/`pyyaml`/`jq` helper deps and
  `bioconductor-rtracklayer`/`mmseqs2`), so each build re-solved them to whatever
  the channels served that day. They are now pinned to the exact versions that
  the validated 1.0.7 container shipped (a known co-solvable set), so the CARP
  envs build reproducibly. The r-* helper pins deliberately differ between the
  two TideCluster envs because they track each env's `r-base` (4.1.3 vs 4.3.3).
  The upstream tools' transitive deps are still unconstrained — that is the
  subject of `docs/dante_tidecluster_dependency_pinning_request.md`.

## 1.0.7
- **Pin the Singularity base image and harden `%post` against base-image drift
  (fixes the SIF build).** The `continuumio/miniconda3` base was unpinned, so
  every build pulled `:latest`; between the working 1.0.4 build and the next one
  it drifted twice and broke the release build (the `%post` and `envs/` were
  byte-identical the whole time — `git diff 1.0.4 1.0.5 -- Singularity envs/` is
  empty). First a newer conda began enforcing Anaconda's channel Terms of
  Service on the default `pkgs/main`/`pkgs/r` channels
  (`CondaToSNonInteractiveError`) — and that conda build enforced it but didn't
  ship the `conda tos` CLI to accept it. Then the base python jumped to 3.14 with
  an Anaconda helper package (`anaconda-channel-guide`) that hard-pinned it,
  blocking the downgrade to 3.11 (`LibMambaUnsatisfiableError`). Fixes:
  (1) **pin** `From: continuumio/miniconda3:24.9.2-0` (a pre-drift conda 24.9.2
  build — the same conda the project's sandbox uses); (2) build only from
  **conda-forge + bioconda**, never Anaconda's default channels (drops `defaults`
  from the config + `--override-channels` on the bootstrap installs) — avoids the
  ToS gate on any conda and keeps CARP off Anaconda's commercially-licensed
  channels; (3) stop downgrading the base env's python — snakemake 8.12.0 is a
  noarch package, so it installs onto the base python (3.12) as-is; the old
  `conda install python=3.11` rewrote the base env and broke the base `conda`
  package, whose next operation then failed plugin discovery
  (`PluginError: Conflicting post_solves plugins: signature-verification`);
  (4) remove the Anaconda telemetry helper packages (`anaconda-anon-usage`,
  `anaconda-channel-guide`) so `conda info --json` stays clean for Snakemake.
  (The earlier 1.0.4 `jq` system-install fix for the Bioconductor
  `genomeinfodbdata` post-link is retained.)
- **HTML report: fix large-genome crash and make report generation non-fatal.**
  `make_repeat_report.R` aborted on a large genome while building the density
  panels: a bin midpoint on a chromosome > ~1.07 Gbp overflowed 32-bit integer
  arithmetic (`starts + ends` → `NA`), and those `NA`s drove the recursive
  `rdp_simplify` line-simplifier into O(n) recursion → `Error: C stack usage is
  too close to the limit`. Fixes: (1) do the midpoint sum in double precision so
  it never overflows; (2) rewrite `rdp_simplify` with an explicit stack (no
  recursion → cannot overflow the C stack), verified identical to the old output
  on 300 random cases. Robustness: (3) each chart section is now built via
  `safe_build()` so a failure in one panel leaves a visible "not generated"
  placeholder and the rest of the report still renders; (4) the `make_repeat_report`
  rule writes a minimal placeholder HTML if the script fails entirely, so a
  report failure can never abort the pipeline (the annotation outputs are
  unaffected).

## 1.0.4
- **Fix Singularity/SIF build failure (bioconda `genomeinfodbdata` post-link),
  properly this time.** `bioconductor-genomeinfodbdata`'s post-link script runs
  `installBiocDataPackage.sh` → `yq` → `jq`; with `jq` absent the link errored
  and the whole SIF build failed (`yq: Error starting jq`). 1.0.3 tried pinning
  `jq` inside the two Bioconductor env yamls, but that did **not** fix it: conda
  gives no ordering guarantee that the env-level `jq` links before
  `genomeinfodbdata`'s post-link runs (no dependency between them), so `jq` can
  still be missing at post-link time. The authoritative fix is to install `jq`
  at the **system level** in the `Singularity` `%post` (`apt-get install jq`,
  before the conda env creation), so it is on `PATH` for every post-link
  regardless of link order. The env-yaml `jq` is kept as redundancy. (1.0.2 and
  1.0.3 published no artifacts — their builds never completed — so this is the
  first release to ship the flank-alignment grouping change below.)
- **Bounded, deterministic grouping for the DANTE_TIR_FALLBACK / DANTE_LINE flank
  alignment (OOM fix on large genomes).** The all-vs-all flank aligner
  (`scripts/global_local_aln.py`) is O(N²) in memory and compute; on a big
  genome a single high-copy TIR subtype (e.g. CACTA/Mutator) has thousands of
  near-identical TPase anchors, so even after the mmseqs prefilter the candidate
  pairs stay ~O(N²) and the run OOMs. It now mirrors DANTE_TIR's `--max_class_size`
  strategy: when a subtype (fallback) or the LINE pattern set exceeds
  `dante_tir_fallback_max_group_size` / `dante_line_max_group_size` (default
  `1000`), the sequences are split into deterministic, mmseqs-clustering-based
  groups of at most that size and aligned per group, dropping cost to O(N·group).
  Clustering keeps similar flanks together, so the split is near-lossless
  (verified byte-identical `Selected_Length` vs the ungrouped path on distinct
  families); at or below the threshold the code path is unchanged. New unit test
  `tests/test_aln_grouping.py`.

## 1.0.1
- **TideCluster 1.16.1 → 1.16.3.** 1.16.2 fixed the TRC-superfamily
  empty-fallback naming (canonical `<prefix>_trc_superfamilies.csv` with header,
  even when empty), giving the map consumed by
  `tidecluster_reannotate_superfamily_merge` a stable name/schema. 1.16.3 fixes
  the `tc_reannotate` chunk-pool RepeatMasker library-build race that could
  silently truncate the reannotation on a fresh env's first run (RM libraries
  are now built once before the pool, and a failed chunk fails the run instead
  of being silently dropped).
- **New `tidecluster_reannotate_superfamily_merge` option (default `True`).**
  Superfamily-aware array recovery for the RM-on-TideCluster reannotation.
  `tc_reannotate`'s array-length filter is strictly per-TRC, so a real tandem
  array tiled by ≥2 TRCs of the **same superfamily** is fragmented below
  threshold and lost — the underlying TE then prevails, and enabling
  `tidecluster_reannotate_culling_limit` only masks this by collapsing each locus
  to one TRC (making the annotation culling-dependent). The pipeline now runs
  `tc_reannotate --debug` and re-filters its raw hits
  (`scripts/tc_reannotate_sf_filter.py`) grouping sibling TRCs by superfamily, so
  those arrays are recovered **deterministically, independent of culling**, while
  every feature keeps its bare `TRC_<n>` `Name`.
- **New `rm_tc_tandem_gate` option (default `True`).** A Tier-4 RM-on-TideCluster
  satellite may override a Tier-5 TE call only where it has independent tandem
  evidence (raw TideHunter); an unsupported RM_TC array — a short AT-rich
  consensus tiling a genuinely non-tandem TE (e.g. a Tekay LTR-RT) — is demoted
  below the TE, preventing spurious TE→satellite over-masking. Genuine satellites
  and RM_TC over non-TE sequence are unaffected. Measured on a Merodon genome:
  recovers ~7 Mb of TE that culling had re-labelled satellite, while keeping the
  ~22 Mb of real TE-derived tandem arrays.

## 1.0.0
- **New `repeatmasker_culling_limit` option (default `0` = off).** Passes rmblastn
  `-culling_limit N` to the RepeatMasker step to cap the redundant per-locus HSP
  explosion a de-novo library causes (each genomic TE copy aligns to hundreds of
  near-identical consensi). `2` → ~3× faster RepeatMasker at ~−0.7 % masked bp,
  classification preserved. Injected via an `RMBLAST_DIR` shim
  (`scripts/rmblast_culling_shim.py`) — no RepeatMasker/container patch, works in
  the existing image.
- **New `tidecluster_reannotate_culling_limit` option (default `0` = off).** The
  same culling for TideCluster's internal RepeatMasker in `tidecluster_reannotate`
  — the strongest culling target (a consensus dimer matches a satellite array at
  every rotational phase). `2` → ~3.7× faster reannotation at ~−0.8 % masked bp.
  Independent of `repeatmasker_culling_limit`; validate masked bp on
  large-satellite genomes before enabling (may under-mask megasatellite arrays).
- **`dante_tir_fallback` performance fix.** Reverted a subtype-parallelisation
  that regressed the rule ~3.6× on genomes with size-imbalanced TIR subtypes (the
  dominant subtype was starved to `threads/N` while small subtypes left cores
  idle). Subtypes now run serially with the full thread budget; output unchanged.
- **Tool bumps:** TideCluster 1.16.0 → 1.16.1, dante_tir 0.2.5 → 0.2.6.
- **rDNA classification restructured (breaking; output schema v3).** The flat
  top-level classes `rDNA_45S` / `rDNA_5S` are now nested under a single `rDNA`
  parent — `rDNA/45S_rDNA`, `rDNA/45S_rDNA/{18S,25S,5.8S,IGS,ITS1,ITS2}`,
  `rDNA/5S_rDNA`, `rDNA/5S_rDNA/5S` (`classification_vocabulary.yaml` +
  `data/rdna_library.fasta` headers + the TideCluster rDNA map in
  `make_unified_annotation.R`). This **renames the per-class rDNA GFF3/BigWig
  files** (`rDNA_45S.* → rDNA.45S_rDNA.*`, etc.); the `rDNA` **rollups**
  (`rDNA_RepeatMasker.gff3`, `rDNA_{10k,100k}.bw`) and the REQUIRED outputs keep
  their names. No backward-compat symlinks; the manifest `schema_version` bump to
  `3` is the migration signal. Consumers that resolve by manifest logical name
  are unaffected; those globbing `rDNA_45S*` filenames must update. See
  [`docs/output_schema.md`](docs/output_schema.md).
- **HTML report:** classification table and pie now use a **fixed category
  order** (Class_I, Class_II, rDNA, Tandem_repeats, Simple_repeat,
  Low_complexity, Unknown) instead of ordering by content; the pie is
  **genome-relative** (inner ring = Repeats vs grey Non-repetitive, so every
  percentage matches the table); overview cards show **LTR-RT / DNA transposons /
  LINE as % of genome** and the "Sequences in charts" card was removed. The
  table+pie show **rDNA at the array level** (45S / 5S; subunit detail kept in
  the GFF3/density outputs but not broken out). The "Unspecified" row tag was
  renamed **"Unclassified"**. The per-sequence repeat-content panel now draws
  **two reference lines** (shown-sequence average + whole-genome average) with a
  note on how many contigs/bp are omitted from the bars. Fixed an off-by-one that
  dropped the last chromosome separator line in `summary_plots.pdf`.
- **Per-sequence repeat-content panel reworked + Level-1-only statistics.** The
  per-sequence composition bars are now computed at base precision from the
  unified GFF3 Level-1 union (not the smoothed density BigWigs, which inflated
  the totals), so each bar sums to the sequence's true repeat fraction. All
  sub-threshold contigs are aggregated into one prominent "Other (M contigs,
  X Mb)" bar, so the bars cover the whole genome; the redundant "shown avg" line
  is dropped (one "genome avg" line remains). `summary_statistics.csv` and the
  per-class split now count **Level-1 features only** — Level-2 nested children
  (LTR_RT_TR member copies, simple repeats nested in satellites) are no longer
  double-counted, so per-class numbers and the total equal the true union and
  the per-class partition is exactly disjoint. The Class_II bar/legend label is
  now just "Class II" (it includes Helitron, not only TIR).
- **Fixed tier-4 rDNA mislabelling in the unified annotation:** rDNA arrays
  annotated only via the RM-on-TideCluster pass were labelled
  `Satellite/TideCluster/TRC_<n>` instead of rDNA; the authoritative TideCluster
  `<prefix>_rdna.tsv` now drives rDNA labelling across all TideCluster tiers.

## 1.0.0rc3
- Density BigWig tracks no longer exceed 1.0: overlapping features are merged into a **strand-agnostic union** before coverage (`calculate_density.R` / `calculate_density_batch.R`). The Unified annotation tolerates overlap (L1 `Simple_repeat`/`Low_complexity` over a TE; nested L2 children), which previously stacked the total track to ~3.5× and per-class tracks to ~2.2×; every density track is now a true union fraction in `[0, 1]`. Validated across three full assemblies (Boechera, Dunaliella, *A. thaliana* Col-CC T2T).
- **Tandem LTR-RT (`LTR_RT_TR`)**: head-to-tail, same-lineage LTR-RT arrays that share boundary LTRs are collapsed to one Level-1 container with the member copies as Level-2 children (`scripts/resolve_ltr_tandems.py`), so a shared-LTR array is annotated once instead of double-counting the overlapping LTRs. Containers are also written to a new top-level `DANTE_LTR_tandems.gff3` (header-only when none; logical name `dante_ltr_tandems_gff3`). See [`docs/dante_ltr_tandem_feature_request.md`](dante_ltr_tandem_feature_request.md).
- **TE-derived satellite conflict** resolved: where a TideCluster satellite is a tandem of complete LTR-RTs, the satellite wins the region and is tagged `TE_origin` (plus `TE_origin_structure=tandem_LTR_RT` for the full-LTR-RT case), and the underlying structural elements/members are trimmed out of the unified file — eliminating the double annotation. The unified-GFF3 drift guard and spec were extended to cover the new attributes.

## 1.0.0rc2
- Output contract for `Repeat_Annotation_Unified.gff3` written down ([`docs/unified_annotation_gff3_spec.md`](unified_annotation_gff3_spec.md)) and enforced by an executable drift guard (`scripts/validate_unified_gff3.py` + `tests/test_unified_gff3_spec.py`, run in CI and the release gate).
- Fixed a satellite `Name` regression that silently emptied the per-family BigWig outputs — every TideCluster satellite `Name` stays the bare `TRC_<n>` (downstream apps and `split_gff_by_name.R --name-prefix TRC_` key on it); rDNA is routed by `classification` instead.
- New `/release` skill: one-command version bump + cheap-CI gate + tag.

## 1.0.0rc1
- **TideCluster upgraded to 1.16.0**, adding three default-on behaviours the unified annotation consumes: array-level **rDNA identification** (`rDNA_45S`/`rDNA_5S`, labelled clearly as rDNA while still counted as a tandem family), cross-TRC overlap resolution, and a chunked/pooled `tc_reannotate` RepeatMasker (fixes `-pa`-with-custom-`-lib` under-parallelism). New config knobs `tidecluster_detect_rdna` / `tidecluster_rdna_library` / `tidecluster_keep_trc_overlaps` / `tidecluster_chunk_size`.
- Unified annotation now **resolves TR-from-TE overlaps**: a tandem array built from multiple same-family structural TEs is reported once as the satellite with a `TE_origin` tag rather than double-annotated; a non-fatal `Repeat_Annotation_Unified.overlaps.tsv` reports any residual Level-1 overlaps.
- Major performance/memory work: fixed a `reduce_library` merge OOM (~296 GB RSS → bounded), parallelised the density batch / DANTE_TIR fallback / `reduce_library` BLAST+CAP3, and packed small scaffolds into shared RepeatMasker chunks.

## 0.9.0rc10
- New **second-round containment reduction** of the RepeatMasker library (`reduce_library_containment`, default True, independent of `reduce_library`). After the per-class CAP3/mmseqs reduction, a `reduce_library_containment` rule drops short repeat fragments fully contained in a longer element of the **same class** (greedy blastn self-comparison, `scripts/containment_reduce_library.py`); RepeatMasker masks their genomic copies via the container, so masking **and** classification are preserved. Validated masked-bp-lossless with RepeatMasker on the Pisum pangenome over two genome regions (−0.09…−0.12%) while cutting ~22% of library bp and ~30% of RepeatMasker wall-time. Thresholds via `containment_min_identity` (80) / `containment_min_coverage` (0.90); set `reduce_library_containment: False` to feed RepeatMasker the full per-class-reduced library.

## 0.9.0rc9
- Fixed a hard pipeline failure in `tidecluster_reannotate`: the TideCluster dimer-library reducer (`reduce_dimer_library.py`) ran `mmseqs easy-search` with the default nucleotide k-mer (15), which **segfaults on short tandem monomers** (the query is the first-half monomer, down to ~21 bp). Pinned `-k 7`; groups with monomers below the prefilter floor are skipped and any residual `mmseqs` error retries single-threaded then keeps the group unreduced — so the reducer can no longer abort a run.
- Dimer reduction made lossless on more families: replaced single-linkage clustering with greedy/star clustering (each dropped dimer aligns directly to a kept representative) plus a length guard (nested/harmonic periods no longer collapse onto the longest rep). Validated masked-bp-lossless on tiny_pea (−0.02%) and GCA_041296365.1. Known limitation: large satellite families (dimer ≳ a few kb) can still under-mask because RepeatMasker needs the full consensus diversity — set `reduce_tidecluster_library: False` to mask with the full dimer library.

## 0.9.0rc8
- Per-family / by-class density BigWigs build ~140× faster on assemblies with many scaffolds: `calculate_density_batch.R` now parallelises across families (`-t`, wired to `workflow.cores`) and computes density only on the scaffolds each family occupies instead of re-binning the whole genome per family. Track values are byte-identical to rc7; empty scaffolds are no longer written to the BigWig header (occupied-only). On a 2 Gbp / 1888-scaffold assembly with 477 tandem families the step drops from ~6.5 h to ~3 min.
- CI now asserts `carp_manifest.json` matches the produced output tree — every declared logical-name→path must exist after a run (`scripts/assert_manifest_outputs.py`, run on the release-gate fixture) — so the manifest can no longer drift out of sync with the pipeline's outputs.

## 0.9.0rc7
- Emit `carp_manifest.json` at the output root (machine-readable output contract: `schema_version` + logical-name→path map; written by `run_pipeline.py` on success and failure). See [`docs/output_schema.md`](output_schema.md).
- Faster `tidecluster_reannotate`: rotation-invariant reduction of the TideCluster consensus dimer library (aligns each monomer against the already-doubled dimers, collapsing redundant phase variants per TRC). Decoupled from `reduce_library` via the new `reduce_tidecluster_library` option (default True). Lossless for masking (validated; masked bp unchanged within ±0.15%), library ~5–20× smaller, remasking ~1.5–5× faster.
- `summary_plots.pdf` hardened: empty repeat categories render a placeholder panel instead of crashing; any render failure (e.g. data too large) falls back to a one-page placeholder PDF so the rule never fails.
- Snakefile-derived workflow schematic ([`figs/workflow_overview.svg`](../figs/workflow_overview.svg), regenerated by `scripts/make_workflow_diagram.py`).

## 0.9.0rc6
- Density BigWigs reworked: all annotation-derived tracks now sourced from `Repeat_Annotation_Unified.gff3`; BigWig outputs renamed/relocated (breaking — see "Migration" in [`docs/output_reference.md`](output_reference.md)); genome-wide total now includes tandems (`Repeat_density/Repeat_density_total_*.bw`); per-family tandem tracks come in two flavours per `TRC_<n>` (structural TideCluster + Unified union); fixed crash on empty `RM_on_TideCluster_Library.gff3`.

## 0.9.0rc5
- Density BigWigs written run-length-merged (sparse, lossless, much smaller); new `Class_II.Subclass_1.TIR` density rollup; new per-family BigWigs for the RepeatMasker tandem pass (superseded by rc6).

## 0.8.0
- Bug fix in DANTE_LINE, filtering of tandem repeats from DANTE_LINE added.

## 0.7.4
- TideCluster updated to v1.8.0 with `--long` option added to detect tandem repeats with monomer up to 25 kb. Bugfix in subtracting tandem repeats from dispersed repeats.

## 0.7.2
- DANTE_LINE added.

## 0.7.1
- DANTE_TIR added.

## 0.6.7
- More efficient calculation of bigwig files.

## 0.6.6
- DANTE_LTR updated to 0.6.0.4, TideCluster updated to 1.6.

## 0.6.5
- Bugfix in gff3 merging.

## 0.6.4
- dante_ltr runs on smaller chunks (50000000) → better memory usage.

## 0.6.3
- DANTE update to 0.2.5 — bugfix.

## 0.6.2
- Bugfix in bigwig calculation.

## 0.6.1
- DANTE_LTR update to 0.4.0.3 (bugfix).

## 0.6.0
- REXdb Viridiplantae v4.0, library size reduction added, RepeatMasker parallelization added, missing full LTR-RT handling added.

## 0.5.2
- RepeatMasker sensitivity can be set.

## 0.5.1
- Graphical output to PDF added.
