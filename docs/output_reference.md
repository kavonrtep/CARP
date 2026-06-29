# CARP output reference

Complete reference for the CARP output directory. For the everyday view (main
output + key files) see the [README](../README.md); for the
`Repeat_Annotation_Unified.gff3` field contract see
[unified_annotation_gff3_spec.md](unified_annotation_gff3_spec.md); for the
machine-readable output map see [output_schema.md](output_schema.md).

## Top-level files

Main annotation (symlinks to the per-tool files for easy access):

- `Repeat_Annotation_Unified.gff3` — the authoritative non-overlapping
  annotation (see [the spec](unified_annotation_gff3_spec.md))
- `DANTE_filtered.gff3` — filtered DANTE protein-domain annotations
- `DANTE_LTR.gff3` — complete LTR retrotransposons
- `DANTE_LTR_tandems.gff3` — head-to-tail LTR-RT arrays (`LTR_RT_TR` containers;
  header-only when none detected)
- `DANTE_TIR.gff3` — TIR DNA transposons (→ `DANTE_TIR_combined.gff3`)
- `DANTE_LINE.gff3` — LINE elements (experimental)
- `Tandem_repeats_TideCluster.gff3` — tandem repeats from TideCluster
- `Tandem_repeats_TideCluster_annotated.gff3` — annotated tandems (if a
  `tandem_repeat_library` was provided)
- `Tandem_repeats_RepeatMasker.gff3` — RepeatMasker annotation of the tandem library
- `Mobile_elements_RepeatMasker.gff3` — Class_I and Class_II transposons
- `Simple_repeats_RepeatMasker.gff3` — simple repeats
- `Low_complexity_RepeatMasker.gff3` — low-complexity regions
- `rDNA_RepeatMasker.gff3` — ribosomal DNA
- `All_Ty1_Copia_RepeatMasker.gff3`, `All_Ty3_Gypsy_RepeatMasker.gff3` — LTR-RT rollups

Summary and masking:

- `summary_statistics.csv` — genome-wide repeat statistics by classification
  (computed from `Repeat_Annotation_Unified.gff3`; the `Tandem_repeats` row
  aggregates TideCluster + TideHunter + RepeatMasker `Satellite/*`)
- `summary_plots.pdf` — repeat-content visualisations
- `all_repeats_for_masking.bed` — merged coordinates of all repeats
- `gaps_10plus.bed` — assembly gaps (N runs ≥ 10 bp)

Metadata (written by `run_pipeline.py`):

- `carp_manifest.json` — machine-readable output contract: `schema_version` plus
  a map of stable logical names → paths (see [output_schema.md](output_schema.md))
- `run_provenance.json` — pipeline version, git SHA, config, and per-environment
  tool versions

HTML reports:

- `TideCluster_report.html` — tandem-repeat analysis
- `DANTE_LTR_report.html` — LTR retrotransposon findings

## Subdirectories

**`DANTE/`** — protein-domain detection
- `DANTE.gff3` — raw domain annotations
- `DANTE_filtered.gff3` — filtered domains
- `DANTE_filtered.fasta` — protein sequences for filtered domains

**`DANTE_LTR/`** — complete LTR retrotransposons
- `DANTE_LTR.gff3` — LTR-RT annotations (both LTRs present)
- `DANTE_LTR_summary.html` — summary report
- `LTR_RTs_library.fasta` — representative LTR-RT sequences for RepeatMasker
- `library/` — clustering and library-construction detail

**`DANTE_TIR/`** — TIR DNA transposons
- `DANTE_TIR_final.gff3` / `.fasta` — final primary TIR annotations / sequences
- `DANTE_TIR_combined.gff3` — primary + non-overlapping fallback (top-level
  `DANTE_TIR.gff3` points here)
- `DANTE_TIR_fallback_filtered.fasta` — non-overlapping fallback sequences
- `TIR_classification_summary.txt` — classification statistics
- `all_representative_elements_combined.fasta` — primary-only library used by
  RepeatMasker
- `fallback_library.fasta` — optional fallback-derived library (empty unless
  `include_dante_tir_fallback_in_library: true`)
- `fallback_library_dropped.tsv` — audit log of fallback-library decisions

**`DANTE_LINE/`** — LINE elements (experimental; may not work well on all genomes)
- `DANTE_LINE.gff3` — LINE annotations
- `LINE_rep_lib.fasta` — representative library for RepeatMasker
- `LINE_regions.fasta`, `LINE_regions_extended.fasta` — extracted / extended sequences

**`TideCluster/`** — tandem repeats
- `default/` — standard run (40–5000 bp monomers)
  - `TideCluster_clustering.gff3` — clustered tandem families
  - `TideCluster_tidehunter.gff3` — raw TideHunter detections
  - `TideCluster_annotation.gff3` — annotated against the custom library (if any)
  - `TideCluster_index.html` — interactive report
  - `TideCluster_consensus_dimer_library.fasta` — consensus sequences
  - `TideCluster_clustering_split_files/` — one GFF3 per cluster (`TRC_<n>`),
    feeding the structural per-family BigWigs
  - `RM_on_TideCluster_Library.gff3` — RepeatMasker annotation using the tandem library
- `short_monomer/` — short-monomer run (10–39 bp); same structure as `default/`
- `TideCluster_clustering_default_and_short_merged.gff3` — merged results

**`Libraries/`** — repeat libraries
- `combined_library.fasta` — full library: `LTR_RTs_library_clean.fasta` +
  custom library (if any) + primary TIR + fallback (only when enabled) + LINE +
  rDNA
- `combined_library_short_names.fasta` — shortened IDs
- `combined_library_reduced.fasta` — size-reduced library (if `reduce_library: True`)
- `LTR_RTs_library_clean.fasta` — LTR library after Class_II-contamination filtering
- `class_ii_library.fasta` — Class_II/Subclass_1 elements used for that filtering

**`RepeatMasker/`** — similarity annotation
- `RM_on_combined_library.out` / `.gff3` — RepeatMasker output / GFF3
- `RM_on_combined_library_plus_DANTE.gff3` — RepeatMasker merged with DANTE
- `Repeat_Annotation_NoSat.gff3` — RepeatMasker annotation excluding tandems;
  used **only** for `all_repeats_for_masking.bed` (not for statistics or density)

**`Repeat_density/`** — genome-wide total density
- `Repeat_density_total_10k.bw`, `Repeat_density_total_100k.bw` — density of all
  repeats from `Repeat_Annotation_Unified.gff3` (includes tandems)

**`Repeat_Annotation_NoSat_split_by_class_gff3/`** — per-class GFF3s
- One GFF3 per major class, derived from `Repeat_Annotation_Unified.gff3`
  (directory name kept for backward compatibility)
- `Tandem_repeats.gff3` aggregates TideCluster default + short + TideHunter +
  RepeatMasker `Satellite/*`
- `Class_II.Subclass_1.TIR.gff3` aggregates all TIR superfamilies

**`Repeat_density_by_class_bigwig/`** — per-class density (see below)
- One BigWig per class/superfamily at 10 kb and 100 kb, derived from
  `Repeat_Annotation_Unified.gff3`, plus the rollups `All_Ty1_Copia_*`,
  `All_Ty3_Gypsy_*`, `Class_II.Subclass_1.TIR_*`

**Tandem-repeat density (top level)**
- `Tandem_repeats_TideCluster_{10k,100k}.bw` — structural TideCluster total
- `Tandem_repeats_TideCluster_split_by_family_bigwig/{10k,100k}/TRC_<n>_*.bw` —
  structural per-family
- `Tandem_repeats_unified_split_by_family_bigwig/{10k,100k}/TRC_<n>_*.bw` —
  unified per-family (structural ∪ similarity, tier-resolved)
- `Tandem_repeats_unified_split_by_family_gff3/TRC_<n>.gff3` — the GFF3s those
  BigWigs are built from

## Annotation categories — full hierarchy

The complete leaf list (and tool-native aliases) is the single source of truth
in [`classification_vocabulary.yaml`](../classification_vocabulary.yaml). The
mobile-element hierarchy below is REXdb Viridiplantae v4.0; satellite, rDNA,
simple/low-complexity, and unknown classes come from RepeatMasker / TideCluster.

```text
Class_I/DIRS
Class_I/LINE
Class_I/LTR/Ty1_copia
Class_I/LTR/Ty1_copia/Ale
Class_I/LTR/Ty1_copia/Alesia
Class_I/LTR/Ty1_copia/Alexandra
Class_I/LTR/Ty1_copia/Angela
Class_I/LTR/Ty1_copia/Bianca
Class_I/LTR/Ty1_copia/Bryana
Class_I/LTR/Ty1_copia/Bryco
Class_I/LTR/Ty1_copia/Ferco
Class_I/LTR/Ty1_copia/Gymco-I
Class_I/LTR/Ty1_copia/Gymco-II
Class_I/LTR/Ty1_copia/Gymco-III
Class_I/LTR/Ty1_copia/Gymco-IV
Class_I/LTR/Ty1_copia/Ikeros
Class_I/LTR/Ty1_copia/Ivana
Class_I/LTR/Ty1_copia/Lyco
Class_I/LTR/Ty1_copia/Osser
Class_I/LTR/Ty1_copia/SIRE
Class_I/LTR/Ty1_copia/TAR
Class_I/LTR/Ty1_copia/Tork
Class_I/LTR/Ty1_copia/Ty1-outgroup
Class_I/LTR/Ty3_gypsy/chromovirus
Class_I/LTR/Ty3_gypsy/chromovirus/Chlamyvir
Class_I/LTR/Ty3_gypsy/chromovirus/chromo-outgroup
Class_I/LTR/Ty3_gypsy/chromovirus/chromo-unclass
Class_I/LTR/Ty3_gypsy/chromovirus/CRM
Class_I/LTR/Ty3_gypsy/chromovirus/Ferney
Class_I/LTR/Ty3_gypsy/chromovirus/Galadriel
Class_I/LTR/Ty3_gypsy/chromovirus/Reina
Class_I/LTR/Ty3_gypsy/chromovirus/Tcn1
Class_I/LTR/Ty3_gypsy/chromovirus/Tekay
Class_I/LTR/Ty3_gypsy/non-chromovirus/non-chromo-outgroup
Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Athila
Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Tat/Ogre
Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Tat/Retand
Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Tat/TatI
Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Tat/TatII
Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Tat/TatIII
Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Tatius
Class_I/LTR/Ty3_gypsy/non-chromovirus/Phygy
Class_I/LTR/Ty3_gypsy/non-chromovirus/Selgy
Class_I/pararetrovirus
Class_I/Penelope
Class_I/SINE
Class_II/Subclass_1/TIR/EnSpm_CACTA
Class_II/Subclass_1/TIR/hAT
Class_II/Subclass_1/TIR/Kolobok
Class_II/Subclass_1/TIR/Merlin
Class_II/Subclass_1/TIR/MITE
Class_II/Subclass_1/TIR/MITE/Stowaway
Class_II/Subclass_1/TIR/MuDR_Mutator
Class_II/Subclass_1/TIR/Novosib
Class_II/Subclass_1/TIR/P
Class_II/Subclass_1/TIR/PIF_Harbinger
Class_II/Subclass_1/TIR/PiggyBac
Class_II/Subclass_1/TIR/Sola1
Class_II/Subclass_1/TIR/Sola2
Class_II/Subclass_1/TIR/Tc1_Mariner
Class_II/Subclass_2/Helitron
rDNA
rDNA/45S_rDNA
rDNA/45S_rDNA/18S
rDNA/45S_rDNA/25S
rDNA/45S_rDNA/5.8S
rDNA/45S_rDNA/IGS
rDNA/45S_rDNA/ITS1
rDNA/45S_rDNA/ITS2
rDNA/5S_rDNA
rDNA/5S_rDNA/5S
```

## Density BigWig tracks

All `*.bw` outputs are per-window **density** tracks with identical value
semantics, ready for genome browsers (JBrowse, IGV, UCSC). **Every
annotation-derived track is built from `Repeat_Annotation_Unified.gff3`**; the
only exception is the structural TideCluster tracks, built from the TideCluster
clustering alone.

- **Value** = fraction of the window covered by that track's features (`0`–`1`),
  computed as mean per-base coverage per bin and smoothed with a 10-bin moving
  average. Overlapping features are merged into a **strand-agnostic union**
  before coverage, so the value is **always within `[0, 1]`** even though the
  Unified annotation tolerates some overlap (a Level-1 `Simple_repeat`/
  `Low_complexity` over a TE; nested Level-2 children).
- **Window size:** `_10k` files use 1 kb bins (~10 kb effective window); `_100k`
  files use 10 kb bins (~100 kb effective window).
- **Encoding (sparse):** tracks are run-length-merged — consecutive windows with
  an exactly equal value are collapsed into one interval. Lossless and needs no
  consumer change; makes per-family tracks much smaller.

### Track families

- **Total** — `Repeat_density/Repeat_density_total_{10k,100k}.bw`: all repeats
  genome-wide (Unified).
- **Per class / superfamily** — `Repeat_density_by_class_bigwig/{10k,100k}/…`:
  one track per class, plus the `All_Ty1_Copia`, `All_Ty3_Gypsy`, and
  `Class_II.Subclass_1.TIR` rollups (Unified).

> **Partition vs. roll-up (do not double-count).** The **exact-classification
> tracks** (`Class_I.LTR.Ty1_copia.Ale`, `Tandem_repeats`, `rDNA.45S_rDNA`,
> `Simple_repeat`, `Low_complexity`, `Unknown`, …) are a **disjoint partition** —
> each feature appears in exactly one, at its leaf classification — so they sum
> (approximately, modulo smoothing) to the **Total** track. The **roll-up
> tracks** — `All_Ty1_Copia`, `All_Ty3_Gypsy`, `Class_II.Subclass_1.TIR`, and
> `Mobile_elements` — are cumulative supersets that deliberately **overlap** the
> partition. For a stacked/summed view use the partition tracks; never sum a
> roll-up together with the leaf tracks it contains.

### Tandem-repeat tracks — which one to use

Tandem repeats get an aggregate track plus **two per-family views** (structural
and unified). They are **not** interchangeable:

| Track | Detection basis | Scope | Per-family? |
|---|---|---|---|
| `Repeat_density_by_class_bigwig/{10k,100k}/Tandem_repeats_*.bw` | Unified (structural ∪ similarity, tier-resolved) | aggregate of all families | no |
| `Tandem_repeats_unified_split_by_family_bigwig/{10k,100k}/TRC_<n>_*.bw` | Unified `Name=TRC_<n>` | one track per family | yes |
| `Tandem_repeats_TideCluster_split_by_family_bigwig/{10k,100k}/TRC_<n>_*.bw` | structural only (TideCluster clustering) | one track per family | yes |
| `Tandem_repeats_TideCluster_{10k,100k}.bw` | structural only (TideCluster clustering) | aggregate of all families | no |

Everything Unified-derived is **tier-resolved**: structure-based mobile elements
(LTR/TIR/LINE/DANTE) outrank tandem layers, and within tandem the TideCluster
structural calls outrank the RepeatMasker remasking. Any tandem territory
overlapping a higher-priority element is reassigned to that element.

Guidance:

- For the authoritative "where are the tandem repeats" view, use the **Unified**
  tracks — the aggregate `Tandem_repeats_*.bw`, or the unified per-family tracks
  for one family.
- For each `TRC_<n>` you get a *unified* track (authoritative) and a *structural*
  one (TideCluster only). Unified adds the RepeatMasker-remasking territory and
  removes any stretch reassigned to a higher-priority element — usually larger
  than structural, but smaller where a family overlaps a mobile element. Pair
  them to see structural vs. similarity capture. (When the RM-on-TideCluster pass
  is empty and nothing overlaps, the two are identical.)
- Per-family tandem tracks cover only `TRC_<n>`-named families; RM
  `Satellite/Unknown` and TideHunter residuals contribute to the aggregate but
  have no family track.

### Migration — BigWig paths changed (0.9.0rc6)

BigWig outputs were renamed/relocated and re-sourced from Unified. GFF3 paths are
unchanged. Repoint any genome-browser tracks:

| Old (≤ 0.9.0rc5) | New (≥ 0.9.0rc6) | Note |
|---|---|---|
| `RepeatMasker/Repeat_Annotation_NoSat_{10k,100k}.bw` | `Repeat_density/Repeat_density_total_{10k,100k}.bw` | source: RM-only-no-tandems → **Unified (all repeats)** |
| `Repeat_Annotation_NoSat_split_by_class_bigwig/` | `Repeat_density_by_class_bigwig/` | values unchanged (already Unified) |
| `TideCluster/default/TideCluster_clustering_{10k,100k}.bw` | `Tandem_repeats_TideCluster_{10k,100k}.bw` | values unchanged (structural TC) |
| `TideCluster/default/TideCluster_clustering_split_files_bigwig/` | `Tandem_repeats_TideCluster_split_by_family_bigwig/` | values unchanged (structural TC) |
| `Tandem_repeats_RepeatMasker_split_files_bigwig/` (rc-only) | `Tandem_repeats_unified_split_by_family_bigwig/` | source: raw RM-on-TC → **Unified `TRC_<n>`** |
