# Annotation Rules for Unified Repeat Annotation Output

## 1. Purpose

This document defines the rules for producing a single unified repeat annotation GFF3 from all
pipeline annotation layers. The goal is to minimise redundant annotation while preserving
biological signal and clearly communicating annotation confidence.

---

## 2. Annotation Layers (in order of decreasing reliability)

### Tier 1 — Structure-based element annotation (highest confidence)

These annotations are derived from detecting structural hallmarks of transposable elements and
are considered the most reliable. A region covered by a Tier 1 annotation should generally
not be overwritten by lower tiers.

| Source | Pipeline output | Description |
|--------|-----------------|-------------|
| DANTE_LTR | `DANTE_LTR/DANTE_LTR.gff3` | Complete LTR retrotransposons (both LTRs + internal domain). Also includes partial elements where LTRs were not identified but the internal domain structure was detected. |
| DANTE_TIR | `DANTE_TIR/DANTE_TIR_final.gff3` | Complete DNA transposons with Terminal Inverted Repeats identified by DANTE_TIR. |
| DANTE_LINE | `DANTE_LINE/DANTE_LINE.gff3` | LINE retrotransposons identified via protein domain clustering and structural analysis. |

### Tier 2 — Filtered protein domain annotation (high confidence)

Protein domain hits that passed DANTE quality filters but were **not** assembled into a complete
or partial element by Tier 1 tools. These represent real TE-related sequences (solo LTRs,
truncated / degenerate elements, orphan domains) and are reliable at the domain level even
though the full element boundary is unknown.

| Source | Pipeline output | Description |
|--------|-----------------|-------------|
| DANTE filtered | `DANTE/DANTE_filtered.gff3` | Filtered TE protein domains not already covered by a Tier 1 element annotation. |

### Tier 3 — Structure-based tandem/satellite repeat annotation (high confidence for satellites)

TideCluster identifies tandem repeats by clustering raw TideHunter arrays; clusters represent
distinct satellite families with well-defined consensus monomers. The annotation boundaries are
based on TideHunter array detection (structural), not sequence similarity.

Within this tier, the **default** run (monomer ≥ 40 bp) takes priority over the **short
monomer** run (10–39 bp) when the two overlap.

| Source | Pipeline output | Priority |
|--------|-----------------|----------|
| TideCluster (default) | `TideCluster/default/TideCluster_clustering.gff3` | Higher — standard monomer length (≥ 40 bp) |
| TideCluster (short monomer) | `TideCluster/short_monomer/TideCluster_clustering.gff3` | Lower — short monomers (10–39 bp) |

### Tier 4 — Similarity-based tandem/satellite annotation (moderate confidence)

RepeatMasker run using the TideCluster-derived consensus dimer library. Extends satellite
coverage to dispersed/degenerate copies that did not trigger TideCluster clustering.

| Source | Pipeline output | Description |
|--------|-----------------|-------------|
| RM on TideCluster library | `TideCluster/default/RM_on_TideCluster_Library.gff3` | RepeatMasker similarity hits to TideCluster consensus sequences. |

### Tier 5 — Similarity-based TE annotation (moderate confidence)

RepeatMasker run using the combined TE library (LTR-RT representatives, DANTE_TIR library,
LINE library, custom library, rDNA library). Covers degenerate/truncated elements and element
copies whose structures were not detected by Tier 1 tools. Lower confidence than structure-based
because similarity can be misleading near diverged or chimeric elements.

| Source | Pipeline output | Description |
|--------|-----------------|-------------|
| RepeatMasker + DANTE merged | `RepeatMasker/RM_on_combined_library_plus_DANTE.gff3` | RM similarity hits merged with Tier 2 DANTE domains. |

### Tier 6 — Residual tandem repeat annotation (lowest confidence)

Raw TideHunter arrays that did **not** reach TideCluster clustering thresholds. These are
genuine tandem repeats but their classification is uncertain (no matching satellite family).

| Source | Pipeline output | Description |
|--------|-----------------|-------------|
| TideHunter (default) | `TideCluster/default/TideCluster_tidehunter_short.gff3` | Sub-threshold tandem arrays, standard run. |
| TideHunter (short monomer) | `TideCluster/short_monomer/TideCluster_tidehunter_short.gff3` | Sub-threshold tandem arrays, short monomer run. |

---

## 3. Conflict Resolution Rules

### 3.1 General principle — higher tier wins, lower tier is trimmed

When a genomic region is covered by annotations from different tiers, the highest-tier
annotation takes priority. The preferred resolution is **trimming** the lower-tier feature
to its non-overlapping portion rather than removing it entirely. This is important because
similarity-based annotation boundaries are imprecise: a Tier 5 RM hit often extends beyond the
boundaries of a Tier 1 structural element, and that non-overlapping flanking portion may
represent a real (degenerate) copy fragment worth retaining.

Removal (full suppression) applies only when the lower-tier feature is **fully contained**
within the higher-tier feature and is of the same classification — in that case it is
redundant and adds no information.

A minimum overlap threshold (suggested: ≥ 80 % of the lower-tier feature length) should be
applied before deciding to suppress rather than trim.

### 3.2 Intra-tier conflicts

Within the same tier, prefer the annotation with the higher classification resolution (more
specific taxonomic/family assignment). When two same-tier annotations of different families
overlap significantly (≥ 50 % of the shorter one), use a **hierarchical resolution**: the
common ancestor in the classification hierarchy becomes the annotation for the overlapping
region (e.g., two hits classified as `Class_I/LTR/Ty1_Copia` and `Class_I/LTR/Ty3_Gypsy`
resolve to `Class_I/LTR`). This approach is already implemented in the RepeatMasker output
parser (`merge_repeat_annotations.R`).

If no common ancestor can be determined above the root, flag with `conflict=true` and retain
both names.

### 3.3 Acceptable overlaps (do NOT suppress)

The following overlap patterns are biologically expected and must be preserved as nested
(parent/child) annotations rather than flattened:

| Inner annotation | Outer annotation | Rationale |
|-----------------|-----------------|-----------|
| Tier 3–4 satellite | Tier 1 LTR-RT (internal region) | Satellite arrays frequently originate from internal duplications of LTR retrotransposons. A tandem array inside an LTR-RT element is a valid nested annotation. |
| Simple repeat / Low complexity (Tier 5) | Any TE (Tier 1–5) | Simple sequence composition can arise within TE bodies and is not a competing annotation. |
| Tier 6 TideHunter residual | Tier 1 LTR-RT | Unclassified short tandem arrays within an LTR-RT body — these typically do not reach the TideCluster size threshold and are better understood as internal TE structure rather than standalone satellites. Retain as nested. |

**Note on Tier 2 domains inside Tier 1 elements**: a DANTE protein domain that falls inside a
structurally annotated Tier 1 element should be suppressed at that region — it is already
represented by the more complete Tier 1 annotation.

### 3.4 Tandem repeat vs TE conflict (satellite priority rule)

A common and important conflict arises when the TE repeat library contains sequences with
regions similar to known satellite families. RepeatMasker will then annotate those satellite
regions as LTR-RT elements via similarity (Tier 5). The correct resolution is:

- **Tier 3 TideCluster annotation always takes priority over Tier 5 RM annotation** for the
  same region, even when the Tier 5 hit nominally has a "higher" TE classification.
- When a Tier 5 RM hit overlaps a Tier 3 satellite cluster, the RM annotation is trimmed or
  suppressed for the overlapping region. The satellite (Tier 3) is the Level 1 feature there.
- The nested case (§3.3 first row) — a short tandem array *inside* a Tier 1 LTR element body
  that does not reach TideCluster thresholds — is distinct from this conflict and is handled
  as a nested feature, not a conflict.

### 3.5 Unacceptable overlaps (must be resolved)

| Situation | Resolution |
|-----------|------------|
| Same region annotated as two different TE families at Tier 5 | Suppress the hit with lower bitscore / lower percent identity; use hierarchical name resolution for the overlap zone. |
| Tier 5 RM hit fully inside a Tier 1 element of the **same** family | Suppress the Tier 5 record (redundant). |
| Tier 5 RM hit fully inside a Tier 1 element of a **different** family | Flag as potential mis-annotation or chimeric element; retain both with `conflict=true`. |
| Tier 3 short-monomer cluster fully inside a Tier 3 default cluster | Suppress the short-monomer record (Tier 3 default has priority within tier). |

---

## 4. Tandem Repeat ↔ Mobile Element Interaction

This is the most complex overlap case in the pipeline and deserves its own section.

### Background

Long satellite arrays (e.g., CentO, FabTR-type) often derive evolutionarily from internal
regions of LTR retrotransposons. In practice this creates two common annotation problems:

1. **RepeatMasker annotates satellite regions as LTR-RT** (via similarity to the TE library):
   this is a frequent and well-known problem. The satellite sequence is similar enough to the
   LTR-RT consensus (because it evolved from it) that RM assigns an LTR-RT classification.
   TideCluster Tier 3 annotation has priority here (see §3.4).

2. **Short tandem arrays inside LTR-RT elements**: TideCluster has a minimum array size
   threshold, so short internal tandem regions of LTR-RT elements do not produce Tier 3 clusters.
   These may appear as Tier 6 TideHunter residuals. They are retained as nested features inside
   the Tier 1 element (not as standalone satellite annotations).

Note: a TideCluster cluster overlapping a DANTE_LTR element boundary does not occur in
practice — the TideCluster size threshold ensures that only large, independently organized
satellite arrays pass clustering.

### Proposed refined handling in unified output

1. **Tier 1 LTR element with internal satellite array (Tier 3/4)**: retain both annotations.
   The outer element boundary comes from DANTE_LTR (Tier 1); the inner tandem array is a
   nested Tier 3 feature. Both present with a `Parent` relationship.

2. **Tier 3/4 satellite cluster NOT inside a Tier 1 element**: treat as standalone satellite.
   Suppress any overlapping Tier 5 RM hit to the LTR-RT library for the same region
   (§3.4 satellite priority rule).

3. **Tier 5 RM hit to LTR-RT library overlapping a Tier 3/4 satellite region**: suppress or
   trim the RM hit. This is the core satellite-priority rule and is already partially
   implemented via the `subtract_satellites_from_rm` rule in the current pipeline.

4. **Tier 6 TideHunter residual inside Tier 3/4 satellite cluster**: suppress the Tier 6
   record (redundant with the more reliable Tier 3/4 annotation).

5. **Tier 6 TideHunter residual NOT covered by any higher tier**: retain as unclassified
   tandem repeat.

---

## 5. GFF3 Format and Normalisation

Before merging the annotation layers, each input GFF3 must be normalised to a consistent
format. The tools produce structurally different GFF3 files that require specific handling.

### 5.1 Input GFF3 structure per tool

| Tool | Top-level feature type | Child features | Notes |
|------|----------------------|----------------|-------|
| DANTE_LTR | `transposable_element` | `long_terminal_repeat`, `target_site_duplication`, `primer_binding_site`, `protein_domain` / `protein_match` | Rich attribute set; uppercase custom attributes must be lowercased; child features often lack IDs. See `dante_ltr_gff3_to_canonical`. |
| DANTE_TIR | `transposable_element` | Currently none (single flat feature per element; may gain children in future versions) | Simpler structure; `Final_Classification` attribute carries family assignment. |
| DANTE_LINE | `repeat_region` or `transposable_element` | `protein_domain` children included | Domain children carry the LINE family evidence. |
| DANTE filtered | `protein_domain` | None | Flat features; `Name` attribute carries classification in `Class/Family` form. |
| TideCluster | `repeat_region` | `satellite` child arrays per cluster | `Name` carries cluster ID; `Satellite` attribute carries monomer length and repeat unit. |
| RM on TideCluster library | `repeat_region` | None | Flat; `Name` = matched TideCluster consensus ID. |
| RepeatMasker (TE library) | `repeat_region` | None | Flat; `Name` carries RM classification string (may use `|` as separator, converted to `/`). |
| TideHunter residuals | `repeat_region` | None | Flat; no family classification. |

### 5.2 Required normalisation steps

The following transformations must be applied before merging:

1. **Attribute case normalisation** (DANTE_LTR): uppercase custom attributes renamed to
   lowercase per `dante_ltr_gff3_to_canonical` (e.g., `Final_Classification` →
   `final_classification`, `LTR_Identity` → `ltr_identity`).

2. **Feature type standardisation** using Sequence Ontology terms where possible:

   | Current type | Canonical SO type | SO accession |
   |-------------|-------------------|--------------|
   | `protein_domain` | `protein_match` | SO:0000349 |
   | `transposable_element` | `transposable_element` | SO:0000101 |
   | `long_terminal_repeat` | `long_terminal_repeat` | SO:0000286 |
   | `target_site_duplication` | `target_site_duplication` | SO:0000434 |
   | `primer_binding_site` | `primer_binding_site` | SO:0005850 |
   | `repeat_region` | `repeat_region` | SO:0000657 |
   | `satellite_DNA` / `satellite` | `satellite_DNA` | SO:0000005 |

3. **Unique IDs**: every feature must have a unique `ID` attribute. Child features in
   DANTE_LTR output frequently lack IDs — generate them from parent ID + feature type +
   counter (as implemented in `dante_ltr_gff3_to_canonical`).

4. **Classification attribute**: add a uniform `classification` attribute to every top-level
   feature using the `Class_I/LTR/Ty3_Gypsy` style path. Source mapping:
   - DANTE_LTR: from `Final_Classification` (after lowercasing)
   - DANTE_TIR: from `Final_Classification`
   - DANTE_LINE: from `Name` or domain-level classification
   - DANTE filtered: from `Name`
   - TideCluster / TideHunter: `Satellite/TideCluster/<cluster_id>` or `Satellite/Unknown`
   - RepeatMasker: from `Name` (separator `|` → `/`)

5. **Separator normalisation**: ensure classification path uses `/` throughout. The RM
   output and some library headers use `|` — convert before merging.

6. **Source field**: set the `source` column (column 2) to the originating tool:
   `DANTE_LTR`, `DANTE_TIR`, `DANTE_LINE`, `DANTE`, `TideCluster`, `TideCluster_RM`,
   `RepeatMasker`, `TideHunter`.

7. **Strip internal/verbose attributes**: remove DANTE_LTR alignment attributes
   (`DB_Seq`, `Region_Seq`, `Query_Seq`, `Best_Hit_DB_Pos`) to reduce output file size,
   unless explicitly requested.

### 5.3 Unified output attributes

Every top-level (Level 1) feature in `Repeat_Annotation_Unified.gff3` carries:

```
ID=<unique_id>
Name=<classification_path>
classification=<classification_path>
source_tier=<1..6>
source_tool=<tool_name>
```

Child (Level 2) features additionally carry:
```
Parent=<Level1_ID>
```

Retain tool-specific informative attributes (e.g., `ltr_identity`, `ndomains`, `tsd`,
`ltr5_length`, `ltr3_length` from DANTE_LTR) as optional attributes on the relevant
features.

---

## 6. Implementation Design

### 6.1 Language and libraries

**R with GenomicRanges / rtracklayer** — all required packages (`GenomicRanges`, `IRanges`,
`rtracklayer`, `BiocParallel`) are already present in the `tidecluster` conda environment and
are used by existing pipeline scripts. No new dependencies are needed.

Key primitives used:

| Operation | GenomicRanges function | Used for |
|-----------|----------------------|----------|
| Load GFF3 with attributes | `rtracklayer::import.gff3` | Read each input tier |
| Non-overlapping decomposition + back-map | `disjoin(gr, with.revmap=TRUE)` | Build Level 1 output |
| Trim lower-tier to non-overlapping part | `setdiff(lower, higher)` | Preserve flanking fragments |
| Detect overlapping pairs | `findOverlaps(query, subject)` | Tier conflict detection |
| Merge same-classification neighbours | `reduce(split(gr, gr$Name))` | Consolidate same-family hits |
| LCA name resolution | `resolve_name()` (existing) | Hierarchical conflict resolution |
| Export valid GFF3 | `rtracklayer::export(format="gff3")` | Write output |
| Per-sequence chunking | `split(gr, seqnames(gr))` | Memory management |
| Parallel processing | `BiocParallel::bplapply` | Speed on multi-core |

The existing `resolve_name()` and `gff_cleanup()` functions from `clean_rm_output.R` and
`merge_repeat_annotations.R` can be reused directly.

### 6.2 Memory management — adaptive chunking strategy

Loading all annotation tiers for a whole genome simultaneously can exceed available RAM for
large or highly repetitive genomes. Processing is therefore split by sequence (chromosome /
contig), but naive per-sequence splitting creates excessive overhead when the assembly
contains many thousands of small contigs.

**Chunking is triggered only when total annotation size exceeds a threshold** (suggested:
total bp of all input GFF3 intervals > 500 Mb, or configurable). Below the threshold,
the entire genome is processed as a single batch.

When chunking is required, sequences are partitioned into **balanced batches** rather than
processed one-by-one:

```
1. Sort sequences by length (descending).
2. Any sequence whose length alone exceeds the per-batch threshold is assigned its own batch.
3. Remaining sequences are packed into batches greedily (next-fit decreasing bin packing):
   accumulate sequences into the current batch until adding the next would exceed the
   target batch size, then start a new batch.
4. Each batch is processed independently; results are concatenated.
```

This ensures:
- A handful of large chromosomes each get a dedicated batch (no RAM spike from combining them).
- Hundreds of small contigs are grouped into a small number of batches (no process-launch
  overhead per contig).
- Total number of batches ≈ `ceil(total_annotation_size / batch_target_size)`, independent
  of contig count.

Suggested defaults: per-batch target = 200 Mb of annotation; minimum batch count = 1;
maximum batch count = number of available cores (so batches map 1:1 to parallel workers).

```r
make_batches <- function(seqlengths, batch_target_bp = 200e6) {
  # seqlengths: named integer vector, sequences sorted largest first
  sl <- sort(seqlengths, decreasing = TRUE)
  batches <- list()
  current <- c()
  current_size <- 0
  for (i in seq_along(sl)) {
    if (sl[i] >= batch_target_bp) {
      # Large sequence: own batch
      batches <- c(batches, list(names(sl)[i]))
    } else if (current_size + sl[i] > batch_target_bp && length(current) > 0) {
      # Current batch full: flush and start new
      batches <- c(batches, list(current))
      current <- names(sl)[i]
      current_size <- sl[i]
    } else {
      current <- c(current, names(sl)[i])
      current_size <- current_size + sl[i]
    }
  }
  if (length(current) > 0) batches <- c(batches, list(current))
  batches
}
```

### 6.3 Script structure

The implementation will be a single R script `scripts/make_unified_annotation.R` with the
following internal organisation:

```
make_unified_annotation.R
├── parse_args()          — optparse: input paths per tier, output path, thresholds
├── load_tier()           — rtracklayer::import.gff3 + normalise attributes (§5.2)
├── make_batches()        — adaptive chunking (§6.2)
├── process_batch()       — core logic for one batch of sequences:
│   ├── apply_tier_priority()   — findOverlaps + setdiff per tier pair
│   ├── build_level1()          — disjoin + revmap + resolve_name → non-overlapping output
│   └── build_level2()          — acceptable overlaps as nested children (§3.3)
├── merge_batches()       — c() results + sortSeqlevels + sort
└── export_gff3()         — assign IDs + rtracklayer::export
```

---

## 7. Proposed Output Structure

The unified annotation output `Repeat_Annotation_Unified.gff3` should contain features at two
levels:

### Level 1 — Repeat regions (non-overlapping)

Each genomic position belongs to at most one Level 1 feature. Conflicts are resolved by the
rules above. This track is suitable for masking and genome-wide statistics.

Attributes: `tier=1..6`, `classification=<name>`, `method=<source tool>`

### Level 2 — Nested features (acceptable overlaps)

Features from lower tiers that overlap a Level 1 feature acceptably (§3.3) are stored as
Level 2 nested children with a `Parent=<Level1_ID>` attribute.

Examples:
- Tandem array nested inside an LTR-RT element.
- Simple repeat nested inside a TE body.

### Summary statistics

Statistics should be computed from **Level 1 features only** to avoid double-counting.
Nested Level 2 features are reported separately (e.g., "satellite bp inside LTR-RTs").

---

## 8. Open Questions / Points for Further Discussion

1. **Minimum overlap threshold**: what fraction of a lower-tier feature must be covered before
   it is suppressed rather than trimmed? Suggested 80 % but needs benchmarking.

2. **Satellite inside LTR-RT — boundary definition**: when a satellite array is inside an
   LTR-RT, the Level 1 feature is the LTR-RT (Tier 1) and the satellite is Level 2. This
   avoids classifying the full LTR-RT extent as "satellite" in summary stats, but means the
   satellite bp within LTR-RTs must be tracked separately.

3. **Tier 2 DANTE domains outside any TE body**: retain as Level 1 with classification from
   the domain type (e.g., `Class_I/LTR/Ty3_gypsy`, feature type `protein_domain`). A Tier 5
   RM hit for the same classification and region → suppress Tier 5.

4. **LINE annotation boundary extension**: DANTE_LINE uses protein domains + structural
   clustering but may not cover the full element. Should the DANTE_LINE boundary be extended
   by RM hits to the LINE library? The current pipeline does not do this — worth considering.

5. **rDNA**: added via a fixed rDNA library to the RM library. rDNA is not a transposable
   element and should always be retained regardless of tier conflicts.

6. **Unknown / unclassified RM hits**: hits with no clear classification (e.g., hits to
   custom library sequences annotated as `Unknown`) — retain at Tier 5 with
   `classification=Unknown`.