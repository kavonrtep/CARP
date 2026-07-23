# Unified Repeat Annotation — Detailed Implementation Plan

## 1. Overview

Script: `scripts/make_unified_annotation.R`
Output: `{output_dir}/Repeat_Annotation_Unified.gff3`
Rule: `make_unified_annotation` (to be added to `Snakefile`)

---

## 2. Inputs

All paths are confirmed from `tmp/output/` test data.

| Tier | Variable | Path | Feature types | Size (test) |
|------|----------|------|---------------|-------------|
| 1 | `ltr_gff` | `DANTE_LTR/DANTE_LTR.gff3` | `transposable_element`, `protein_domain`, `long_terminal_repeat`, `target_site_duplication`, `primer_binding_site` | 33 850 lines |
| 1 | `tir_gff` | `DANTE_TIR/DANTE_TIR_final.gff3` | `sequence_feature` (flat) | — |
| 1 | `line_gff` | `DANTE_LINE/DANTE_LINE.gff3` | `LINE_element`, `protein_domain` | — |
| 2 | `dante_gff` | `DANTE/DANTE_filtered.gff3` | `protein_domain` (flat) | — |
| 3 | `tc_default_gff` | `TideCluster/default/TideCluster_clustering.gff3` | `tandem_repeat` | — |
| 3 | `tc_short_gff` | `TideCluster/short_monomer/TideCluster_clustering.gff3` | `tandem_repeat` | — |
| 4 | `tc_rm_gff` | `TideCluster/default/RM_on_TideCluster_Library.gff3` | `repeat` | — |
| 5 | `rm_gff` | `RepeatMasker/RM_on_combined_library_plus_DANTE.gff3` | `repeat_region` | 243 398 lines |
| 6 | `th_default_gff` | `TideCluster/default/TideCluster_tidehunter_short.gff3` | `tandem_repeat` | — |
| 6 | `th_short_gff` | `TideCluster/short_monomer/TideCluster_tidehunter_short.gff3` | `tandem_repeat` | — |

Plus:
- `genome_fai`: `genome_cleaned.fasta.fai` — sequence names and lengths for chunking

---

## 3. CLI Arguments

```r
option_list <- list(
  make_option("--ltr",        type="character", help="DANTE_LTR GFF3"),
  make_option("--tir",        type="character", help="DANTE_TIR GFF3"),
  make_option("--line",       type="character", help="DANTE_LINE GFF3"),
  make_option("--dante",      type="character", help="DANTE filtered GFF3"),
  make_option("--tc_default", type="character", help="TideCluster default GFF3"),
  make_option("--tc_short",   type="character", help="TideCluster short monomer GFF3"),
  make_option("--tc_rm",      type="character", help="RM on TideCluster library GFF3"),
  make_option("--rm",         type="character", help="RepeatMasker+DANTE merged GFF3"),
  make_option("--th_default", type="character", help="TideHunter default residuals GFF3"),
  make_option("--th_short",   type="character", help="TideHunter short residuals GFF3"),
  make_option("--fai",        type="character", help="Genome FAI file"),
  make_option("--output",     type="character", help="Output unified GFF3"),
  make_option("--threads",    type="integer",   default=4),
  make_option("--batch_size", type="double",    default=200e6,
              help="Target annotation bp per processing batch [200000000]"),
  make_option("--chunk_threshold", type="double", default=500e6,
              help="Total annotation bp above which chunking is activated [500000000]"),
  make_option("--min_feature_length", type="integer", default=50,
              help="Minimum bp for a trimmed feature fragment to be retained [50]")
)
```

---

## 4. GFF3 Normalisation (per-tool)

Each tool has a `load_<tool>()` function that returns a GRanges with uniform attributes:
`Name`, `classification`, `source_tier` (integer 1–6), `source_tool` (string).
Top-level features are separated from child features.

### 4.1 Attribute discrepancies observed in test data

| Tool | Feature type as found | Correct SO type | Classification attr | Separator | Notes |
|------|-----------------------|-----------------|--------------------|-----------|-|
| DANTE_LTR | `transposable_element` | keep | `Final_Classification` | `\|` → `/` | Uppercase attrs need lowercasing; `TE_partial_*` IDs = partial elements |
| DANTE_TIR | `sequence_feature` | → `transposable_element` | `Classification` | `_` → `/` e.g. `Class_II_Subclass_1_TIR_hAT` → `Class_II/Subclass_1/TIR/hAT` |
| DANTE_LINE | `LINE_element` | → `transposable_element` | `Final_Classification` | `\|` → `/` | Source column is `DANTE` (uppercase) |
| DANTE filtered | `protein_domain` | keep | `Final_Classification` | `\|` → `/` | `Name` = domain name (GAG/RT/RH etc.) |
| TideCluster | `tandem_repeat` | → `repeat_region` | `Name` = cluster ID (e.g. `TRC_5`) | n/a | `classification = paste0("Satellite/TideCluster/", Name)` |
| RM on TC lib | `repeat` | → `repeat_region` | `Name` = cluster ID | n/a | same mapping as TideCluster |
| RepeatMasker+DANTE | `repeat_region` | keep | `Name` already `/`-separated | `/` | Already processed by `clean_rm_output.R` |
| TideHunter | `tandem_repeat` | → `repeat_region` | none | n/a | `classification = "Satellite/Unknown"`; has `consensus_sequence`, `consensus_length` |

### 4.2 DANTE_TIR classification conversion

`Classification` uses `_` as word separator and as hierarchy separator, e.g.:
`Class_II_Subclass_1_TIR_hAT` → path components: `Class_II`, `Subclass_1`, `TIR`, `hAT`

The conversion rule:
```r
convert_tir_classification <- function(x) {
  # Known prefixes that are single tokens containing underscores:
  x <- sub("^Class_II_Subclass_1_", "Class_II/Subclass_1/", x)
  x <- sub("^Class_II_Subclass_2_", "Class_II/Subclass_2/", x)
  # Remaining underscores in the suffix are word separators within a token,
  # except for known compound superfamily names (EnSpm_CACTA, MuDR_Mutator, PIF_Harbinger)
  # → leave them as-is (they are single token names)
  x
}
# Example: "Class_II_Subclass_1_TIR_hAT" → "Class_II/Subclass_1/TIR/hAT"
# Example: "Class_II_Subclass_1_TIR_EnSpm_CACTA" → "Class_II/Subclass_1/TIR/EnSpm_CACTA"
```

### 4.3 load_tier1_ltr()

```r
load_tier1_ltr <- function(path) {
  raw <- import.gff3(path)
  # Separate top-level (transposable_element) from children
  top <- raw[raw$type == "transposable_element"]
  children <- raw[raw$type %in% c("protein_domain", "long_terminal_repeat",
                                   "target_site_duplication", "primer_binding_site")]
  # Normalise classification: replace | with /
  cls <- gsub("|", "/", top$Final_Classification, fixed=TRUE)
  top$classification <- cls
  top$Name <- cls
  top$source_tier <- 1L
  top$source_tool <- "DANTE_LTR"
  # Mark partial vs complete
  top$element_type <- ifelse(grepl("^TE_partial_", top$ID), "partial", "complete")
  # Strip verbose alignment attrs from children
  STRIP <- c("DB_Seq","Region_Seq","Query_Seq","Best_Hit_DB_Pos",
             "cls_prefilter","neighbors_count","upstream_domain","downstream_domain","domain_order")
  for (a in STRIP) mcols(children)[[a]] <- NULL
  list(top=top, children=children)
}
```

### 4.4 load_tier1_tir()

```r
load_tier1_tir <- function(path) {
  raw <- import.gff3(path)
  raw$type <- "transposable_element"
  raw$classification <- convert_tir_classification(raw$Classification)
  raw$Name <- raw$classification
  raw$source_tier <- 1L
  raw$source_tool <- "DANTE_TIR"
  raw
}
```

### 4.5 load_tier1_line()

```r
load_tier1_line <- function(path) {
  raw <- import.gff3(path)
  top <- raw[raw$type == "LINE_element"]
  children <- raw[raw$type == "protein_domain"]
  top$type <- "transposable_element"
  top$classification <- "Class_I/LINE"   # LINE elements all Class_I/LINE
  top$Name <- "Class_I/LINE"
  top$source_tier <- 1L
  top$source_tool <- "DANTE_LINE"
  list(top=top, children=children)
}
```

### 4.6 load_tier2_dante()

```r
load_tier2_dante <- function(path) {
  raw <- import.gff3(path)
  # Keep only protein_domain features (flat GFF3)
  raw <- raw[raw$type == "protein_domain"]
  raw$classification <- gsub("|", "/", raw$Final_Classification, fixed=TRUE)
  raw$Name <- raw$classification
  raw$source_tier <- 2L
  raw$source_tool <- "DANTE"
  raw
}
```

### 4.7 load_tier3_tidecluster()

```r
load_tier3_tidecluster <- function(path_default, path_short) {
  def <- import.gff3(path_default)
  sho <- import.gff3(path_short)
  def$type <- "repeat_region"
  sho$type <- "repeat_region"
  def$classification <- paste0("Satellite/TideCluster/", def$Name)
  sho$classification <- paste0("Satellite/TideCluster/", sho$Name)
  def$source_tier <- 3L;  def$source_tool <- "TideCluster_default"
  sho$source_tier <- 3L;  sho$source_tool <- "TideCluster_short"
  list(default=def, short=sho)
}
```

### 4.8 load_tier4_tc_rm(), load_tier5_rm(), load_tier6_tidehunter()

Similar pattern: `import.gff3`, normalise `type` to `repeat_region`, set `classification`,
`source_tier`, `source_tool`. Tier 5 already has `/`-separated `Name` from
`clean_rm_output.R`; no separator conversion needed.

---

## 5. Tier Resolution Algorithm

Runs per-batch (see §6 for batching). Input is a list of pre-loaded, normalised GRanges for
the current batch's sequences.

```
INPUTS (GRanges for one batch):
  t1     = all Tier 1 top-level features (LTR + TIR + LINE)
  t2     = Tier 2 DANTE filtered domains
  t3_def = Tier 3 TideCluster default
  t3_sho = Tier 3 TideCluster short monomer
  t4     = Tier 4 RM on TideCluster library
  t5_te  = Tier 5 RM features (excluding Simple_repeat / Low_complexity)
  t5_sc  = Tier 5 Simple_repeat + Low_complexity features
  t6     = Tier 6 TideHunter residuals

OUTPUT:
  level1 = non-overlapping GRanges  (one classification per position)
  level2 = nested GRanges with Parent attribute
```

### Step 1 — Tier 1: add all to Level 1 (no trimming)

```r
level1 <- t1
```

### Step 2 — Tier 2: trim to non-Tier-1 regions

```r
t2_trimmed <- trim_to_nonoverlap(t2, t1, min_len=opt$min_feature_length)
level1 <- c(level1, t2_trimmed)
```

### Step 3 — Tier 3 default: satellite priority

TideCluster default clusters are standalone Level 1 features.
In the rare case a cluster falls entirely within a Tier 1 LTR element (size threshold
makes this uncommon), make it Level 2.

```r
t1_ltr <- t1[t1$source_tool == "DANTE_LTR"]
t3_inside <- subsetByOverlaps(t3_def, t1_ltr, type="within")
t3_level1  <- t3_def[!t3_def$ID %in% t3_inside$ID]
level1  <- c(level1, t3_level1)
# Assign Parent to nested satellite features
t3_inside$Parent <- get_parent_id(t3_inside, t1_ltr)
level2 <- t3_inside
```

### Step 4 — Tier 3 short monomer: trim against default, then against Tier 1

```r
t3s_trimmed <- trim_to_nonoverlap(t3_sho, c(t3_level1, t1),
                                  min_len=opt$min_feature_length)
level1 <- c(level1, t3s_trimmed)
```

### Step 5 — Tier 4 RM on TideCluster: trim against Tiers 1–3

```r
higher_1_to_3 <- reduce(c(t1, t2_trimmed, t3_level1, t3s_trimmed))
t4_trimmed <- trim_to_nonoverlap(t4, higher_1_to_3, min_len=opt$min_feature_length)
level1 <- c(level1, t4_trimmed)
```

### Step 6 — Tier 5 TE hits: trim against Tiers 1–4

```r
higher_1_to_4 <- reduce(c(higher_1_to_3, t4_trimmed))
t5_trimmed <- trim_to_nonoverlap(t5_te, higher_1_to_4, min_len=opt$min_feature_length)
level1 <- c(level1, t5_trimmed)
```

### Step 7 — Tier 5 Simple_repeat / Low_complexity: nested where overlapping TE

```r
current_level1 <- reduce(c(higher_1_to_4, t5_trimmed))
# Simple/Low that do NOT overlap any annotation → Level 1
t5_sc_l1 <- subsetByOverlaps(t5_sc, current_level1, invert=TRUE)
# Simple/Low that DO overlap → Level 2 nested under the overlapping feature
t5_sc_l2 <- subsetByOverlaps(t5_sc, current_level1)
t5_sc_l2$Parent <- get_parent_id(t5_sc_l2, level1)
level1 <- c(level1, t5_sc_l1)
level2 <- c(level2, t5_sc_l2)
```

### Step 8 — Tier 6 TideHunter residuals

```r
all_higher <- reduce(c(current_level1, t5_sc_l1))
# TideHunter inside LTR-RT body → Level 2
t6_in_ltr <- subsetByOverlaps(t6, t1_ltr, type="within")
t6_in_ltr$Parent <- get_parent_id(t6_in_ltr, t1_ltr)
# TideHunter not covered by anything → Level 1
t6_l1 <- subsetByOverlaps(t6, all_higher, invert=TRUE)
level1 <- c(level1, t6_l1)
level2 <- c(level2, t6_in_ltr)
```

---

## 6. Key Helper Functions

### trim_to_nonoverlap()

Trims each feature in `lower` to its non-overlapping portion against `higher`.
Features fully covered by `higher` are dropped.
Resulting fragments shorter than `min_len` are dropped.

```r
trim_to_nonoverlap <- function(lower, higher, min_len=50) {
  if (length(lower) == 0 || length(higher) == 0) return(lower)
  hits <- findOverlaps(lower, reduce(higher))
  overlapping_idx <- unique(queryHits(hits))
  intact <- lower[-overlapping_idx]
  if (length(overlapping_idx) == 0) return(intact)

  higher_reduced <- reduce(higher)
  trimmed_list <- lapply(overlapping_idx, function(i) {
    feat <- lower[i]
    remaining <- setdiff(feat, higher_reduced)
    if (length(remaining) == 0) return(NULL)
    remaining <- remaining[width(remaining) >= min_len]
    if (length(remaining) == 0) return(NULL)
    # Propagate metadata from source feature
    for (col in c("Name","classification","source_tier","source_tool","ID")) {
      if (!is.null(mcols(lower)[[col]]))
        mcols(remaining)[[col]] <- mcols(lower)[[col]][1]
    }
    remaining
  })
  trimmed_list <- trimmed_list[!sapply(trimmed_list, is.null)]
  if (length(trimmed_list) == 0) return(intact)
  c(intact, do.call(c, trimmed_list))
}
```

### get_parent_id()

Finds the Level 1 feature with the largest overlap for each child feature and returns
its ID as the parent.

```r
get_parent_id <- function(children, parents) {
  hits <- findOverlaps(children, parents)
  if (length(hits) == 0) return(rep(NA_character_, length(children)))
  # For each child, pick the parent with maximum overlap
  ov_width <- width(pintersect(children[queryHits(hits)], parents[subjectHits(hits)]))
  best <- tapply(seq_along(hits), queryHits(hits), function(idx) idx[which.max(ov_width[idx])])
  parent_ids <- rep(NA_character_, length(children))
  for (child_i in as.integer(names(best))) {
    parent_ids[child_i] <- parents$ID[subjectHits(hits[best[[as.character(child_i)]]])]
  }
  parent_ids
}
```

### make_batches()

See annotation_rules.md §6.2 for the bin-packing logic.
Input: named integer vector of sequence lengths from FAI.
Output: list of character vectors (sequence name batches).

---

## 7. Output Assembly

After all batches are processed and concatenated:

```r
finalise_output <- function(level1, level2) {
  # 1. Sort
  level1 <- sort(sortSeqlevels(level1))
  level2 <- sort(sortSeqlevels(level2))

  # 2. Assign unique IDs
  level1$ID <- paste0(level1$classification, "_L1_",
                      formatC(seq_along(level1), width=8, flag="0"))
  level2$ID <- paste0(level2$classification, "_L2_",
                      formatC(seq_along(level2), width=8, flag="0"))

  # 3. Set GFF3 columns
  level1$type   <- ifelse(grepl("^Satellite", level1$classification),
                          "repeat_region", level1$type)
  level1$source <- level1$source_tool
  level2$type   <- level2$type  # preserve original
  level2$source <- level2$source_tool

  # 4. Combine and export
  all_features <- c(level1, level2)
  all_features <- sort(all_features)
  export(all_features, opt$output, format="gff3")
}
```

---

## 8. Attributes in Output GFF3

### Level 1 feature (mandatory)

```
ID=<classification>_L1_<zero-padded-integer>
Name=<classification>
classification=<Class_I/LTR/Ty3_Gypsy or Satellite/TideCluster/TRC_5 etc.>
source_tier=<1..6>
source_tool=<DANTE_LTR|DANTE_TIR|DANTE_LINE|DANTE|TideCluster_default|TideCluster_short|TideCluster_RM|RepeatMasker|TideHunter>
```

### Level 1 feature (optional, retained when present)

```
# From DANTE_LTR top-level:
ltr_identity=<float>
ltr5_length=<int>
ltr3_length=<int>
tsd=<string>
ndomains=<int>
element_type=<complete|partial>
rank=<DLTP|DL|D|...>

# From DANTE_TIR:
tir_seq5=<string>
tir_seq3=<string>
tsd=<string>
cluster_id=<int>

# From DANTE_LINE:
pattern_type=<ENDO-RT-RH|ENDO-RT|...>
extension_5prime=<int>
extension_3prime=<int>
```

### Level 2 feature

```
ID=<classification>_L2_<zero-padded-integer>
Name=<classification>
classification=<...>
Parent=<ID of the Level 1 parent>
source_tier=<1..6>
source_tool=<...>
```

---

## 9. Snakemake Rule

```python
rule make_unified_annotation:
    input:
        ltr        = F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3",
        tir        = F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_final.gff3",
        line       = F"{config['output_dir']}/DANTE_LINE/DANTE_LINE.gff3",
        dante      = F"{config['output_dir']}/DANTE/DANTE_filtered.gff3",
        tc_default = F"{config['output_dir']}/TideCluster/default/TideCluster_clustering.gff3",
        tc_short   = F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_clustering.gff3",
        tc_rm      = F"{config['output_dir']}/TideCluster/default/RM_on_TideCluster_Library.gff3",
        rm         = F"{config['output_dir']}/RepeatMasker/RM_on_combined_library_plus_DANTE.gff3",
        th_default = F"{config['output_dir']}/TideCluster/default/TideCluster_tidehunter_short.gff3",
        th_short   = F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_tidehunter_short.gff3",
        fai        = genome_fasta_cleaned + ".fai"
    output:
        F"{config['output_dir']}/Repeat_Annotation_Unified.gff3"
    log:
        stdout = F"{config['output_dir']}/logs/make_unified_annotation.log",
        stderr = F"{config['output_dir']}/logs/make_unified_annotation.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/make_unified_annotation.tsv"
    conda:
        "envs/tidecluster.yaml"
    threads: workflow.cores
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        make_unified_annotation.R \
            --ltr        {input.ltr}        \
            --tir        {input.tir}        \
            --line       {input.line}       \
            --dante      {input.dante}      \
            --tc_default {input.tc_default} \
            --tc_short   {input.tc_short}   \
            --tc_rm      {input.tc_rm}      \
            --rm         {input.rm}         \
            --th_default {input.th_default} \
            --th_short   {input.th_short}   \
            --fai        {input.fai}        \
            --output     {output}           \
            --threads    {threads}
        """
```

Add `Repeat_Annotation_Unified.gff3` to `rule all` and `BENCHMARKED_RULES`.

---

## 10. Testing Plan

### 10.1 Data

Test genome: `tmp/output/genome_cleaned.fasta` (~316 Mb, ~16 chromosomes from the FAI).
All required input GFF3s are present in `tmp/output/`.

### 10.2 Test command

```bash
RSCRIPT=/home/petr/PycharmProjects/assembly_repeat_annotation_pipeline/.snakemake/conda/1fe3a845921ac5cb77d600d3e9b79fa5_/bin/Rscript

conda run -p /home/petr/PycharmProjects/assembly_repeat_annotation_pipeline/.snakemake/conda/1fe3a845921ac5cb77d600d3e9b79fa5_ \
  Rscript scripts/make_unified_annotation.R \
    --ltr        tmp/output/DANTE_LTR/DANTE_LTR.gff3 \
    --tir        tmp/output/DANTE_TIR/DANTE_TIR_final.gff3 \
    --line       tmp/output/DANTE_LINE/DANTE_LINE.gff3 \
    --dante      tmp/output/DANTE/DANTE_filtered.gff3 \
    --tc_default tmp/output/TideCluster/default/TideCluster_clustering.gff3 \
    --tc_short   tmp/output/TideCluster/short_monomer/TideCluster_clustering.gff3 \
    --tc_rm      tmp/output/TideCluster/default/RM_on_TideCluster_Library.gff3 \
    --rm         tmp/output/RepeatMasker/RM_on_combined_library_plus_DANTE.gff3 \
    --th_default tmp/output/TideCluster/default/TideCluster_tidehunter_short.gff3 \
    --th_short   tmp/output/TideCluster/short_monomer/TideCluster_tidehunter_short.gff3 \
    --fai        tmp/output/genome_cleaned.fasta.fai \
    --output     tmp/output/Repeat_Annotation_Unified.gff3 \
    --threads    4
```

### 10.3 Sanity checks

After running, verify:

1. **No overlapping Level 1 features** (same seqname):
   ```r
   gr <- import.gff3("tmp/output/Repeat_Annotation_Unified.gff3")
   l1 <- gr[is.na(gr$Parent)]
   stopifnot(length(reduce(l1)) == length(l1))
   ```

2. **Total coverage ≤ genome size** (Level 1 only):
   ```r
   fai <- read.table("tmp/output/genome_cleaned.fasta.fai")
   genome_bp <- sum(fai$V2)
   l1_bp <- sum(width(reduce(l1)))
   stopifnot(l1_bp <= genome_bp)
   ```

3. **All Level 2 features have a valid Parent**:
   ```r
   l2 <- gr[!is.na(gr$Parent)]
   stopifnot(all(l2$Parent %in% l1$ID))
   ```

4. **Classification coverage comparison** (Level 1 bp per class vs existing
   `summary_statistics.csv`) — values should be in the same ballpark, not identical
   because the unified output uses different priority rules.

5. **GFF3 validity** (if `gt` is available):
   ```bash
   gt gff3 -tidy -sort -retainids tmp/output/Repeat_Annotation_Unified.gff3 > /dev/null
   ```

---

## 11. Implementation Order

1. Write `load_*()` functions and test each individually against the test data,
   printing feature counts and checking attribute availability.
2. Write `trim_to_nonoverlap()` and `get_parent_id()` helpers; unit-test on small
   synthetic GRanges examples.
3. Write `make_batches()` and verify batching logic on the test FAI.
4. Implement `process_batch()` using Steps 1–8 from §5 for a single chromosome first.
5. Wrap in `BiocParallel::bplapply` over batches.
6. Implement `finalise_output()` and GFF3 export.
7. Run full test (§10.2) and apply all sanity checks (§10.3).
8. Add Snakemake rule; dry-run; full pipeline run.
