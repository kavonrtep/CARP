# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is an Assembly Repeat Annotation Pipeline that uses multiple bioinformatics tools to comprehensively identify and classify repetitive elements in genomic sequences. The pipeline is implemented as a Snakemake workflow and is designed to run in a Singularity container environment.

## Upstream tools reference

The pipeline wraps five upstream tools built around the REXdb protein-domain
classification: **DANTE** (domain annotation), **DANTE_LTR** (intact LTR-RTs
from DANTE output), **DANTE_TIR** (TIR DNA transposons from DANTE TPase
domains), **TideCluster** (tandem repeats), and **REXdb** (the reference
database that backs all DANTE-family tools). A consolidated reference
covering each tool's algorithm, CLI, GFF3 attributes, classification
dialect, and the relevant results from the two foundational papers
(Novák et al. 2024, Neumann et al. 2019) lives in
[`context/tools_summary.md`](context/tools_summary.md). Read it before
changing any rule that produces or consumes these tools' outputs.

## Core Pipeline Architecture

The pipeline follows a multi-stage workflow:

1. **DANTE/DANTE_LTR Analysis**: Identifies intact LTR retrotransposons using protein domain-based detection
2. **DANTE_TIR Analysis**: Identifies DNA transposons with Terminal Inverted Repeats (TIRs) from DANTE domain annotations
3. **TideCluster Analysis**: Detects tandem repeats with two separate runs (default and short monomer)
4. **Library Construction**: Creates custom repeat libraries from discovered elements, with filtering to remove composite elements
5. **RepeatMasker Annotation**: Performs similarity-based annotation using the constructed libraries
6. **Output Processing**: Generates GFF3 files, density tracks (BigWig), and summary statistics

## Key Components

### Main Entry Points
- `run_pipeline.py`: Main execution script that wraps Snakemake execution
- `Snakefile`: Core workflow definition with all rules and dependencies
- `Singularity`: Container definition for reproducible execution environment

### Configuration
- Config files follow pattern `config*.yaml` (e.g., `config.yaml`, `config_minimal.yaml`)
- Required parameters: `genome_fasta`, `output_dir`
- Optional parameters: `custom_library`, `tandem_repeat_library`, `repeatmasker_sensitivity`, `reduce_library`, `reduce_tidecluster_library`, `reduce_library_containment` (+ `containment_min_identity`, `containment_min_coverage`), and the TideCluster 1.16.0 knobs `tidecluster_detect_rdna` (default True), `tidecluster_rdna_library` (default empty → bundled), `tidecluster_keep_trc_overlaps` (default False), `tidecluster_chunk_size` (default 50_000_000)
- `reduce_library` (main RepeatMasker/CAP3 library reduction) and `reduce_tidecluster_library` (rotation-invariant reduction of the TideCluster consensus *dimer* library used by `tidecluster_reannotate`, default True) are **independent** flags. The dimer reducer (`scripts/reduce_dimer_library.py`) aligns each monomer against the already-doubled dimers (`mmseqs easy-search`, `--cov-mode 2`) so rotational phase variants collapse per-TRC — ~5–20× smaller dimer library, faster remasking. Clustering is **greedy/star** (each dropped dimer aligns directly to a kept rep, not transitively chained) with a **length guard** (only near-equal-length dimers merge; nested/harmonic periods stay separate) so the kept rep can mask everything it replaces. It runs `mmseqs easy-search` with an explicit small nucleotide k-mer (`-k 7`): the mmseqs **default k=15 segfaults on short tandem monomers** (the query is the first-half monomer, often ≤21 bp), which is the common case here; a small k only raises prefilter sensitivity, with specificity still gated by `--min-seq-id`/coverage. Two safety nets keep a pathological group from failing the run: groups whose shortest monomer is below the prefilter floor (<12 bp) are skipped, and any residual `mmseqs` error triggers a single-threaded retry then an unreduced fallback. **Known limitation (accepted):** reduction is lossless for small/medium tandem families but can under-mask *large satellite families* (dimer ≳ a few kb) — `tidecluster_reannotate` keeps only masked regions ≥ the family's max dimer length, and RepeatMasker needs the full consensus diversity to tile such big arrays contiguously, so thinning them drops some masking (measured ~15% of reannotate masked bp on a 17 kb-dimer megasatellite). The structural TideCluster clustering still annotates those arrays directly, so the net effect on the unified annotation is limited; set `reduce_tidecluster_library: False` to mask with the full dimer library.
- `reduce_library_containment` (default True, **independent** of `reduce_library`) adds a **second-round containment reduction** of the per-class-reduced RepeatMasker library, in the `reduce_library_containment` rule (`reduce_library → reduce_library_containment → repeatmasker`). The per-class CAP3/mmseqs `reduce_library` leaves many short fragments fully contained in a longer element of the same class; `scripts/containment_reduce_library.py` (greedy longest-first blastn self-comparison) drops a fragment when a retained, strictly-longer, **same-class** sequence covers ≥ `containment_min_coverage` (0.90) of it at ≥ `containment_min_identity` (80) % identity. RepeatMasker masks those fragments' genomic copies via the container, so masking **and** classification are preserved. Validated masked-bp-lossless with RepeatMasker on the Pisum pangenome over two genome regions (−0.09…−0.12%, within the ±0.15% bar) while cutting ~22% of library bp and ~30% of RepeatMasker wall-time. blastn-unavailable / failure copies the library through unchanged (the reduction is an optimisation, never a correctness requirement). Note: *consensus modeling* of the abundant families was tested and rejected — at comparable reduction it crosses the losslessness bar (the conserved-core redundancy is partial/local and can't be merged without dropping the unique flanks RepeatMasker needs); whole-sequence containment is the lossless win.
- **Unified annotation: rDNA labelling + TR-from-structural-TE resolution** (`scripts/make_unified_annotation.R`, consumes TideCluster 1.16.0 output). Two behaviours, both judged inside the tier resolver:
  - **rDNA (array-level).** TideCluster 1.16.0 flags rDNA TRCs with `rDNA_type=45S|5S` in the clustering GFF3; `normalise_tc_satellite()` puts the array-level label `rDNA_45S` / `rDNA_5S` in the **`classification`** of those arrays (no internal 18S/ITS/5.8S/IGS/25S substructure — that detail still comes from the RepeatMasker `data/rdna_library.fasta` layer where present). **Invariant: the satellite `Name` is ALWAYS the bare `TRC_<n>`** for every TideCluster satellite (rDNA included) — downstream apps and `split_gff_by_name.R --name-prefix TRC_` (the per-family BigWig rule) key on it, so it must stay stable across releases. A plain-satellite record is therefore identical to the previous release; rDNA only changes `classification`. Routing: `calculate_statistics_and_make_groups.R` relabels `Name` **from `classification`** in-memory (`^rDNA`) so rDNA arrays count as **rDNA** (not folded into `Tandem_repeats`); on disk `Name` stays `TRC_<n>`. Because the rDNA Name stays `TRC_<n>`, rDNA arrays also still appear in the per-family tandem split (rDNA is both a tandem family and rDNA).
  - **TR-from-structural-TE (`TE_origin`).** A length-qualified Tier-3 *clustering* TRC (never TideHunter residuals, never RepeatMasker) that tandemly stacks **multiple same-family structural TEs** (Tier 1: DANTE_LTR/TIR/LINE) is a TE-derived satellite (e.g. OZ408687.1 TRC_13 over an Ale array). Policy (`identify_te_derived_trcs()` pre-pass in `process_batch`): the satellite **wins** the region, is tagged `TE_origin=<LCA class of the covered structural TEs>`, and the covered structural TEs are **trimmed out of the unified file** (they remain in `DANTE_*.gff3`). Trigger thresholds are constants in the script: `TE_ORIGIN_MIN_ELEMENTS=2`, `TE_ORIGIN_MIN_LCA_DEPTH=3` (shared lineage ≥ `Class_I/LTR/Ty1_copia`), `TE_ORIGIN_MIN_COVERAGE=0.25` (fraction of the TRC's array bp covered by those TEs — calibrated on OZ408687.1 TRC_13 at 0.35 vs the genome's genuine large satellites TRC_1/TRC_3 at ~0.00; DANTE annotates only *complete* elements, so a degraded tandem-TE array is only partially covered). This is the only place a satellite outranks a Tier-1 structural call. **Two sub-types** (the satellite↔tandem-LTR-RT conflict): when the TE-derived satellite overlaps a tandem `LTR_RT_TR` container — a tandem of *complete* LTR-RTs — it is tagged `TE_origin_structure=tandem_LTR_RT`, and the container **plus its member copies** (Level-2) are dropped along with the standalone TEs (the pre-pass trims `t1`, `t2`, **and `t1_members`** against `te_sat`; dropping members is essential — otherwise a member whose container is trimmed orphans into a malformed feature). A satellite over scattered complete TEs (no container) gets `TE_origin` without the sub-type. A purely degraded array with no complete elements is not flagged here (structural trigger only) and stays a plain `Satellite`. Conversely, every **non**-TE-derived clustering satellite (default Step 3, short Step 4) is trimmed against everything already placed (`te_sat` + Tier 1 + Tier 2) by standard tier priority, so a single inserted/abutting structural TE wins its contested span and the satellite keeps the rest (within-LTR clusters still nest as Level 2).
  - **Verification (non-fatal).** `make_unified_annotation` also writes `Repeat_Annotation_Unified.overlaps.tsv` — every Level-1 vs Level-1 overlap where neither partner is `Simple_repeat`/`Low_complexity`. With the Step-3/4 trimming above it is **empty** on the OZ408 calibration genome (no satellite/rDNA-vs-structural-TE rows); anything residual (e.g. same-tool partial-TE overlaps on other genomes) is surfaced (warned, not failed).

### Critical Workflow Rules
- `dante` → `dante_ltr` → `make_library_of_ltrs`: LTR retrotransposon discovery and library creation
- `dante` → `dante_tir` → `dante_tir_fallback` → `merge_dante_tir_with_fallback` → `make_tir_combined_library`: DNA transposon identification with TIR fallback for partial elements
- `tidecluster_long` + `tidecluster_short`: Tandem repeat detection with different parameters
- `filter_ltr_rt_library`: Removes Class_II elements from LTR library to prevent misannotation
- `repeatmasker` → `merge_rm_and_dante`: Final repeat annotation combining all approaches
- `make_summary_statistics_and_split_by_class`: Generates classification and statistics

## Common Development Commands

### Running the Pipeline
```bash
# Standard execution with config file
python run_pipeline.py -c config.yaml -t 20

# With Singularity container
singularity run -B /path/to/data -B $PWD assembly_repeat_annotation_pipeline.sif -c config.yaml -t 20

# Dry run to check workflow
python run_pipeline.py -c config.yaml -t 20 -S="--dry-run"
```

### Development and Testing
```bash
# Check Snakemake syntax
snakemake --lint

# Generate workflow visualization
snakemake --dag | dot -Tsvg > workflow.svg
snakemake --rulegraph | dot -Tsvg > rulegraph.svg

# Run specific rules for testing
snakemake --use-conda target_file.gff3 -c config.yaml --cores 4
```

### Container Development
```bash
# Build Singularity container
sudo singularity build assembly_repeat_annotation_pipeline.sif Singularity

# Test container interactively
singularity shell -B $PWD assembly_repeat_annotation_pipeline.sif
```

### Releasing
`version.py` is the single source of truth; a release is a version bump + an
**unprefixed** tag equal to it (`1.0.0rc1`, not `v1.0.0rc1`). Use the `/release`
skill or run its helper directly — it validates, bumps `version.py`, runs the
cheap CI gates, commits and tags (no push):
```bash
.claude/skills/release/cut-release.sh <VERSION>   # e.g. 1.0.0rc1
```
Then push from the **host** (the sandbox has no ssh): `git push origin main && git push origin <VERSION>`.
The tag push drives `.github/workflows/release.yml` (SIF build → in-container
fixture → GHCR → GitHub Release → Zenodo). See `.claude/skills/release/SKILL.md`.

## Important Architecture Details

### Conda Environment Management
- Environments defined in `envs/` directory (tidecluster.yaml, tidecluster_run.yaml, bedtools.yaml, dante_tir.yaml, etc.)
- Pipeline automatically manages conda environments through Snakemake
- Main computational environment: `tidecluster.yaml` (contains DANTE, DANTE_LTR, RepeatMasker, and the R/report stack — used by ~18 rules). Despite the file name it no longer ships TideCluster.
- TideCluster environment: `tidecluster_run.yaml` (TideCluster 1.16.0 + bundled TideHunter/mmseqs2/RepeatMasker/rtracklayer) — used only by the `tidecluster_long`, `tidecluster_short`, and `tidecluster_reannotate` rules. Split out because TideCluster ≥1.15 needs r-base ≥4.2 (via r-igraph) while DANTE pins r-base 4.1.3, so the two cannot share one conda env. **1.16.0 adds three default-on behaviours the pipeline consumes** (see the unified-annotation bullet below): rDNA TRC identification (`rDNA_type=45S|5S` in the clustering GFF3), cross-TRC overlap resolution (clustering GFF3 made non-overlapping across satellite TRCs), and a chunked/pooled `tc_reannotate` RepeatMasker that fixes the `-pa`-with-custom-`-lib` under-parallelism (Dfam #274; `--chunk_size` default 50 Mb). All four `tidecluster_*` config knobs default to TideCluster's own defaults, so a default run adds no extra flags.
- Specialized environment: `dante_tir.yaml` (contains DANTE_TIR tool for TIR transposon detection)

### DANTE_TIR Fallback Workflow
- `dante_tir_fallback` (`scripts/dante_tir_fallback.py`): Identifies partial TIR elements from TPase domain flanking-region analysis that the primary `dante_tir` tool may have missed
- `merge_dante_tir_with_fallback` (`scripts/merge_tir_fallback.py`): Removes fallback elements that overlap any primary DANTE_TIR element (any coordinate overlap = discard), labels surviving elements with `source=DANTE_TIR_fallback`, `Status=partial`, and `_partial` ID suffix
- `make_tir_combined_library`: Builds the per-genome TIR library from **primary DANTE_TIR elements only**. Fallback elements remain visible in `DANTE_TIR_combined.gff3` (and the unified annotation) as low-confidence partials but are not trusted enough to seed RepeatMasker by default
- `build_fallback_tir_library` (`scripts/build_fallback_tir_library.py`, **optional**, default OFF): When `include_dante_tir_fallback_in_library: true` is set in the config, re-clusters the post-overlap fallback survivors, applies a Multiplicity floor (`dante_tir_fallback_library_min_multiplicity`, defaults to inheriting `dante_tir_min_multiplicity` = 3), and runs a strict class-aware blastn filter against the union of LTR / DANTE_TIR primary / LINE / custom libraries. A fallback rep is dropped if any blast hit (e-value 1e-19) lands on a subject whose canonical classification is not on the rep's lineage chain — siblings such as CACTA-vs-hAT count as incompatible. The kept reps are appended to `combined_library.fasta`. The fallback library is intentionally NOT routed through `make_subclass_2_library` (so it does not contribute to filtering the LTR library); the policy is that the fallback layer is treated as less reliable than the primary library and a misclassified entry must not damage the LTR side of the annotation
- The combined GFF3 (`DANTE_TIR_combined.gff3`) and library (`all_representative_elements_combined.fasta`, plus the optional `fallback_library.fasta`) feed downstream rules (unified annotation, library construction, RepeatMasker)
- Audit log of fallback-library decisions: `DANTE_TIR/fallback_library_dropped.tsv` (rep_id, rep_class, status, conflicting_subject_id, conflicting_subject_class). Read this when chasing an annotation shift introduced by the fallback library
- Conda environment: `dante_line.yaml` (seqkit, mmseqs2, parasail-python, blast)

### Singularity container sync

The pipeline ships as a Singularity container; files that the pipeline reads at run time must be listed in the `%files` section of `Singularity`, otherwise the container builds cleanly but fails at run time.

This rule is enforced by `.githooks/pre-commit`. Enable it once per checkout:

```bash
git config core.hooksPath .githooks
```

With the hook active, `git commit` fails if a new file under a runtime-relevant path is staged without a matching `%files` update. Dev-only paths (tests, docs, CI config, figures) are exempt via `EXCLUDE_PATTERNS` in the hook. Bypass with `git commit --no-verify` if you must.

Current `%files` coverage: `envs/`, `Snakefile`, `config.yaml`, `classification_vocabulary.yaml`, `data/rdna_library.fasta`, `run_pipeline.py`, `scripts/`.

### Calling helper scripts from rules

Rules in the Snakefile must invoke pipeline helper scripts using the
**dual-context PATH-prepend pattern**, so the same rule works whether
snakemake runs from the container (PBS / HPC deployments) or directly
from a repo checkout (CI, local development).

```
shell:
    """
    exec > {log.stdout} 2> {log.stderr}
    set -euo pipefail
    set -x
    scripts_dir=$(realpath scripts)
    export PATH="$scripts_dir:$PATH"
    my_helper.py --foo {input.x}
    """
```

Why this works in both contexts:

- **Local / CI:** `realpath scripts` resolves to `<repo>/scripts/` (cwd is
  the repo root); the prepend puts the local checkout's scripts on PATH;
  bare-name lookup hits them.
- **Container (`singularity run`):** the launcher cd's to `$SCRATCHDIR`,
  which has no `scripts/` directory. `realpath scripts` returns a
  non-existent path. The prepend is harmless — `bash` walks past the
  missing dir and falls through to `/opt/pipeline/scripts`, already on
  PATH from `Singularity` `%environment`. Bare-name lookup finds the
  script there.

Two anti-patterns. Both broke production runs in April 2026:

- ❌ `python3 "$scripts_dir/foo.py"` — explicit path. When `scripts_dir`
  resolves to a non-existent scratch dir (Metacentrum case), the open
  fails. Use bare name and let PATH resolve.
- ❌ `scripts_dir=$(realpath {workflow.basedir}/scripts)` — works today
  because `workflow.basedir` substitutes to the Snakefile's directory in
  both contexts, but it's a different convention from every other rule
  in the file. A future contributor copying that idiom into a new rule
  may pair it with the explicit-path anti-pattern and reproduce the
  Metacentrum failure. Stay on the canonical `realpath scripts` form.

For a helper script to be reachable by bare name it must (a) live under
`scripts/`, (b) carry a `#!/usr/bin/env python3` (or `Rscript`) shebang,
and (c) have the executable bit set. The `Singularity` `%files` copy
preserves permissions, so `chmod +x` in the repo flows through to the
container. New helper scripts must be `chmod +x`'d before the first
commit.

### Classification handling

- Canonical classification form is slash-separated (e.g. `Class_I/LTR/Ty1_copia/Ale`). Underscores inside a path component are always part of the name, never separators (`Ty1_copia`, `Tc1_Mariner`, `EnSpm_CACTA`, `Class_II`, `Subclass_1`).
- `classification_vocabulary.yaml` is the single source of truth. Every valid classification is listed there, together with tool-native aliases (e.g. REXdb `Ty1/copia` → `Ty1_copia`) and the DANTE_TIR underscore-prefix → slash-prefix mapping.
- Normalisation is centralised in `scripts/classification.py` and the R mirror `scripts/classification.R`. Both modules read the same YAML and pass the same shared test vector `tests/classification_cases.tsv`.
- **Never sed / awk / regex-convert classifications in Snakefile or scripts.** Call `canonicalise(x, source="DANTE_TIR" | "DANTE_LTR" | "DANTE" | "RepeatMasker" | ...)` instead. Shell rules invoke the module via `classification.py canonicalise-fasta-headers --source <TOOL>` (or the `canonicalise-gff3-attribute` / `validate` subcommands), using the bare-name + PATH-prepend pattern from "Calling helper scripts from rules" above.
- The `validate_classifications` Snakemake rule runs after every upstream tool produces its GFF3; it fails the pipeline before RepeatMasker / unified annotation start if any classification is not in the vocabulary. A new tool release emitting an unknown leaf fails loudly here — the fix is a one-line edit to `classification_vocabulary.yaml`.
- Test both parsers: `python3 tests/test_classification.py && Rscript tests/test_classification.R`.

### Library Filtering Logic
- Critical step: `filter_ltr_rt_library` rule prevents misannotation by removing LTR sequences similar to Class_II elements
- Class_II/Subclass_1 elements from custom library are used to screen LTR library
- This addresses composite retrotransposons containing inserted DNA transposons

### Output Structure
- Main outputs organized in subdirectories: `DANTE/`, `DANTE_TIR/`, `DANTE_TIR_FALLBACK/`, `DANTE_LTR/`, `TideCluster/`, `Libraries/`, `RepeatMasker/`
- Top-level symlinks created for easy access to key results (including `DANTE_TIR.gff3` → `DANTE_TIR_combined.gff3`)
- BigWig density tracks generated at 10kb and 100kb windows
- Classification splits: Mobile_elements, Simple_repeats, Tandem_repeats, Low_complexity, rDNA, etc.
- `summary_statistics.csv` and the per-class split GFFs in `Repeat_Annotation_NoSat_split_by_class_gff3/` are computed from `Repeat_Annotation_Unified.gff3` (not `Repeat_Annotation_NoSat.gff3`). Unified is non-overlapping by tier-priority, so DANTE-direct calls outside any RepeatMasker hit count toward the user-visible totals. The directory name keeps the historical `_NoSat_` prefix for backward compat
- The `Tandem_repeats` aggregation row replaces the older `Satellites` row; it spans TideCluster default + short clusters, TideHunter short tandems, and any RepeatMasker `Satellite/*` calls. Defined in `classification_vocabulary.yaml::aggregation_buckets`

### Parallel Processing Strategy
- TideCluster runs with default and short monomer parameters separately then merges
- RepeatMasker uses the ncbi engine
- BigWig calculation optimized for batch processing of classification groups

### Custom Library Integration
- Supports user-provided repeat libraries with specific naming convention: `>unique_id#family/subfamily`
- Enforces standard classification scheme (Class_I/Class_II hierarchy)
- Optional tandem repeat library for TideCluster similarity-based annotation

## File Organization Patterns

- Scripts in `scripts/` directory (mix of R and Python)
- Data files expected in `data/` directory
- All outputs directed to configurable `output_dir`
- Intermediate files organized by tool (DANTE, TideCluster, etc.)

## Testing make_repeat_report.R

Test data lives in `tmp/output/`. The R script requires packages (`rtracklayer`, `IRanges`, `GenomicRanges`, `jsonlite`, `optparse`) that are only available inside the Snakemake-managed conda environments. The correct env (hash may differ) is identified by having all required packages:

```bash
# Find the tidecluster env with rtracklayer and jsonlite
RSCRIPT=/home/petr/PycharmProjects/assembly_repeat_annotation_pipeline/.snakemake/conda/1fe3a845921ac5cb77d600d3e9b79fa5_/bin/Rscript

# Run on test data
conda run -p /home/petr/PycharmProjects/assembly_repeat_annotation_pipeline/.snakemake/conda/1fe3a845921ac5cb77d600d3e9b79fa5_ \
  Rscript scripts/make_repeat_report.R --output_dir tmp/output
```

If the env hash changed, find the right one with:
```bash
for d in .snakemake/conda/*/bin/; do
  [ -f "${d}Rscript" ] && ${d}Rscript -e "requireNamespace('jsonlite',quietly=TRUE) && requireNamespace('rtracklayer',quietly=TRUE) && cat('OK\n')" 2>/dev/null
done
```

Note: `jsonlite` may need to be manually installed into the env if it's missing (`r-jsonlite` was added to `envs/tidecluster.yaml` for future builds).


Note: `hermit` directory is not part of the project - it is setup for agent managing this project from the container - sandbox.
