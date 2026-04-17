# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is an Assembly Repeat Annotation Pipeline that uses multiple bioinformatics tools to comprehensively identify and classify repetitive elements in genomic sequences. The pipeline is implemented as a Snakemake workflow and is designed to run in a Singularity container environment.

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
- Optional parameters: `custom_library`, `tandem_repeat_library`, `repeatmasker_sensitivity`, `reduce_library`

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

## Important Architecture Details

### Conda Environment Management
- Environments defined in `envs/` directory (tidecluster.yaml, bedtools.yaml, dante_tir.yaml, etc.)
- Pipeline automatically manages conda environments through Snakemake
- Main computational environment: `tidecluster.yaml` (contains DANTE, DANTE_LTR, TideCluster tools)
- Specialized environment: `dante_tir.yaml` (contains DANTE_TIR tool for TIR transposon detection)

### DANTE_TIR Fallback Workflow
- `dante_tir_fallback` (`scripts/dante_tir_fallback.py`): Identifies partial TIR elements from TPase domain flanking-region analysis that the primary `dante_tir` tool may have missed
- `merge_dante_tir_with_fallback` (`scripts/merge_tir_fallback.py`): Removes fallback elements that overlap any primary DANTE_TIR element (any coordinate overlap = discard), labels surviving elements with `source=DANTE_TIR_fallback`, `Status=partial`, and `_partial` ID suffix
- `make_tir_combined_library`: Concatenates primary DANTE_TIR element sequences with non-overlapping fallback sequences and re-clusters everything with mmseqs2 to produce a unified representative library
- The combined GFF3 (`DANTE_TIR_combined.gff3`) and library (`all_representative_elements_combined.fasta`) replace the original DANTE_TIR outputs in all downstream rules (unified annotation, library construction, RepeatMasker)
- Conda environment: `dante_line.yaml` (seqkit, mmseqs2, parasail-python)

### Singularity container sync

The pipeline ships as a Singularity container; files that the pipeline reads at run time must be listed in the `%files` section of `Singularity`, otherwise the container builds cleanly but fails at run time.

This rule is enforced by `.githooks/pre-commit`. Enable it once per checkout:

```bash
git config core.hooksPath .githooks
```

With the hook active, `git commit` fails if a new file under a runtime-relevant path is staged without a matching `%files` update. Dev-only paths (tests, docs, CI config, figures) are exempt via `EXCLUDE_PATTERNS` in the hook. Bypass with `git commit --no-verify` if you must.

Current `%files` coverage: `envs/`, `Snakefile`, `config.yaml`, `classification_vocabulary.yaml`, `data/rdna_library.fasta`, `run_pipeline.py`, `scripts/`.

### Classification handling

- Canonical classification form is slash-separated (e.g. `Class_I/LTR/Ty1_copia/Ale`). Underscores inside a path component are always part of the name, never separators (`Ty1_copia`, `Tc1_Mariner`, `EnSpm_CACTA`, `Class_II`, `Subclass_1`).
- `classification_vocabulary.yaml` is the single source of truth. Every valid classification is listed there, together with tool-native aliases (e.g. REXdb `Ty1/copia` → `Ty1_copia`) and the DANTE_TIR underscore-prefix → slash-prefix mapping.
- Normalisation is centralised in `scripts/classification.py` and the R mirror `scripts/classification.R`. Both modules read the same YAML and pass the same shared test vector `tests/classification_cases.tsv`.
- **Never sed / awk / regex-convert classifications in Snakefile or scripts.** Call `canonicalise(x, source="DANTE_TIR" | "DANTE_LTR" | "DANTE" | "RepeatMasker" | ...)` instead. Shell rules invoke the module via `python3 scripts/classification.py canonicalise-fasta-headers --source <TOOL>` (or the `canonicalise-gff3-attribute` / `validate` subcommands).
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
- Classification splits: Mobile_elements, Simple_repeats, Low_complexity, rDNA, etc.

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
