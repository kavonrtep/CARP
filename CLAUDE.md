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
