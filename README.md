<img src="figs/CARP.png" alt="CARP logo" width="400"/>

The **Comprehensive Annotation of Repeats Pipeline (CARP)** integrates several
specialised tools to identify and classify the repetitive elements of a genome
assembly, then merges everything into a single, non-overlapping annotation.

- **DANTE / DANTE_LTR** — intact LTR retrotransposons
- **DANTE_TIR** — TIR DNA transposons
- **DANTE_LINE** — LINE elements (experimental)
- **TideCluster** — tandem repeats / satellites and rDNA arrays
- A custom repeat **library** is built from the discovered elements (optionally
  extended with user databases) and used by **RepeatMasker** for genome-wide
  similarity masking
- Outputs: a unified GFF3 annotation, per-class GFF3 files, density tracks
  (BigWig), summary statistics, and HTML reports

> **Limitation:** CARP was developed for **plant** genomes. Do not use it on
> animal genomes.

## Tools

| Tool | Role in CARP |
|------|--------------|
| [**DANTE**](https://github.com/kavonrtep/dante) | Annotates conserved transposon protein domains against REXdb; foundation for the DANTE_* tools |
| [**DANTE_LTR**](https://github.com/kavonrtep/dante_ltr) | Reconstructs intact LTR retrotransposons from DANTE domains |
| [**DANTE_TIR**](https://github.com/kavonrtep/dante_tir) | Calls TIR DNA transposons from DANTE transposase domains |
| [**DANTE_LINE**](https://github.com/kavonrtep/dante_line) | Calls LINE elements (experimental) |
| [**TideCluster**](https://github.com/kavonrtep/TideCluster) | Detects tandem repeats in two passes (default + short monomer) and flags 45S/5S rDNA arrays |
| [**RepeatMasker**](https://www.repeatmasker.org/) | Similarity masking of the genome with the CARP-built repeat library |
| [**REXdb**](https://github.com/repeatexplorer/rexdb) | Viridiplantae protein-domain reference database backing the DANTE family |

## Example output

Explore a live example of CARP's output — the
[**interactive HTML report**](https://kavonrtep.github.io/carp-example-report/),
produced on the *Boechera stricta* genome
([GCA_018361405.1](https://www.ebi.ac.uk/ena/browser/view/GCA_018361405.1)) — to
see what the pipeline produces (genome-wide repeat composition, the LTR / TIR
structural summaries, and the tandem-repeat / satellite family drill-down)
without running it yourself. The page hosts the HTML reports only; a full run
also emits the complete GFF3 annotation, density tracks (BigWig) and masked
sequence.

## Requirements

- **Container runtime:** [Apptainer](https://apptainer.org/) ≥ 1.0 (recommended)
  or [Singularity](https://sylabs.io/singularity/) ≥ 3.7 (needed for `oras://`
  pulls). The two share the same CLI; either works at run time.
- **Host OS:** Linux x86_64. The image ships GNU/Linux binaries; macOS and
  Windows can run it only inside a Linux VM.

If no runtime is installed system-wide, install one with conda:

```bash
conda create -n apptainer -c conda-forge "apptainer>=1.0"
conda activate apptainer
```

## Install

Each tagged release publishes the same image to GHCR (as an ORAS/OCI artefact)
and to a GitHub Release. Pull whichever is convenient:

```bash
# From GHCR (recommended):
apptainer pull oras://ghcr.io/kavonrtep/carp/sif:1.0.0rc3
# or always-latest:
apptainer pull oras://ghcr.io/kavonrtep/carp/sif:latest

# Alternative: the .sif attached to a GitHub Release
# https://github.com/kavonrtep/assembly_repeat_annotation_pipeline/releases
```

Tags follow PEP 440 unprefixed (`1.0.0rc3`, `1.0.0`). Both sources receive the
*same* image; publication is gated by an in-container fixture run, so any
published tag has passed the test fixture end-to-end.

## Usage

### 1. Write a `config.yaml`

```yaml
genome_fasta: data/genome.fasta        # required: assembly to annotate
output_dir: output                     # required: results directory
custom_library: data/custom_lib.fasta  # optional: extra repeat library
tandem_repeat_library: data/TR_lib.fasta  # optional: tandem-repeat reference
repeatmasker_sensitivity: default      # rush | default | quick
reduce_library: True                   # deduplicate the RepeatMasker library
```

See **[Configuration parameters](#configuration-parameters)** for the full list.

### 2. Library formats (optional inputs)

Both libraries are FASTA. Sequence IDs encode the classification:

```text
# custom_library — RepeatMasker convention  >name#class/subclass
>myCopia_1#Class_I/LTR/Ty1_copia/Ale

# tandem_repeat_library — used by TideCluster to name tandem families
>PisTR-B#Satellite
```

Use the canonical slash-separated classification scheme listed under
[Annotation categories](#annotation-categories).

### 3. Run

```bash
apptainer run -B /path/to/data -B $PWD carp_1.0.0rc3.sif -c config.yaml -t 20
# equivalently with singularity:
singularity run -B /path/to/data -B $PWD carp_1.0.0rc3.sif -c config.yaml -t 20
```

- `-t` sets the thread count.
- `-B` binds host directories into the container so it can read inputs and write
  outputs. The `config.yaml` and every path it references must live under a
  bound path (above, `$PWD` covers the current directory).

On Metacentrum, use
[`scripts/annotate_repeats_metacentrum.sh`](scripts/annotate_repeats_metacentrum.sh)
and adjust the input/output/image paths.

## Configuration parameters

Everyday knobs (full reference, including library-reduction, DANTE_TIR-fallback,
and TideCluster tuning, in **[docs/configuration.md](docs/configuration.md)**):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `genome_fasta` | — (required) | Genome assembly FASTA to annotate |
| `output_dir` | — (required) | Output directory |
| `custom_library` | none | Extra repeat library merged into the RepeatMasker library |
| `tandem_repeat_library` | none | Reference used by TideCluster to name tandem families |
| `repeatmasker_sensitivity` | `default` | RepeatMasker mode: `rush`, `default`, or `quick` |
| `repeatmasker_culling_limit` | `0` (off) | rmblastn `-culling_limit` for RepeatMasker — caps redundant per-locus HSPs; `2` ≈ 3× faster at ~−0.7 % masked bp |
| `tidecluster_reannotate_culling_limit` | `0` (off) | Same culling for TideCluster reannotation; `2` ≈ 3.7× faster. With superfamily-merge on (below) the result is culling-independent, so this becomes a pure speed knob |
| `tidecluster_reannotate_superfamily_merge` | `True` | Group sibling TRCs by superfamily when applying the RM-on-TideCluster array-length filter, so a real tandem array tiled by several near-identical satellite TRCs is recovered instead of fragmented and lost (and the result no longer depends on culling) |
| `rm_tc_tandem_gate` | `True` | A TideCluster-RM satellite may override a TE call only where independent tandem evidence (TideHunter) supports it; an unsupported satellite tiling a non-tandem TE is demoted below the TE, preventing spurious TE→satellite over-masking |
| `reduce_library` | `True` | Deduplicate the RepeatMasker library (smaller, faster) |
| `cleanup_intermediates` | `minimal` | Delete per-tool scratch after a successful run: `minimal` (default), `maximal` (also the big TideCluster trees + tool workdirs), or `none`. `run_pipeline.py --keep-all` forces `none`; manifest outputs are never touched |

The pipeline also screens the LTR library against `Class_II/Subclass_1` elements
from the custom library: a `Class_I` sequence that resembles a DNA transposon is
likely a composite element and is removed to avoid mis-annotation.

## Main output — `Repeat_Annotation_Unified.gff3`

The primary deliverable: an **authoritative, non-overlapping** GFF3 that merges
every tool's calls and resolves conflicts by **tier priority** — a higher-tier
call wins contested territory, so each base is annotated once (the only
tolerated overlaps are nested children and `Simple_repeat`/`Low_complexity` over
a TE).

**Feature levels.** Features are **Level 1** (top-level) or **Level 2** (nested
children carrying a `Parent`, e.g. member copies inside a tandem LTR-RT array).

**Column 3 (`type`)** is `repeat_region` when the classification starts with
`Satellite`, `Simple_repeat`, `Low_complexity`, `rDNA`, or `Unknown`; otherwise
`transposable_element`.

**Key attributes (column 9):**

| Attribute | Meaning |
|-----------|---------|
| `ID` | Unique id; `UA_L1_…` (Level 1) or `UA_L2_…` (Level 2) |
| `Name` | Classification path, or the bare `TRC_<n>` for a tandem family |
| `classification` | Slash-separated label (drives `type`) |
| `source_tool` / `source_tier` | Which tool produced the call and its priority |
| `element_type` | `complete` or `partial` (DANTE_LTR elements) |
| `structure` / `copy_number` | `LTR_RT_TR` tandem LTR-RT container and its copy count |
| `TE_origin` / `TE_origin_structure` | A satellite that is actually a tandem of TEs |
| `Parent` | Links a Level-2 child to its Level-1 container |

**Tier priority (high → low):**

1. Structural TEs — DANTE_LTR, DANTE_TIR, DANTE_LINE
2. DANTE protein domains
3. TideCluster tandem clusters (default + short)
4. RepeatMasker on the TideCluster library
5. RepeatMasker similarity (TEs, simple/low complexity, rDNA subunits)
6. TideHunter residual tandems

Full field-by-field contract:
[docs/unified_annotation_gff3_spec.md](docs/unified_annotation_gff3_spec.md).

### Annotation categories

Top-level buckets in the annotation:

- **Mobile elements** (`transposable_element`)
  - **Class I** (retrotransposons): `LTR/Ty1_copia`, `LTR/Ty3_gypsy`, `LINE`,
    `SINE`, `DIRS`, `Penelope`, `pararetrovirus`
  - **Class II** (DNA transposons): `Subclass_1/TIR` (`EnSpm_CACTA`, `hAT`,
    `MuDR_Mutator`, `PIF_Harbinger`, `Tc1_Mariner`, `MITE`, …),
    `Subclass_2/Helitron`
- **Tandem repeats / Satellite** — TideCluster families (`Name=TRC_<n>`)
- **rDNA** — `rDNA/45S_rDNA`, `rDNA/5S_rDNA` (with nested subunits)
- **Simple_repeat**, **Low_complexity**, **Unknown**

The mobile-element hierarchy comes from
[REXdb Viridiplantae v4.0](https://github.com/repeatexplorer/rexdb); satellite,
rDNA, simple/low-complexity, and unknown classes come from RepeatMasker /
TideCluster. The complete, machine-readable leaf list (and tool-native aliases)
is the single source of truth in
[`classification_vocabulary.yaml`](classification_vocabulary.yaml); the full
expanded hierarchy is also in
[docs/output_reference.md](docs/output_reference.md).

## Other top-level outputs

A short selection (full directory tree in
**[docs/output_reference.md](docs/output_reference.md)**):

| Output | Description |
|--------|-------------|
| `Repeat_Annotation_Unified.gff3` | Main annotation (above) |
| `DANTE_filtered.gff3`, `DANTE_LTR.gff3`, `DANTE_TIR.gff3`, `DANTE_LINE.gff3` | Per-tool annotations |
| `Mobile_elements_*.gff3`, `Tandem_repeats_*.gff3`, `rDNA_*.gff3`, `Simple_repeats_*.gff3`, `Low_complexity_*.gff3` | Per-class splits |
| `summary_statistics.csv`, `summary_plots.pdf` | Genome-wide statistics and plots |
| `Repeat_density/`, `Repeat_density_by_class_bigwig/` | Density tracks (BigWig) for genome browsers |
| `all_repeats_for_masking.bed`, `gaps_10plus.bed` | Masking coordinates and assembly gaps |
| `carp_manifest.json`, `run_provenance.json` | Machine-readable output map and run provenance |
| `TideCluster_report.html`, `DANTE_LTR_report.html` | Interactive reports |

A high-level workflow schematic is in
[`figs/workflow_overview.svg`](figs/workflow_overview.svg).

## Documentation

- [docs/output_reference.md](docs/output_reference.md) — complete output tree and
  density-track reference
- [docs/configuration.md](docs/configuration.md) — every configuration parameter
- [docs/unified_annotation_gff3_spec.md](docs/unified_annotation_gff3_spec.md) —
  `Repeat_Annotation_Unified.gff3` field contract
- [docs/development.md](docs/development.md) — building the container, releasing,
  regenerating the workflow diagram
- [CHANGELOG.md](CHANGELOG.md) — release history

## License

See [LICENSE](LICENSE).
