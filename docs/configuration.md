# CARP configuration reference

Every CARP parameter, grouped by area. Pass a YAML config with
`-c config.yaml`. Only `genome_fasta` and `output_dir` are required; every other
knob has a default that reproduces the standard run. A complete annotated
template is in [`config_full.yaml`](../config_full.yaml).

## Required

| Parameter | Description |
|-----------|-------------|
| `genome_fasta` | Path to the genome assembly FASTA to annotate. |
| `output_dir` | Output directory (created if missing). |

## Core options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `custom_library` | none | Extra repeat library (FASTA, `>name#class/subclass`) merged into the RepeatMasker library and used to screen the LTR library for Class_II contamination. |
| `tandem_repeat_library` | none | Reference library used by TideCluster to name discovered tandem families by similarity. |
| `repeatmasker_sensitivity` | `default` | RepeatMasker mode: `rush`, `default`, or `quick`. |
| `repeatmasker_culling_limit` | `0` (off) | rmblastn `-culling_limit` for the RepeatMasker step. A de-novo library carries many near-identical consensi per family, so each genomic TE copy aligns to hundreds of them — a per-locus HSP explosion that dominates rmblastn extension and RepeatMasker's post-processing. `-culling_limit N` discards a hit enveloped by ≥ N higher-scoring hits, collapsing that redundancy at search time (the best hit per locus is never culled, so masking and classification are preserved). `2` → ~3× faster RepeatMasker at ~−0.7 % masked bp; combine with `repeatmasker_sensitivity: rush` for more. Injected via an `RMBLAST_DIR` shim — no RepeatMasker/container patch, works in the existing image. Validate on your data before adopting as default. |
| `tidecluster_reannotate_culling_limit` | `0` (off) | The same culling for TideCluster's internal RepeatMasker in the `tidecluster_reannotate` step, **independent** of `repeatmasker_culling_limit`. This is the strongest culling target in the pipeline: a consensus dimer aligns to a satellite array at every rotational phase, so each locus accrues many overlapping HSPs. `2` → ~3.7× faster reannotation at ~−0.8 % masked bp (measured on tiny_pea). **Caution:** may under-mask *large* satellite families (megasatellites) the same way `reduce_tidecluster_library` can — validate masked bp on a big-satellite genome before enabling. |

## Library reduction

These three flags are **independent** of each other.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `reduce_library` | `True` | Per-class CAP3/mmseqs deduplication of the RepeatMasker library — smaller library, faster RepeatMasker. |
| `reduce_tidecluster_library` | `True` | Rotation-invariant reduction of the TideCluster consensus *dimer* library used by the genome-wide tandem remasking (`tidecluster_reannotate`). Lossless for masking (validated); ~1.5–5× faster remasking. Set `False` to remask with the full dimer library (recommended only when masking very large satellite families, dimer ≳ a few kb). |
| `reduce_library_containment` | `True` | Second-round containment reduction: drops short repeat fragments fully contained in a longer element of the **same class** (RepeatMasker masks their copies via the container, so masking and classification are preserved). ~−22 % library bp, ~−30 % RepeatMasker wall-time on the Pisum pangenome. |
| `containment_min_identity` | `80` | Minimum % identity for the containment pass to drop a fragment. |
| `containment_min_coverage` | `0.90` | Minimum fraction of a fragment's length that a longer same-class sequence must cover for it to be dropped. |
| `reduce_library_max_parallel_bp` | `50000000` | Per-class input-FASTA size (bytes) above which a `reduce_library` CAP3 class runs **sequentially** (Phase 1) instead of in parallel — caps physical memory peak. Raise on machines with plenty of RAM to trade memory for parallelism; lower on tight setups. |
| `reduce_library_max_big_cap3_parallel` | `4` | How many *big* (≥ `reduce_library_max_parallel_bp`) Class_I/LTR CAP3 classes `reduce_library` runs concurrently in its Phase 1b. CAP3 is single-threaded, so this bounds the peak of several big CAP3 jobs at once. Set to `1` for strictly-sequential Phase 1 on memory-constrained machines. |

`80` / `0.90` is the validated masked-bp-lossless default; raise either threshold
for a more conservative reduction. The two `reduce_library_max_*` knobs only
affect scheduling/memory, never the reduced library's contents.

## DANTE_TIR

| Parameter | Default | Description |
|-----------|---------|-------------|
| `dante_tir_min_multiplicity` | `3` | Primary TIR elements with `Multiplicity` below this are dropped **from the library only** (they remain in `DANTE_TIR_combined.gff3` and the unified annotation). Lower to `1` to disable. |
| `dante_tir_fallback_min_alignments` | `3` | How many anchors must support an extension length before the fallback applies it (larger = more conservative extensions). |
| `dante_tir_fallback_min_cluster_size` | `3` | Minimum mmseqs cluster size whose representative is kept in the fallback repeat library. |

### Optional fallback library (off by default)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `include_dante_tir_fallback_in_library` | `false` | When `true`, re-cluster the post-overlap fallback survivors, apply a multiplicity floor, run a strict class-aware blastn filter against the LTR / DANTE_TIR primary / LINE / custom libraries, and append the kept reps to `combined_library.fasta`. A rep is dropped if any hit lands on a subject of incompatible classification (siblings such as CACTA vs hAT count as incompatible). Audit log: `DANTE_TIR/fallback_library_dropped.tsv`. The fallback library is intentionally **not** used to filter the LTR library. |
| `dante_tir_fallback_library_min_multiplicity` | `null` | Multiplicity floor applied to fallback cluster sizes when the flag above is `true`. `null` inherits `dante_tir_min_multiplicity` (3). |

## TideCluster (1.16.0)

All four default to TideCluster's own defaults, so a default CARP run adds no
extra flags.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `tidecluster_detect_rdna` | `true` | Flag 45S/5S rDNA TRCs; the unified annotation labels them `rDNA_45S` / `rDNA_5S` (array level) instead of a generic Satellite TRC. `false` = `--no_rdna`. |
| `tidecluster_rdna_library` | `""` (bundled) | Override for TideCluster's rDNA reference library. Empty uses the bundled `data/rdna_library.fasta`. |
| `tidecluster_keep_trc_overlaps` | `false` | When `false`, TideCluster makes the clustering GFF3 non-overlapping across satellite TRCs (dominant-TRC-wins). `true` keeps raw overlapping regions (`--keep_overlaps`). |
| `tidecluster_chunk_size` | `50000000` | Genome chunk size (bp) for the parallel, pooled RepeatMasker in `tc_reannotate`. Sequences below 2× this are only packed (byte-identical); larger ones split with < 0.15 % masked-bp drift. |
