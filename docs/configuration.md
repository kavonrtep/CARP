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
| `repeatmasker_culling_limit` | `0` (off) | rmblastn `-culling_limit` for the RepeatMasker step. A de-novo library carries many near-identical consensi per family, so each genomic TE copy aligns to hundreds of them — a per-locus HSP explosion that dominates rmblastn extension and RepeatMasker's post-processing. `-culling_limit N` discards a hit enveloped by ≥ N higher-scoring hits, collapsing that redundancy at search time (the best hit per locus is never culled, so masking and classification are preserved). `2` → ~3× faster RepeatMasker at ~−0.7 % masked bp; combine with `repeatmasker_sensitivity: rush` for more. Injected via an `RMBLAST_DIR` shim — no RepeatMasker/container patch, works in the existing image. The shim also pins rmblastn to `-num_threads 1`, which makes culling **bit-reproducible**: BLAST `-culling_limit` is otherwise non-deterministic across threads (RepeatMasker runs rmblastn with `-num_threads 4` even at `-pa 1`), so culling ON previously jittered ~0.03–0.07 % of hits run-to-run. Throughput is preserved because parallelism runs at the chunk level (the RepeatMasker-wrapper pool), which also removes the old thread-over-subscription. Validate on your data before adopting as default. |
| `tidecluster_reannotate_culling_limit` | `0` (off) | The same culling for TideCluster's internal RepeatMasker in the `tidecluster_reannotate` step, **independent** of `repeatmasker_culling_limit`. This is the strongest culling target in the pipeline: a consensus dimer aligns to a satellite array at every rotational phase, so each locus accrues many overlapping HSPs. `2` → ~3.7× faster reannotation at ~−0.8 % masked bp (measured on tiny_pea). **Caution:** may under-mask *large* satellite families (megasatellites) the same way `reduce_tidecluster_library` can — validate masked bp on a big-satellite genome before enabling. With `tidecluster_reannotate_superfamily_merge` on (default), culling no longer changes which regions are masked — it becomes a pure speed knob. Uses the same `RMBLAST_DIR` shim as `repeatmasker_culling_limit`, so it inherits the `-num_threads 1` pin that makes culling **bit-reproducible** (TideCluster's internal RepeatMasker is made deterministic with no TideCluster patch). |
| `tidecluster_reannotate_superfamily_merge` | `True` | Superfamily-aware array recovery for the RM-on-TideCluster reannotation. TideCluster's `tc_reannotate` keeps an array only where a **single TRC** reaches the monomer-length threshold; when one real tandem array is tiled by ≥2 TRCs of the same superfamily (common for short, near-identical satellites) each per-TRC piece is sub-threshold and the whole array is discarded — so a TE-derived tandem array is missed and the underlying TE prevails. (Enabling `tidecluster_reannotate_culling_limit` masks this by collapsing each locus to one TRC, making the annotation culling-dependent.) When `True`, the pipeline runs `tc_reannotate --debug` and re-filters its raw hits (`scripts/tc_reannotate_sf_filter.py`) grouping sibling TRCs by superfamily, so real arrays are recovered **deterministically (independent of culling)** while every feature keeps its bare `TRC_<n>` `Name`. |
| `rm_tc_tandem_gate` | `True` | Tandem gate for the unified annotation. A Tier-4 RM-on-TideCluster satellite may override a Tier-5 TE call **only where it has independent tandem evidence** (raw TideHunter). An RM_TC array with no tandem support — a short AT-rich consensus tiling a genuinely non-tandem TE (e.g. a Tekay LTR-RT) — is demoted below the TE, so the TE is not spuriously re-labelled satellite. RM_TC over non-TE sequence and genuine satellites are unaffected. |

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
| `dante_tir_fallback_max_group_size` | `1000` | Guards the fallback's O(N²) all-vs-all flank alignment against OOM on large genomes. When a TIR subtype has **more** than this many TPase anchors, they are split into deterministic, mmseqs-clustering-based groups of at most this size and aligned per group (compute/memory drop from O(N²) to O(N·group)). At or below the threshold the run is byte-identical to no grouping. Clustering keeps similar flanks together, so the split is near-lossless. Must be ≥ 2; raise it to split only enormous families (closer to ungrouped output, higher peak RAM). |

### Optional fallback library (off by default)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `include_dante_tir_fallback_in_library` | `false` | When `true`, re-cluster the post-overlap fallback survivors, apply a multiplicity floor, run a strict class-aware blastn filter against the LTR / DANTE_TIR primary / LINE / custom libraries, and append the kept reps to `combined_library.fasta`. A rep is dropped if any hit lands on a subject of incompatible classification (siblings such as CACTA vs hAT count as incompatible). Audit log: `DANTE_TIR/fallback_library_dropped.tsv`. The fallback library is intentionally **not** used to filter the LTR library. |
| `dante_tir_fallback_library_min_multiplicity` | `null` | Multiplicity floor applied to fallback cluster sizes when the flag above is `true`. `null` inherits `dante_tir_min_multiplicity` (3). |

## DANTE_LINE

| Parameter | Default | Description |
|-----------|---------|-------------|
| `dante_line_max_group_size` | `1000` | Same OOM guard as `dante_tir_fallback_max_group_size`, applied to the LINE all-vs-all flank alignment. LINE families are usually small (few copies), so this rarely fires; it is a safety limit for LINE-rich large genomes. When more than this many LINE patterns are found, they are split into deterministic, mmseqs-clustering-based groups of at most this size and aligned per group; at or below the threshold the run is byte-identical to no grouping. Must be ≥ 2. |

## TideCluster (1.16.0)

All four default to TideCluster's own defaults, so a default CARP run adds no
extra flags.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `tidecluster_detect_rdna` | `true` | Flag 45S/5S rDNA TRCs; the unified annotation labels them `rDNA_45S` / `rDNA_5S` (array level) instead of a generic Satellite TRC. `false` = `--no_rdna`. |
| `tidecluster_rdna_library` | `""` (bundled) | Override for TideCluster's rDNA reference library. Empty uses the bundled `data/rdna_library.fasta`. |
| `tidecluster_keep_trc_overlaps` | `false` | When `false`, TideCluster makes the clustering GFF3 non-overlapping across satellite TRCs (dominant-TRC-wins). `true` keeps raw overlapping regions (`--keep_overlaps`). |
| `tidecluster_chunk_size` | `50000000` | Genome chunk size (bp) for the parallel, pooled RepeatMasker in `tc_reannotate`. Sequences below 2× this are only packed (byte-identical); larger ones split with < 0.15 % masked-bp drift. |

### Output cleanup

Handled by `run_pipeline.py` (not a Snakemake rule), so it applies to full-pipeline runs launched through the wrapper / container.

| Parameter | Default | Description |
|---|---|---|
| `cleanup_intermediates` | `minimal` | Post-run deletion of per-tool intermediate scratch, done after a **successful** run only (never on failure or `--dry-run`, since `--keep-incomplete` preserves partials). `minimal` (default) removes clearly-unconsumed scratch — staged genome copies (`RepeatMasker/genome_cleaned.fasta`, `TideCluster/genome_cleaned.fasta`), `DANTE_TIR/DANTE_TIR.RData`, `DANTE_LTR/library` mmseqs/`TE*.fasta` scratch, `filter_ltr_rt_library` side-files, the DANTE tmp GFF3. `maximal` additionally purges the large `TideCluster_tarean` / `TideCluster_kite` / `TideCluster_consensus` trees (can be multiple GB) and the tool workdirs (`Libraries/workdir`, `RepeatMasker/workdir`, `DANTE_TIR/mmseqs_combined`, …). `none` keeps everything. Files listed in `carp_manifest.json` (and their symlink targets) and the CI/count-checked outputs are **never** deleted. The CLI flag `run_pipeline.py --keep-all` forces `none`, overriding this key. |
