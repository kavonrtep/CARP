# Pipeline Output Cleanup Audit

## Scope

This note summarizes cleanup candidates before implementing a `debug` switch.
The review is based on:

- `Snakefile` outputs, especially `rule all`, `directory(...)` outputs, and ad hoc scratch paths.
- A real run in `tmp/output` on 2026-04-09.

Observed sample output size:

- `tmp/output`: 1803 files, about `2.0G`

## Keep In Normal Mode

These are primary deliverables or are directly exposed via top-level symlinks and reports:

- Top-level symlinks such as `DANTE_LTR.gff3`, `DANTE_TIR.gff3`, `TideCluster_report.html`, `Mobile_elements_RepeatMasker.gff3`
- `summary_statistics.csv`, `summary_plots.pdf`, `repeat_annotation_report.html`
- `RepeatMasker/*.gff3`, `RepeatMasker/*.bw`, `all_repeats_for_masking.bed`, `gaps_10plus.bed`
- `DANTE/DANTE.gff3`, `DANTE/DANTE_filtered.gff3`
- `DANTE_LTR/DANTE_LTR.gff3`, `DANTE_LTR/LTR_RTs_library.fasta`
- `DANTE_TIR/DANTE_TIR_final.gff3`, `DANTE_TIR/all_representative_elements_min3.fasta`
- `DANTE_LINE/DANTE_LINE.gff3`, `DANTE_LINE/LINE_rep_lib.fasta`
- `TideCluster/default/TideCluster_clustering.gff3`, `TideCluster/default/RM_on_TideCluster_Library.gff3`, `TideCluster/default/TideCluster_clustering_10k.bw`, `TideCluster/default/TideCluster_clustering_100k.bw`
- `Repeat_Annotation_NoSat_split_by_class_gff3/` and `Repeat_Annotation_NoSat_split_by_class_bigwig/`

The last two directories are large (`82M` and `105M`) but are part of the current published output set and should not be removed by default.

## Strong Cleanup Candidates

These look safe to remove automatically in normal mode after a successful run.

| Path | Files | Size | Why it is a cleanup candidate |
| --- | ---: | ---: | --- |
| `.cache/` | 3 | 152K | Runtime cache only |
| `.snakemake/` | 71 | 456K | Snakemake metadata, useful only for debugging/restart state |
| `logs/` | 18 | 48K | Helpful for troubleshooting, not a deliverable |
| `benchmarks/` | 28 | 116K | Needed only for performance inspection |
| `Libraries/workdir/` | 323 | 23M | CAP3, BLAST, and per-class scratch created by `reduce_library` |
| `RepeatMasker/workdir/` | 0 | 4K | Scratch directory created by `repeatmasker_wrapper.py` |
| `RepeatMasker/genome_cleaned.fasta` | 1 | 306M | Copy of input genome staged only for RepeatMasker |
| `RepeatMasker/combined_library_reduced.fasta` | 1 | 2.2M | Copy of reduced library staged only for RepeatMasker |
| `DANTE/DANTE_filtered.gff3.tmp.gff3` | 1 | 18M | Temporary cleaned DANTE file used only during merge |
| `DANTE_LTR/LTR_RTs_library.fasta.reformatted` | 1 | 2.1M | Intermediate reformatted library |

Approximate reclaimable space from this set: `348.5 MiB`

## Debug-Only Candidates

These are not required for the current top-level outputs, but removing them may reduce drill-down detail or make component-specific debugging harder.

| Path | Files | Size | Notes |
| --- | ---: | ---: | --- |
| `DANTE_LINE/mmseqs/` | 3 | 1.2M | Internal clustering detail behind `LINE_rep_lib.fasta` |
| `DANTE_LTR/library/` | 13 | 96M | Internal LTR library construction, including MMseqs intermediates |
| `DANTE_TIR/DANTE_TIR.RData` | 1 | 77M | Large analysis workspace, likely debug/report support only |
| `DANTE_TIR/mmseqs2/` | 8 | 632K | Internal clustering detail |
| `DANTE_TIR/log/` | 2 | 12K | Internal component logs |
| `DANTE_TIR/img/` | 6 | 156K | Per-family images; useful for debugging/report detail |
| `DANTE_TIR/report.html` | 1 | 2.5K | Internal report, not linked from top-level outputs |
| `TideCluster/default/TideCluster_consensus/` | 114 | 1.9M | Per-TRC consensus FASTA detail |
| `TideCluster/default/TideCluster_kite/` | 176 | 14M | Deep per-TRC HTML and RData detail |
| `TideCluster/default/TideCluster_tarean/` | 274 | 206M | Largest TideCluster internal working tree |
| `TideCluster/default/TideCluster_clustering_split_files/` | 57 | 280K | Per-cluster split GFFs |
| `TideCluster/default/TideCluster_clustering_split_files_bigwig/` | 114 | 133M | Per-cluster bigWig tracks |

Approximate reclaimable space from this set: `526.0 MiB`

Important caveat:

- `TideCluster/default/TideCluster_kite/` and `TideCluster/default/TideCluster_tarean/` may support deep links or secondary content behind `TideCluster_index.html`. They are good `debug: true` retainers, but not the first targets for default cleanup unless report behavior is verified afterward.

## Wrapper-Specific Extras

These were present in `tmp/output`, but they are not produced by `Snakefile`. They come from `scripts/annotate_repeats_metacentrum.sh`.

| Path | Size | Notes |
| --- | ---: | --- |
| `genome.fasta` | 306M | Copy of original input genome |
| `config.yaml` | 77B | Runtime config snapshot |
| `env.sh` | 6.1K | Environment dump |
| `pbs_script.sh` | 5.0K | Submitted PBS script |
| `pbs_log.txt` | 37K | Cluster job log |

Approximate reclaimable space from this set: `305.2 MiB`

## Recommended First Implementation

For a first `debug` implementation:

- `debug: true`: keep current behavior, preserve everything.
- `debug: false`: remove the Strong Cleanup Candidates after the workflow succeeds.

Do not remove these in the first pass:

- `Repeat_Annotation_NoSat_split_by_class_gff3/`
- `Repeat_Annotation_NoSat_split_by_class_bigwig/`
- `TideCluster/default/TideCluster_kite/`
- `TideCluster/default/TideCluster_tarean/`
- `DANTE_LTR/library/`

Those can be revisited later as an `aggressive_cleanup` mode after verifying report integrity and downstream usability.
