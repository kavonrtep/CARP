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
| `bigwig_max_workers` | `0` (no ceiling) | Hard ceiling on concurrent workers in the BigWig density rules (`make_bigwig_density`, `make_unified_tandem_per_family_bigwig`, `make_tidecluster_tandem_per_family_bigwig`), independent of `threads`. Caps **concurrency only** — every task writes its own independent `.bw`, so this cannot change any track. The default is `0` because those rules size the pool from a **measurement** rather than a guess: the heaviest task (largest input, finest resolution) runs first on its own, and the pool is then `budget × 0.8 / its peak RSS`, where the budget is `resources.mem_mb` when a profile sets it, else the tightest of the cgroup limit (walking **up** the hierarchy, so a PBS/Slurm job-scope limit is honoured rather than the node's free memory) and `/proc/meminfo MemAvailable`. **Why it exists:** on a 94.26 Gbp assembly `make_tidecluster_tandem_per_family_bigwig` peaked at **568.5 GB** — 74 % of a 768 GB host — producing 1,603 small per-family tracks, because each of the 3,206 concurrent tasks built the whole-genome 94.26 M-bin tile grid before discarding all but the bins it occupied. That grid is gone (the tiling is now streamed one sequence at a time), and this gates what remains: GFF3 import and the merged output track. Set > 0 to pin the worker count directly, e.g. on a scheduler where no memory budget can be detected. Integer ≥ 0. |
| `rm_tc_tandem_gate` | `True` | Tandem gate for the unified annotation. A Tier-4 RM-on-TideCluster satellite may override a Tier-5 TE call **only where it has independent tandem evidence** (raw TideHunter). An RM_TC array with no tandem support — a short AT-rich consensus tiling a genuinely non-tandem TE (e.g. a Tekay LTR-RT) — is demoted below the TE, so the TE is not spuriously re-labelled satellite. RM_TC over non-TE sequence and genuine satellites are unaffected. |

## Resources

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_memory_gb` | `0` (auto-detect) | Memory available to this run, in **GB** — the allocation the scheduler granted (`qsub -l mem=128gb`, `--mem=128G`), not a per-rule figure. Sets the default `resources.mem_mb` of every memory-gated rule, so one value sizes the `make_unified_annotation` worker pool, the three BigWig density pools, DANTE_TIR's concurrent CAP3 assemblies and TideCluster's TideHunter/TAREAN pools. `run_pipeline.py -m 128` is the same thing from the command line (and wins over the config file). An HPC profile or `--set-resources` still overrides it per rule. |

**Why you should set this on a cluster.** With `0` the budget is detected, and
detection can be wrong in exactly the situation CARP is usually run in — a `.sif`
under PBS/Slurm. `/proc/meminfo MemAvailable` reports the **host**, because the
kernel does not namespace it; the cgroup limit that will actually kill the job is
normally set on an *ancestor* job scope, which may not be reachable from inside
the container (controllers not delegated, `/sys/fs/cgroup` not mounted, or a
scheduler that enforces by polling instead). Detection then cannot distinguish
"no limit" from "a limit I cannot see", every pool sizes itself against the whole
node, and the run is OOM-killed hours in with no earlier symptom. This is
[TideCluster issue #6](https://github.com/kavonrtep/TideCluster/issues/6), where a
128 GB job believed it had 1.6 TB.

**Detection order** when `max_memory_gb` is `0` (first hit wins,
`scripts/mem_utils.py` / `.R`, matching TideCluster 1.20.0's chain):

| # | Source | Notes |
|---|--------|-------|
| 1 | `max_memory_gb` / `-m` | the explicit allocation |
| 2 | `AGENT_MEMORY` | GB |
| 3 | `PBS_RESC_MEM` (bytes), `SLURM_MEM_PER_NODE` (MB), `LSB_MAX_MEM_RUSAGE` (KB), `SLURM_MEM_PER_CPU` × `SLURM_CPUS_ON_NODE` | scheduler variables cross the container boundary intact, so they work where cgroups do not — unless stripped by `--cleanenv` / `--containall` |
| 4/5 | the tighter of the cgroup limit (walked leaf → root, v1 + v2, minus current usage) and `/proc/meminfo MemAvailable` | exceeding either kills the job |
| 6 | nothing readable | pools are not memory-gated |

Sources 1–3 name an *allocation*, so 80 % of it is offered to the pools, leaving
room for the parent process and page cache; 4–5 are already availability
readings and are used as-is. The resolved figure and its source are printed as
the second line of every run (`[mem] budget 104858 MB from PBS_RESC_MEM …`) and
recorded in `run_provenance.json` under `resources`. If the budget ends up
coming from host memory while a scheduler job id is set, the run says so with a
warning at startup — the case where you want `-m`.

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
| `dante_tir_cap3_max_memory_gb` | `0` (auto) | CAP3 memory budget (GB) for DANTE_TIR 0.3.0's concurrent assemblies (a split high-copy class can need tens of GB per assembly). `0` = let DANTE_TIR auto-detect: it reads the cgroup memory limit (v2 `memory.max` / v1 `limit_in_bytes`), falls back to node RAM, and uses 60 %. The `dante_tir` rule overrides this by deriving an explicit budget = 60 % of `resources.mem_mb` when a scheduler/profile allocates memory, so CAP3 respects the job's allocation even on schedulers that enforce **no** cgroup limit; this value is the fallback used when `mem_mb` is unset. Set > 0 to pin an explicit budget (non-cgroup schedulers, or very large genomes where an unbounded CAP3 pool would over-allocate on a big shared node). Integer, ≥ 0. |
| `make_unified_max_workers` | `8` | Hard ceiling on concurrent `mclapply` workers inside `make_unified_annotation`, independent of `threads`. Caps **concurrency only** — batch composition comes from `threads`/`batch_size`, so this cannot change the annotation, only how many batches run at once (verified: byte-identical output at 1 vs 5 workers). **Why it exists:** `mclapply` forks, and R's GC dirties the inherited parent heap, so each worker's peak RSS converges on the **parent's** regardless of its own batch size — measured on a 94 Gbp genome at 48.3 GB per worker for *every* one of 55 batches (min = median = max), against a 48.4 GB parent, for batches spanning 143 Mb to 2.15 Gb. At `threads: workflow.cores` = 96 that demanded ~2.66 TB on a 768 GB host: the run thrashed (388 M minor page faults ≈ 1.5 TB of copy-on-write copying at only 4 workers), 15 workers were OOM-killed, and the rule died after 5 h 50 m. The identical work at 4 workers took **885 s**. Tier resolution is only ~13 % of this rule's wall time (loading and GFF3 export are serial), so a low ceiling costs almost nothing. The rule *also* applies a memory gate on top of this: budget × 0.8 / parent RSS, where the budget is `resources.mem_mb` when a profile sets it, else the tightest of the cgroup limit (walking **up** the hierarchy, so a PBS/Slurm job-scope limit is honoured rather than the node's free memory) and `/proc/meminfo MemAvailable`. Integer ≥ 0; `0` = no ceiling (memory gate still applies). |
| `dante_tir_fallback_max_group_size` | `1000` | Guards the fallback's O(N²) all-vs-all flank alignment against OOM on large genomes. When a TIR subtype has **more** than this many TPase anchors, they are split into deterministic, mmseqs-clustering-based groups of at most this size and aligned per group (compute/memory drop from O(N²) to O(N·group)). At or below the threshold the run is byte-identical to no grouping. Clustering keeps similar flanks together, so the split is near-lossless. Must be ≥ 2; raise it to split only enormous families (closer to ungrouped output, higher peak RAM). |

**Fallback flank bounds are not config keys.** How far `dante_tir_fallback` may
extend a TPase anchor on each side is read directly from `max_extension_per_side`
in `classification_vocabulary.yaml`, **per superfamily**: `EnSpm_CACTA` 5800,
`MuDR_Mutator` 5700, `Tc1_Mariner` 3200, `hAT` 2300, `PIF_Harbinger` 1400, and
6000 for a TIR superfamily not listed. These numbers are **measured, not
guessed** — the p95 of the true core→boundary distance over 13,760 element-sides
taken from DANTE_TIR *primary* successes (elements with `tir_seq5` + `tir_seq3` +
`tsd`, so the boundary is exact) across wheat, Lycopus and Boechera. The medians
span 6× between superfamilies (`PIF_Harbinger` 951 bp vs `EnSpm_CACTA` 3544 bp),
which is why one shared cap cannot work: the previous shared `1500` cut into the
**median** real element for four of the five superfamilies. p95 rather than p99 is
deliberate — these guards must fail toward truncation, because an over-extended
consensus poisons the whole RepeatMasker search while a truncated one only costs
a little masking. Edit the table in `classification_vocabulary.yaml` to change
them; `dante_line_max_extension` does **not** apply here.

### Optional fallback library (off by default)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `include_dante_tir_fallback_in_library` | `false` | When `true`, re-cluster the post-overlap fallback survivors, apply a multiplicity floor, run a strict class-aware blastn filter against the LTR / DANTE_TIR primary / LINE / custom libraries, and append the kept reps to `combined_library.fasta`. A rep is dropped if any hit lands on a subject of incompatible classification (siblings such as CACTA vs hAT count as incompatible). Audit log: `DANTE_TIR/fallback_library_dropped.tsv`. The fallback library is intentionally **not** used to filter the LTR library. |
| `dante_tir_fallback_library_min_multiplicity` | `null` | Multiplicity floor applied to fallback cluster sizes when the flag above is `true`. `null` inherits `dante_tir_min_multiplicity` (3). |

## DANTE_LINE

| Parameter | Default | Description |
|-----------|---------|-------------|
| `dante_line_max_group_size` | `1000` | Same OOM guard as `dante_tir_fallback_max_group_size`, applied to the LINE all-vs-all flank alignment. LINE families are usually small (few copies), so this rarely fires; it is a safety limit for LINE-rich large genomes. When more than this many LINE patterns are found, they are split into deterministic, mmseqs-clustering-based groups of at most this size and aligned per group; at or below the threshold the run is byte-identical to no grouping. Must be ≥ 2. |
| `dante_tir_fallback_support_fraction` | `0.5` | Same group-scaled support rule as `dante_line_support_fraction`, applied to the DANTE_TIR fallback — it shares `dante_line`'s flank-extension engine (it imports `analyze_alignment_lengths`) and shows the same failure. Measured on a wheat run: **85 %** of all fallback base pairs come from inferred flanks rather than the TPase anchor, 3,681 elements hit the 10,000 bp flank ceiling, and 13–35 % of anchor groups have fewer than 5 alignments, so the old fixed "3rd largest" filtered nothing. `0` restores the previous behaviour **and takes a single-pass code path** — worth knowing because the fallback's alignment tables are large (43 GB for EnSpm_CACTA 3' on a 14.5 Gb genome; each row carries its aligned sequences). |
| `dante_tir_fallback_min_group_alignments` | `5` | Fallback anchor groups with fewer flank alignments than this get **no** extension — the element keeps its TPase anchor only. `0` disables. |
| `screen_library_cross_class` | `True` | Run the cross-class contamination screen (`screen_library_cross_class.py`) between `reduce_library` and `reduce_library_containment`. Every other reduction CARP runs is **within** a class — `reduce_library_size.py` clusters per classification and `containment_reduce_library.py` only drops a fragment into a same-class container — so a consensus that is part one class and part another is invisible to both, and RepeatMasker masks the foreign part genome-wide under the wrong label. This blastn-screens every consensus against the rest of the library and **truncates it to its longest span carrying no foreign material at all**; a consensus is dropped only when too little survives. (An earlier design peeled conflict blocks inward from the two ends and left internal ones alone — on real data that stranded 54 % of the foreign base pairs it had already identified, because the foreign material in a chimera is a mosaic rather than one clean block: the region the boundary ran out into is itself a patchwork of decayed insertions with unclassifiable rubble between them.) Every decision lands in `Libraries/cross_class_screen.tsv`. Measured on a wheat library already carrying the `dante_line` boundary guards: 156 consensi truncated, 0.93 % of library bp, leaving **zero** identified foreign base pairs; together they project wheat `Class_I/LINE` at **2.51 %** against 15.93 % as shipped, inside the 1.41–3.64 % core-anchored band. blastn-unavailable or failing copies the library through unchanged. |
| `cross_class_min_identity` | `80` | Minimum percent identity for a blastn hit to count as cross-class contamination. |
| `cross_class_min_length` | `200` | Minimum alignment length (bp) for a blastn hit to count as cross-class contamination. |
| `cross_class_max_shared_depth` | `1` | Two classifications conflict only when their lowest common ancestor is at this depth or shallower. At `1`: `Class_I/LINE` vs `Class_I/LTR/…/Retand` (LCA `Class_I`, depth 1) conflicts, as does anything `Class_I` vs `Class_II` (depth 0) or TE vs `rDNA`; but sibling LTR lineages such as `Ty1_copia/Ale` vs `Ty1_copia/Angela` (LCA depth 3) and even `Ty1_copia` vs `Ty3_gypsy` (LCA `Class_I/LTR`, depth 2) do **not** — they share real sequence homology and trimming them against each other would be destructive. An ancestor/descendant pair is never a conflict. |
| `cross_class_ownership_margin` | `0.10` | A shared region is trimmed out of the *query* only when it covers this much more of the *subject* than of the query — i.e. the region is more of what the subject **is**. A blastn hit is symmetric and does not by itself say which consensus is the chimera; without this rule the screen trims both, and on a wheat library that meant **1,114 `Ty1_copia/Angela` consensi trimmed against 86 LINE consensi** — damaging the correct library because of the broken one. Ties and sub-margin differences are left alone. |
| `dante_line_library_source` | `core` | Which span of each DANTE_LINE element seeds the **RepeatMasker library**. `core` uses only the ENDO..RT/RH span DANTE actually annotated; the inferred flanks stay in `DANTE_LINE.gff3` (the structural layer still reports the element extent it found) but never enter the similarity search. `element` restores the previous behaviour of clustering the extended regions. **Why the asymmetry:** DANTE_LINE exists to give partial coverage where no structural tool can, so its boundaries are inferred rather than observed — and a wrong flank in the GFF3 costs one element, while a wrong flank in the library is amplified by every RepeatMasker search that follows. Measured on a wheat run by re-attributing the existing hits against each library span: `Class_I/LINE` **15.91 %** of the genome with the full extension, **5.82 %** capped, **3.95 %** core-only, against a published wheat LINE content of ~1–2 %. The territory is not lost — 73.5 % of it is claimed by Ty1_copia / Ty3_gypsy consensi once the chimeras stop out-scoring them, taking `Class_I/LTR` from 53.4 % back to ~67 %, in line with the literature. |
| `dante_tir_fallback_library_source` | `core` | The same choice for the DANTE_TIR fallback: `core` seeds the library from the TPase anchor span only. The fallback exists to cover TIR elements where DANTE_TIR found no complete structural element but DANTE still reports a TPase, so it is the least reliable layer in the pipeline and 85 % of its base pairs are inferred flank. `element` restores the previous behaviour. |
| `dante_line_support_fraction` | `0.5` | Fraction of a LINE group's flank alignments that must reach the extension kept for that group. `dante_line` sets an element's boundary to the *k*-th largest of its pairwise flank alignments; `k` used to be `--min-num-alignments` (3) regardless of group size, which is satisfiable by construction — a group with exactly three alignments has its third-largest **equal to its minimum**, so nothing is filtered. `k` is now `max(min_num_alignments, ceil(support_fraction × n_alignments))`, so a handful of near-identical outlier pairs cannot set the boundary for a large group. Measured on a wheat run: `LINE_group_2147` took a 5,968 bp extension from 3 of its 51 partners, where the 25th percentile of the same group was 751 bp. `0` restores the previous behaviour. Must be in [0, 1]. |
| `dante_line_min_group_alignments` | `5` | LINE groups with fewer flank alignments than this get **no** extension — too few partners to place a boundary, so the element keeps its ENDO/RT domain core. Most runaway extensions came from groups of 3–4 members. `0` disables. |
| `dante_line_max_extension` | `0` (use the measured bounds) | Symmetric **escape hatch** on the DANTE_LINE flank extension: any value > 0 caps **both** sides at that many bp. It is a **complete manual override** — it disables the per-side bounds, the core→3' end span bound, *and* the unconverged-inference gate, so a value you set is never quietly tightened by a vocabulary number. At the default `0` the bounds are read from `classification_vocabulary.yaml` and are **measured, not guessed** — 3400 bp 5', a 2500 bp 3' backstop, a 4500 bp cap on (domain core + 3' extension), and a tighter 2000 bp 5' / 800 bp 3' fallback for an element whose flank inference never converged; see **LINE flank bounds are not config keys** below the table. Why a bound at all: the uncapped inference can reach `--flank` (10,000) per side, producing 16–22 kb "LINE" consensi; no plant LINE is that long. The bounds are deliberately set to fail toward truncation: DANTE_LINE returns elements whose boundaries are inferred, not observed, so a **truncated** consensus is the safe failure — it still masks its own family — whereas an over-extended one mislabels whatever its tail matches across the entire genome. Extensions are trimmed, not dropped; the domain core is always kept. This key reaches **DANTE_LINE only**; the DANTE_TIR fallback has no equivalent config key and always reads its per-superfamily bounds from `max_extension_per_side` (see [DANTE_TIR](#dante_tir)). |
| `dante_line_max_element_length` | `8000` | Cap in bp on the whole LINE element (core + both extensions). The equivalent bound for the DANTE_TIR fallback is not a config key: it is read per TIR superfamily from `max_consensus_length` in `classification_vocabulary.yaml`, because CACTA genuinely reaches ~20 kb where `Tc1_Mariner` cannot. Applied last, after the per-side bounds (`dante_line_max_extension` / `max_extension_per_side`) and the `max_core_to_3prime_end` span bound; when the two sides together exceed the remaining budget it is split evenly, except that a side asking for less than its half keeps what it asked for and yields the rest. A core already at or over this length gets no extension and is **never itself trimmed**. Full-length plant LINEs run 4–7 kb around an ENDO+RT core of ~2.1 kb. `0` disables. |
| `line_complete_elements` | `True` | Mark LINE elements whose **both ends were directly observed** — a poly-A tail plus a target-site duplication (TSD) — as `Status=complete` in `DANTE_LINE.gff3`, **replace their inferred span with the measured one**, and record the evidence (`TSD=`, `PolyA_length=`). Every other element gets `Status=inferred`. Sequences also go to `DANTE_LINE/LINE_complete_elements.fasta` (RepeatMasker header convention, so it is usable as a library directly). Because the span changes this **does affect the annotation**, for the 0.7–5 % of elements that qualify — see *Confirmed LINE boundaries* below. ~4 s on a 5.4 Gb genome. `False` skips it. |

**LINE flank bounds are not config keys.** At the default `dante_line_max_extension:
0`, how far `dante_line` may extend an ENDO/RT domain core is read from three
`Class_I/LINE` tables in `classification_vocabulary.yaml`:

- `max_extension_per_side` — **3400 bp 5'**, **2500 bp 3'**, the allowance for an
  extension the flank alignment converged on, bounding each flank separately.
- `max_core_to_3prime_end` — **4500 bp**, bounding (domain core length + 3'
  extension), i.e. the distance from the **start** of the domain core to the
  element's 3' end. Applied *after* the per-side caps. **LINE only** — TIR is
  deliberately not covered by this table and stays on per-side bounds.
- `unconverged_max_extension_per_side` — **2000 bp 5'**, **800 bp 3'**, used
  *instead* of the allowance above for a side whose flank inference never found a
  boundary. See **Converged vs unconverged extensions** below. **LINE only** —
  TIR is deliberately not covered by this table either.

Setting `dante_line_max_extension` above `0` disables all three; it is a complete
manual override.

The 3' per-side number is a **backstop only**: the observed maximum 3' extension is
2338 bp (p99 2180), so it never binds in practice; it exists so a pathological
inference is still bounded. The span bound is what actually constrains the 3' side.

**Converged vs unconverged extensions.** `dante_line` infers how far an element
extends by aligning the flanking sequence of many copies. Sometimes that alignment
stops on its own, well inside the search window — the evidence chose a length.
Sometimes it runs until it hits the window edge, and the "length" is just wherever
we stopped looking. These are treated as two confidence classes:

- Raw inferred extension **within** its allowance (`max_extension_per_side`, and
  for the 3' side also `max_core_to_3prime_end`) → **converged**: the alignment
  stopped on evidence, so the full inferred length is kept.
- Raw inferred extension **exceeding** that allowance → **not converged**: no
  boundary was found, so the element falls back to
  `unconverged_max_extension_per_side` (2000 bp 5' / 800 bp 3') rather than being
  handed the generous bound. Those are deliberately the values shipped in
  1.7.0rc1.

The fallback can never **loosen** a bound: the clamp takes the tighter of the
fallback and the allowance, because the span rule can squeeze the 3' allowance
below 800 bp (a core already at or past the 4500 bp span has an allowance of `0`).

**Why the gate exists.** An end-to-end check on GCA_973357735.1 (550 LINE
elements) found its raw 5' inference has a median of **7,613 bp against a
10,000 bp search window** — for **63 %** of elements the alignment never converged
at all, so the number was the window edge, not a boundary. Granting those the
converged allowance raised total extension by **74.7 %**, moved the median element
from 4,960 bp to 7,429 bp, and put **322 of 550** elements above 7 kb — over the
documented 4–7 kb range for plant LINEs. With the gate the median element is
6,033 bp and only 22 elements exceed 7 kb. On a genome whose inference does
converge the gate changes nothing at all (byte-identical result), so it costs
nothing where it is not needed. Like every other bound here it fails toward
truncation.

**Why the 3' side is bounded as a span, not as a per-side distance.** The 3'
extension is not a biological distance — it is whatever is left over after DANTE
stops annotating the core. Measured over 148 confirmed loci, Pearson *r*(core
length, 3' extension) = **−0.936**: the two are very nearly a constant sum. The
extension alone varies **13.8×** (p90/p10), while the span from the core start
varies only **1.16×**, with per-genome medians of 3779–4247 bp across seven
genomes. Concretely, two Boechera elements both end ~3750 bp after their ENDO
domain starts, but one has a 2164 bp core with a 1568 bp extension and the other a
3893 bp core with a 155 bp extension — no single per-side number fits both.

**Measurement basis.** 218 LINE loci across 7 genomes where **both** ends were
observed (a poly-A tail found and a TSD clearing a measured chance floor), after
discarding saturated 20 bp matches (not TSDs) and elements over 7500 bp (not plant
LINEs). 5' p95 = 3360 → **3400**. Span p95 = 4274 → **4500**, chosen just above p95
and below the observed maximum of 4507 so the guard fails toward truncation — the
same p95 rule used for the TIR per-superfamily caps. On those 218 loci the new
rules truncate **9 loci (4.1 %)** on the 5' side and **1 (0.5 %)** on the 3'; the
previously shipped flat 2000 / 800 caps truncated **125 (57.3 %)** and **157
(72.0 %)**, losing 89,196 bp and 149,884 bp respectively. Edit the tables in
`classification_vocabulary.yaml` to change them.

**Confirmed LINE boundaries — `Status=complete`.** `dante_line` normally *infers*
where an element ends by aligning the flanks of many copies. A small minority need
no inference: where a recent insertion still carries both a **poly-A tail** and a
**target-site duplication (TSD)**, the extent was measured directly.

Those elements are marked `Status=complete`, **their span is replaced with the
measured one**, and the evidence is recorded as attributes. Every LINE element
ends up with one of three values, so the confidence classes are distinguishable
downstream:

| `Status` | meaning | span |
|---|---|---|
| `complete` | both ends observed — poly-A tail and TSD, given in `TSD=` and `PolyA_length=` | measured |
| `inferred` | the flank comparison extended the core; `Extension_5prime`/`Extension_3prime` say by how much | estimated |
| `core` | the comparison could not extend it — too few comparable copies | the domain core alone |

On *Boechera stricta* the split is 17 / 89 / 237 of 343 elements, so `core` is
the majority case: most LINEs never get an extension at all.

Because the span changes, this **does affect the annotation** — for the 0.7–5 % of
elements that qualify. Measured on two genomes: `Class_I/LINE` +55,889 bp on one
(32 elements complete; total annotation +0.028 %, essentially all from previously
unannotated sequence) and no change at all on the other (nothing qualified). The
measured boundary is much longer than the inferred one — median +1,740 bp at the
5′ end, +1,558 bp at the 3′ end — which is why the span is replaced rather than
merely labelled.

A confirmed call is right about **97 %** of the time. The poly-A and TSD tests are
statistically independent, so their error rates multiply: 2.8 % spurious poly-A ×
2.6 % chance TSD = 0.073 % predicted, 0.073 % observed over 2,743 loci in four
genomes. Two guards reject a candidate TSD that is itself a tandem repeat, or that
matches more than `--max-tsd-hits` places in its window — a real duplication is
unique, one taken from inside a tandem array is not. Without them the plant
telomere repeat was read as a 16 bp TSD on one genome, displacing ~5 kb of genuine
telomere annotation.

Background: [docs/line_boundaries.md](line_boundaries.md).

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
| `cleanup_intermediates` | `minimal` | Post-run deletion of per-tool intermediate scratch, done after a **successful** run only (never on failure or `--dry-run`, since `--keep-incomplete` preserves partials). `minimal` (default) removes clearly-unconsumed scratch — staged genome copies (`RepeatMasker/genome_cleaned.fasta`, `TideCluster/genome_cleaned.fasta`), `DANTE_TIR/DANTE_TIR.RData`, `DANTE_LTR/library` mmseqs/`TE*.fasta` scratch, `filter_ltr_rt_library` side-files, the DANTE tmp GFF3. `maximal` additionally purges TideCluster's scratch **files** — the k-mer counts, the write-only `ggmin`/`monomers.RData`, the kite periodogram and the duplicated per-array FASTAs inside `TideCluster_tarean` / `TideCluster_kite` / `TideCluster_consensus` — and the tool workdirs (`Libraries/workdir`, `RepeatMasker/workdir`, `DANTE_TIR/mmseqs_combined`, …). `none` keeps everything. Files listed in `carp_manifest.json` (and their symlink targets) and the CI/count-checked outputs are **never** deleted. Those three trees are no longer deleted whole: that also removed the inputs for TideCluster's per-array consensus (`tc_per_tra_consensus.py`) and report re-render (`tc_rerender_report.py`), and the kite monomer-size CSV the CARP report reads to label a family with its monomer length. Measured on a 94 Gbp run the file-level set still frees **40.2 GB of the trees' 44.6 GB (90 %)**. The CLI flag `run_pipeline.py --keep-all` forces `none`, overriding this key. |
