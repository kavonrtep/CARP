import os
import re
def create_dirs(*dirs):
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)
print(config)
subdirs = [config['output_dir']+"/"+i for i in ['DANTE', 'DANTE_TIR', 'DANTE_TIR_FALLBACK',
                                                'DANTE_LINE', 'DANTE_LTR',
                                                'TideCluster/default',
                                                'TideCluster/short_monomer',
                                                'Libraries', 'RepeatMasker',
                                                'logs', 'benchmarks']]
create_dirs(*subdirs)
snakemake_dir = os.path.dirname(workflow.snakefile)

BENCHMARKED_RULES = [
    "clean_genome_fasta", "dante", "dante_tir", "dante_tir_fallback",
    "merge_dante_tir_with_fallback", "make_tir_combined_library", "filter_dante",
    "dante_line", "dante_ltr", "make_library_of_ltrs",
    "tidecluster_long", "tidecluster_short", "tidecluster_reannotate",
    "merge_tidecluster_default_and_short", "build_fallback_tir_library",
    "make_subclass_2_library",
    "filter_ltr_rt_library", "concatenate_libraries", "reduce_library",
    "repeatmasker", "subtract_satellites_from_rm", "merge_rm_and_dante",
    "make_track_for_masking", "make_track_for_Ns",
    "make_summary_statistics_and_split_by_class", "make_bigwig_density",
    "add_top_level_outputs", "calculate_bigwig_density", "add_html_outputs",
    "make_unified_tandem_per_family_bigwig",
    "make_tidecluster_tandem_per_family_bigwig",
    "calculate_seqlengths", "make_summary_plots", "make_repeat_report",
    "make_unified_annotation",
]
print(snakemake_dir)
def filter_fasta(input_file, output_file, filter_string):
    """Filter FASTA files based on a filter_string, which is a regular expression."""
    with open(input_file, "r") as f:
        with open(output_file, "w") as o:
            write_sequence = False
            for line in f:
                if line.startswith(">"):
                    if re.search(filter_string, line):
                        write_sequence = True
                        o.write(line)
                    else:
                        write_sequence = False
                elif write_sequence:
                    o.write(line)

# repeatmasker sensitivity:
if "repeatmasker_sensitivity" not in config:
    config["repeatmasker_sensitivity"] = "default"

rm_sensitivity_option = {
    "default": "",
    "rush":    "-qq",
    "quick":   "-q",
    "" :       ""
    }[config["repeatmasker_sensitivity"]]

# RepeatMasker rmblastn culling (-culling_limit). 0 = off (default; behaviour
# unchanged). N>0 caps redundant HSPs per genomic locus: on a redundant de-novo
# library (each genomic TE copy otherwise hits hundreds of near-identical
# consensi) this collapses the per-locus HSP explosion. Measured on a Pisum-type
# genome through the real wrapper path: 2 -> ~3x faster RepeatMasker at -0.7%
# masked bp, classification preserved; combine with repeatmasker_sensitivity:
# rush for ~7x at -1.8%. Injected into rmblastn via an RMBLAST_DIR shim in
# repeatmasker_wrapper.py (no container patch required).
if "repeatmasker_culling_limit" not in config:
    config["repeatmasker_culling_limit"] = 0

# Same rmblastn culling for TideCluster's internal RepeatMasker in the
# tidecluster_reannotate rule (independent knob, default 0 = off). Tandem-array
# reannotation is the strongest culling target: a dimer consensus matches a
# satellite array at every rotational phase, producing many overlapping HSPs per
# locus that culling collapses. Measured on tiny_pea: ~8x faster tc_reannotate at
# -0.8% masked bp. CAUTION: tc_reannotate already risks under-masking LARGE
# satellite families (RepeatMasker needs full consensus diversity to tile
# megasatellite arrays); culling may worsen that on big-satellite genomes, so
# validate masked bp there before enabling. Default 0 keeps behaviour unchanged.
if "tidecluster_reannotate_culling_limit" not in config:
    config["tidecluster_reannotate_culling_limit"] = 0

# Superfamily-aware array recovery for the RM-on-TideCluster reannotation.
# TideCluster's tc_reannotate keeps an array only where a *single TRC* reaches
# the monomer-length threshold; when one real tandem array is tiled by >=2 TRCs
# of the same superfamily (common for short near-identical satellites) each
# per-TRC piece is sub-threshold and the array is lost -> a TE-derived tandem
# array is missed and the underlying TE prevails. (Culling masks this by
# collapsing each locus to one TRC, making the annotation culling-dependent.)
# When True the pipeline runs tc_reannotate in --debug mode and re-filters its
# raw hits (scripts/tc_reannotate_sf_filter.py) grouping sibling TRCs by
# superfamily, so real arrays are recovered deterministically (culling-
# independent) while every feature keeps its bare TRC_<n> Name. Default True.
if "tidecluster_reannotate_superfamily_merge" not in config:
    config["tidecluster_reannotate_superfamily_merge"] = True
elif config["tidecluster_reannotate_superfamily_merge"] not in [True, False]:
    raise ValueError("Invalid value for tidecluster_reannotate_superfamily_merge. "
                     "Must be either True or False.")

# Tandem gate for the unified annotation. When True, a Tier-4 RM-on-TideCluster
# satellite may override a Tier-5 TE call only where it has independent tandem
# evidence (raw TideHunter); an RM_TC array with no tandem support (a short
# AT-rich consensus tiling a non-tandem TE) is demoted below the TE, so the TE
# is not spuriously re-labelled satellite. RM_TC over non-TE sequence and
# genuine satellites are unaffected. Default True.
if "rm_tc_tandem_gate" not in config:
    config["rm_tc_tandem_gate"] = True
elif config["rm_tc_tandem_gate"] not in [True, False]:
    raise ValueError("Invalid value for rm_tc_tandem_gate. Must be either True or False.")

# TideCluster sensitivity preset (--sensitivity {quick,default,rush}).
# Kept in sync with the RepeatMasker sensitivity setting: TideCluster uses
# RepeatMasker internally for its reannotation step.
tc_sensitivity = {
    "default": "default",
    "rush":    "rush",
    "quick":   "quick",
    "":        "default",
}.get(config["repeatmasker_sensitivity"], "default")

if "reduce_library" not in config:
    config["reduce_library"] = True
else:
    # check validity of the value
    if config["reduce_library"] not in [True, False]:
        raise ValueError("Invalid value for reduce_library_size. Must be either True or False.")

# Rotation-invariant reduction of the TideCluster consensus *dimer* library
# used by tidecluster_reannotate. DECOUPLED from `reduce_library` (which
# governs the main RepeatMasker / CAP3 library) and defaults to True: it is
# lossless for masking (validated: masked bp unchanged within ±0.15% on
# tiny_pea and an 800 Mbp genome) while shrinking the dimer library ~5–20×,
# making the genome-wide tandem remasking ~1.5–5× faster. See
# scripts/reduce_dimer_library.py.
if "reduce_tidecluster_library" not in config:
    config["reduce_tidecluster_library"] = True
elif config["reduce_tidecluster_library"] not in [True, False]:
    raise ValueError("Invalid value for reduce_tidecluster_library. Must be either True or False.")

# Second-round CONTAINMENT reduction of the per-class-reduced library before
# RepeatMasker (independent of reduce_library). Removes short fragments fully
# contained in a longer same-class element; RepeatMasker masks their copies via
# the container, so masking + classification are preserved. See
# scripts/containment_reduce_library.py. Validated masked-bp-lossless on the
# Pisum pangenome at 80% identity / 0.90 coverage (~-22% library bp, ~-30%
# RepeatMasker wall-time).
if "reduce_library_containment" not in config:
    config["reduce_library_containment"] = True
elif config["reduce_library_containment"] not in [True, False]:
    raise ValueError("Invalid value for reduce_library_containment. Must be either True or False.")
if "containment_min_identity" not in config:
    config["containment_min_identity"] = 80
if "containment_min_coverage" not in config:
    config["containment_min_coverage"] = 0.90

# DANTE_TIR_FALLBACK stringency knobs. Both default to 3; both must be
# positive integers. See dante_tir_fallback.py for semantics.
for _key in ("dante_tir_fallback_min_alignments", "dante_tir_fallback_min_cluster_size"):
    if _key not in config:
        config[_key] = 3
    if not isinstance(config[_key], int) or config[_key] < 1:
        raise ValueError(f"Invalid value for {_key}: must be a positive integer.")

# DANTE_TIR CAP3 memory budget (GB) for concurrent assemblies. DANTE_TIR 0.3.0
# sizes its CAP3 pool to a memory budget (a split high-copy class can want tens
# of GB per assembly). Default 0 = let DANTE_TIR auto-detect: it reads the cgroup
# limit (v2 memory.max / v1 limit_in_bytes), falls back to node RAM, and uses
# 60%. That is safe under cgroup-enforced schedulers but over-detects on
# schedulers that do NOT set a cgroup mem limit. The dante_tir rule therefore
# also derives an explicit budget = 60% of resources.mem_mb when a profile
# allocates that, and honours this knob as the fallback. Integer GB, >= 0.
if "dante_tir_cap3_max_memory_gb" not in config:
    config["dante_tir_cap3_max_memory_gb"] = 0
if (not isinstance(config["dante_tir_cap3_max_memory_gb"], int)
        or isinstance(config["dante_tir_cap3_max_memory_gb"], bool)
        or config["dante_tir_cap3_max_memory_gb"] < 0):
    raise ValueError("Invalid value for dante_tir_cap3_max_memory_gb: must be an integer >= 0.")

# All-vs-all flank-alignment grouping cap for the fallback (per TIR subtype)
# and for DANTE_LINE (LINE pattern set). When the number of anchors/patterns
# exceeds this value, the O(N^2) parasail all-vs-all is split into
# deterministic, clustering-based groups of at most this many sequences and run
# per group, bounding the peak memory that OOMs on large families (big
# genomes). At or below the threshold the result is byte-identical to no
# grouping, so LINE and small TIR subtypes are unaffected. Grouping keeps
# similar flanks together (mmseqs easy-cluster), so it is near-lossless. Must be
# an integer >= 2. See group_sequences_for_alignment() in global_local_aln.py.
for _key in ("dante_tir_fallback_max_group_size", "dante_line_max_group_size"):
    if _key not in config:
        config[_key] = 1000
    if not isinstance(config[_key], int) or isinstance(config[_key], bool) or config[_key] < 2:
        raise ValueError(f"Invalid value for {_key}: must be an integer >= 2.")

# The DANTE_TIR fallback shares dante_line's flank-extension engine, and shows the
# same failure: on a wheat run 85% of all fallback base pairs come from inferred
# flanks rather than the TPase anchor, 3,681 elements hit the 10 kb flank ceiling,
# and 13-35% of anchor groups have fewer than 5 alignments (so a fixed "3rd
# largest" filters nothing). It gets the same support rule. The length bound is
# per TIR superfamily, from max_consensus_length in classification_vocabulary.yaml,
# because CACTA genuinely reaches ~20 kb where Tc1_Mariner cannot.
if "dante_tir_fallback_support_fraction" not in config:
    config["dante_tir_fallback_support_fraction"] = 0.5
if (isinstance(config["dante_tir_fallback_support_fraction"], bool)
        or not isinstance(config["dante_tir_fallback_support_fraction"], (int, float))
        or not 0.0 <= config["dante_tir_fallback_support_fraction"] <= 1.0):
    raise ValueError("Invalid value for dante_tir_fallback_support_fraction: "
                     "must be a number in [0, 1].")
if "dante_tir_fallback_min_group_alignments" not in config:
    config["dante_tir_fallback_min_group_alignments"] = 5
if (not isinstance(config["dante_tir_fallback_min_group_alignments"], int)
        or isinstance(config["dante_tir_fallback_min_group_alignments"], bool)
        or config["dante_tir_fallback_min_group_alignments"] < 0):
    raise ValueError("Invalid value for dante_tir_fallback_min_group_alignments: "
                     "must be an integer >= 0.")

# Which span of an inferred element seeds the RepeatMasker library. DANTE_LINE
# and the DANTE_TIR fallback exist to give partial coverage where the structural
# tools find nothing, so their output is the least reliable in the pipeline: the
# domain span is what DANTE observed, the flanks around it are inferred. A wrong
# flank in the GFF3 costs one element; a wrong flank in the LIBRARY is amplified
# by every RepeatMasker search that follows. 'core' (default) therefore keeps the
# inferred flanks out of the library while leaving them in the GFF3, so the
# structural layer still reports the element extent it found. 'element' restores
# the previous behaviour.
#
# Measured on a wheat run by re-attributing the existing RepeatMasker hits
# against each library span: Class_I/LINE 15.91% of the genome with the full
# extension, 5.82% with the extension capped, 3.95% core-only -- against a
# published wheat LINE content of ~1-2%. The territory does not vanish: 73.5% of
# it is claimed by Ty1_copia / Ty3_gypsy consensi once the chimeras stop
# out-scoring them, taking Class_I/LTR from 53.4% back to ~67%.
for _key in ("dante_line_library_source", "dante_tir_fallback_library_source"):
    if _key not in config:
        config[_key] = "core"
    if config[_key] not in ("core", "element"):
        raise ValueError(f"Invalid value for {_key}: must be 'core' or 'element'.")

# Cross-class library screen. Every other reduction CARP runs is WITHIN a class
# (reduce_library_size.py clusters per classification; containment_reduce_library.py
# only drops a fragment into a SAME-class container), so a consensus that is part
# one class and part another is invisible to both and RepeatMasker masks the
# foreign part genome-wide under the wrong label. This screens every de-novo
# consensus against the other classes and TRIMS the contaminated span rather than
# dropping the sequence.
#
#   screen_library_cross_class (default True): master switch.
#   cross_class_min_identity (80) / cross_class_min_length (200): what counts as a
#     conflicting hit.
#   cross_class_max_shared_depth (default 1): two classes conflict only when their
#     lowest common ancestor is at this depth or shallower. 1 catches LINE-vs-LTR
#     (LCA Class_I) and Class_I-vs-Class_II, while leaving sibling LTR lineages
#     (Ale vs Angela, LCA Class_I/LTR/Ty1_copia) alone -- they share real homology
#     and trimming them against each other would be destructive.
#   cross_class_ownership_margin (default 0.10): a shared region is trimmed out of
#     the query only when it covers this much more of the subject than of the
#     query. Without it the screen is symmetric and a chimera drags its
#     correctly-built neighbour down with it (measured on a wheat library: 1,114
#     Ty1_copia/Angela consensi trimmed against 86 LINE).
if "screen_library_cross_class" not in config:
    config["screen_library_cross_class"] = True
if not isinstance(config["screen_library_cross_class"], bool):
    raise ValueError("Invalid value for screen_library_cross_class: must be True or False.")

for _key, _default in (("cross_class_min_length", 200),
                       ("cross_class_max_shared_depth", 1)):
    if _key not in config:
        config[_key] = _default
    if (not isinstance(config[_key], int) or isinstance(config[_key], bool)
            or config[_key] < 0):
        raise ValueError(f"Invalid value for {_key}: must be an integer >= 0.")

for _key, _default, _lo, _hi in (("cross_class_min_identity", 80.0, 0.0, 100.0),
                                 ("cross_class_ownership_margin", 0.10, 0.0, 1.0)):
    if _key not in config:
        config[_key] = _default
    if (isinstance(config[_key], bool)
            or not isinstance(config[_key], (int, float))
            or not _lo <= config[_key] <= _hi):
        raise ValueError(f"Invalid value for {_key}: must be a number in [{_lo}, {_hi}].")

# DANTE_LINE consensus-boundary guards. dante_line infers where a LINE ends
# from an all-vs-all alignment of its flanks; in a repeat-dense genome the flank
# is frequently another, far more abundant repeat, so the extension runs to the
# --flank ceiling (10 kb per side) and the "LINE consensus" becomes mostly
# foreign sequence that RepeatMasker then masks genome-wide as Class_I/LINE.
# Measured over 87 assemblies: the median share of Class_I/LINE base pairs
# coming from consensus extensions rather than the element's own ENDO/RT core
# rises from 16% below 0.5 Gb to 61% above 4 Gb.
#
#   dante_line_support_fraction (default 0.5): fraction of a group's flank
#     alignments that must reach the kept extension. Replaces a fixed "3rd
#     largest", which is satisfiable by construction -- a group with three
#     alignments has its third largest equal to its minimum. 0 restores the old
#     behaviour.
#   dante_line_min_group_alignments (default 5): groups with fewer flank
#     alignments get no extension; too few partners to place a boundary. 0 off.
#   dante_line_max_extension (default 0 = use the vocabulary): symmetric
#     per-side cap in bp, an escape hatch that overrides both sides. At the
#     default 0 the bounds come from max_extension_per_side in
#     classification_vocabulary.yaml, which is ASYMMETRIC for LINE -- 2000 bp
#     5' (inter-ORF + ORF1 + 5'UTR sit outside the ENDO+RT core) but only 800 bp
#     3' (ORF2 C-terminus + 3'UTR + polyA). The former shared 1500 was too tight
#     5' and too loose 3' at once. Deliberately tighter than a full-length LINE
#     needs either way: a truncated consensus still masks its family; an
#     over-extended one mislabels its neighbour genome-wide, and a false
#     positive here is amplified across the whole RepeatMasker search.
#     dante_tir_fallback takes no such knob -- its per-side bounds are per
#     superfamily from the same table (measured p95, medians span 6x).
#   dante_line_max_element_length (default 8000): whole-element cap in bp. The
#     domain core is never trimmed, only the appended flanks. No plant LINE is
#     longer than this; full-length elements run 4-7 kb. 0 off.
if "dante_line_support_fraction" not in config:
    config["dante_line_support_fraction"] = 0.5
if (isinstance(config["dante_line_support_fraction"], bool)
        or not isinstance(config["dante_line_support_fraction"], (int, float))
        or not 0.0 <= config["dante_line_support_fraction"] <= 1.0):
    raise ValueError("Invalid value for dante_line_support_fraction: must be a number in [0, 1].")

# Emit the subset of LINE elements whose BOTH ends were directly observed --
# a poly-A tail plus a target-site duplication, the two marks target-primed
# reverse transcription leaves behind. Their extent is measured rather than
# inferred, so they are the only LINEs in the run that need no boundary guard.
#
# EVALUATION OUTPUT: it feeds nothing. Not the RepeatMasker library, not the
# unified annotation, not any downstream rule. It exists so the idea can be
# judged across many genomes before anything depends on it -- measured on two,
# adding these to the LINE library moved masking +28.2% and +0.7%, clean in both
# (0.7-0.9% contamination vs the cores' own 8.8-11.4%) but far too
# genome-dependent to enable blindly. Costs ~4 s on a 5.4 Gb genome.
if "line_complete_elements" not in config:
    config["line_complete_elements"] = True
if not isinstance(config["line_complete_elements"], bool):
    raise ValueError("Invalid value for line_complete_elements: must be true or false.")

for _key, _default in (("dante_line_min_group_alignments", 5),
                       ("dante_line_max_extension", 0),
                       ("dante_line_max_element_length", 8000)):
    if _key not in config:
        config[_key] = _default
    if (not isinstance(config[_key], int) or isinstance(config[_key], bool)
            or config[_key] < 0):
        raise ValueError(f"Invalid value for {_key}: must be an integer >= 0.")

# DANTE_TIR primary-element library filter (Multiplicity floor). Default 3
# matches the pre-fallback behaviour where the library was sourced from
# `all_representative_elements_min3.fasta` produced by `dante_tir_summary.R`
# with `--min_cluster_size 3`. Set to 1 to disable. Affects
# make_tir_combined_library; partials and low-Multiplicity primaries remain
# in DANTE_TIR_combined.gff3 regardless of this value.
if "dante_tir_min_multiplicity" not in config:
    config["dante_tir_min_multiplicity"] = 3
if not isinstance(config["dante_tir_min_multiplicity"], int) or config["dante_tir_min_multiplicity"] < 1:
    raise ValueError("Invalid value for dante_tir_min_multiplicity: must be a positive integer.")

# Hard ceiling on concurrent mclapply workers inside make_unified_annotation.
# Caps CONCURRENCY only -- batch composition is set from `threads`, so this can
# never change the annotation, only how many batches run at once. Needed because
# mclapply forks and R's GC dirties the inherited parent heap, so each worker's
# peak RSS converges on the PARENT's regardless of its batch size (measured on a
# 94 Gbp genome: 48.3 GB per worker for every one of 55 batches, against a 48.4 GB
# parent). At `threads: workflow.cores` = 96 that demanded ~2.66 TB on a 768 GB
# host: the run thrashed, 15 workers were OOM-killed and the rule died after
# 5h50m; the same work at 4 workers took 885 s. Tier resolution is only ~13% of
# this rule's wall time (loading and GFF3 export are serial), so a low ceiling
# costs almost nothing. 0 = no ceiling.
if "make_unified_max_workers" not in config:
    config["make_unified_max_workers"] = 8
if not isinstance(config["make_unified_max_workers"], int) or config["make_unified_max_workers"] < 0:
    raise ValueError("Invalid value for make_unified_max_workers: must be a non-negative integer (0 = no ceiling).")

# Hard ceiling on concurrent workers in the BigWig density rules
# (calculate_density_batch.R: make_bigwig_density and the two per-family rules).
# Default 0 = no hard ceiling, because those rules size their pool from a
# MEASUREMENT: the heaviest task runs first on its own and the memory budget is
# divided by its peak RSS. That gate is what the rules actually needed -- on the
# 94 Gbp run-000156 make_tidecluster_tandem_per_family_bigwig peaked at 568.5 GB
# of a 768 GB host, because every one of 3,206 concurrent tasks built the whole
# 94.26 M-bin genome grid. Both halves are fixed: density_track() now streams the
# grid one sequence at a time, and the pool is gated. Set > 0 to pin the worker
# count directly (e.g. non-cgroup HPC where the budget cannot be detected).
if "bigwig_max_workers" not in config:
    config["bigwig_max_workers"] = 0
if not isinstance(config["bigwig_max_workers"], int) or config["bigwig_max_workers"] < 0:
    raise ValueError("Invalid value for bigwig_max_workers: must be a non-negative integer (0 = no ceiling).")

# Memory allocated to this run, in GB — the number the scheduler granted
# (`qsub -l mem=128gb`, `--mem=128G`), not a per-rule figure. Default 0 = detect
# it (scripts/mem_utils.py resolution chain: AGENT_MEMORY, the scheduler
# environment, the cgroup limit, then MemAvailable).
#
# Why an explicit knob is worth having: CARP normally runs from a .sif, and
# inside a container /proc/meminfo reports the HOST — the kernel does not
# namespace it — while the cgroup limit that will actually kill the job usually
# sits on an ancestor job scope that may not be reachable (controllers not
# delegated, /sys/fs/cgroup not mounted, or a scheduler enforcing by polling).
# Detection then cannot distinguish "no limit" from "a limit I cannot see", and
# every memory-gated pool sizes itself against the whole node. This is
# TideCluster issue #6. `run_pipeline.py -m <GB>` sets this.
#
# It reaches the tools through `resources.mem_mb`, so an HPC profile or
# --set-resources still overrides it per rule. When unset MAX_MEMORY_MB is 0 and
# every rule below is byte-identical to before this knob existed.
if "max_memory_gb" not in config:
    config["max_memory_gb"] = 0
if (isinstance(config["max_memory_gb"], bool)
        or not isinstance(config["max_memory_gb"], (int, float))
        or config["max_memory_gb"] < 0):
    raise ValueError("Invalid value for max_memory_gb: must be a number >= 0 (0 = auto-detect).")
MAX_MEMORY_MB = int(config["max_memory_gb"] * 1024)

# Optional inclusion of DANTE_TIR_FALLBACK reps in the RepeatMasker
# library. Default OFF to preserve previous behaviour byte-for-byte.
# When ON, build_fallback_tir_library re-clusters fallback survivors,
# applies a Multiplicity floor, and drops any rep whose blast hits
# include a default-library entry of an incompatible classification
# (strict path-prefix: same path or one ancestor of the other; siblings
# like CACTA-vs-hAT count as incompatible). The fallback library is
# *not* used to filter the LTR library — it is treated as less
# reliable than the primary library on purpose.
if "include_dante_tir_fallback_in_library" not in config:
    config["include_dante_tir_fallback_in_library"] = False
if not isinstance(config["include_dante_tir_fallback_in_library"], bool):
    raise ValueError(
        "Invalid value for include_dante_tir_fallback_in_library: "
        "must be a boolean."
    )

# Multiplicity floor for fallback reps. None inherits
# dante_tir_min_multiplicity so the user-facing default matches the
# primary library exactly (3 unless overridden).
if "dante_tir_fallback_library_min_multiplicity" not in config:
    config["dante_tir_fallback_library_min_multiplicity"] = None
_v = config["dante_tir_fallback_library_min_multiplicity"]
if _v is not None and (not isinstance(_v, int) or _v < 1):
    raise ValueError(
        "Invalid value for dante_tir_fallback_library_min_multiplicity: "
        "must be null or a positive integer."
    )

# Memory cap for the parallel-CAP3 phase of `reduce_library`. LTR
# classes whose per-class input FASTA is at least this many bytes run
# sequentially in the rule's Phase 1; smaller classes run in parallel
# in Phase 2. Default 50 MB matches the threshold the Python rewrite
# was profiled against (peak resident drops from ~6 GB physical /
# ~30 GB snakemake-reported on 4-way parallel down to ~2 GB physical
# on this setting). Raise on machines with > 32 GB RAM available to
# trade memory for parallelism; lower on tighter setups.
if "reduce_library_max_parallel_bp" not in config:
    config["reduce_library_max_parallel_bp"] = 50_000_000
if (not isinstance(config["reduce_library_max_parallel_bp"], int)
        or config["reduce_library_max_parallel_bp"] < 1):
    raise ValueError(
        "Invalid value for reduce_library_max_parallel_bp: "
        "must be a positive integer (bytes)."
    )

# How many big (>= reduce_library_max_parallel_bp) Class_I/LTR CAP3 classes
# `reduce_library` runs CONCURRENTLY in its Phase 1b. CAP3 is single-threaded,
# so >1 overlaps the big-class assemblies that dominate Phase-1 wall time, at
# the cost of that many CAP3 working sets resident at once (a few hundred MB
# each). mmseqs2 classes stay strictly sequential regardless (8 GB floor).
# Output is byte-identical for any value. Default 4 suits the large-RAM HPC
# nodes this pipeline targets; set to 1 to restore strictly-sequential Phase 1
# on memory-constrained machines.
if "reduce_library_max_big_cap3_parallel" not in config:
    config["reduce_library_max_big_cap3_parallel"] = 4
if (not isinstance(config["reduce_library_max_big_cap3_parallel"], int)
        or config["reduce_library_max_big_cap3_parallel"] < 1):
    raise ValueError(
        "Invalid value for reduce_library_max_big_cap3_parallel: "
        "must be a positive integer."
    )

# ── TideCluster 1.16.0 knobs ─────────────────────────────────────────────────
# rDNA identification (run_all/tarean, default-on in TideCluster). When True the
# clustering GFF3 gains rDNA_type=45S|5S on rDNA TRCs, which
# make_unified_annotation.R surfaces as array-level rDNA_45S / rDNA_5S. Set
# False to pass --no_rdna and keep those arrays as generic Satellite TRCs.
if "tidecluster_detect_rdna" not in config:
    config["tidecluster_detect_rdna"] = True
if not isinstance(config["tidecluster_detect_rdna"], bool):
    raise ValueError("Invalid value for tidecluster_detect_rdna: must be a boolean.")

# Optional override for TideCluster's rDNA reference library (blastn subject for
# 45S/5S calling). Empty string → TideCluster's bundled data/rdna_library.fasta
# (validated on the OZ408684.1 calibration array). Point at the pipeline's richer
# data/rdna_library.fasta only if the bundled one under-calls on your taxon.
if "tidecluster_rdna_library" not in config:
    config["tidecluster_rdna_library"] = ""
if not isinstance(config["tidecluster_rdna_library"], str):
    raise ValueError("Invalid value for tidecluster_rdna_library: must be a string path.")

# Cross-TRC overlap resolution (default-on in TideCluster 1.16.0): the clustering
# GFF3 is made non-overlapping across satellite TRCs (dominant-TRC-wins). Set
# True here to pass --keep_overlaps and retain the raw overlapping regions.
if "tidecluster_keep_trc_overlaps" not in config:
    config["tidecluster_keep_trc_overlaps"] = False
if not isinstance(config["tidecluster_keep_trc_overlaps"], bool):
    raise ValueError("Invalid value for tidecluster_keep_trc_overlaps: must be a boolean.")

# Genome chunk size (bp) for the parallel, pooled RepeatMasker in tc_reannotate
# (TideCluster 1.16.0). Sequences below 2*chunk_size are only packed
# (byte-identical); larger ones split, with <0.15% masked-bp drift at the cuts.
# Default 50 Mb matches TideCluster's default and the pipeline's own RM wrapper.
if "tidecluster_chunk_size" not in config:
    config["tidecluster_chunk_size"] = 50_000_000
if (not isinstance(config["tidecluster_chunk_size"], int)
        or config["tidecluster_chunk_size"] < 1):
    raise ValueError("Invalid value for tidecluster_chunk_size: must be a positive integer (bytes).")


# Define path to cleaned genome (will be created by clean_genome_fasta rule)
genome_fasta_cleaned = F"{config['output_dir']}/genome_cleaned.fasta"

rule all:
    input:
        genome_fasta_cleaned,
        F"{config['output_dir']}/DANTE/DANTE.gff3",
        F"{config['output_dir']}/DANTE/DANTE_filtered.gff3",
        F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_final.gff3",
        F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_combined.gff3",
        F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3",
        F"{config['output_dir']}/DANTE_LTR/LTR_RTs_library.fasta",
        F"{config['output_dir']}/TideCluster/default/TideCluster_clustering.gff3",
        F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_clustering.gff3",
        F"{config['output_dir']}/TideCluster/default/RM_on_TideCluster_Library.gff3",
        F"{config['output_dir']}/TideCluster/TideCluster_clustering_default_and_short_merged.gff3",
        F"{config['output_dir']}/Tandem_repeats_TideCluster_split_by_family_bigwig/.done",
        F"{config['output_dir']}/Libraries/class_ii_library.fasta",
        F"{config['output_dir']}/Libraries/LTR_RTs_library_clean.fasta",
        F"{config['output_dir']}/Libraries/combined_library.fasta",
        F"{config['output_dir']}/RepeatMasker/RM_on_combined_library.out",
        F"{config['output_dir']}/RepeatMasker/RM_on_combined_library_plus_DANTE.gff3",
        F"{config['output_dir']}/RepeatMasker/Repeat_Annotation_NoSat.gff3",
        F"{config['output_dir']}/all_repeats_for_masking.bed",
        F"{config['output_dir']}/DANTE_LTR.gff3",
        F"{config['output_dir']}/TideCluster_report.html",
        F"{config['output_dir']}/DANTE_LTR_report.html",
        F"{config['output_dir']}/gaps_10plus.bed",
        F"{config['output_dir']}/summary_statistics.csv",
        F"{config['output_dir']}/Repeat_density/Repeat_density_total_10k.bw",
        F"{config['output_dir']}/Repeat_density/Repeat_density_total_100k.bw",
        F"{config['output_dir']}/Tandem_repeats_TideCluster_10k.bw",
        F"{config['output_dir']}/Repeat_density_by_class_bigwig/.done",
        F"{config['output_dir']}/Tandem_repeats_unified_split_by_family_bigwig/.done",
        F"{config['output_dir']}/summary_plots.pdf",
        F"{config['output_dir']}/benchmark_report.html",
        F"{config['output_dir']}/repeat_annotation_report.html",
        F"{config['output_dir']}/Libraries/library_health.tsv",
        F"{config['output_dir']}/Repeat_Annotation_Unified.gff3",
        F"{config['output_dir']}/.classifications_validated"

rule clean_genome_fasta:
    """
    Clean FASTA headers by removing everything after the first whitespace character.
    Accepts either plain FASTA or gzip-compressed input. Detection is by the
    file's leading bytes (gzip magic 1f 8b), not by filename extension — so a
    gzipped file named `.fa` or `.fasta` is still handled correctly. Output is
    always plain `genome_cleaned.fasta` because every downstream tool (DANTE,
    RepeatMasker, TideCluster, …) wants an unzipped reference.
    Cleaning headers ensures consistent sequence IDs across all downstream
    analyses and prevents issues with tools that handle whitespace differently
    in FASTA headers.
    """
    input:
        config["genome_fasta"]
    output:
        genome_fasta_cleaned
    log:
        stdout=F"{config['output_dir']}/logs/clean_genome_fasta.log",
        stderr=F"{config['output_dir']}/logs/clean_genome_fasta.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/clean_genome_fasta.tsv"
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        # Sniff the first two bytes — gzip files start with 0x1f 0x8b
        # regardless of extension. `gzip -t` then fails fast on corruption
        # so we don't silently produce truncated output.
        magic=$(head -c2 {input} | od -An -t x1 | tr -d ' \n')
        if [ "$magic" = "1f8b" ]; then
            gzip -t {input}
            READER="gzip -dc {input}"
        else
            READER="cat {input}"
        fi
        # Clean FASTA headers - keep only ID before first whitespace
        eval "$READER" | awk '/^>/ {{split($1, a, " "); print a[1]; next}} {{print}}' > {output}

        # Sanity check: output must start with the FASTA header byte '>' and
        # MUST NOT carry the gzip magic 1f 8b. Catches the case where a
        # pre-fix container or a stale leftover slipped a binary stream
        # through awk into a file named .fasta — DANTE / RepeatMasker /
        # everything downstream then fails far away from the cause.
        first_byte=$(head -c1 {output} | od -An -t x1 | tr -d ' \n')
        gz_magic=$(head -c2 {output} | od -An -t x1 | tr -d ' \n')
        if [ "$gz_magic" = "1f8b" ]; then
            echo "ERROR: {output} is gzipped (magic 1f 8b). The clean_genome_fasta rule should always emit plain FASTA." >&2
            rm -f {output}
            exit 1
        fi
        if [ "$first_byte" != "3e" ]; then
            echo "ERROR: {output} does not start with a FASTA header ('>'); first byte is 0x$first_byte." >&2
            rm -f {output}
            exit 1
        fi
        """

rule index_genome_fasta:
    """Create FASTA index (.fai) for the cleaned genome."""
    input:
        genome_fasta_cleaned
    output:
        F"{config['output_dir']}/genome_cleaned.fasta.fai"
    conda:
        "envs/tidecluster.yaml"
    shell:
        """
        seqkit faidx {input}
        """

rule dante:
    input:
        genome_fasta_cleaned,
    output:
        F"{config['output_dir']}/DANTE/DANTE.gff3"
    log:
        stdout=F"{config['output_dir']}/DANTE/dante.log",
        stderr=F"{config['output_dir']}/DANTE/dante.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/dante.tsv"
    conda:
        # dante has its own env since 0.2.11 (r-base 4.2.3) — it can no longer
        # share envs/tidecluster.yaml with dante_ltr (r-base <4.2). See envs/dante.yaml.
        "envs/dante.yaml"
    threads: workflow.cores
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        dante -q {input} -o {output} -c {threads}
        """

rule dante_tir:
    input:
        gff=F"{config['output_dir']}/DANTE/DANTE.gff3",
        fasta=genome_fasta_cleaned
    output:
        checkpoint=F"{config['output_dir']}/DANTE_TIR/.done",
        gff=F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_final.gff3",
        fasta=F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_final.fasta",
        summary=F"{config['output_dir']}/DANTE_TIR/TIR_classification_summary.txt",
        dante_tir_lib=F"{config['output_dir']}/DANTE_TIR/all_representative_elements_min3.fasta"
    params:
        output_dir=F"{config['output_dir']}/DANTE_TIR",
        cap3_max_memory_gb=config["dante_tir_cap3_max_memory_gb"]
    log:
        stdout=F"{config['output_dir']}/DANTE_TIR/dante_tir.log",
        stderr=F"{config['output_dir']}/DANTE_TIR/dante_tir.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/dante_tir.tsv"
    conda:
        "envs/dante_tir.yaml"
    threads: workflow.cores
    # mem_mb comes from config max_memory_gb (0 = unset); an HPC profile /
    # --set-resources overrides it per rule. When set, the shell caps DANTE_TIR's
    # CAP3 memory budget at 60% of it so CAP3 respects the job allocation even on
    # schedulers that enforce no cgroup limit.
    resources:
        mem_mb=MAX_MEMORY_MB
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        # Resolve the CAP3 memory budget passed to DANTE_TIR 0.3.0's --cap3_max_memory
        # (GB). Precedence: the explicit per-tool knob wins, then 60% of the memory
        # allocation (resources.mem_mb, from config max_memory_gb or a profile),
        # else pass nothing (0.3.0 then auto-detects the cgroup/node memory itself,
        # which over-detects inside a container). 0 => omit the flag.
        # The knob is checked FIRST so the global max_memory_gb cannot silently
        # override a budget the user pinned for this tool specifically.
        cap3_mem=0
        if [ {params.cap3_max_memory_gb} -gt 0 ]; then
            cap3_mem={params.cap3_max_memory_gb}
        elif [ {resources.mem_mb} -gt 0 ]; then
            cap3_mem=$(( {resources.mem_mb} * 6 / 10 / 1024 ))
        fi
        cap3_arg=""
        [ "$cap3_mem" -gt 0 ] && cap3_arg="--cap3_max_memory $cap3_mem"
        echo "DANTE_TIR CAP3 memory budget: ${{cap3_mem}} GB (0 = tool auto-detects)"
        # Run dante_tir.py and check exit status
        if dante_tir.py -g {input.gff} -f {input.fasta} -o {params.output_dir} -c {threads} $cap3_arg; then
            # dante_tir.py succeeded
            echo "DANTE_TIR completed successfully"

            # Check if DANTE_TIR_final.fasta exists and is not empty
            if [ -s {params.output_dir}/DANTE_TIR_final.fasta ]; then
                echo "Running dante_tir_summary.R on non-empty results"
                dante_tir_summary.R -g {params.output_dir}/DANTE_TIR_final.gff3 -f {input.fasta} -o {params.output_dir}
            else
                echo "No TIR elements found - skipping summary step"
            fi
        else
            echo "DANTE_TIR failed with non-zero exit status"
        fi

        # Ensure all expected output files exist (create empty ones if needed).
        # All five declared outputs must be touched — when DANTE_TIR finds no
        # elements, it skips producing the representative library, and omitting
        # it here trips snakemake's MissingOutputException on small inputs.
        touch {output.gff} {output.fasta} {output.summary} {output.dante_tir_lib}

        # Inject the DANTE TPase protein_domain row as a child of each
        # DANTE_TIR sequence_feature parent. Idempotent: a future dante_tir
        # release that emits TPase children already will leave them in
        # place; this step then becomes a no-op and can be retired.
        # Bare-name invocation via PATH — see CLAUDE.md "Calling helper
        # scripts from rules" for the dual-context contract.
        scripts_dir=$(realpath scripts)
        export PATH="$scripts_dir:$PATH"
        enrich_dante_tir_with_tpase.py \
            --dante-gff {input.gff} \
            --dante-tir-gff {output.gff}

        # Create checkpoint file to indicate completion
        touch {output.checkpoint}
        """


rule dante_tir_fallback:
    """
    Fallback TIR detection using TPase domain flanking-region analysis.
    Identifies partial TIR elements that DANTE_TIR may have missed.
    """
    input:
        gff=F"{config['output_dir']}/DANTE/DANTE.gff3",
        genome=genome_fasta_cleaned,
        dante_tir_checkpoint=F"{config['output_dir']}/DANTE_TIR/.done",
        mask_gff=F"{config['output_dir']}/TideCluster/default/TideCluster_tidehunter.gff3",
        # Same flank masks as dante_line: a TPase flank must not run through an
        # element another tool has already delimited structurally. The primary
        # DANTE_TIR set is the right mask here -- masking on this rule's OWN
        # (fallback) output would be circular.
        mask_dante_ltr=F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3",
        mask_dante_tir=F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_final.gff3"
    output:
        gff=F"{config['output_dir']}/DANTE_TIR_FALLBACK/DANTE_TIR_FALLBACK.gff3",
        rep_lib=F"{config['output_dir']}/DANTE_TIR_FALLBACK/TIR_fallback_rep_lib.fasta",
        extended_fasta=F"{config['output_dir']}/DANTE_TIR_FALLBACK/TIR_fallback_extended.fasta"
    params:
        output_dir=F"{config['output_dir']}/DANTE_TIR_FALLBACK",
        min_alignments=config["dante_tir_fallback_min_alignments"],
        min_cluster_size=config["dante_tir_fallback_min_cluster_size"],
        max_group_size=config["dante_tir_fallback_max_group_size"],
        support_fraction=config["dante_tir_fallback_support_fraction"],
        min_group_alignments=config["dante_tir_fallback_min_group_alignments"],
        library_source=config["dante_tir_fallback_library_source"]
    log:
        stdout=F"{config['output_dir']}/DANTE_TIR_FALLBACK/dante_tir_fallback.log",
        stderr=F"{config['output_dir']}/DANTE_TIR_FALLBACK/dante_tir_fallback.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/dante_tir_fallback.tsv"
    conda:
        "envs/dante_line.yaml"
    threads: workflow.cores
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x

        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH

        dante_tir_fallback.py \
            -g {input.genome} \
            -a {input.gff} \
            -o {params.output_dir} \
            -t {threads} \
            --min-num-alignments {params.min_alignments} \
            --min-cluster-size {params.min_cluster_size} \
            --max-group-size {params.max_group_size} \
            --support-fraction {params.support_fraction} \
            --min-group-alignments {params.min_group_alignments} \
            --library-source {params.library_source} \
            --mask-gff3 {input.mask_gff} \
            --mask-gff3 {input.mask_dante_ltr} \
            --mask-gff3 {input.mask_dante_tir}

        # Ensure outputs exist even if no TIR elements were found
        touch {output.gff} {output.rep_lib} {output.extended_fasta}
        """


rule merge_dante_tir_with_fallback:
    """
    Merge primary DANTE_TIR and fallback annotations.
    Fallback elements overlapping any primary element are discarded.
    Surviving fallback elements are labeled as partial.
    """
    input:
        primary_gff=F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_final.gff3",
        fallback_gff=F"{config['output_dir']}/DANTE_TIR_FALLBACK/DANTE_TIR_FALLBACK.gff3",
        fallback_fasta=F"{config['output_dir']}/DANTE_TIR_FALLBACK/TIR_fallback_extended.fasta"
    output:
        combined_gff=F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_combined.gff3",
        filtered_fallback_fasta=F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_fallback_filtered.fasta"
    log:
        stdout=F"{config['output_dir']}/DANTE_TIR/merge_tir_fallback.log",
        stderr=F"{config['output_dir']}/DANTE_TIR/merge_tir_fallback.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/merge_dante_tir_with_fallback.tsv"
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x

        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH

        merge_tir_fallback.py \
            --primary-gff {input.primary_gff} \
            --fallback-gff {input.fallback_gff} \
            --fallback-fasta {input.fallback_fasta} \
            --output-gff {output.combined_gff} \
            --output-fasta {output.filtered_fallback_fasta}
        """


rule make_tir_combined_library:
    """
    Build the TIR repeat library from primary DANTE_TIR elements only.

    DANTE_TIR_FALLBACK partial elements are intentionally excluded — they
    remain visible in DANTE_TIR_combined.gff3 (and the unified annotation)
    as low-confidence partials, but are not trusted enough to seed
    RepeatMasker. The library is built by mmseqs2-clustering the primary
    DANTE_TIR_final.fasta after canonicalising headers and (optionally)
    dropping primaries below the configured Multiplicity floor.
    """
    input:
        primary_fasta=F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_final.fasta",
        primary_gff=F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_final.gff3"
    output:
        combined_lib=F"{config['output_dir']}/DANTE_TIR/all_representative_elements_combined.fasta"
    params:
        mmseqs_dir=F"{config['output_dir']}/DANTE_TIR/mmseqs_combined",
        min_multiplicity=config["dante_tir_min_multiplicity"]
    log:
        stdout=F"{config['output_dir']}/DANTE_TIR/make_tir_combined_library.log",
        stderr=F"{config['output_dir']}/DANTE_TIR/make_tir_combined_library.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/make_tir_combined_library.tsv"
    conda:
        "envs/dante_line.yaml"
    threads: workflow.cores
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x

        mkdir -p {params.mmseqs_dir}

        # Bare-name invocation via PATH — see CLAUDE.md "Calling helper
        # scripts from rules" for the dual-context contract.
        scripts_dir=$(realpath scripts)
        export PATH="$scripts_dir:$PATH"

        # Canonicalise primary DANTE_TIR FASTA headers into canonical slash
        # form (#Class_II/Subclass_1/TIR/hAT) via classification_vocabulary.yaml
        # so any unknown leaf fails loudly here. Fallback sequences are
        # deliberately not included — see rule docstring.
        CANON_INPUT={params.mmseqs_dir}/primary_canonical.fasta
        classification.py canonicalise-fasta-headers \
            --source DANTE_TIR {input.primary_fasta} "$CANON_INPUT"

        # Optional Multiplicity floor: drop primary elements whose parent
        # row in DANTE_TIR_final.gff3 has Multiplicity < threshold. Default
        # threshold is 1 (no filter). This affects only the library; the
        # GFF still carries every primary element.
        FILTERED_INPUT={params.mmseqs_dir}/primary_filtered.fasta
        if [ "{params.min_multiplicity}" -gt 1 ] && [ -s "$CANON_INPUT" ] && [ -s {input.primary_gff} ]; then
            filter_dante_tir_by_multiplicity.py \
                --gff {input.primary_gff} \
                --fasta-in "$CANON_INPUT" \
                --fasta-out "$FILTERED_INPUT" \
                --min-multiplicity {params.min_multiplicity}
        else
            cp "$CANON_INPUT" "$FILTERED_INPUT"
        fi

        # If filtered input is empty, create empty output
        if [ ! -s "$FILTERED_INPUT" ]; then
            touch {output.combined_lib}
            exit 0
        fi

        # Deterministic clustering: mmseqs easy-cluster is order-sensitive
        # (same sequences in a different order -> different representatives and
        # a different cluster count), so canonically sort the input by sequence
        # first. See scripts/canonical_sort_fasta.py.
        SORTED_INPUT={params.mmseqs_dir}/primary_sorted.fasta
        canonical_sort_fasta.py "$FILTERED_INPUT" "$SORTED_INPUT" \
            --threads {threads} --tmpdir {params.mmseqs_dir}

        # Run mmseqs2 clustering on the primary-only set
        mmseqs easy-cluster \
            "$SORTED_INPUT" \
            {params.mmseqs_dir}/cluster \
            {params.mmseqs_dir}/tmp \
            --threads {threads}

        # Use cluster representatives as the library
        if [ -s {params.mmseqs_dir}/cluster_rep_seq.fasta ]; then
            cp {params.mmseqs_dir}/cluster_rep_seq.fasta {output.combined_lib}
        else
            touch {output.combined_lib}
        fi
        """


rule filter_dante:
    input:
        F"{config['output_dir']}/DANTE/DANTE.gff3"
    output:
        gff=F"{config['output_dir']}/DANTE/DANTE_filtered.gff3",
        fasta=F"{config['output_dir']}/DANTE/DANTE_filtered.fasta"
    log:
        stdout=F"{config['output_dir']}/DANTE/filter_dante.log",
        stderr=F"{config['output_dir']}/DANTE/filter_dante.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/filter_dante.tsv"
    conda:
        # dante_gff_output_filtering.py ships with the DANTE package, which now
        # lives in its own env (dante 0.2.11, r-base 4.2.3). See envs/dante.yaml.
        "envs/dante.yaml"
    threads: 1
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        dante_gff_output_filtering.py --dom_gff {input} --domains_filtered {output.gff} --domains_prot_seq {output.fasta}
        """


rule validate_classifications:
    """
    Fail fast if any upstream tool emits a classification not listed in
    classification_vocabulary.yaml. Runs after every producer of a
    classification-bearing file; consumers (library construction,
    RepeatMasker, unified annotation) depend on this rule's marker.

    Also validates the OPTIONAL user-supplied libraries (`custom_library`,
    `tandem_repeat_library`). Those are inputs, not pipeline products, so
    nothing else canonicalises them — `concatenate_libraries` cats the custom
    library into the RepeatMasker library verbatim. A non-vocabulary class
    there flows straight through RepeatMasker into
    Repeat_Annotation_Unified.gff3 and summary_statistics.csv, splitting one
    class into parallel rows that do not sum (observed: a curated library
    carrying `rDNA_45S/18S` put 9 such records beside the canonical
    `rDNA/45S_rDNA/*` ones). Caught here, before the library is built.
    """
    input:
        dante_filtered=F"{config['output_dir']}/DANTE/DANTE_filtered.gff3",
        dante_ltr=F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3",
        dante_tir_final=F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_final.gff3",
        dante_tir_combined=F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_combined.gff3",
        dante_line=F"{config['output_dir']}/DANTE_LINE/DANTE_LINE.gff3"
    params:
        custom_library=config.get("custom_library", ""),
        tandem_repeat_library=config.get("tandem_repeat_library", "")
    output:
        marker=F"{config['output_dir']}/.classifications_validated"
    log:
        stdout=F"{config['output_dir']}/logs/validate_classifications.log",
        stderr=F"{config['output_dir']}/logs/validate_classifications.err"
    conda:
        "envs/tidecluster.yaml"
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        # Bare-name invocation via PATH — see CLAUDE.md "Calling helper
        # scripts from rules" for the dual-context contract.
        scripts_dir=$(realpath scripts)
        export PATH="$scripts_dir:$PATH"
        CLS="classification.py validate --mode gff3"
        $CLS --source DANTE       --attribute Final_Classification {input.dante_filtered}
        $CLS --source DANTE_LTR   --attribute Final_Classification {input.dante_ltr}
        $CLS --source DANTE_TIR   --attribute Classification       {input.dante_tir_final}
        $CLS --source DANTE_TIR   --attribute Classification       {input.dante_tir_combined}
        $CLS --source DANTE_LINE  --attribute Final_Classification {input.dante_line}

        # User-supplied libraries: RepeatMasker `name#class` headers. Optional,
        # so guard on set-and-non-empty (same pattern as the other rules that
        # take these). --mode fasta reads the class from each header.
        CLSF="classification.py validate --mode fasta --source RepeatMasker"
        if [ -n "{params.custom_library}" ] && [ -s "{params.custom_library}" ]; then
            $CLSF "{params.custom_library}"
        fi
        if [ -n "{params.tandem_repeat_library}" ] && [ -s "{params.tandem_repeat_library}" ]; then
            $CLSF "{params.tandem_repeat_library}"
        fi
        touch {output.marker}
        """

rule dante_line:
    input:
        gff=F"{config['output_dir']}/DANTE/DANTE_filtered.gff3",
        gff3_tidehunter=F"{config['output_dir']}/TideCluster/default/TideCluster_tidehunter.gff3",
        # Flank masks. A LINE flank must not run through an element another tool
        # has already delimited structurally, so DANTE_LTR elements and primary
        # DANTE_TIR elements terminate it alongside the TideHunter arrays.
        # DANTE_TIR_final.gff3 (primary) rather than DANTE_TIR_combined.gff3:
        # the fallback partials come from this same flank-extension engine, so
        # masking on them would propagate a boundary of the kind being fixed.
        gff3_dante_ltr=F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3",
        gff3_dante_tir=F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_final.gff3",
        genome=genome_fasta_cleaned
    output:
        line_rep_lib=F"{config['output_dir']}/DANTE_LINE/LINE_rep_lib.fasta",
        gff_out=F"{config['output_dir']}/DANTE_LINE/DANTE_LINE.gff3",
        line_regions=F"{config['output_dir']}/DANTE_LINE/LINE_regions.fasta",
        line_regions_extended=F"{config['output_dir']}/DANTE_LINE/LINE_regions_extended.fasta",
        # Sequences of the elements marked Status=complete below. Always created
        # (empty when the feature is off or nothing qualifies) so the DAG is
        # independent of the knob.
        complete_fasta=F"{config['output_dir']}/DANTE_LINE/LINE_complete_elements.fasta"
    params:
        output_dir=F"{config['output_dir']}/DANTE_LINE",
        max_group_size=config["dante_line_max_group_size"],
        support_fraction=config["dante_line_support_fraction"],
        min_group_alignments=config["dante_line_min_group_alignments"],
        max_extension=config["dante_line_max_extension"],
        max_element_length=config["dante_line_max_element_length"],
        library_source=config["dante_line_library_source"],
        mark_complete=config["line_complete_elements"]
    log:
        stdout=F"{config['output_dir']}/DANTE_LINE/dante_line.log",
        stderr=F"{config['output_dir']}/DANTE_LINE/dante_line.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/dante_line.tsv"
    priority: 50
    conda:
        "envs/dante_line.yaml"
    threads: workflow.cores
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x

        # Add scripts directory to PATH
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH

        # dante_line.py exits 3 when the genome legitimately has too few LINE
        # features / no valid domain patterns (common on small CI fixtures). Only
        # that benign case is turned into empty outputs so downstream continues;
        # ANY OTHER non-zero exit is a real failure and must stop the pipeline
        # (previously every failure was swallowed into empty LINE outputs, so a
        # crash silently dropped the entire LINE layer genome-wide).
        mkdir -p {params.output_dir}
        rc=0
        dante_line.py -g {input.genome} -a {input.gff} \
                -o {params.output_dir} -t {threads} \
                --max-group-size {params.max_group_size} \
                --support-fraction {params.support_fraction} \
                --min-group-alignments {params.min_group_alignments} \
                --max-extension {params.max_extension} \
                --max-element-length {params.max_element_length} \
                --library-source {params.library_source} \
                --mask-gff3 {input.gff3_tidehunter} \
                --mask-gff3 {input.gff3_dante_ltr} \
                --mask-gff3 {input.gff3_dante_tir} || rc=$?
        if [ "$rc" -eq 3 ]; then
            echo "dante_line.py: no LINE content (too few features); creating empty outputs"
        elif [ "$rc" -ne 0 ]; then
            echo "dante_line.py FAILED with exit code $rc" >&2
            exit "$rc"
        fi
        [ -f {output.line_rep_lib} ]          || : > {output.line_rep_lib}
        [ -f {output.gff_out} ]               || echo "##gff-version 3" > {output.gff_out}
        [ -f {output.line_regions} ]          || : > {output.line_regions}
        [ -f {output.line_regions_extended} ] || : > {output.line_regions_extended}

        # Mark the elements whose BOTH ends were directly observed -- a poly-A
        # tail plus a target-site duplication -- and replace their inferred span
        # with the measured one. Runs inside this rule because it rewrites this
        # rule's own GFF3; a separate rule cannot write another rule's output.
        # Costs ~4 s on a 5.4 Gb genome.
        if [ "{params.mark_complete}" = "True" ] && [ -s {output.gff_out} ]; then
            line_complete_elements.py \
                -g {input.genome} \
                -a {output.gff_out} \
                -o {output.complete_fasta} \
                --annotate-gff3
        fi
        [ -f {output.complete_fasta} ] || : > {output.complete_fasta}
        """

rule dante_ltr:
    input:
        fasta=genome_fasta_cleaned,
        gff=F"{config['output_dir']}/DANTE/DANTE.gff3"

    output:
        gff = F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3",
        html = F"{config['output_dir']}/DANTE_LTR/DANTE_LTR_summary.html"
    params:
        prefix = lambda wildcards, output: output.gff.replace(".gff3", "")
    log:
        stdout=F"{config['output_dir']}/DANTE_LTR/dante_ltr.log",
        stderr=F"{config['output_dir']}/DANTE_LTR/dante_ltr.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/dante_ltr.tsv"
    conda:
        "envs/tidecluster.yaml"
    threads: workflow.cores
    # dante_ltr 0.5.4.0 gates its chunk pool at 0.8 * budget / per-chunk peak RSS.
    # Left to itself it resolves the budget (scheduler env -> cgroup ->
    # MemAvailable), but MemAvailable is the NODE's inside a .sif — upstream issue
    # #13 — and pool_size is linear in it, so pass the allocation when we know it.
    resources:
        mem_mb=MAX_MEMORY_MB
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        # STOPGAP (large genomes): dante_ltr opens one temp-file handle per chunk
        # simultaneously (int(total_size / -S) handles, ~1800 at 90 Gbp / 50 Mb),
        # which exceeds the default soft open-file limit of 1024 and aborts with
        # "OSError: [Errno 24] Too many open files". Raise the soft limit to the
        # hard limit for this rule. Durable fix is upstream —
        # see docs/archive/dante_ltr_too_many_open_files_request.md.
        ulimit -n "$(ulimit -Hn)" || true
        # --max_memory takes GB; below 1 GB pass nothing and let dante_ltr resolve
        # its own budget (its chain matches scripts/mem_utils.py).
        max_mem_arg=""
        if [ {resources.mem_mb} -ge 1024 ]; then
            max_mem_arg="--max_memory $(( {resources.mem_mb} / 1024 ))"
        fi
        dante_ltr -o {params.prefix} -s {input.fasta} -g {input.gff} -c {threads} -M 1 -S 50000000 $max_mem_arg
        # if exit status is 0 and gff3 file was created but html is missing, create an empty file
        echo "DANTE LTR-RTs finished"
        if [ -f {output.gff} ]; then
            if [ ! -f {output.html} ]; then
                echo "Creating an empty html file"
                echo "No complete LTR-RTs found" > {output.html}
            fi
        fi
        """


rule resolve_ltr_tandems:
    """
    Detect tandem LTR-RT (LTR_RT_TR) arrays — head-to-tail same-lineage LTR-RTs
    sharing a boundary LTR (Macko-Podgorni et al., Mobile DNA 2025) — and write a
    SMALL side file of just the derived containers (one per array, with the member
    element IDs). DANTE_LTR.gff3 is read-only and UNTOUCHED, so the LTR library,
    masking track, repeat report and dante_line keep seeing every individual
    element. make_unified_annotation reads DANTE_LTR.gff3 + this side file and
    emits one Level-1 container per array with the member copies nested as Level-2
    children. A genome with no tandems yields a header-only (empty) file.
    See scripts/resolve_ltr_tandems.py.
    """
    input:
        gff=F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3"
    output:
        gff=F"{config['output_dir']}/DANTE_LTR/DANTE_LTR_tandems.gff3"
    log:
        stdout=F"{config['output_dir']}/DANTE_LTR/resolve_ltr_tandems.log",
        stderr=F"{config['output_dir']}/DANTE_LTR/resolve_ltr_tandems.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/resolve_ltr_tandems.tsv"
    conda:
        "envs/tidecluster.yaml"
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        scripts_dir=$(realpath scripts)
        export PATH="$scripts_dir:$PATH"
        resolve_ltr_tandems.py -i {input.gff} -o {output.gff}
        """


rule make_library_of_ltrs:
    input:
        gff3=F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3",
        genome_fasta=genome_fasta_cleaned
    output:
        dir=directory(F"{config['output_dir']}/DANTE_LTR/library"),
        fasta=F"{config['output_dir']}/DANTE_LTR/LTR_RTs_library.fasta"
    log:
        stdout=F"{config['output_dir']}/DANTE_LTR/make_library_of_ltrs.log",
        stderr=F"{config['output_dir']}/DANTE_LTR/make_library_of_ltrs.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/make_library_of_ltrs.tsv"
    conda:
        "envs/tidecluster.yaml"
    threads: workflow.cores
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        # run only if gff3 contains some records
        # (number of lines not starting with # is greater than 1)
        # but check just first 30 lines
        mkdir -p {output.dir}
        if [ $(head -n 30 {input.gff3} | grep -v "^#" | wc -l) -gt 1 ]; then
            # dante_ltr_to_library can fail on small inputs (mmseq_clustering.R
            # raises "attempt to set too many names (1) on GroupedIRanges
            # object of length 0" when a cluster has 1 member). In that case
            # we still want the pipeline to continue — create an empty library
            # and carry on; RepeatMasker will simply see no LTR consensi.
            if dante_ltr_to_library --gff {input.gff3} --output_dir {output.dir} -s {input.genome_fasta} -c {threads}; then
                ln -sf library/mmseqs2/mmseqs_representative_seq_clean.fasta {output.fasta}
            else
                echo "dante_ltr_to_library failed (too few LTRs to cluster); creating empty library"
                : > {output.fasta}
            fi
        else
            echo "No LTR-RTs found, creating an empty file"
            : > {output.fasta}
        fi
        """


rule tidecluster_long:
    input:
        genome_fasta=genome_fasta_cleaned,
        library= config.get("tandem_repeat_library", []),
        genome_seqlengths=F"{config['output_dir']}/genome_seqlengths.rds"
    output:
        gff3_clust=F"{config['output_dir']}/TideCluster/default/TideCluster_clustering.gff3",
        gff3_tidehunter=F"{config['output_dir']}/TideCluster/default/TideCluster_tidehunter.gff3",
        dimer_library_default=F"{config['output_dir']}/TideCluster/default/TideCluster_consensus_dimer_library.fasta",
        tr_default_short=F"{config['output_dir']}/TideCluster/default/TideCluster_tidehunter_short.gff3",
        html=F"{config['output_dir']}/TideCluster/default/TideCluster_index.html",
        split_files=directory(F"{config['output_dir']}/TideCluster/default/TideCluster_clustering_split_files")
    params:
        prefix = lambda wildcards, output: output.gff3_clust.replace("_clustering.gff3", ""),
        library = config.get("tandem_repeat_library", ""),
        rdna_flag = "" if config["tidecluster_detect_rdna"] else "--no_rdna",
        overlaps_flag = "--keep_overlaps" if config["tidecluster_keep_trc_overlaps"] else "",
        rdna_library = config["tidecluster_rdna_library"]
    log:
        stdout=F"{config['output_dir']}/TideCluster/default/tidecluster_long.log",
        stderr=F"{config['output_dir']}/TideCluster/default/tidecluster_long.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/tidecluster_long.tsv"
    conda:
        "envs/tidecluster_run.yaml"
    threads: workflow.cores
    # TideCluster 1.20.0 sizes its TideHunter pool and TAREAN thread cap from a
    # memory budget. Left to itself it resolves one (scheduler env -> cgroup ->
    # MemAvailable), but MemAvailable is the NODE's inside a .sif — issue #6 —
    # so pass the job's allocation explicitly whenever we know it.
    resources:
        mem_mb=MAX_MEMORY_MB

    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        # STOPGAP (large genomes): TideCluster's chunked run opens one handle per
        # chunk (~genome_bp/chunk_size, ~1800 at 90 Gbp) and, with >10k contigs,
        # more — exceeding the default 1024 soft limit ("Too many open files").
        # Same class as dante_ltr/tc_reannotate. Raise soft→hard for this rule.
        ulimit -n "$(ulimit -Hn)" || true
        wd=$(dirname {output.gff3_clust})
        prefix=$(basename {params.prefix})
        original_dir=$PWD
        genome_absolute_path=$(realpath {input.genome_fasta})
        genome_seqlengths=$(realpath {input.genome_seqlengths})
        # TideCluster 1.16.0 knobs (empty by default → TideCluster's own defaults).
        rdna_lib_arg=""
        if [ -n "{params.rdna_library}" ]; then
            rdna_lib_arg="--rdna_library $(realpath {params.rdna_library})"
        fi
        # --max_memory takes GB; below 1 GB pass nothing and let TideCluster
        # resolve its own budget (its chain matches scripts/mem_utils.py).
        max_mem_arg=""
        if [ {resources.mem_mb} -ge 1024 ]; then
            max_mem_arg="--max_memory $(( {resources.mem_mb} / 1024 ))"
        fi
        tc_extra="{params.rdna_flag} {params.overlaps_flag} $rdna_lib_arg $max_mem_arg"
        # NOTE - there is a bug in tidecluster - it does not correctly format html links; solution is
        # to run it in the directory where the output will be created
        echo "Library: {input.library}"
        if [ -z "{params.library}" ]; then
            cd $wd
            echo "Running TideCluster without a library"
            TideCluster.py run_all -pr $prefix -c {threads}  -f $genome_absolute_path --long $tc_extra
        else
            library_absolute_path=$(realpath {params.library})
            echo "Running TideCluster with a custom library"
            cd $wd
            TideCluster.py run_all -pr $prefix -c {threads} -f $genome_absolute_path -l $library_absolute_path --long $tc_extra
        fi
        # TideCluster may not create any of its outputs when TideHunter finds
        # zero candidates on a low-satellite genome (common for small CI
        # fixtures). Create empty stubs for every declared output so the
        # rule satisfies snakemake's missing-output check in all cases.
        cd $original_dir
        [ -f {output.gff3_clust} ]            || echo "##gff-version 3" > {output.gff3_clust}
        [ -f {output.gff3_tidehunter} ]       || echo "##gff-version 3" > {output.gff3_tidehunter}
        [ -f {output.tr_default_short} ]      || echo "##gff-version 3" > {output.tr_default_short}
        [ -f {output.dimer_library_default} ] || : > {output.dimer_library_default}
        [ -f {output.html} ]                  || : > {output.html}
        cd $wd
        # Ensure the per-cluster split-files directory exists so the declared
        # directory() output is satisfied even when TideCluster produced no
        # clusters. Per-family density BigWigs are built downstream by the
        # dedicated make_tidecluster_tandem_per_family_bigwig rule.
        mkdir -p TideCluster_clustering_split_files
        """


rule tidecluster_short:
    input:
        genome_fasta=genome_fasta_cleaned,
    output:
        gff3_clust=F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_clustering.gff3",
        gff3_tidehunter=F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_tidehunter.gff3",
        dimer_library_short=F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_consensus_dimer_library.fasta",
        tr_short_short=F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_tidehunter_short.gff3"
    params:
        prefix = lambda wildcards, output: output.gff3_clust.replace("_clustering.gff3", ""),
        library = config.get("tandem_repeat_library", ""),
        rdna_flag = "" if config["tidecluster_detect_rdna"] else "--no_rdna",
        overlaps_flag = "--keep_overlaps" if config["tidecluster_keep_trc_overlaps"] else "",
        rdna_library = config["tidecluster_rdna_library"]
    log:
        stdout=F"{config['output_dir']}/TideCluster/short_monomer/tidecluster_short.log",
        stderr=F"{config['output_dir']}/TideCluster/short_monomer/tidecluster_short.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/tidecluster_short.tsv"
    conda:
        "envs/tidecluster_run.yaml"
    threads: workflow.cores
    # TideCluster 1.20.0 sizes its TideHunter pool and TAREAN thread cap from a
    # memory budget. Left to itself it resolves one (scheduler env -> cgroup ->
    # MemAvailable), but MemAvailable is the NODE's inside a .sif — issue #6 —
    # so pass the job's allocation explicitly whenever we know it.
    resources:
        mem_mb=MAX_MEMORY_MB

    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        # STOPGAP (large genomes): see tidecluster_long — chunked run + >10k contigs
        # can exceed the default 1024 open-file soft limit. Raise soft→hard.
        ulimit -n "$(ulimit -Hn)" || true
        wd=$(dirname {output.gff3_clust})
        prefix=$(basename {params.prefix})
        original_dir=$PWD
        genome_absolute_path=$(realpath {input.genome_fasta})
        # TideCluster 1.16.0 knobs (empty by default → TideCluster's own defaults).
        rdna_lib_arg=""
        if [ -n "{params.rdna_library}" ]; then
            rdna_lib_arg="--rdna_library $(realpath {params.rdna_library})"
        fi
        # --max_memory takes GB; below 1 GB pass nothing and let TideCluster
        # resolve its own budget (its chain matches scripts/mem_utils.py).
        max_mem_arg=""
        if [ {resources.mem_mb} -ge 1024 ]; then
            max_mem_arg="--max_memory $(( {resources.mem_mb} / 1024 ))"
        fi
        tc_extra="{params.rdna_flag} {params.overlaps_flag} $rdna_lib_arg $max_mem_arg"
        # NOTE - there is a bug in tidecluster - it does not correctly format html links; solution is
        # to run it in the directory where the output will be created
        if [ -z "{params.library}" ]; then
            cd $wd
            echo "Running TideCluster without a library"
            TideCluster.py run_all -pr $prefix -c {threads} -f $genome_absolute_path -T "-p 10 -P 39 -c 5 -e 0.25" -m 5000 $tc_extra
        else
            echo "Running TideCluster with a custom library"
            library_absolute_path=$(realpath {params.library})
            cd $wd
            TideCluster.py run_all -pr $prefix -c {threads} -f $genome_absolute_path -l $library_absolute_path -T "-p 10 -P 39 -c 5 -e 0.25" -m 5000 $tc_extra
        fi
        # Same defensive stubs as tidecluster_long (see comment there).
        cd $original_dir
        [ -f {output.gff3_clust} ]          || echo "##gff-version 3" > {output.gff3_clust}
        [ -f {output.gff3_tidehunter} ]     || echo "##gff-version 3" > {output.gff3_tidehunter}
        [ -f {output.tr_short_short} ]      || echo "##gff-version 3" > {output.tr_short_short}
        [ -f {output.dimer_library_short} ] || : > {output.dimer_library_short}
        """

rule tidecluster_reannotate:
    input:
        genome_fasta=genome_fasta_cleaned,
        dimer_library_default=F"{config['output_dir']}/TideCluster/default/TideCluster_consensus_dimer_library.fasta",
    output:
        gff3=F"{config['output_dir']}/TideCluster/default/RM_on_TideCluster_Library.gff3"
    params:
        outdir=directory(F"{config['output_dir']}/TideCluster"),
        tc_sensitivity=tc_sensitivity,
        reduce_dimer=config["reduce_tidecluster_library"],
        chunk_size=config["tidecluster_chunk_size"],
        culling_limit=config["tidecluster_reannotate_culling_limit"],
        superfamily_merge=config["tidecluster_reannotate_superfamily_merge"]
    log:
        stdout=F"{config['output_dir']}/TideCluster/default/tidecluster_reannotate.log",
        stderr=F"{config['output_dir']}/TideCluster/default/tidecluster_reannotate.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/tidecluster_reannotate.tsv"
    conda:
        "envs/tidecluster_run.yaml"
    threads: workflow.cores
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        # NOTE: tc_reannotate.py has no --max_memory upstream (unlike run_all /
        # tidehunter / tarean in TideCluster 1.20.0) — its chunked RepeatMasker
        # pool is sized by -c alone — so there is no memory budget to pass here.
        # STOPGAP (large genomes): tc_reannotate's chunked RepeatMasker opens one
        # temp-file handle per genome chunk simultaneously (~genome_bp/chunk_size,
        # ~1800 at 90 Gbp / 50 Mb) in TideCluster's split_fasta_to_chunk_files,
        # exceeding the default soft open-file limit of 1024 ("Too many open
        # files"). Same bug class as dante_ltr. Raise the soft limit to the hard
        # limit for this rule. Durable fix is upstream —
        # see docs/archive/tidecluster_large_genome_request.md.
        ulimit -n "$(ulimit -Hn)" || true
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH

        # skip entirely if the dimer library is empty
        if [ ! -s {input.dimer_library_default} ]; then
            echo "No dimer library found, skipping reannotation"
            : > {output}
            exit 0
        fi

        if [ "{params.reduce_dimer}" = "True" ]; then
            reduced={params.outdir}/default/TideCluster_consensus_dimer_library_reduced.fasta
            reduce_dimer_library.py \
                -i {input.dimer_library_default} \
                -o $reduced \
                -t {threads}
            gf_absolute_path=$(realpath $reduced)
        else
            gf_absolute_path=$(realpath {input.dimer_library_default})
        fi

        dl_absolute_path=$(realpath {input.genome_fasta})
        dl_basename=$(basename {input.genome_fasta})
        gff_absolute_path=$(realpath {output.gff3})
        # Locate TideCluster's TRC->superfamily map. Its basename varies across
        # TideCluster builds (TideCluster_trc_superfamilies.csv vs
        # TideCluster_superfamilies.csv) and it is empty when clustering found no
        # multi-TRC superfamilies — take the first non-empty candidate.
        sf_map=""
        for cand in TideCluster_trc_superfamilies.csv TideCluster_superfamilies.csv; do
            c="{params.outdir}/default/$cand"
            if [ -s "$c" ]; then sf_map=$(realpath "$c"); break; fi
        done
        cd {params.outdir}
        cp $dl_absolute_path .
        # Optional rmblastn culling for TideCluster's internal RepeatMasker.
        # The shim makes rmblastn append -culling_limit N (no RM patch); tandem
        # arrays produce many phase-redundant HSPs per locus that culling
        # collapses. 0 -> shim helper prints nothing, RMBLAST_DIR stays unset.
        shim=$(rmblast_culling_shim.py -n {params.culling_limit} -d "$PWD" || true)
        if [ -n "$shim" ]; then
            export RMBLAST_DIR="$shim"
            echo "tc_reannotate culling: -culling_limit {params.culling_limit} (RMBLAST_DIR=$shim)"
        fi
        if [ "{params.superfamily_merge}" = "True" ] && [ -z "$sf_map" ]; then
            echo "NOTE: tidecluster_reannotate_superfamily_merge is on but no populated" \
                 "TRC superfamily map was found (TideCluster found no multi-TRC" \
                 "superfamilies, or the file is named/located unexpectedly) — using the" \
                 "per-TRC filter (nothing to merge)."
        fi
        if [ "{params.superfamily_merge}" = "True" ] && [ -n "$sf_map" ]; then
            # Superfamily-aware array recovery (see config comment). Run
            # tc_reannotate in --debug mode so its raw per-TRC hits (rm.gff3)
            # survive, then re-filter grouping sibling TRCs by superfamily so a
            # real array tiled by >1 same-family TRC is not lost to the per-TRC
            # length filter. Makes the reannotation culling-independent.
            tc_log=$(mktemp)
            tc_reannotate.py -s $dl_basename -f $gf_absolute_path -o $gff_absolute_path -c {threads} --sensitivity {params.tc_sensitivity} --chunk_size {params.chunk_size} --debug > "$tc_log" 2>&1
            cat "$tc_log"
            tc_tmp=$(grep -m1 'Temp directory:' "$tc_log" | sed 's/.*Temp directory: //')
            if [ -n "$tc_tmp" ] && [ -s "$tc_tmp/rm.gff3" ]; then
                tc_reannotate_sf_filter.py \
                    --rm-gff3 "$tc_tmp/rm.gff3" \
                    --superfamily-map "$sf_map" \
                    --dimer-lib $gf_absolute_path \
                    --output $gff_absolute_path
                rm -rf "$tc_tmp"
            else
                echo "WARNING: tc_reannotate rm.gff3 not found; keeping per-TRC filtered output"
            fi
            rm -f "$tc_log"
        else
            tc_reannotate.py -s $dl_basename -f $gf_absolute_path -o $gff_absolute_path -c {threads} --sensitivity {params.tc_sensitivity} --chunk_size {params.chunk_size}
        fi
        rm $dl_basename
        [ -n "${{shim:-}}" ] && rm -rf "$shim" || true
        """

rule merge_tidecluster_default_and_short:
    input:
        gff3_default=F"{config['output_dir']}/TideCluster/default/TideCluster_clustering.gff3",
        gff3_short=F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_clustering.gff3"
    output:
        F"{config['output_dir']}/TideCluster/TideCluster_clustering_default_and_short_merged.gff3"
    log:
        stdout=F"{config['output_dir']}/TideCluster/merge_tidecluster.log",
        stderr=F"{config['output_dir']}/TideCluster/merge_tidecluster.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/merge_tidecluster_default_and_short.tsv"
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        # Concatenate both clustering GFF3s. Strip the ##gff-version (and any
        # other comment lines) from the short-monomer file so the merged
        # output doesn't carry a duplicate header — bedtools subtract refuses
        # to parse a GFF3 with a directive after a data row.
        cat {input.gff3_default} > {output}
        # replace TRC_ with TRC_S_ in the short monomer clusters to avoid name conflicts
        grep -v '^#' {input.gff3_short} | sed 's/TRC_/TRC_S_/g' >> {output} || true
        """


rule build_fallback_tir_library:
    """
    Build the optional DANTE_TIR_FALLBACK-derived library.

    Default OFF (include_dante_tir_fallback_in_library=false): emits an
    empty FASTA so concatenate_libraries can append unconditionally and
    behaviour is byte-identical to a pre-flag run.

    When ON, the helper script:

      1. Re-clusters the post-overlap fallback survivors with mmseqs2.
         Re-clustering is essential because the cluster sizes recorded
         earlier in dante_tir_fallback.py predate the primary-overlap
         filter (some members have since been dropped).
      2. Keeps cluster reps with size >= the configured Multiplicity
         floor (dante_tir_fallback_library_min_multiplicity, default
         inherits from dante_tir_min_multiplicity = 3).
      3. Canonicalises the default-library headers (LTR / DANTE_TIR
         primary / LINE / optional custom) and concatenates them into a
         BLAST DB.
      4. blastn the surviving reps against that DB
         (-evalue 1e-19, -max_target_seqs 10 — same stringency as
         filter_ltr_rt_library).
      5. Strict class-aware filter: drop any rep whose hits include a
         subject of incompatible classification (siblings like CACTA
         vs hAT count as incompatible — only same-path or
         ancestor/descendant pairs are kept).

    Audit log (fallback_library_dropped.tsv) records one row per kept
    rep and one row per (dropped rep × conflicting subject) pair.

    The fallback library is intentionally NOT fed into
    make_subclass_2_library. The DANTE_TIR_FALLBACK layer is treated
    as less reliable than the primary library, so it must not be used
    to filter the LTR library — a misclassified fallback would
    disproportionately damage the LTR side of the annotation.
    """
    input:
        fallback_fasta=F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_fallback_filtered.fasta",
        ltr_lib=F"{config['output_dir']}/DANTE_LTR/LTR_RTs_library.fasta",
        tir_primary_lib=F"{config['output_dir']}/DANTE_TIR/all_representative_elements_combined.fasta",
        line_lib=F"{config['output_dir']}/DANTE_LINE/LINE_rep_lib.fasta"
    output:
        library=F"{config['output_dir']}/DANTE_TIR/fallback_library.fasta",
        dropped_tsv=F"{config['output_dir']}/DANTE_TIR/fallback_library_dropped.tsv"
    params:
        enabled=config["include_dante_tir_fallback_in_library"],
        min_multiplicity=(
            config["dante_tir_fallback_library_min_multiplicity"]
            if config["dante_tir_fallback_library_min_multiplicity"] is not None
            else config["dante_tir_min_multiplicity"]
        ),
        custom_library=config.get("custom_library", ""),
        workdir=F"{config['output_dir']}/DANTE_TIR/fallback_library_workdir"
    log:
        stdout=F"{config['output_dir']}/DANTE_TIR/build_fallback_tir_library.log",
        stderr=F"{config['output_dir']}/DANTE_TIR/build_fallback_tir_library.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/build_fallback_tir_library.tsv"
    conda:
        "envs/dante_line.yaml"
    threads: workflow.cores
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x

        scripts_dir=$(realpath scripts)
        export PATH="$scripts_dir:$PATH"

        ENABLED_FLAG=""
        if [ "{params.enabled}" = "True" ]; then
            ENABLED_FLAG="--enabled"
        fi

        CUSTOM_FLAG=""
        if [ -n "{params.custom_library}" ] && [ -s "{params.custom_library}" ]; then
            CUSTOM_FLAG="--custom-library {params.custom_library}"
        fi

        build_fallback_tir_library.py \
            $ENABLED_FLAG \
            --fallback-fasta {input.fallback_fasta} \
            --ltr-library {input.ltr_lib} \
            --tir-primary-library {input.tir_primary_lib} \
            --line-library {input.line_lib} \
            $CUSTOM_FLAG \
            --workdir {params.workdir} \
            --min-multiplicity {params.min_multiplicity} \
            --threads {threads} \
            --output-fasta {output.library} \
            --output-dropped-tsv {output.dropped_tsv}
        """


rule make_subclass_2_library:
    params:
        library=config.get("custom_library", "")
    input:
        dante_tir_lib=F"{config['output_dir']}/DANTE_TIR/all_representative_elements_combined.fasta",
    output:
        library=F"{config['output_dir']}/Libraries/class_ii_library.fasta"
    log:
        stdout=F"{config['output_dir']}/Libraries/make_subclass_2_library.log",
        stderr=F"{config['output_dir']}/Libraries/make_subclass_2_library.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/make_subclass_2_library.tsv"
    run:
        import sys
        with open(log.stdout, "w") as _stdout, open(log.stderr, "w") as _stderr:
            sys.stdout = _stdout
            sys.stderr = _stderr
            if params.library:
                print("Custom library provided, filtering FASTA.")
                filter_fasta(params.library, output.library, "Class_II/Subclass_1")
                # add dante_tir sequences to the library
                with open(output.library, "a") as f_out:
                    with open(input.dante_tir_lib, "r") as f_in:
                        f_out.write(f_in.read())
            else:
                print("No custom library provided, using only DANTE_TIR sequences.")
                with open(output.library, "w") as f:
                    f.write("")
                    with open(input.dante_tir_lib, "r") as f_in:
                        f.write(f_in.read())



rule filter_ltr_rt_library:
    input:
        dante_library=F"{config['output_dir']}/DANTE_LTR/LTR_RTs_library.fasta",
        subclass_2_library=F"{config['output_dir']}/Libraries/class_ii_library.fasta"
    output:
        library=F"{config['output_dir']}/Libraries/LTR_RTs_library_clean.fasta"
    log:
        stdout=F"{config['output_dir']}/Libraries/filter_ltr_rt_library.log",
        stderr=F"{config['output_dir']}/Libraries/filter_ltr_rt_library.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/filter_ltr_rt_library.tsv"
    conda:
        "envs/tidecluster.yaml"
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        # Canonicalise DANTE_LTR library headers (#Class_I|LTR|Ty1/copia|Ale
        # → #Class_I/LTR/Ty1_copia/Ale) via classification_vocabulary.yaml.
        # Replaces the previous two-step sed that relied on the accident that
        # legitimate leaf underscores never sit next to a pipe.
        # Bare-name invocation via PATH — see CLAUDE.md "Calling helper
        # scripts from rules" for the dual-context contract.
        scripts_dir=$(realpath scripts)
        export PATH="$scripts_dir:$PATH"
        classification.py canonicalise-fasta-headers \
            --source DANTE_LTR {input.dante_library} {input.dante_library}.reformatted
        # if the input.subclass_2_library is empty, just copy the reformatted library
        if [ ! -s {input.subclass_2_library} ]; then
            cp {input.dante_library}.reformatted {output.library}
        else
            # if the input.subclass_2_library is not empty, filter the reformatted library using blast
            makeblastdb -in {input.subclass_2_library} -dbtype nucl
            blastn -task blastn -query {input.dante_library}.reformatted -db {input.subclass_2_library} -outfmt 6 -evalue 1e-19 -max_target_seqs 10 -out {output.library}.blast.csv
            # get the list of sequences that passed the filter
            cut -f1 {output.library}.blast.csv | sort | uniq > {output.library}.filtered_ids
            # filter the library
            seqkit grep -v -f {output.library}.filtered_ids {input.dante_library}.reformatted > {output.library}
        fi
        """


rule concatenate_libraries:
    input:
        ltr_rt_library=F"{config['output_dir']}/Libraries/LTR_RTs_library_clean.fasta",
        dante_tir_lib=F"{config['output_dir']}/DANTE_TIR/all_representative_elements_combined.fasta",
        line_rep_lib=F"{config['output_dir']}/DANTE_LINE/LINE_rep_lib.fasta",
        # Optional DANTE_TIR_FALLBACK library. Empty when the
        # include_dante_tir_fallback_in_library flag is off, in which
        # case the append below is a no-op.
        fallback_tir_lib=F"{config['output_dir']}/DANTE_TIR/fallback_library.fasta",
        # Barrier: don't build the combined library (or run the expensive
        # RepeatMasker that depends on it) until every upstream tool's
        # classifications are known to resolve against the vocabulary.
        validation_marker=F"{config['output_dir']}/.classifications_validated"
    output:
        full_names=F"{config['output_dir']}/Libraries/combined_library.fasta",
        short_names=F"{config['output_dir']}/Libraries/combined_library_short_names.fasta",
    params:
        custom_library = config.get("custom_library", ""),
        rdna_library = os.path.join(snakemake_dir, "data/rdna_library.fasta")
    log:
        stdout=F"{config['output_dir']}/Libraries/concatenate_libraries.log",
        stderr=F"{config['output_dir']}/Libraries/concatenate_libraries.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/concatenate_libraries.tsv"
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        # Start with LTR library or custom library
        if [ -z "{params.custom_library}" ]; then
            cp {input.ltr_rt_library} {output.full_names}
        else
            cat {input.ltr_rt_library} {params.custom_library} > {output.full_names}
        fi
        # Append DANTE_TIR library if not empty
        if [ -s {input.dante_tir_lib} ]; then
            cat {input.dante_tir_lib} >> {output.full_names}
        fi
        # Append DANTE_TIR_FALLBACK-derived library if not empty.
        # Empty when include_dante_tir_fallback_in_library is off.
        if [ -s {input.fallback_tir_lib} ]; then
            cat {input.fallback_tir_lib} >> {output.full_names}
        fi
        # Append LINE library if not empty
        if [ -s {input.line_rep_lib} ]; then
            cat {input.line_rep_lib} >> {output.full_names}
        fi
        # Append rDNA library
        cat {params.rdna_library} >> {output.full_names}

        # Deterministic downstream clustering: reduce_library re-clusters this
        # library per class with CAP3 / mmseqs easy-cluster, both order-sensitive
        # (same sequences in a different order -> different representatives). Sort
        # the combined library canonically by sequence here — a single choke
        # point that makes every per-class clustering input order-invariant, so
        # the reduced library (and the RepeatMasker annotation) is reproducible
        # run-to-run. Bare-name via PATH (see CLAUDE.md). canonical_sort_fasta.py
        # is out-of-core (GNU sort), so this stays cheap on large libraries.
        scripts_dir=$(realpath scripts)
        export PATH="$scripts_dir:$PATH"
        canonical_sort_fasta.py {output.full_names} {output.full_names}.sorted \
            --tmpdir "$(dirname {output.full_names})"
        mv {output.full_names}.sorted {output.full_names}

        # Create short names version
        awk '/^>/{{count++; split($0,a,"#"); print ">" count "#" a[2]; next}} {{print}}' {output.full_names} > {output.short_names}
        """

rule reduce_library:
    """
    Per-classification CAP3 / mmseqs2 dedup of the combined library.

    The Python implementation (scripts/reduce_library_size.py) replaces
    the R version (scripts/reduce_library_size.R, retained for parity
    testing only). Both produce byte-identical output; the Python
    rewrite drops actual peak memory from ~6 GB to ~2 GB on the medium
    fixture and eliminates the snakemake benchmark over-count caused
    by COW-shared pages from forked R workers (which previously
    reported 30+ GB max_rss on inputs of a few hundred KB and triggered
    GHA-runner OOM kills in release.yml's test-in-container step).

    Two-pass scheduling caps physical peak: LTR classes whose input
    FASTA size meets `reduce_library_max_parallel_bp` (default 50 MB)
    run sequentially in Phase 1; everything else runs in parallel in
    Phase 2. tests/test_reduce_library_parity.sh asserts byte-identity
    against the R reference on any input.
    """
    input:
        library=F"{config['output_dir']}/Libraries/combined_library.fasta"
    output:
        library_reduced=F"{config['output_dir']}/Libraries/combined_library_reduced.fasta"
    params:
        reduce_library_size = config["reduce_library"],
        max_parallel_bp = config["reduce_library_max_parallel_bp"],
        max_big_cap3_parallel = config["reduce_library_max_big_cap3_parallel"]
    log:
        stdout=F"{config['output_dir']}/Libraries/reduce_library.log",
        stderr=F"{config['output_dir']}/Libraries/reduce_library.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/reduce_library.tsv"
    conda: "envs/tidecluster.yaml"
    threads: workflow.cores
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        workdir=$(dirname {output.library_reduced})/workdir
        echo "Reduce library size: {params.reduce_library_size}"
        if [ "{params.reduce_library_size}" = "False" ]; then
            cp  {input.library} {output.library_reduced}
            exit 0
        fi
        reduce_library_size.py -i {input.library} -o {output.library_reduced} \
            -t {threads} -d $workdir \
            --max-parallel-bp {params.max_parallel_bp} \
            --max-big-cap3-parallel {params.max_big_cap3_parallel} \
            --max-blast-parallel {threads}
        """

rule screen_library_cross_class:
    """
    Cross-class contamination screen, between the per-class reduction and the
    containment pass.

    reduce_library_size.py clusters WITHIN a classification and
    containment_reduce_library.py only drops a fragment into a SAME-class
    container, so neither can see a consensus that is part one class and part
    another — a LINE whose inferred boundary ran out into the neighbouring
    Retand array looks like an ordinary, slightly long LINE. This blastn-screens
    every consensus against the rest of the library and TRUNCATES it to its
    longest span carrying no foreign material at all.

    Placed BEFORE reduce_library_containment on purpose: the containment pass
    can then collapse fragments that only become redundant once the foreign
    tails are gone, and it — not this rule — remains the single place that
    canonically sorts the final RepeatMasker library.

    Measured on a wheat library already carrying the dante_line boundary guards:
    156 consensi truncated (0.93% of library bp), leaving ZERO identified foreign
    base pairs in the library. With the guards this projects wheat Class_I/LINE at
    2.51% against 15.93% as shipped, inside the 1.41-3.64% core-anchored band.
    Every decision is recorded in
    Libraries/cross_class_screen.tsv. blastn-unavailable / failure copies the
    library through unchanged.
    """
    input:
        library_reduced=F"{config['output_dir']}/Libraries/combined_library_reduced.fasta"
    output:
        library_screened=F"{config['output_dir']}/Libraries/combined_library_screened.fasta",
        audit=F"{config['output_dir']}/Libraries/cross_class_screen.tsv"
    params:
        enabled=config["screen_library_cross_class"],
        min_identity=config["cross_class_min_identity"],
        min_length=config["cross_class_min_length"],
        max_shared_depth=config["cross_class_max_shared_depth"],
        ownership_margin=config["cross_class_ownership_margin"]
    log:
        stdout=F"{config['output_dir']}/Libraries/screen_library_cross_class.log",
        stderr=F"{config['output_dir']}/Libraries/screen_library_cross_class.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/screen_library_cross_class.tsv"
    conda: "envs/tidecluster.yaml"
    threads: workflow.cores
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        workdir=$(dirname {output.library_screened})/cross_class_workdir
        disabled_flag=""
        if [ "{params.enabled}" = "False" ]; then disabled_flag="--disabled"; fi
        screen_library_cross_class.py -i {input.library_reduced} \
            -o {output.library_screened} -a {output.audit} \
            -t {threads} -d $workdir \
            --min-identity {params.min_identity} \
            --min-length {params.min_length} \
            --max-shared-depth {params.max_shared_depth} \
            --ownership-margin {params.ownership_margin} \
            $disabled_flag
        """

rule library_health:
    """
    Advisory summary of the repeat library and the inferred element boundaries.

    Two defects shipped for several releases and were invisible in every output
    the pipeline produced: the DANTE_TIR library was empty on every run (so the
    RepeatMasker library held no Class_II sequences at all), and DANTE_LINE built
    16-22 kb "LINE" consensi whose flanks were a different, far more abundant
    repeat. Neither was hard to see once someone looked at the library; nothing
    in the run invited anyone to look.

    This writes the numbers that would have made both obvious -- per-class
    consensus counts and length distributions against the class bounds in
    classification_vocabulary.yaml, the share of each element builder's output
    that is inferred flank rather than anchoring domain, and how many elements
    had a flank alignment reach the --flank ceiling (a count that tracked the
    scale of the problem across 87 assemblies). Warnings go to the rule's log;
    nothing here can fail a run.
    """
    input:
        library=F"{config['output_dir']}/Libraries/combined_library_reduced_containment.fasta",
        dante_line=F"{config['output_dir']}/DANTE_LINE/DANTE_LINE.gff3",
        fallback=F"{config['output_dir']}/DANTE_TIR_FALLBACK/DANTE_TIR_FALLBACK.gff3",
        screen_audit=F"{config['output_dir']}/Libraries/cross_class_screen.tsv"
    output:
        tsv=F"{config['output_dir']}/Libraries/library_health.tsv"
    log:
        stdout=F"{config['output_dir']}/Libraries/library_health.log",
        stderr=F"{config['output_dir']}/Libraries/library_health.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/library_health.tsv"
    conda: "envs/tidecluster.yaml"
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        library_health.py -o {output.tsv} \
            --library {input.library} \
            --dante-line-gff3 {input.dante_line} \
            --fallback-gff3 {input.fallback} \
            --screen-audit {input.screen_audit}
        """

rule reduce_library_containment:
    """
    Second-round CONTAINMENT reduction of the per-class-reduced library, run
    just before RepeatMasker. The per-class CAP3/mmseqs `reduce_library` step
    leaves many short fragments fully contained in a longer element of the same
    class; this pass (`containment_reduce_library.py`, blastn greedy) drops a
    fragment when a retained, strictly longer, SAME-class sequence covers
    >= containment_min_coverage of it at >= containment_min_identity identity.
    RepeatMasker masks those fragments' genomic copies via the container, so
    masking AND classification are preserved — validated masked-bp-lossless on
    the Pisum pangenome (~-22% library bp, ~-30% RepeatMasker wall-time at
    80 / 0.90). Toggle with `reduce_library_containment` (default True);
    blastn-unavailable / failure copies the library through unchanged.
    """
    input:
        library_reduced=F"{config['output_dir']}/Libraries/combined_library_screened.fasta"
    output:
        library_containment=F"{config['output_dir']}/Libraries/combined_library_reduced_containment.fasta"
    params:
        enabled=config["reduce_library_containment"],
        min_identity=config["containment_min_identity"],
        min_coverage=config["containment_min_coverage"]
    log:
        stdout=F"{config['output_dir']}/Libraries/reduce_library_containment.log",
        stderr=F"{config['output_dir']}/Libraries/reduce_library_containment.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/reduce_library_containment.tsv"
    conda: "envs/tidecluster.yaml"
    threads: workflow.cores
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        if [ "{params.enabled}" = "False" ]; then
            cp {input.library_reduced} {output.library_containment}
        else
            workdir=$(dirname {output.library_containment})/containment_workdir
            containment_reduce_library.py -i {input.library_reduced} -o {output.library_containment} \
                -t {threads} -d $workdir \
                --min-identity {params.min_identity} --min-coverage {params.min_coverage}
        fi

        # Determinism: this file is the RepeatMasker library. RepeatMasker builds
        # its search DB with makeblastdb, which assigns OIDs by input order, and
        # tie-breaks equal-scoring HSPs by OID — so a non-canonical library order
        # makes masking non-reproducible run-to-run (measured: ~600 bp of Class_I/
        # LINE flipping between two full-genome runs of identical input, because
        # abundant near-identical LINE consensi are exactly where score ties
        # occur). The order is unstable because reduce_library_size.py streams
        # mmseqs cluster representatives in mmseqs' native, thread-scheduling-
        # dependent order (deliberately, to stay byte-identical to its R
        # reference), so combined_library_reduced.fasta — and this containment
        # output derived from it — has stable CONTENT but an unstable ORDER.
        # Canonically sort by sequence here, the single choke point immediately
        # before RepeatMasker, so masking is reproducible. (concatenate_libraries
        # does the same for the clustering input; the dimer library that feeds
        # TideCluster's internal RepeatMasker is already canonically ordered by
        # reduce_dimer_library.py, so no analogous fix is needed there.)
        canonical_sort_fasta.py {output.library_containment} {output.library_containment}.sorted \
            --tmpdir "$(dirname {output.library_containment})"
        mv {output.library_containment}.sorted {output.library_containment}
        """

rule repeatmasker:
    input:
        genome_fasta=genome_fasta_cleaned,
        library=F"{config['output_dir']}/Libraries/combined_library.fasta",
        library_short=F"{config['output_dir']}/Libraries/combined_library_short_names.fasta",
        library_reduced=F"{config['output_dir']}/Libraries/combined_library_reduced_containment.fasta"


    output:
        out=F"{config['output_dir']}/RepeatMasker/RM_on_combined_library.out",
        gff=F"{config['output_dir']}/RepeatMasker/RM_on_combined_library.gff3"
    params:
        rm_dir=directory(F"{config['output_dir']}/RepeatMasker"),
        rm_sensitivity_option=rm_sensitivity_option,
        rm_sensitivity=config["repeatmasker_sensitivity"],
        culling_limit=config["repeatmasker_culling_limit"]
    log:
        stdout=F"{config['output_dir']}/RepeatMasker/repeatmasker.log",
        stderr=F"{config['output_dir']}/RepeatMasker/repeatmasker.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/repeatmasker.tsv"
    conda:
        "envs/tidecluster.yaml"
    threads: workflow.cores
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        # STOPGAP (large genomes): RepeatMasker with -pa {threads} over >10k contigs
        # opens many temp files at once; raise the soft open-file limit to the hard
        # limit so it does not hit "Too many open files".
        ulimit -n "$(ulimit -Hn)" || true
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        # clean_rm_output.R reads CPU_COUNT for its mclapply worker count;
        # export the rule's thread budget so it no longer hard-codes 8.
        export CPU_COUNT={threads}
        library_absolute_path=$(realpath {input.library_reduced})
        genome_absolute_path=$(realpath {input.genome_fasta})
        out_absolute_path=$(realpath {output.out})
        gff_absolute_path=$(realpath {output.gff})
        cd {params.rm_dir}
        cp $library_absolute_path .
        cp $genome_absolute_path .
        lib_name=$(basename $library_absolute_path)
        gen_name=$(basename $genome_absolute_path)
        repeatmasker_wrapper.py -f $gen_name -l $lib_name -o $out_absolute_path  -s {params.rm_sensitivity} -p {threads} -d workdir --culling-limit {params.culling_limit}
        clean_rm_output.R $out_absolute_path $gff_absolute_path
        """

rule subtract_satellites_from_rm:
    input:
        rm_gff=F"{config['output_dir']}/RepeatMasker/RM_on_combined_library_plus_DANTE.gff3",
        satellite_annotation=F"{config['output_dir']}/TideCluster/TideCluster_clustering_default_and_short_merged.gff3"
    output:
        F"{config['output_dir']}/RepeatMasker/Repeat_Annotation_NoSat.gff3"
    log:
        stdout=F"{config['output_dir']}/RepeatMasker/subtract_satellites_from_rm.log",
        stderr=F"{config['output_dir']}/RepeatMasker/subtract_satellites_from_rm.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/subtract_satellites_from_rm.tsv"
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        bedtools subtract -a {input.rm_gff} -b {input.satellite_annotation} > {output}
        """

rule merge_rm_and_dante:
    input:
        rm_gff=F"{config['output_dir']}/RepeatMasker/RM_on_combined_library.gff3",
        dante_gff=F"{config['output_dir']}/DANTE/DANTE_filtered.gff3"
    output:
        gff=F"{config['output_dir']}/RepeatMasker/RM_on_combined_library_plus_DANTE.gff3"
    log:
        stdout=F"{config['output_dir']}/RepeatMasker/merge_rm_and_dante.log",
        stderr=F"{config['output_dir']}/RepeatMasker/merge_rm_and_dante.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/merge_rm_and_dante.tsv"
    conda:
        "envs/tidecluster.yaml"
        # dante_ltr is already used and it contains the necessary tools (rtracklayer and optparse)
    threads: workflow.cores
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        export CPU_COUNT={threads}
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        # Names in the merged gff3 must be consistent with the names used for
        # RepeatMasker, so the DANTE side is renamed from Final_Classification.
        # merge_repeat_annotations.R does that on the data it already imports —
        # this used to be a separate clean_DANTE_names.R pass whose only product
        # was a rewritten copy of the DANTE gff3 (10.8 GB on a 94 Gbp genome)
        # that was read straight back in and deleted.
        merge_repeat_annotations.R {input.rm_gff} {input.dante_gff} {output.gff}
        """


rule make_unified_annotation:
    """
    Produce a single, tier-prioritised repeat annotation GFF3 from all pipeline layers.
    Structure-based annotations (DANTE_LTR, DANTE_TIR, DANTE_LINE) take priority over
    similarity-based ones (RepeatMasker). See docs/archive/annotation_rules.md for full tier hierarchy.
    """
    input:
        ltr=F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3",
        ltr_tandems=F"{config['output_dir']}/DANTE_LTR/DANTE_LTR_tandems.gff3",
        tir=F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_combined.gff3",
        line=F"{config['output_dir']}/DANTE_LINE/DANTE_LINE.gff3",
        dante=F"{config['output_dir']}/DANTE/DANTE_filtered.gff3",
        tc_default=F"{config['output_dir']}/TideCluster/default/TideCluster_clustering.gff3",
        tc_short=F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_clustering.gff3",
        tc_rm=F"{config['output_dir']}/TideCluster/default/RM_on_TideCluster_Library.gff3",
        rm=F"{config['output_dir']}/RepeatMasker/RM_on_combined_library_plus_DANTE.gff3",
        th_default=F"{config['output_dir']}/TideCluster/default/TideCluster_tidehunter_short.gff3",
        th_short=F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_tidehunter_short.gff3",
        th_raw_default=F"{config['output_dir']}/TideCluster/default/TideCluster_tidehunter.gff3",
        th_raw_short=F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_tidehunter.gff3",
        fai=F"{config['output_dir']}/genome_cleaned.fasta.fai",
        validation_marker=F"{config['output_dir']}/.classifications_validated"
    output:
        gff=F"{config['output_dir']}/Repeat_Annotation_Unified.gff3",
        overlaps=F"{config['output_dir']}/Repeat_Annotation_Unified.overlaps.tsv",
        te_derived_trc=F"{config['output_dir']}/Repeat_Annotation_Unified.te_derived_trc.csv"
    params:
        # TideCluster's authoritative per-TRC rDNA calls (45S/5S). Passed as
        # params, not inputs: they are written alongside the clustering GFF3 by
        # the same upstream rule (so ordering is guaranteed via tc_default/
        # tc_short), but are absent when tidecluster_detect_rdna is false — the
        # R script guards with file.exists and falls back to the rDNA_type attr.
        tc_rdna_default=F"{config['output_dir']}/TideCluster/default/TideCluster_rdna.tsv",
        tc_rdna_short=F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_rdna.tsv",
        # Authoritative per-TRC tandem monomer period, from TideCluster's report
        # table (monomer_tarean -> monomer_kite founder fallback). Used both for
        # the TE-derived TRC table's monomer column AND the domain-rhythm gate in
        # identify_te_derived_trcs. NOT the kite `monomer_size` CSV: that is the
        # top k-mer peak and collapses to SSR sub-periods (e.g. 79 bp for a
        # 13 kb TIR-derived monomer), which would mis-tile the occupancy test.
        # trc_table.tsv lives under TideCluster_report/, which — unlike the
        # TideCluster_kite/ tree — SURVIVES `cleanup_intermediates: maximal`.
        # Params, not inputs (written by the tc_default/tc_short rules; guarded
        # with file.exists in the R script, which falls back gracefully).
        tc_trc_table_default=F"{config['output_dir']}/TideCluster/default/TideCluster_report/data/trc_table.tsv",
        tc_trc_table_short=F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_report/data/trc_table.tsv",
        rm_tc_tandem_gate=config["rm_tc_tandem_gate"],
        max_workers=config["make_unified_max_workers"]
    log:
        stdout=F"{config['output_dir']}/Repeat_Annotation_Unified.log",
        stderr=F"{config['output_dir']}/Repeat_Annotation_Unified.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/make_unified_annotation.tsv"
    conda:
        "envs/tidecluster.yaml"
    threads: workflow.cores
    # mem_mb comes from config max_memory_gb (0 = unset); an HPC profile /
    # --set-resources overrides it per rule.
    # When set it is the authoritative budget for sizing the worker pool — under
    # PBS/Slurm the scheduler's allocation, not the node's free memory, is what
    # gets enforced. When unset the script detects the budget itself (scheduler
    # environment, then the tightest of the cgroup limit — walking up to the job
    # scope — and /proc/meminfo).
    resources:
        mem_mb=MAX_MEMORY_MB
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH

        # 0 => let the script auto-detect (cgroup-aware); otherwise pass the
        # scheduler allocation through as the budget.
        mem_budget_gb=0
        if [ {resources.mem_mb} -gt 0 ]; then
            mem_budget_gb=$(( {resources.mem_mb} / 1024 ))
        fi

        make_unified_annotation.R \
            --max_workers   {params.max_workers} \
            --mem_budget_gb "$mem_budget_gb" \
            --ltr      {input.ltr} \
            --ltr_tandems {input.ltr_tandems} \
            --tir      {input.tir} \
            --line     {input.line} \
            --dante    {input.dante} \
            --tc_default {input.tc_default} \
            --tc_short   {input.tc_short} \
            --tc_rm      {input.tc_rm} \
            --tc_rdna_default {params.tc_rdna_default} \
            --tc_rdna_short   {params.tc_rdna_short} \
            --tc_trc_table_default {params.tc_trc_table_default} \
            --tc_trc_table_short   {params.tc_trc_table_short} \
            --rm       {input.rm} \
            --th_default {input.th_default} \
            --th_short   {input.th_short} \
            --th_raw_default {input.th_raw_default} \
            --th_raw_short   {input.th_raw_short} \
            --rm_tc_tandem_gate {params.rm_tc_tandem_gate} \
            --fai      {input.fai} \
            --output   {output.gff} \
            --threads  {threads}
        """


rule make_track_for_masking:
    input:
        rm=F"{config['output_dir']}/RepeatMasker/Repeat_Annotation_NoSat.gff3",
        tr_main=F"{config['output_dir']}/TideCluster/TideCluster_clustering_default_and_short_merged.gff3",
        tr_default_short=F"{config['output_dir']}/TideCluster/default/TideCluster_tidehunter_short.gff3",
        tr_short_short=F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_tidehunter_short.gff3",
        tr_rm=F"{config['output_dir']}/TideCluster/default/RM_on_TideCluster_Library.gff3",
        dante_ltr=F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3"
    output:
        F"{config['output_dir']}/all_repeats_for_masking.bed"
    log:
        stdout=F"{config['output_dir']}/logs/make_track_for_masking.log",
        stderr=F"{config['output_dir']}/logs/make_track_for_masking.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/make_track_for_masking.tsv"
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        export LC_COLLATE=C
        cat {input.rm} {input.tr_main} {input.tr_default_short} {input.tr_short_short} {input.tr_rm} {input.dante_ltr} | sort -k1,1 -k4,4n > {output}.tmp.gff3
        bedtools merge -i {output}.tmp.gff3 > {output}
        rm {output}.tmp.gff3
        """

rule make_track_for_Ns:
    input:
        genome_fasta=genome_fasta_cleaned
    output:
        F"{config['output_dir']}/gaps_10plus.bed"
    log:
        stdout=F"{config['output_dir']}/logs/make_track_for_Ns.log",
        stderr=F"{config['output_dir']}/logs/make_track_for_Ns.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/make_track_for_Ns.tsv"
    conda:
        "envs/seqtk.yaml"
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        seqtk cutN -n 10 -g {input.genome_fasta} > {output}
        """

rule make_summary_statistics_and_split_by_class:
    """
    Compute summary_statistics.csv and the per-class split GFFs from
    Repeat_Annotation_Unified.gff3.

    Switched from Repeat_Annotation_NoSat.gff3 (RepeatMasker-only) to
    the Unified annotation so DANTE-direct calls — DANTE_LTR /
    DANTE_TIR (incl. fallback partials) / DANTE_LINE structural hits
    that have no matching RepeatMasker call — count toward genome-wide
    totals. Unified is non-overlapping by construction (tier-priority
    resolves multi-source overlaps to one label per bp), so the bp
    totals do not need additional dedup here.

    No separate satellite subtraction: TideCluster wins tier-priority
    in Unified wherever it has a structural call, so the LTR-RT-library
    misannotation that originally motivated subtract_satellites_from_rm
    cannot bleed into Unified. The R splitter collapses TC TRC_*,
    RM Satellite/*, and TideHunter-named records into a single
    `Tandem_repeats` aggregation row.

    The output directory keeps its historical name
    `Repeat_Annotation_NoSat_split_by_class_gff3/` so downstream
    consumers and existing user scripts that hard-code the path
    continue to work — the contained GFFs are now sourced from
    Unified, but the directory layout is unchanged.
    """
    input:
        unified=F"{config['output_dir']}/Repeat_Annotation_Unified.gff3",
        genome_fasta=genome_fasta_cleaned
    output:
        csv=F"{config['output_dir']}/summary_statistics.csv",
        dir=directory(F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3"),
        mobile_elements=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3/Mobile_elements.gff3",
        simple_repeats=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3/Simple_repeats.gff3",
        low_complexity=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3/Low_complexity.gff3",
        rdna=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3/rDNA.gff3",
        all_copia=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3/All_Ty1_Copia.gff3",
        all_gypsy=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3/All_Ty3_Gypsy.gff3",
    log:
        stdout=F"{config['output_dir']}/logs/make_summary_statistics_and_split_by_class.log",
        stderr=F"{config['output_dir']}/logs/make_summary_statistics_and_split_by_class.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/make_summary_statistics_and_split_by_class.tsv"
    conda:
        "envs/tidecluster.yaml"
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        calculate_statistics_and_make_groups.R -r {input.unified} -o {output.csv} -g {input.genome_fasta} -d {output.dir} \
        -M {output.mobile_elements} -S {output.simple_repeats} -L {output.low_complexity} -R {output.rdna} -C {output.all_copia} -G {output.all_gypsy}
        """

rule make_bigwig_density:
    input:
        cvs=F"{config['output_dir']}/summary_statistics.csv",  # this file is available if gffs were created
        genome_seqlengths=F"{config['output_dir']}/genome_seqlengths.rds"
    output:
        checkpoint=F"{config['output_dir']}/Repeat_density_by_class_bigwig/.done"
    params:
        bwdir=F"{config['output_dir']}/Repeat_density_by_class_bigwig",
        gffdir=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3",
        genome_fasta=genome_fasta_cleaned,
        max_workers=config["bigwig_max_workers"]
    log:
        stdout=F"{config['output_dir']}/logs/make_bigwig_density.log",
        stderr=F"{config['output_dir']}/logs/make_bigwig_density.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/make_bigwig_density.tsv"
    threads: workflow.cores
    conda:
        "envs/tidecluster.yaml"
    # See make_unified_annotation: unset (0) lets the script auto-detect the
    # budget (scheduler environment, then the cgroup limit walking up to the job
    # scope, else /proc/meminfo). config max_memory_gb, an HPC profile or
    # --set-resources supply the real allocation instead.
    resources:
        mem_mb=MAX_MEMORY_MB
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        mkdir -p {params.bwdir}
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        ls_absolute_path=$(realpath {input.genome_seqlengths})
        mem_budget_gb=0
        if [ {resources.mem_mb} -gt 0 ]; then
            mem_budget_gb=$(( {resources.mem_mb} / 1024 ))
        fi
        calculate_density_batch.R -d {params.gffdir} -o {params.bwdir} -g $ls_absolute_path -t {threads} \
            --max_workers {params.max_workers} --mem_budget_gb "$mem_budget_gb"
        touch {output.checkpoint}
        """

rule add_top_level_outputs:
    input:
        dante=F"{config['output_dir']}/DANTE/DANTE_filtered.gff3",
        dante_ltr=F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3",
        dante_ltr_tandems=F"{config['output_dir']}/DANTE_LTR/DANTE_LTR_tandems.gff3",
        dante_tir=F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_combined.gff3",
        sat_tc=F"{config['output_dir']}/TideCluster/default/TideCluster_clustering.gff3",
        sat_rm=F"{config['output_dir']}/TideCluster/default/RM_on_TideCluster_Library.gff3",
        simple_repeats=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3/Simple_repeats.gff3",
        low_complexity=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3/Low_complexity.gff3",
        rdna=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3/rDNA.gff3",
        mobile_elements=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3/Mobile_elements.gff3",
        all_copia=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3/All_Ty1_Copia.gff3",
        all_gypsy=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3/All_Ty3_Gypsy.gff3"

    output:
        dante=F"{config['output_dir']}/DANTE_filtered.gff3",
        dante_ltr=F"{config['output_dir']}/DANTE_LTR.gff3",
        dante_ltr_tandems=F"{config['output_dir']}/DANTE_LTR_tandems.gff3",
        dante_tir=F"{config['output_dir']}/DANTE_TIR.gff3",
        sat_tc=F"{config['output_dir']}/Tandem_repeats_TideCluster.gff3",
        sat_rm=F"{config['output_dir']}/Tandem_repeats_RepeatMasker.gff3",
        simple_repeats=F"{config['output_dir']}/Simple_repeats_RepeatMasker.gff3",
        low_complexity=F"{config['output_dir']}/Low_complexity_RepeatMasker.gff3",
        rdna=F"{config['output_dir']}/rDNA_RepeatMasker.gff3",
        mobile_elements=F"{config['output_dir']}/Mobile_elements_RepeatMasker.gff3",
        all_copia=F"{config['output_dir']}/All_Ty1_Copia_RepeatMasker.gff3",
        all_gypsy=F"{config['output_dir']}/All_Ty3_Gypsy_RepeatMasker.gff3"
    params:
        sat_tc_annot_in=F"{config['output_dir']}/TideCluster/default/TideCluster_annotation.gff3",
        sat_tc_annot_out=F"{config['output_dir']}/Tandem_repeats_TideCluster_annotated.gff3"
    benchmark:
        F"{config['output_dir']}/benchmarks/add_top_level_outputs.tsv"
    shell:
        """
        # make symbolic links to all the outputs
        ln -fs -r {input.dante} {output.dante}
        ln -fs -r {input.dante_ltr} {output.dante_ltr}
        ln -fs -r {input.dante_ltr_tandems} {output.dante_ltr_tandems}
        ln -fs -r {input.dante_tir} {output.dante_tir}
        ln -fs -r {input.sat_tc} {output.sat_tc}
        ln -fs -r {input.sat_rm} {output.sat_rm}
        ln -fs -r {input.simple_repeats} {output.simple_repeats}
        ln -fs -r {input.low_complexity} {output.low_complexity}
        ln -fs -r {input.rdna} {output.rdna}
        ln -fs -r {input.mobile_elements} {output.mobile_elements}
        ln -fs -r {input.all_copia} {output.all_copia}
        ln -fs -r {input.all_gypsy} {output.all_gypsy}
        if [ -f {params.sat_tc_annot_in} ]; then
             ln -fs -r {params.sat_tc_annot_in} {params.sat_tc_annot_out}
        fi
        """

rule calculate_bigwig_density:
    """
    Genome-wide total density (from the Unified annotation — all repeats,
    incl. tandems) and the structural TideCluster total density (a subset:
    TideCluster clustering only). Per-class and per-family tracks are built
    by make_bigwig_density / the per-family rules.
    """
    input:
        unified=F"{config['output_dir']}/Repeat_Annotation_Unified.gff3",
        tc=F"{config['output_dir']}/TideCluster/default/TideCluster_clustering.gff3",
        genome_seqlengths=F"{config['output_dir']}/genome_seqlengths.rds"
    output:
        total_10k=F"{config['output_dir']}/Repeat_density/Repeat_density_total_10k.bw",
        total_100k=F"{config['output_dir']}/Repeat_density/Repeat_density_total_100k.bw",
        tc_10k=F"{config['output_dir']}/Tandem_repeats_TideCluster_10k.bw",
        tc_100k=F"{config['output_dir']}/Tandem_repeats_TideCluster_100k.bw"
    log:
        stdout=F"{config['output_dir']}/logs/calculate_bigwig_density.log",
        stderr=F"{config['output_dir']}/logs/calculate_bigwig_density.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/calculate_bigwig_density.tsv"
    conda:
        "envs/tidecluster.yaml"
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        mkdir -p $(dirname {output.total_10k})
        calculate_density.R -b {input.unified} -o {output.total_10k}  -f gff3 --window 10000  -g {input.genome_seqlengths}
        calculate_density.R -b {input.unified} -o {output.total_100k} -f gff3 --window 100000 -g {input.genome_seqlengths}
        calculate_density.R -b {input.tc}      -o {output.tc_10k}     -f gff3 --window 10000  -g {input.genome_seqlengths}
        calculate_density.R -b {input.tc}      -o {output.tc_100k}    -f gff3 --window 100000 -g {input.genome_seqlengths}
        """


rule make_unified_tandem_per_family_bigwig:
    """
    Per-family tandem density BigWigs from the Unified annotation.

    Splits Repeat_Annotation_Unified.gff3 to the tandem families
    (Name=TRC_<n>) and emits one density BigWig per family. Because Unified
    is tier-resolved, each family's track is the structural (TideCluster) ∪
    similarity (RM-on-TideCluster) territory for that family, minus any
    region won by a higher-priority mobile element. This is the authoritative
    per-family view; the structural-only counterpart is produced by
    make_tidecluster_tandem_per_family_bigwig. Empty input is a no-op.
    """
    input:
        unified=F"{config['output_dir']}/Repeat_Annotation_Unified.gff3",
        genome_seqlengths=F"{config['output_dir']}/genome_seqlengths.rds"
    output:
        done=F"{config['output_dir']}/Tandem_repeats_unified_split_by_family_bigwig/.done"
    params:
        splitdir=F"{config['output_dir']}/Tandem_repeats_unified_split_by_family_gff3",
        bwdir=F"{config['output_dir']}/Tandem_repeats_unified_split_by_family_bigwig",
        max_workers=config["bigwig_max_workers"]
    log:
        stdout=F"{config['output_dir']}/logs/make_unified_tandem_per_family_bigwig.log",
        stderr=F"{config['output_dir']}/logs/make_unified_tandem_per_family_bigwig.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/make_unified_tandem_per_family_bigwig.tsv"
    threads: workflow.cores
    conda:
        "envs/tidecluster.yaml"
    resources:
        mem_mb=MAX_MEMORY_MB
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        ls_absolute_path=$(realpath {input.genome_seqlengths})
        mkdir -p {params.splitdir} {params.bwdir}
        mem_budget_gb=0
        if [ {resources.mem_mb} -gt 0 ]; then
            mem_budget_gb=$(( {resources.mem_mb} / 1024 ))
        fi
        split_gff_by_name.R -i {input.unified} -o {params.splitdir} --name-prefix TRC_
        if ls {params.splitdir}/*.gff3 >/dev/null 2>&1; then
            calculate_density_batch.R -d {params.splitdir} -o {params.bwdir} -g $ls_absolute_path -t {threads} \
                --max_workers {params.max_workers} --mem_budget_gb "$mem_budget_gb"
        else
            echo "No tandem families to split — nothing to do"
        fi
        touch {output.done}
        """


rule make_tidecluster_tandem_per_family_bigwig:
    """
    Per-family STRUCTURAL tandem density BigWigs (TideCluster clustering
    only — a subset of the Unified per-family union). Reads the per-cluster
    split GFF3s that tidecluster_long emits (default run only) and writes one
    density BigWig per family. Empty / no-cluster input is a no-op.
    """
    input:
        tc_clust=F"{config['output_dir']}/TideCluster/default/TideCluster_clustering.gff3",
        split_files=F"{config['output_dir']}/TideCluster/default/TideCluster_clustering_split_files",
        genome_seqlengths=F"{config['output_dir']}/genome_seqlengths.rds"
    output:
        done=F"{config['output_dir']}/Tandem_repeats_TideCluster_split_by_family_bigwig/.done"
    params:
        bwdir=F"{config['output_dir']}/Tandem_repeats_TideCluster_split_by_family_bigwig",
        max_workers=config["bigwig_max_workers"]
    log:
        stdout=F"{config['output_dir']}/logs/make_tidecluster_tandem_per_family_bigwig.log",
        stderr=F"{config['output_dir']}/logs/make_tidecluster_tandem_per_family_bigwig.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/make_tidecluster_tandem_per_family_bigwig.tsv"
    threads: workflow.cores
    conda:
        "envs/tidecluster.yaml"
    resources:
        mem_mb=MAX_MEMORY_MB
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        ls_absolute_path=$(realpath {input.genome_seqlengths})
        mkdir -p {params.bwdir}
        mem_budget_gb=0
        if [ {resources.mem_mb} -gt 0 ]; then
            mem_budget_gb=$(( {resources.mem_mb} / 1024 ))
        fi
        if ls {input.split_files}/*.gff3 >/dev/null 2>&1; then
            calculate_density_batch.R -d {input.split_files} -o {params.bwdir} -g $ls_absolute_path -t {threads} \
                --max_workers {params.max_workers} --mem_budget_gb "$mem_budget_gb"
        else
            echo "No TideCluster per-cluster split files — nothing to do"
        fi
        touch {output.done}
        """


rule add_html_outputs:
    input:
        tc_index=F"{config['output_dir']}/TideCluster/default/TideCluster_index.html",
        dante_ltr_index=F"{config['output_dir']}/DANTE_LTR/DANTE_LTR_summary.html"
    output:
        tc_index=F"{config['output_dir']}/TideCluster_report.html",
        dante_ltr_index=F"{config['output_dir']}/DANTE_LTR_report.html"
    benchmark:
        F"{config['output_dir']}/benchmarks/add_html_outputs.tsv"
    shell:
        """
        # TideCluster's v2 report (TideCluster_index.html) is NOT self-contained:
        # it loads styling/JS from a sibling TideCluster_report/ directory by
        # relative path. A plain symlink at the output root makes the browser
        # resolve those asset URLs against the root (where TideCluster_report/
        # does not exist), so the page renders unstyled. Write a tiny redirect
        # page instead — the browser navigates to the real path, where the
        # report's relative assets resolve correctly.
        cat > {output.tc_index} <<'HTML'
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta http-equiv="refresh" content="0; url=TideCluster/default/TideCluster_index.html">
<title>TideCluster report</title>
</head>
<body>
<script>location.replace("TideCluster/default/TideCluster_index.html");</script>
<p>Redirecting to the <a href="TideCluster/default/TideCluster_index.html">TideCluster report</a>…</p>
</body>
</html>
HTML
        # DANTE_LTR_summary.html is self-contained (no sibling assets), so a
        # relative symlink renders fine.
        ln -s -r {input.dante_ltr_index} {output.dante_ltr_index}
        """


rule calculate_seqlengths:
    input:
        genome_fasta=genome_fasta_cleaned
    output:
        F"{config['output_dir']}/genome_seqlengths.rds"
    log:
        stdout=F"{config['output_dir']}/logs/calculate_seqlengths.log",
        stderr=F"{config['output_dir']}/logs/calculate_seqlengths.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/calculate_seqlengths.tsv"
    conda:
        "envs/tidecluster.yaml"
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        calculate_seqlengths.R  {input.genome_fasta} {output}
        """


rule make_summary_plots:
    input:
        SL = F"{config['output_dir']}/genome_seqlengths.rds",
        bw1 = F"{config['output_dir']}/Repeat_density/Repeat_density_total_10k.bw",
        bw2_info = F"{config['output_dir']}/Repeat_density_by_class_bigwig/.done",
        bw_tc = F"{config['output_dir']}/Tandem_repeats_TideCluster_split_by_family_bigwig/.done"
        # these inputs gate ordering; the plot script reads the BigWig dirs directly
    output:
        F"{config['output_dir']}/summary_plots.pdf"
    params:
        output_dir = config["output_dir"]
    log:
        stdout=F"{config['output_dir']}/logs/make_summary_plots.log",
        stderr=F"{config['output_dir']}/logs/make_summary_plots.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/make_summary_plots.tsv"
    conda:
        "envs/tidecluster.yaml"
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -uo pipefail   # NOT -e: failures are handled explicitly below
        set -x
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        # Plotting must never fail the workflow. The R script already falls back
        # to a placeholder PDF on internal errors; this is the outer net for the
        # cases R cannot catch (e.g. the process is OOM-killed): if the script
        # exits non-zero or leaves no/empty output, write a placeholder so the
        # rule still produces a valid, non-empty PDF. (No bare `touch`, which
        # would create a 0-byte "PDF".)
        rc=0
        make_summary_plots.R {params.output_dir} {output} || rc=$?
        if [ "$rc" -ne 0 ] || [ ! -s {output} ]; then
            echo "make_summary_plots.R failed (rc=$rc) or produced empty output; writing placeholder PDF"
            make_placeholder_pdf.R {output} "Summary plots not rendered"
        fi
        """


rule make_benchmark_report:
    input:
        expand(F"{config['output_dir']}/benchmarks/{{rule_name}}.tsv",
               rule_name=BENCHMARKED_RULES)
    output:
        F"{config['output_dir']}/benchmark_report.html"
    params:
        benchmark_dir=F"{config['output_dir']}/benchmarks"
    log:
        stdout=F"{config['output_dir']}/logs/make_benchmark_report.log",
        stderr=F"{config['output_dir']}/logs/make_benchmark_report.err"
    conda:
        "envs/tidecluster.yaml"
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        make_benchmark_report.R {params.benchmark_dir} {output}
        """

rule make_repeat_report:
    input:
        fai   = genome_fasta_cleaned,
        stats = F"{config['output_dir']}/summary_statistics.csv",
        bw_rm = F"{config['output_dir']}/Repeat_density_by_class_bigwig/.done",
        bw_tc = F"{config['output_dir']}/Tandem_repeats_TideCluster_split_by_family_bigwig/.done",
        ltr   = F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3",
        tir   = F"{config['output_dir']}/DANTE_TIR/TIR_classification_summary.txt",
        line  = F"{config['output_dir']}/DANTE_LINE/DANTE_LINE.gff3",
        te_derived_trc = F"{config['output_dir']}/Repeat_Annotation_Unified.te_derived_trc.csv"
    output:
        F"{config['output_dir']}/repeat_annotation_report.html"
    params:
        output_dir       = config['output_dir'],
        bin_width        = config.get("report_bin_width", 100000),
        min_len_chart    = config.get("report_min_len_chart", 500000),
        min_len_tracks   = config.get("report_min_len_tracks", 1000000),
        max_tracks       = config.get("report_max_tracks", 50),
        top_sat_clusters = config.get("report_top_sat_clusters", 10)
    log:
        stdout = F"{config['output_dir']}/logs/make_repeat_report.log",
        stderr = F"{config['output_dir']}/logs/make_repeat_report.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/make_repeat_report.tsv"
    conda:
        "envs/tidecluster.yaml"
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH

        # The HTML report is a best-effort, non-essential output: it must never
        # fail the whole pipeline. Per-section failures are already handled
        # inside the script (safe_build -> "not generated" placeholders). This
        # guard is the last resort for a TOTAL failure (e.g. an unforeseen
        # large-genome edge case or OOM before any HTML is written): log it and
        # emit a minimal placeholder so the declared output exists and the
        # `all` target still completes. The annotation outputs are unaffected.
        if ! make_repeat_report.R \
                --output_dir {params.output_dir} \
                --bin_width {params.bin_width} \
                --min_len_chart {params.min_len_chart} \
                --min_len_tracks {params.min_len_tracks} \
                --max_tracks {params.max_tracks} \
                --top_sat_clusters {params.top_sat_clusters}; then
            echo "WARNING: make_repeat_report.R failed; writing placeholder report so the pipeline continues" >&2
        fi

        if [ ! -s {output} ]; then
            cat > {output} <<'HTML'
<!DOCTYPE html><html><head><meta charset="utf-8">
<title>Repeat Annotation Report — not generated</title></head>
<body style="font-family:sans-serif;max-width:820px;margin:40px auto;color:#000">
<h1>Repeat annotation report was not generated</h1>
<p>The HTML report could not be built for this run (for example, a very large
genome hit a resource or edge-case limit during report rendering). This does
<b>not</b> affect the annotation results: the GFF3 files, repeat libraries,
BigWig density tracks and <code>summary_statistics.csv</code> are complete and
valid.</p>
<p>See <code>logs/make_repeat_report.err</code> for the cause.</p>
</body></html>
HTML
        fi
        """
