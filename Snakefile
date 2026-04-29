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

# DANTE_TIR_FALLBACK stringency knobs. Both default to 3; both must be
# positive integers. See dante_tir_fallback.py for semantics.
for _key in ("dante_tir_fallback_min_alignments", "dante_tir_fallback_min_cluster_size"):
    if _key not in config:
        config[_key] = 3
    if not isinstance(config[_key], int) or config[_key] < 1:
        raise ValueError(f"Invalid value for {_key}: must be a positive integer.")

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
        F"{config['output_dir']}/TideCluster/default/.bigwig_done",
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
        F"{config['output_dir']}/RepeatMasker/Repeat_Annotation_NoSat_10k.bw",
        F"{config['output_dir']}/RepeatMasker/Repeat_Annotation_NoSat_100k.bw",
        F"{config['output_dir']}/TideCluster/default/TideCluster_clustering_10k.bw",
        F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_bigwig/.done",
        F"{config['output_dir']}/summary_plots.pdf",
        F"{config['output_dir']}/benchmark_report.html",
        F"{config['output_dir']}/repeat_annotation_report.html",
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
        "envs/tidecluster.yaml"
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
        output_dir=F"{config['output_dir']}/DANTE_TIR"
    log:
        stdout=F"{config['output_dir']}/DANTE_TIR/dante_tir.log",
        stderr=F"{config['output_dir']}/DANTE_TIR/dante_tir.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/dante_tir.tsv"
    conda:
        "envs/dante_tir.yaml"
    threads: workflow.cores
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        # Run dante_tir.py and check exit status
        if dante_tir.py -g {input.gff} -f {input.fasta} -o {params.output_dir} -c {threads}; then
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
        mask_gff=F"{config['output_dir']}/TideCluster/default/TideCluster_tidehunter.gff3"
    output:
        gff=F"{config['output_dir']}/DANTE_TIR_FALLBACK/DANTE_TIR_FALLBACK.gff3",
        rep_lib=F"{config['output_dir']}/DANTE_TIR_FALLBACK/TIR_fallback_rep_lib.fasta",
        extended_fasta=F"{config['output_dir']}/DANTE_TIR_FALLBACK/TIR_fallback_extended.fasta"
    params:
        output_dir=F"{config['output_dir']}/DANTE_TIR_FALLBACK",
        min_alignments=config["dante_tir_fallback_min_alignments"],
        min_cluster_size=config["dante_tir_fallback_min_cluster_size"]
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
            --mask-gff3 {input.mask_gff}

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

        # Run mmseqs2 clustering on the primary-only set
        mmseqs easy-cluster \
            "$FILTERED_INPUT" \
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
        "envs/tidecluster.yaml"
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
    """
    input:
        dante_filtered=F"{config['output_dir']}/DANTE/DANTE_filtered.gff3",
        dante_ltr=F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3",
        dante_tir_final=F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_final.gff3",
        dante_tir_combined=F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_combined.gff3",
        dante_line=F"{config['output_dir']}/DANTE_LINE/DANTE_LINE.gff3"
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
        touch {output.marker}
        """

rule dante_line:
    input:
        gff=F"{config['output_dir']}/DANTE/DANTE_filtered.gff3",
        gff3_tidehunter=F"{config['output_dir']}/TideCluster/default/TideCluster_tidehunter.gff3",
        genome=genome_fasta_cleaned
    output:
        line_rep_lib=F"{config['output_dir']}/DANTE_LINE/LINE_rep_lib.fasta",
        gff_out=F"{config['output_dir']}/DANTE_LINE/DANTE_LINE.gff3",
        line_regions=F"{config['output_dir']}/DANTE_LINE/LINE_regions.fasta",
        line_regions_extended=F"{config['output_dir']}/DANTE_LINE/LINE_regions_extended.fasta"
    params:
        output_dir=F"{config['output_dir']}/DANTE_LINE"
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

        # dante_line.py exits non-zero when it can't find any valid domain
        # patterns (happens on small CI fixtures with few LINE hits). Catch
        # the failure and ensure all four declared outputs exist as empty
        # files so downstream RepeatMasker / library construction continues.
        mkdir -p {params.output_dir}
        if ! dante_line.py -g {input.genome} -a {input.gff} \
                -o {params.output_dir} -t {threads} \
                --mask-gff3 {input.gff3_tidehunter}; then
            echo "dante_line.py failed (likely too few LINE features); creating empty outputs"
        fi
        [ -f {output.line_rep_lib} ]          || : > {output.line_rep_lib}
        [ -f {output.gff_out} ]               || echo "##gff-version 3" > {output.gff_out}
        [ -f {output.line_regions} ]          || : > {output.line_regions}
        [ -f {output.line_regions_extended} ] || : > {output.line_regions_extended}
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
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        dante_ltr -o {params.prefix} -s {input.fasta} -g {input.gff} -c {threads} -M 1 -S 50000000
        # if exit status is 0 and gff3 file was created but html is missing, create an empty file
        echo "DANTE LTR-RTs finished"
        if [ -f {output.gff} ]; then
            if [ ! -f {output.html} ]; then
                echo "Creating an empty html file"
                echo "No complete LTR-RTs found" > {output.html}
            fi
        fi
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
        bigwig_done=F"{config['output_dir']}/TideCluster/default/.bigwig_done"
    params:
        prefix = lambda wildcards, output: output.gff3_clust.replace("_clustering.gff3", ""),
        library = config.get("tandem_repeat_library", "")
    log:
        stdout=F"{config['output_dir']}/TideCluster/default/tidecluster_long.log",
        stderr=F"{config['output_dir']}/TideCluster/default/tidecluster_long.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/tidecluster_long.tsv"
    conda:
        "envs/tidecluster.yaml"
    threads: workflow.cores

    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        wd=$(dirname {output.gff3_clust})
        prefix=$(basename {params.prefix})
        original_dir=$PWD
        genome_absolute_path=$(realpath {input.genome_fasta})
        genome_seqlengths=$(realpath {input.genome_seqlengths})
        # NOTE - there is a bug in tidecluster - it does not correctly format html links; solution is
        # to run it in the directory where the output will be created
        echo "Library: {input.library}"
        if [ -z "{params.library}" ]; then
            cd $wd
            echo "Running TideCluster without a library"
            TideCluster.py run_all -pr $prefix -c {threads}  -f $genome_absolute_path --long
        else
            library_absolute_path=$(realpath {params.library})
            echo "Running TideCluster with a custom library"
            cd $wd
            TideCluster.py run_all -pr $prefix -c {threads} -f $genome_absolute_path -l $library_absolute_path --long
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
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        cd $wd
        if [ -d TideCluster_clustering_split_files ]; then
            mkdir -p TideCluster_clustering_split_files_bigwig
            calculate_density_batch.R -d TideCluster_clustering_split_files -o TideCluster_clustering_split_files_bigwig -g $genome_seqlengths
            touch .bigwig_done
        else
            echo "No split files found"
            touch .bigwig_done
        fi
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
        library = config.get("tandem_repeat_library", "")
    log:
        stdout=F"{config['output_dir']}/TideCluster/short_monomer/tidecluster_short.log",
        stderr=F"{config['output_dir']}/TideCluster/short_monomer/tidecluster_short.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/tidecluster_short.tsv"
    conda:
        "envs/tidecluster.yaml"
    threads: workflow.cores

    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        wd=$(dirname {output.gff3_clust})
        prefix=$(basename {params.prefix})
        original_dir=$PWD
        genome_absolute_path=$(realpath {input.genome_fasta})
        # NOTE - there is a bug in tidecluster - it does not correctly format html links; solution is
        # to run it in the directory where the output will be created
        if [ -z "{params.library}" ]; then
            cd $wd
            echo "Running TideCluster without a library"
            TideCluster.py run_all -pr $prefix -c {threads} -f $genome_absolute_path -T "-p 10 -P 39 -c 5 -e 0.25" -m 5000
        else
            echo "Running TideCluster with a custom library"
            library_absolute_path=$(realpath {params.library})
            cd $wd
            TideCluster.py run_all -pr $prefix -c {threads} -f $genome_absolute_path -l $library_absolute_path -T "-p 10 -P 39 -c 5 -e 0.25" -m 5000
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
        reduce_library=config["reduce_library"]
    log:
        stdout=F"{config['output_dir']}/TideCluster/default/tidecluster_reannotate.log",
        stderr=F"{config['output_dir']}/TideCluster/default/tidecluster_reannotate.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/tidecluster_reannotate.tsv"
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

        # skip entirely if the dimer library is empty
        if [ ! -s {input.dimer_library_default} ]; then
            echo "No dimer library found, skipping reannotation"
            : > {output}
            exit 0
        fi

        if [ "{params.reduce_library}" = "True" ]; then
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
        cd {params.outdir}
        cp $dl_absolute_path .
        tc_reannotate.py -s $dl_basename -f $gf_absolute_path -o $gff_absolute_path -c {threads} --sensitivity {params.tc_sensitivity}
        rm $dl_basename
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
        # Create short names version
        awk '/^>/{{count++; split($0,a,"#"); print ">" count "#" a[2]; next}} {{print}}' {output.full_names} > {output.short_names}
        """

rule reduce_library:
    input:
        library=F"{config['output_dir']}/Libraries/combined_library.fasta"
    output:
        library_reduced=F"{config['output_dir']}/Libraries/combined_library_reduced.fasta"
    params:
        reduce_library_size = config["reduce_library"]
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
        reduce_library_size.R -i {input.library} -o {output.library_reduced} -t {threads} -d $workdir
        """

rule repeatmasker:
    input:
        genome_fasta=genome_fasta_cleaned,
        library=F"{config['output_dir']}/Libraries/combined_library.fasta",
        library_short=F"{config['output_dir']}/Libraries/combined_library_short_names.fasta",
        library_reduced=F"{config['output_dir']}/Libraries/combined_library_reduced.fasta"


    output:
        out=F"{config['output_dir']}/RepeatMasker/RM_on_combined_library.out",
        gff=F"{config['output_dir']}/RepeatMasker/RM_on_combined_library.gff3"
    params:
        rm_dir=directory(F"{config['output_dir']}/RepeatMasker"),
        rm_sensitivity_option=rm_sensitivity_option,
        rm_sensitivity=config["repeatmasker_sensitivity"]
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
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        library_absolute_path=$(realpath {input.library_reduced})
        genome_absolute_path=$(realpath {input.genome_fasta})
        out_absolute_path=$(realpath {output.out})
        gff_absolute_path=$(realpath {output.gff})
        cd {params.rm_dir}
        cp $library_absolute_path .
        cp $genome_absolute_path .
        lib_name=$(basename $library_absolute_path)
        gen_name=$(basename $genome_absolute_path)
        repeatmasker_wrapper.py -f $gen_name -l $lib_name -o $out_absolute_path  -s {params.rm_sensitivity} -p {threads} -d workdir
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
        # names in gff3 must be consistent with names used for RepeatMasker
        clean_DANTE_names.R {input.dante_gff}  {input.dante_gff}.tmp.gff3
        merge_repeat_annotations.R {input.rm_gff} {input.dante_gff}.tmp.gff3 {output.gff}
        """


rule make_unified_annotation:
    """
    Produce a single, tier-prioritised repeat annotation GFF3 from all pipeline layers.
    Structure-based annotations (DANTE_LTR, DANTE_TIR, DANTE_LINE) take priority over
    similarity-based ones (RepeatMasker). See annotation_rules.md for full tier hierarchy.
    """
    input:
        ltr=F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3",
        tir=F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_combined.gff3",
        line=F"{config['output_dir']}/DANTE_LINE/DANTE_LINE.gff3",
        dante=F"{config['output_dir']}/DANTE/DANTE_filtered.gff3",
        tc_default=F"{config['output_dir']}/TideCluster/default/TideCluster_clustering.gff3",
        tc_short=F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_clustering.gff3",
        tc_rm=F"{config['output_dir']}/TideCluster/default/RM_on_TideCluster_Library.gff3",
        rm=F"{config['output_dir']}/RepeatMasker/RM_on_combined_library_plus_DANTE.gff3",
        th_default=F"{config['output_dir']}/TideCluster/default/TideCluster_tidehunter_short.gff3",
        th_short=F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_tidehunter_short.gff3",
        fai=F"{config['output_dir']}/genome_cleaned.fasta.fai",
        validation_marker=F"{config['output_dir']}/.classifications_validated"
    output:
        gff=F"{config['output_dir']}/Repeat_Annotation_Unified.gff3"
    log:
        stdout=F"{config['output_dir']}/Repeat_Annotation_Unified.log",
        stderr=F"{config['output_dir']}/Repeat_Annotation_Unified.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/make_unified_annotation.tsv"
    conda:
        "envs/tidecluster.yaml"
    threads: workflow.cores
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        make_unified_annotation.R \
            --ltr      {input.ltr} \
            --tir      {input.tir} \
            --line     {input.line} \
            --dante    {input.dante} \
            --tc_default {input.tc_default} \
            --tc_short   {input.tc_short} \
            --tc_rm      {input.tc_rm} \
            --rm       {input.rm} \
            --th_default {input.th_default} \
            --th_short   {input.th_short} \
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
    input:
        rm=F"{config['output_dir']}/RepeatMasker/Repeat_Annotation_NoSat.gff3",
        sat_tc=F"{config['output_dir']}/TideCluster/default/TideCluster_clustering.gff3",
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
        calculate_statistics_and_make_groups.R -r {input.rm} -s {input.sat_tc} -o {output.csv} -g {input.genome_fasta} -d {output.dir} \
        -M {output.mobile_elements} -S {output.simple_repeats} -L {output.low_complexity} -R {output.rdna} -C {output.all_copia} -G {output.all_gypsy}
        """

rule make_bigwig_density:
    input:
        cvs=F"{config['output_dir']}/summary_statistics.csv",  # this file is available if gffs were created
        genome_seqlengths=F"{config['output_dir']}/genome_seqlengths.rds"
    output:
        checkpoint=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_bigwig/.done"
    params:
        bwdir=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_bigwig",
        gffdir=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3",
        genome_fasta=genome_fasta_cleaned
    log:
        stdout=F"{config['output_dir']}/logs/make_bigwig_density.log",
        stderr=F"{config['output_dir']}/logs/make_bigwig_density.err"
    benchmark:
        F"{config['output_dir']}/benchmarks/make_bigwig_density.tsv"
    conda:
        "envs/tidecluster.yaml"
    shell:
        """
        exec > {log.stdout} 2> {log.stderr}
        set -euo pipefail
        set -x
        mkdir -p {params.bwdir}
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        ls_absolute_path=$(realpath {input.genome_seqlengths})
        calculate_density_batch.R -d {params.gffdir} -o {params.bwdir} -g $ls_absolute_path
        touch {output.checkpoint}
        """

rule add_top_level_outputs:
    input:
        dante=F"{config['output_dir']}/DANTE/DANTE_filtered.gff3",
        dante_ltr=F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3",
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
    input:
        F"{config['output_dir']}/RepeatMasker/Repeat_Annotation_NoSat.gff3",
        F"{config['output_dir']}/TideCluster/default/TideCluster_clustering.gff3",
        genome_seqlengths=F"{config['output_dir']}/genome_seqlengths.rds"
    output:
        F"{config['output_dir']}/RepeatMasker/Repeat_Annotation_NoSat_10k.bw",
        F"{config['output_dir']}/RepeatMasker/Repeat_Annotation_NoSat_100k.bw",
        F"{config['output_dir']}/TideCluster/default/TideCluster_clustering_10k.bw",
        F"{config['output_dir']}/TideCluster/default/TideCluster_clustering_100k.bw"
    log:
        stdout=F"{config['output_dir']}/RepeatMasker/calculate_bigwig_density.log",
        stderr=F"{config['output_dir']}/RepeatMasker/calculate_bigwig_density.err"
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
        calculate_density.R -b {input[0]} -o {output[0]} -f gff3 --window 10000 -g {input.genome_seqlengths}
        calculate_density.R -b {input[0]} -o {output[1]} -f gff3 --window 100000 -g {input.genome_seqlengths}
        calculate_density.R -b {input[1]} -o {output[2]} -f gff3 --window 10000 -g {input.genome_seqlengths}
        calculate_density.R -b {input[1]} -o {output[3]} -f gff3 --window 100000 -g {input.genome_seqlengths}
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
        ln -s -r {input.tc_index} {output.tc_index}
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
        bw1 = F"{config['output_dir']}/RepeatMasker/Repeat_Annotation_NoSat_10k.bw",
        bw2_info = F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_bigwig/.done"
        # these inputs are not used explicitly, but they are necessary for the rule to run
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
        set -euo pipefail
        set -x
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        # command can fail but it should not stop the workflow
        make_summary_plots.R {params.output_dir} {output} || true
        touch {output}
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
        bw_rm = F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_bigwig/.done",
        bw_tc = F"{config['output_dir']}/TideCluster/default/.bigwig_done",
        ltr   = F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3",
        tir   = F"{config['output_dir']}/DANTE_TIR/TIR_classification_summary.txt",
        line  = F"{config['output_dir']}/DANTE_LINE/DANTE_LINE.gff3"
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
        make_repeat_report.R \
            --output_dir {params.output_dir} \
            --bin_width {params.bin_width} \
            --min_len_chart {params.min_len_chart} \
            --min_len_tracks {params.min_len_tracks} \
            --max_tracks {params.max_tracks} \
            --top_sat_clusters {params.top_sat_clusters}
        """
