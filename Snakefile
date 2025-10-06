import os
import re
def create_dirs(*dirs):
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)
print(config)
subdirs = [config['output_dir']+"/"+i for i in ['DANTE', 'DANTE_TIR', 'DANTE_LINE', 'DANTE_LTR',
                                                'TideCluster/default',
                                                'TideCluster/short_monomer',
                                                'Libraries', 'RepeatMasker']]
create_dirs(*subdirs)
snakemake_dir = os.path.dirname(workflow.snakefile)
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
    "sensitive": "-s",
    "quick": "-q",
    "" : ""
    }[config["repeatmasker_sensitivity"]]

if "reduce_library" not in config:
    config["reduce_library"] = True
else:
    # check validity of the value
    if config["reduce_library"] not in [True, False]:
        raise ValueError("Invalid value for reduce_library_size. Must be either True or False.")

# repeatmasker engine:
if "repeatmasker_engine" not in config:
    config["repeatmasker_engine"] = "ncbi"
rm_engine = config["repeatmasker_engine"]
# check that the engine is valid - must be either ncbi or abblast
if rm_engine not in ["ncbi", "abblast"]:
    raise ValueError("Invalid RepeatMasker engine. Must be either 'ncbi' or 'abblast'.")

rule all:
    input:
        F"{config['output_dir']}/DANTE/DANTE.gff3",
        F"{config['output_dir']}/DANTE/DANTE_filtered.gff3",
        F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_final.gff3",
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
        F"{config['output_dir']}/summary_plots.pdf"

rule dante:
    input:
        config["genome_fasta"],
    output:
        F"{config['output_dir']}/DANTE/DANTE.gff3"
    conda:
        "envs/tidecluster.yaml"
    threads: workflow.cores
    shell:
        """
        dir_out=$(dirname {output})
        gff_out=$(basename {output})
        dante -q {input} -o {output} -c {threads}
        """

rule dante_tir:
    input:
        gff=F"{config['output_dir']}/DANTE/DANTE.gff3",
        fasta=config["genome_fasta"]
    output:
        checkpoint=F"{config['output_dir']}/DANTE_TIR/.done",
        gff=F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_final.gff3",
        fasta=F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_final.fasta",
        summary=F"{config['output_dir']}/DANTE_TIR/TIR_classification_summary.txt",
        dante_tir_lib=F"{config['output_dir']}/DANTE_TIR/all_representative_elements_min3.fasta"
    params:
        output_dir=F"{config['output_dir']}/DANTE_TIR"
    conda:
        "envs/dante_tir.yaml"
    threads: workflow.cores
    shell:
        """
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

        # Ensure all expected output files exist (create empty ones if needed)
        touch {output.gff} {output.fasta} {output.summary}

        # Create checkpoint file to indicate completion
        touch {output.checkpoint}
        """




rule filter_dante:
    input:
        F"{config['output_dir']}/DANTE/DANTE.gff3"
    output:
        gff=F"{config['output_dir']}/DANTE/DANTE_filtered.gff3",
        fasta=F"{config['output_dir']}/DANTE/DANTE_filtered.fasta"

    conda:
        "envs/tidecluster.yaml"
    threads: 1
    shell:
        """
        dante_gff_output_filtering.py --dom_gff {input} --domains_filtered {output.gff} --domains_prot_seq {output.fasta}
        """

rule dante_line:
    input:
        gff=F"{config['output_dir']}/DANTE/DANTE_filtered.gff3",
        genome=config["genome_fasta"]
    output:
        line_rep_lib=F"{config['output_dir']}/DANTE_LINE/LINE_rep_lib.fasta",
        gff_out=F"{config['output_dir']}/DANTE_LINE/DANTE_LINE.gff3",
        line_regions=F"{config['output_dir']}/DANTE_LINE/LINE_regions.fasta",
        line_regions_extended=F"{config['output_dir']}/DANTE_LINE/LINE_regions_extended.fasta"
    params:
        output_dir=F"{config['output_dir']}/DANTE_LINE"
    conda:
        "envs/dante_line.yaml"
    threads: workflow.cores
    shell:
        """
        # Add scripts directory to PATH
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH

        # Run dante_line.py - it may fail if no LINE elements are found
        dante_line.py -g {input.genome} -a {input.gff} -o {params.output_dir} -t {threads} || true

        # Ensure all output files exist (create empty ones if needed)
        touch {output.gff_out} {output.line_regions} {output.line_regions_extended} {output.line_rep_lib}
        """

rule dante_ltr:
    input:
        fasta=config["genome_fasta"],
        gff=F"{config['output_dir']}/DANTE/DANTE.gff3"

    output:
        gff = F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3",
        html = F"{config['output_dir']}/DANTE_LTR/DANTE_LTR_summary.html"
    params:
        prefix = lambda wildcards, output: output.gff.replace(".gff3", "")
    conda:
        "envs/tidecluster.yaml"
    threads: workflow.cores
    shell:
        """
        
        dante_ltr -o {params.prefix} -s {input.fasta} -g {input.gff} -c {threads} -M 1 -S 50000000
        # if exit status is 0 and gff3 file was created but html is missing
        #create an empty file
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
        genome_fasta=F"{config['genome_fasta']}"
    output:
        dir=directory(F"{config['output_dir']}/DANTE_LTR/library"),
        fasta=F"{config['output_dir']}/DANTE_LTR/LTR_RTs_library.fasta"
    conda:
        "envs/tidecluster.yaml"
    threads: workflow.cores
    shell:
        """
        # run only gff3 contains some records
        # (number of lines not starting with # is greater than 1)
        # but check just first 30 lines
        if [ $(head -n 30 {input.gff3} | grep -v "^#" | wc -l) -gt 1 ]; then
            dante_ltr_to_library --gff {input.gff3} --output_dir {output.dir} -s {input.genome_fasta} -c {threads}
            ln -s library/mmseqs2/mmseqs_representative_seq_clean.fasta {output.fasta}
        else
            echo "No LTR-RTs found, creating an empty file"
            : > {output.fasta}
            mkdir -p {output.dir}
        fi
        
        """


rule tidecluster_long:
    input:
        genome_fasta=config["genome_fasta"],
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
    conda:
        "envs/tidecluster.yaml"
    threads: workflow.cores

    shell:
        """
        wd=$(dirname {output.gff3_clust})
        prefix=$(basename {params.prefix})
        original_dir=$PWD
        genome_absolute_path=$(realpath {input.genome_fasta})
        genome_seqlengths=$(realpath {input.genome_seqlengths})
        # define library_absolute_path only if it is not empty
        
        # NOTE - there is a bug in tidecluster - it does not correctly formal html links, soluton for now is 
        # to run it in the directory where the output will be created
        echo "Library: {input.library}"
        if [ -z "{params.library}" ]; then
            cd $wd
            echo "Running TideCluster without a library"
            TideCluster.py run_all -pr $prefix -c {threads}  -f $genome_absolute_path
        else
            library_absolute_path=$(realpath {params.library})
            echo "Running TideCluster with a custom library"
            cd $wd
            TideCluster.py run_all -pr $prefix -c {threads} -f $genome_absolute_path -l $library_absolute_path
        fi
        # if gff3_annot was not created but exit code is 0, it means that there were no clusters found, create an empty file
        # but check if gff3_tidehunter was created
        cd $original_dir
        if [ ! -f {output.gff3_clust} ]; then
            # check if gff3_tidehunter was created
            if [ -f {output.gff3_tidehunter} ]; then
                echo "##gff-version 3" > {output.gff3_clust}
                echo "# no clusters found" >> {output.gff3_clust}
            fi
        fi
        if [ ! -f {output.dimer_library_default} ]; then
            # check if gff3_tidehunter was created
            if [ -f {output.gff3_tidehunter} ]; then
                # make empty fasta file
                : > {output.dimer_library_default}
            fi
        fi
        
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        cd $wd
        # check if directory TideCluster_clustering_split_files exists

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
        genome_fasta=config["genome_fasta"],
    output:
        gff3_clust=F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_clustering.gff3",
        gff3_tidehunter=F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_tidehunter.gff3",
        dimer_library_short=F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_consensus_dimer_library.fasta",
        tr_short_short=F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_tidehunter_short.gff3"
    params:
        prefix = lambda wildcards, output: output.gff3_clust.replace("_clustering.gff3", ""),
        library = config.get("tandem_repeat_library", "")
    conda:
        "envs/tidecluster.yaml"
    threads: workflow.cores

    shell:
        """
        wd=$(dirname {output.gff3_clust})
        prefix=$(basename {params.prefix})
        original_dir=$PWD
        genome_absolute_path=$(realpath {input.genome_fasta})
        # NOTE - there is a bug in tidecluster - it does not correctly formal html links, soluton for now is 
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
        # if gff3_annot was not created but exit code is 0, it means that there were no clusters found, create an empty file
        # but check if gff3_tidehunter was created
        cd $original_dir
        # if gff3_clust was not created but exit code is 0, it means that there were no clusters found, create an empty file
        # but check if gff3_tidehunter was created
        if [ ! -f {output.gff3_clust} ]; then
            # check if gff3_tidehunter was created
            if [ -f {output.gff3_tidehunter} ]; then
                echo "##gff-version 3" > {output.gff3_clust}
                echo "# no clusters found" >> {output.gff3_clust}
            fi
        fi
        if [ ! -f {output.dimer_library_short} ]; then
            # check if gff3_tidehunter was created
            if [ -f {output.gff3_tidehunter} ]; then
                # make empty fasta file
                : > {output.dimer_library_short}
            fi
        fi
        """

rule tidecluster_reannotate:
    input:
        genome_fasta=config["genome_fasta"],
        dimer_library_default=F"{config['output_dir']}/TideCluster/default/TideCluster_consensus_dimer_library.fasta",
    output:
        gff3=F"{config['output_dir']}/TideCluster/default/RM_on_TideCluster_Library.gff3"
    params:
        outdir=directory(F"{config['output_dir']}/TideCluster"),
    conda:
        "envs/tidecluster.yaml"
    threads: workflow.cores
    shell:
        """
        # run only if dimer_library_default is not empty
        if [ ! -s {input.dimer_library_default} ]; then
            echo "No dimer library found, skipping reannotation"
            : > {output}
            exit 0
        fi
        gf_absolute_path=$(realpath {input.dimer_library_default})
        dl_absolute_path=$(realpath {input.genome_fasta})
        dl_basename=$(basename {input.genome_fasta})
        gff_absolute_path=$(realpath {output.gff3})
        
        cd {params.outdir}
        cp $dl_absolute_path .
        # tc_reannotate.py -s {input.genome_fasta} -f {input.dimer_library_default} -o {output.gff3}
        tc_reannotate.py -s $dl_basename -f $gf_absolute_path -o $gff_absolute_path  -c {threads}
        rm $dl_basename
        """

rule merge_tidecluster_default_and_short:
    input:
        gff3_default=F"{config['output_dir']}/TideCluster/default/TideCluster_clustering.gff3",
        gff3_short=F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_clustering.gff3"
    output:
        F"{config['output_dir']}/TideCluster/TideCluster_clustering_default_and_short_merged.gff3"
    shell:
        """
        cat {input.gff3_default} > {output}
        # replace TRC_ with TRC_S_ in the short monomer clusters to avoid name conflicts
        sed 's/TRC_/TRC_S_/g' {input.gff3_short} >> {output}
        """



rule make_subclass_2_library:
    params:
        library=config.get("custom_library", "")
    input:
        dante_tir_lib=F"{config['output_dir']}/DANTE_TIR/all_representative_elements_min3.fasta",
    output:
        library=F"{config['output_dir']}/Libraries/class_ii_library.fasta"
    run:
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
    conda:
        "envs/tidecluster.yaml"


    shell:
        """
        # first reformat header( convert | to / and / to _) and then filter using blast 
        # use sed to replace | with / and / with _ in the header
        
        sed  's/\//_/g' {input.dante_library} |  sed 's/|/\//g'  > {input.dante_library}.reformatted
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
        dante_tir_lib=F"{config['output_dir']}/DANTE_TIR/all_representative_elements_min3.fasta",
        line_rep_lib=F"{config['output_dir']}/DANTE_LINE/LINE_rep_lib.fasta"
    output:
        full_names=F"{config['output_dir']}/Libraries/combined_library.fasta",
        short_names=F"{config['output_dir']}/Libraries/combined_library_short_names.fasta",
    params:
        custom_library = config.get("custom_library", ""),
        rdna_library = os.path.join(snakemake_dir, "data/rdna_library.fasta")
    shell:
        """
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
    conda: "envs/tidecluster.yaml"
    threads: workflow.cores
    shell:
        """
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        workdir=$(dirname {output.library_reduced})/workdir
        # if reduce_library_size is set to False, just copy the input library
        echo "Reduce library size: {params.reduce_library_size}"
        if [ "{params.reduce_library_size}" = "False" ]; then
            cp  {input.library} {output.library_reduced}
            exit 0
        fi
        reduce_library_size.R -i {input.library} -o {output.library_reduced} -t {threads} -d $workdir 
        """

rule repeatmasker:
    input:
        genome_fasta=config["genome_fasta"],
        library=F"{config['output_dir']}/Libraries/combined_library.fasta",
        library_short=F"{config['output_dir']}/Libraries/combined_library_short_names.fasta",
        library_reduced=F"{config['output_dir']}/Libraries/combined_library_reduced.fasta"


    output:
        out=F"{config['output_dir']}/RepeatMasker/RM_on_combined_library.out",
        gff=F"{config['output_dir']}/RepeatMasker/RM_on_combined_library.gff3"
    params:
        rm_dir=directory(F"{config['output_dir']}/RepeatMasker"),
        rm_sensitivity_option=rm_sensitivity_option,
        rm_engine=rm_engine,
        rm_sensitivity=config["repeatmasker_sensitivity"]
    conda:
        "envs/tidecluster.yaml"
    threads: workflow.cores
    shell:
        """
        echo $PWD
        # get absolute path of scripts directory
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
        # user repeatmasker wrapper
        repeatmasker_wrapper.py -f $gen_name -l $lib_name -o $out_absolute_path  -s {params.rm_sensitivity} -p {threads} -d workdir 
        clean_rm_output.R $out_absolute_path $gff_absolute_path
        """

rule subtract_satellites_from_rm:
    input:
        rm_gff=F"{config['output_dir']}/RepeatMasker/RM_on_combined_library_plus_DANTE.gff3",
        satellite_annotation=F"{config['output_dir']}/TideCluster/TideCluster_clustering_default_and_short_merged.gff3"
    output:
        F"{config['output_dir']}/RepeatMasker/Repeat_Annotation_NoSat.gff3"
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        bedtools subtract -a {input.rm_gff} -b {input.satellite_annotation} -A > {output}
        """

rule merge_rm_and_dante:
    input:
        rm_gff=F"{config['output_dir']}/RepeatMasker/RM_on_combined_library.gff3",
        dante_gff=F"{config['output_dir']}/DANTE/DANTE_filtered.gff3"
    output:
        gff=F"{config['output_dir']}/RepeatMasker/RM_on_combined_library_plus_DANTE.gff3"
    conda:
        "envs/tidecluster.yaml"
        # dante_ltr is already used and it contains the necessary tools (rtracklayer and optparse)
    threads: workflow.cores
    shell:
        """
        export CPU_COUNT={threads}
        # get absolute path of scripts directory
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        # make temp file from date - names in gff3 must be consistent with names used for RepeatMasker
        clean_DANTE_names.R {input.dante_gff}  {input.dante_gff}.tmp.gff3
        merge_repeat_annotations.R {input.rm_gff} {input.dante_gff}.tmp.gff3 {output.gff}
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
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        # concatenate and sort the gff files
        export LC_COLLATE=C
        cat {input.rm} {input.tr_main} {input.tr_default_short} {input.tr_short_short} {input.tr_rm} {input.dante_ltr} | sort -k1,1 -k4,4n > {output}.tmp.gff3
        bedtools merge -i {output}.tmp.gff3 > {output}
        rm {output}.tmp.gff3
        """

rule make_track_for_Ns:
    input:
        genome_fasta=config["genome_fasta"]
    output:
        F"{config['output_dir']}/gaps_10plus.bed"
    conda:
        "envs/seqtk.yaml"
    shell:
        """
        seqtk cutN -n 10 -g {input.genome_fasta} >  {output}
        """

rule make_summary_statistics_and_split_by_class:
    input:
        rm=F"{config['output_dir']}/RepeatMasker/Repeat_Annotation_NoSat.gff3",
        sat_tc=F"{config['output_dir']}/TideCluster/default/TideCluster_clustering.gff3",
        genome_fasta=config["genome_fasta"]
    output:
        csv=F"{config['output_dir']}/summary_statistics.csv",
        dir=directory(F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3"),
        mobile_elements=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3/Mobile_elements.gff3",
        simple_repeats=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3/Simple_repeats.gff3",
        low_complexity=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3/Low_complexity.gff3",
        rdna=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3/rDNA.gff3",
        all_copia=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3/All_Ty1_Copia.gff3",
        all_gypsy=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3/All_Ty3_Gypsy.gff3",

    conda:
        "envs/tidecluster.yaml"
    shell:
        """
        # get absolute path of scripts directory
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
        genome_fasta=config["genome_fasta"]
    conda:
        "envs/tidecluster.yaml"
    shell:
        """
        mkdir -p {params.bwdir}
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        echo "Calculating bigwig densities"
        ls_absolute_path=$(realpath {input.genome_seqlengths})
        calculate_density_batch.R -d {params.gffdir} -o {params.bwdir} -g $ls_absolute_path 
        touch {output.checkpoint}
        """

rule add_top_level_outputs:
    input:
        dante=F"{config['output_dir']}/DANTE/DANTE_filtered.gff3",
        dante_ltr=F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3",
        dante_tir=F"{config['output_dir']}/DANTE_TIR/DANTE_TIR_final.gff3",
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
        sat_tc={output.sat_tc}

        # check if the file exists
        outdir=$(dirname {output.sat_tc})
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
    conda:
        "envs/tidecluster.yaml"
    shell:
        """
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
    shell:
        """
        ln -s -r {input.tc_index} {output.tc_index}
        ln -s -r {input.dante_ltr_index} {output.dante_ltr_index}
        """


rule calculate_seqlengths:
    input:
        genome_fasta=config["genome_fasta"]
    output:
        F"{config['output_dir']}/genome_seqlengths.rds"
    conda:
        "envs/tidecluster.yaml"
    shell:
        """
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
    conda:
        "envs/tidecluster.yaml"
    shell:
        """
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        # command can fail but it should not stop the workflow
        make_summary_plots.R {params.output_dir} {output} || true
        touch output
        """
