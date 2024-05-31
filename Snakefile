import os
import re
def create_dirs(*dirs):
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)
print(config)
subdirs = [config['output_dir']+"/"+i for i in ['DANTE', 'DANTE_LTR',
                                                'TideCluster/default',
                                                'TideCluster/short_monomer',
                                                'Libraries', 'RepeatMasker']]
create_dirs(*subdirs)

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

rule all:
    input:
        F"{config['output_dir']}/DANTE/DANTE.gff3",
        F"{config['output_dir']}/DANTE/DANTE_filtered.gff3",
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
        F"{config['output_dir']}/RepeatMasker/RM_on_combined_library_minus_satellites.gff3",
        F"{config['output_dir']}/Repeat_Annotation_NoSat.gff3",
        F"{config['output_dir']}/all_repeats_for_masking.bed",
        F"{config['output_dir']}/DANTE_LTR.gff3",
        F"{config['output_dir']}/TideCluster_report.html",
        F"{config['output_dir']}/DANTE_LTR_report.html",
        F"{config['output_dir']}/gaps_10plus.bed",
        F"{config['output_dir']}/summary_statistics.csv",
        F"{config['output_dir']}/Repeat_Annotation_NoSat_10k.bw",
        F"{config['output_dir']}/Repeat_Annotation_NoSat_100k.bw",
        F"{config['output_dir']}/TideCluster/default/TideCluster_clustering_10k.bw",
        F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_bigwig/.done"

rule dante:
    input:
        config["genome_fasta"],
    output:
        F"{config['output_dir']}/DANTE/DANTE.gff3"
    conda:
        "envs/dante.yaml"
    threads: workflow.cores
    shell:
        """
        dir_out=$(dirname {output})
        gff_out=$(basename {output})
        dante -q {input} -o {output} -c {threads} 
        """

rule filter_dante:
    input:
        F"{config['output_dir']}/DANTE/DANTE.gff3"
    output:
        gff=F"{config['output_dir']}/DANTE/DANTE_filtered.gff3",
        fasta=F"{config['output_dir']}/DANTE/DANTE_filtered.fasta"

    conda:
        "envs/dante.yaml"
    threads: 1
    shell:
        """
        dante_gff_output_filtering.py --dom_gff {input} --domains_filtered {output.gff} --domains_prot_seq {output.fasta}
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
        "envs/dante_ltr.yaml"
    threads: workflow.cores
    shell:
        """
        
        dante_ltr -o {params.prefix} -s {input.fasta} -g {input.gff} -c {threads} -M 1 
        """


rule make_library_of_ltrs:
    input:
        gff3=F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3",
        genome_fasta=F"{config['genome_fasta']}"
    output:
        dir=directory(F"{config['output_dir']}/DANTE_LTR/library"),
        fasta=F"{config['output_dir']}/DANTE_LTR/LTR_RTs_library.fasta"
    conda:
        "envs/dante_ltr.yaml"
    threads: workflow.cores
    shell:
        """
        dante_ltr_to_library --gff {input.gff3} --output_dir {output.dir} -s {input.genome_fasta} -c {threads}
        ln -s library/mmseqs2/mmseqs_representative_seq_clean.fasta {output.fasta}
        """


rule tidecluster_long:
    input:
        genome_fasta=config["genome_fasta"],
        library= config.get("tandem_repeat_library", [])
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
        echo "this is a PATH: $PATH"
        
        cd $wd
        # check if directory TideCluster_clustering_split_files exists
        if [ -d TideCluster_clustering_split_files ]; then
            mkdir -p TideCluster_clustering_split_files_bigwig
            for f in TideCluster_clustering_split_files/*.gff3; do
                base=$(basename $f);
                # remove the extension
                base_noext=$(basename "$base" .gff3)
                base_bw10k=TideCluster_clustering_split_files_bigwig/"$base_noext"_10k.bw
                base_bw100k=TideCluster_clustering_split_files_bigwig/"$base_noext"_100k.bw
                calculate_density.R -b $f -o $base_bw10k -f gff3 --window 10000
                calculate_density.R -b $f -o $base_bw100k -f gff3 --window 100000
            done
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
        gff_absolute_path=$(realpath {output.gff3})
        cd {params.outdir}
        # tc_reannotate.py -s {input.genome_fasta} -f {input.dimer_library_default} -o {output.gff3}
        tc_reannotate.py -s $dl_absolute_path -f $gf_absolute_path -o $gff_absolute_path
        """

rule merge_tidecluster_default_and_short:
    input:
        gff3_default=F"{config['output_dir']}/TideCluster/default/TideCluster_clustering.gff3",
        gff3_short=F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_clustering.gff3"
    output:
        F"{config['output_dir']}/TideCluster/TideCluster_clustering_default_and_short_merged.gff3"
    shell:
        """
        echo $PWD
        cat {input.gff3_default} > {output}
        # replace TRC_ with TRC_S_ in the short monomer clusters to avoid name conflicts
        sed 's/TRC_/TRC_S_/g' {input.gff3_short} >> {output}
        """



rule make_subclass_2_library:
    params:
        library=config.get("custom_library", "")
    output:
        library=F"{config['output_dir']}/Libraries/class_ii_library.fasta"
    run:
        if params.library:
            print("Custom library provided, filtering FASTA.")
            filter_fasta(params.library, output.library, "Class_II/Subclass_1")
        else:
            print("No custom library provided, creating an empty file.")
            with open(output.library, "w") as f:
                f.write("")


rule filter_ltr_rt_library:
    input:
        dante_library=F"{config['output_dir']}/DANTE_LTR/LTR_RTs_library.fasta",
        subclass_2_library=F"{config['output_dir']}/Libraries/class_ii_library.fasta"
    output:
        library=F"{config['output_dir']}/Libraries/LTR_RTs_library_clean.fasta"
    conda:
        "envs/blast.yaml"

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
    output:
        F"{config['output_dir']}/Libraries/combined_library.fasta"
    params:
        custom_library = config.get("custom_library", "")
    shell:
        """
        if [ -z "{params.custom_library}" ]; then
            cp {input.ltr_rt_library} {output}
        else
            cat {input.ltr_rt_library} {params.custom_library} > {output}
        fi
        """


rule repeatmasker:
    input:
        genome_fasta=config["genome_fasta"],
        library=F"{config['output_dir']}/Libraries/combined_library.fasta"
    output:
        out=F"{config['output_dir']}/RepeatMasker/RM_on_combined_library.out",
        gff=F"{config['output_dir']}/RepeatMasker/RM_on_combined_library.gff3"
    params:
        rm_dir=directory(F"{config['output_dir']}/RepeatMasker")
    conda:
        "envs/repeatmasker.yaml"
    threads: workflow.cores
    shell:
        """
        echo $PWD
        # get absolute path of scripts directory
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        library_absolute_path=$(realpath {input.library})
        genome_absolute_path=$(realpath {input.genome_fasta})
        out_absolute_path=$(realpath {output.out})
        gff_absolute_path=$(realpath {output.gff})
        cd {params.rm_dir}
        
        RepeatMasker -pa {threads}  $genome_absolute_path  -lib $library_absolute_path -dir $PWD -s -xsmall -e ncbi -no_is
        mv `basename {input.genome_fasta}`.out $out_absolute_path
        mkdir -p RM_files
        mv  `basename {input.genome_fasta}`.* RM_files
        clean_rm_output.R $out_absolute_path $gff_absolute_path
        """

rule subtract_satellites_from_rm:
    input:
        rm_gff=F"{config['output_dir']}/RepeatMasker/RM_on_combined_library.gff3",
        satellite_annotation=F"{config['output_dir']}/TideCluster/TideCluster_clustering_default_and_short_merged.gff3"
    output:
        F"{config['output_dir']}/RepeatMasker/RM_on_combined_library_minus_satellites.gff3"
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        bedtools subtract -a {input.rm_gff} -b {input.satellite_annotation} -A > {output}
        """

rule merge_rm_and_dante:
    input:
        rm_gff=F"{config['output_dir']}/RepeatMasker/RM_on_combined_library_minus_satellites.gff3",
        dante_gff=F"{config['output_dir']}/DANTE/DANTE_filtered.gff3"
    output:
        gff=F"{config['output_dir']}/Repeat_Annotation_NoSat.gff3"
    conda:
        "envs/dante_ltr.yaml"
        # dante_ltr is already used and it contains the necessary tools (rtracklayer and optparse)
    shell:
        """
        echo $PWD
        # get absolute path of scripts directory
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        # make temp file from date - names in gff3 must be consistent with names used for RepeatMasker
        clean_DANTE_names.R {input.dante_gff}  {input.dante_gff}.tmp.gff3
        merge_repeat_annotations.R {input.rm_gff} {input.dante_gff}.tmp.gff3 {output.gff}
        """


rule make_track_for_masking:
    input:
        rm=F"{config['output_dir']}/Repeat_Annotation_NoSat.gff3",
        tr_main=F"{config['output_dir']}/TideCluster/TideCluster_clustering_default_and_short_merged.gff3",
        tr_default_short=F"{config['output_dir']}/TideCluster/default/TideCluster_tidehunter_short.gff3",
        tr_short_short=F"{config['output_dir']}/TideCluster/short_monomer/TideCluster_tidehunter_short.gff3"
    output:
        F"{config['output_dir']}/all_repeats_for_masking.bed"
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        # concatenate and sort the gff files
        export LC_COLLATE=C
        cat {input.rm} {input.tr_main} {input.tr_default_short} {input.tr_short_short} | sort -k1,1 -k4,4n > {output}.tmp.gff3
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
        rm=F"{config['output_dir']}/Repeat_Annotation_NoSat.gff3",
        sat_tc=F"{config['output_dir']}/TideCluster/default/TideCluster_clustering.gff3",
        genome_fasta=config["genome_fasta"]
    output:
        csv=F"{config['output_dir']}/summary_statistics.csv",
        dir=directory(F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3")
    conda:
        "envs/dante_ltr.yaml"
    shell:
        """
        # get absolute path of scripts directory
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH

        calculate_statistics.R -r {input.rm} -s {input.sat_tc} -o {output.csv} -g {input.genome_fasta} -d {output.dir}
        """

rule make_bigwig_density:
    input:
        gffdir=directory(F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_gff3"),
        cvs=F"{config['output_dir']}/summary_statistics.csv"  # this file is available if gffs were created
    output:
        bwdir=directory(F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_bigwig"),
        checkpoint=F"{config['output_dir']}/Repeat_Annotation_NoSat_split_by_class_bigwig/.done"

    conda:
        "envs/dante_ltr.yaml"
    shell:
        """
        mkdir -p {output.bwdir}
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        
        for f in {input.gffdir}/*.gff3; do
            base=$(basename $f);
            # remove the extension
            base_noext=$(basename "$base" .gff3)
            base_bw10k={output.bwdir}/"$base_noext"_10k.bw
            base_bw100k={output.bwdir}/"$base_noext"_100k.bw
            calculate_density.R -b $f -o $base_bw10k -f gff3 --window 100000
            calculate_density.R -b $f -o $base_bw100k -f gff3 --window 1000000
        done
        touch {output.checkpoint}
        """

rule add_top_level_outputs:
    input:
        dante=F"{config['output_dir']}/DANTE/DANTE_filtered.gff3",
        dante_ltr=F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3",
        sat_tc=F"{config['output_dir']}/TideCluster/default/TideCluster_clustering.gff3",
        sat_tc_bw10=F"{config['output_dir']}/TideCluster/default/TideCluster_clustering_10k.bw",
        sat_tc_bw100=F"{config['output_dir']}/TideCluster/default/TideCluster_clustering_100k.bw",
        sat_rm=F"{config['output_dir']}/TideCluster/default/RM_on_TideCluster_Library.gff3"
    output:
        dante=F"{config['output_dir']}/DANTE_filtered.gff3",
        dante_ltr=F"{config['output_dir']}/DANTE_LTR.gff3",
        sat_tc=F"{config['output_dir']}/Tandem_repeats_TideCluster.gff3",
        sat_tc_bw10=F"{config['output_dir']}/Tandem_repeats_TideCluster_10k.bw",
        sat_tc_bw100=F"{config['output_dir']}/Tandem_repeats_TideCluster_100k.bw",
        sat_rm=F"{config['output_dir']}/Tandem_repeats_RepeatMasker.gff3"
    params:
        sat_tc_annot_in=F"{config['output_dir']}/TideCluster/default/TideCluster_annotation.gff3",
        sat_tc_annot_out=F"{config['output_dir']}/Tandem_repeats_TideCluster_annotated.gff3"
    shell:
        """
        # make symbolic links to all the outputs
        ln -fs -r {input.dante} {output.dante}
        ln -fs -r {input.dante_ltr} {output.dante_ltr}
        ln -fs -r {input.sat_tc} {output.sat_tc}
        ln -fs -r {input.sat_rm} {output.sat_rm}
        ln -fs -r {input.sat_tc_bw10} {output.sat_tc_bw10}
        ln -fs -r {input.sat_tc_bw100} {output.sat_tc_bw100}
        sat_tc={output.sat_tc}

        # check if the file exists
        outdir=$(dirname {output.sat_tc})
        if [ -f {params.sat_tc_annot_in} ]; then
             ln -fs -r {params.sat_tc_annot_in} {params.sat_tc_annot_out}
        fi
        """

rule calculate_bigwig_density:
    input:
        F"{config['output_dir']}/Repeat_Annotation_NoSat.gff3",
        F"{config['output_dir']}/TideCluster/default/TideCluster_clustering.gff3",
    output:
        F"{config['output_dir']}/Repeat_Annotation_NoSat_10k.bw",
        F"{config['output_dir']}/Repeat_Annotation_NoSat_100k.bw",
        F"{config['output_dir']}/TideCluster/default/TideCluster_clustering_10k.bw",
        F"{config['output_dir']}/TideCluster/default/TideCluster_clustering_100k.bw"
    conda:
        "envs/dante_ltr.yaml"
    shell:
        """
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        echo "this is a PATH: $PATH"
        calculate_density.R -b {input[0]} -o {output[0]} -f gff3 --window 10000
        calculate_density.R -b {input[0]} -o {output[1]} -f gff3 --window 100000
        calculate_density.R -b {input[1]} -o {output[2]} -f gff3 --window 10000
        calculate_density.R -b {input[1]} -o {output[3]} -f gff3 --window 100000
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
