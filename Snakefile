#configfile: "config.yaml"

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
        F"{config['output_dir']}/Libraries/class_ii_library.fasta",
        F"{config['output_dir']}/Libraries/LTR_RTs_library_clean.fasta",
        F"{config['output_dir']}/Libraries/combined_library.fasta",
        F"{config['output_dir']}/RepeatMasker/RM_on_combined_library.out",
        F"{config['output_dir']}/RepeatMasker/RM_on_combined_library_minus_satellites.gff3",
        F"{config['output_dir']}/Repeat_Annotation_NoSat.gff3",
        F"{config['output_dir']}/all_repeats_for_masking.bed",
        F"{config['output_dir']}/DANTE_LTR.gff3",
        F"{config['output_dir']}/TideCluster_report.html",
        F"{config['output_dir']}/DANTE_LTR_report.html"
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
        html=F"{config['output_dir']}/TideCluster/default/TideCluster_index.html"
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
        library_absolute_path=$(realpath {params.library})
        # NOTE - there is a bug in tidecluster - it does not correctly formal html links, soluton for now is 
        # to run it in the directory where the output will be created
        echo "Library: {input.library}"
        if [ -z "{params.library}" ]; then
            cd $wd
            echo "Running TideCluster without a library"
            TideCluster.py run_all -pr $prefix -c {threads} -M 1 -f $genome_absolute_path
        else
            echo "Running TideCluster with a custom library"
            cd $wd
            TideCluster.py run_all -pr $prefix -c {threads} -M 1 -f $genome_absolute_path -l $library_absolute_path
        fi
        # if gff3_annot was not created but exit code is 0, it means that there were no clusters found, create an empty file
        # but check if gff3_tidehunter was created
        cd $original_dir
        if [ ! -f {output.gff3_clust} ]; then
            # check if gff3_tidehunter was created
            if [ -f {output.gff3_tidehunter} ]; then
                echo "##gff-version 3" > {output.gff3_clust}
                echo "# no clusters found" >> {output.gff3_clust}
            else
                exit 1
            fi
        fi
        """

rule tidecluster_short:
    input:
        genome_fasta=config["genome_fasta"],
        library= config.get("tandem_repeat_library","")
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
        library_absolute_path=$(realpath {params.library})
        # NOTE - there is a bug in tidecluster - it does not correctly formal html links, soluton for now is 
        # to run it in the directory where the output will be created
        echo "Library: {input.library}"
        # TODO - set corretly -M parameter
        if [ -z "{params.library}" ]; then
            cd $wd
            echo "Running TideCluster without a library"
            TideCluster.py run_all -pr $prefix -c {threads} -M 1 -f $genome_absolute_path -T "-p 10 -P 39 -c 5 -e 0.25" -m 5000
        else
            echo "Running TideCluster with a custom library"
            cd $wd
            TideCluster.py run_all -pr $prefix -c {threads} -M 1 -f $genome_absolute_path -l $library_absolute_path -T "-p 10 -P 39 -c 5 -e 0.25" -m 5000
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
                # make empty fasta file
                echo "" > {output.dimer_library_short}
            else
                exit 1
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
    input:
        library=config.get("custom_library", [])
    output:
        library=F"{config['output_dir']}/Libraries/class_ii_library.fasta"
    run:
        # check if the input has attribute library but do not use input.library - it is causing an error is attribute is not present
        if not hasattr(input, "library"):
            print("no library provided")
            print("No custom library provided, creating an empty file")
            with open(output[0], "w") as f:
                f.write("")
        else:
            filter_fasta(input.library, output.library, "Class_II/Subclass_1")


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
        custom_library=F"{config['custom_library']}"
    output:
        F"{config['output_dir']}/Libraries/combined_library.fasta"
    shell:
        """
        cat {input.ltr_rt_library} {input.custom_library} > {output}
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
        "envs/repeatmasker.yaml"
    shell:
        """
        echo $PWD
        # get absolute path of scripts directory
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        merge_repeat_annotations.R {input.rm_gff} {input.dante_gff} {output.gff}
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


rule add_top_level_ouputs:
    input:
        dante=F"{config['output_dir']}/DANTE/DANTE_filtered.gff3",
        dante_ltr=F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3",
        sat_tc=F"{config['output_dir']}/TideCluster/default/TideCluster_clustering.gff3",
        sar_rm=F"{config['output_dir']}/TideCluster/default/RM_on_TideCluster_Library.gff3"
    output:
        dante=F"{config['output_dir']}/DANTE_filtered.gff3",
        dante_ltr=F"{config['output_dir']}/DANTE_LTR.gff3",
        sat_tc=F"{config['output_dir']}/Tandem_repeats_TideCluster.gff3",
        sat_rm=F"{config['output_dir']}/Tandem_repeats_RepeatMasker.gff3"

    shell:
        """
        # make symbolic links to all the outputs
        ln -s -r {input.dante} {output.dante}
        ln -s -r {input.dante_ltr} {output.dante_ltr}
        ln -s -r {input.sat_tc} {output.sat_tc}
        ln -s -r {input.sar_rm} {output.sat_rm}
        sat_tc={output.sat_tc}
        sat_tc_annot=$(echo {output.sat_tc} | sed 's/TideCluster_clustering.gff3/TideCluster_annotation.gff3/')
        echo $sat_tc_annot
        # check if the file exists
        outdir=$(dirname {output.sat_tc})
        if [ -f $sat_tc_annot ]; then
            ln -s -r  $sat_tc_annot $outdir/Tandem_repeats_TideCluster_annotated.gff3
        fi
        
        """


rule add_html_oupouts:
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
