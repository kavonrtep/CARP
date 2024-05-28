configfile: "config.yaml"
import os
def create_dirs(*dirs):
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

subdirs = [config['output_dir']+"/"+i for i in ['DANTE', 'DANTE_LTR', 'TideCluster']]
create_dirs(*subdirs)

rule all:
    input:
        F"{config['output_dir']}/DANTE/DANTE.gff3",
        F"{config['output_dir']}/DANTE/DANTE_filtered.gff3",
        F"{config['output_dir']}/DANTE_LTR/DANTE_LTR.gff3",
        F"{config['output_dir']}/DANTE_LTR/LTR_RTs_library.fasta",
        F"{config['output_dir']}/TideCluster/TideCluster_default_clustering.gff3",
        F"{config['output_dir']}/TideCluster/TideCluster_short_monomer_clustering.gff3",
        F"{config['output_dir']}/TideCluster/TideCluster_default_similarity_based_annotation.gff3",
        F"{config['output_dir']}/TideCluster/TideCluster_default_and_short_clustering_merged.gff3"
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
        genome_fasta=config["genome_fasta"]
    output:
        gff3_clust=F"{config['output_dir']}/TideCluster/TideCluster_default_clustering.gff3",
        gff3_tidehunter=F"{config['output_dir']}/TideCluster/TideCluster_default_tidehunter.gff3",
        dimer_library_default=F"{config['output_dir']}/TideCluster/TideCluster_default_consensus_dimer_library.fasta",
    params:
        prefix = lambda wildcards, output: output.gff3_clust.replace("_clustering.gff3", "")
    conda:
        "envs/tidecluster.yaml"
    threads: workflow.cores

    shell:
        """
        TideCluster.py run_all -pr {params.prefix} -c {threads} -M 1 -f {input.genome_fasta} 
        # if gff3_annot was not created but exit code is 0, it means that there were no clusters found, create an empty file
        # but check if gff3_tidehunter was created
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
        genome_fasta=config["genome_fasta"]
    output:
        gff3_clust=F"{config['output_dir']}/TideCluster/TideCluster_short_monomer_clustering.gff3",
        gff3_tidehunter=F"{config['output_dir']}/TideCluster/TideCluster_short_monomer_tidehunter.gff3",
        dimer_library_short=F"{config['output_dir']}/TideCluster/TideCluster_short_consensus_dimer_library.fasta"

    params:
        prefix = lambda wildcards, output: output.gff3_clust.replace("_clustering.gff3", "")
    conda:
        "envs/tidecluster.yaml"
    threads: workflow.cores

    shell:
        """
        TideCluster.py run_all -pr {params.prefix} -c {threads} -M 1 -f {input.genome_fasta}  -T "-p 10 -P 39 -c 5 -e 0.25" -m 5000
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
        outdir=F"{config['output_dir']}/TideCluster",
        dimer_library_default=F"{config['output_dir']}/TideCluster/TideCluster_default_consensus_dimer_library.fasta",
    output:
        gff3=F"{config['output_dir']}/TideCluster/TideCluster_default_similarity_based_annotation.gff3",
    conda:
        "envs/tidecluster.yaml"
    threads: workflow.cores
    shell:
        """
        gf_absolute_path=$(realpath {input.dimer_library_default})
        dl_absolute_path=$(realpath {input.genome_fasta})
        gff_absolute_path=$(realpath {output.gff3})
        cd {input.outdir}
        # tc_reannotate.py -s {input.genome_fasta} -f {input.dimer_library_default} -o {output.gff3}
        tc_reannotate.py -s $dl_absolute_path -f $gf_absolute_path -o $gff_absolute_path
        """

rule merge_tidecluster_default_and_short:
    input:
        gff3_default=F"{config['output_dir']}/TideCluster/TideCluster_default_clustering.gff3",
        gff3_short=F"{config['output_dir']}/TideCluster/TideCluster_short_monomer_clustering.gff3"
    output:
        F"{config['output_dir']}/TideCluster/TideCluster_default_and_short_clustering_merged.gff3"
    shell:
        """
        echo $PWD
        cat {input.gff3_default} {input.gff3_short} > {output}
        """