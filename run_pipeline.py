#!/usr/bin/env python
""" this script take as argument config.yaml file and run snakemake pipeline """

import argparse
import os
import subprocess
import yaml

def show_singularity_settings(config_object):
    # get absolute paths and show singularity --bind options
    output_dir = os.path.abspath(config_object['output_dir'])
    # get dirname of output_dir
    genome_dir = os.path.dirname(os.path.abspath(config_object['genome_fasta']))
    custom_library = os.path.dirname(os.path.abspath(config_object['custom_library']))
    tandem_repeat_library = os.path.dirname(os.path.abspath(
        config_object['tandem_repeat_library']))
    # get unique directories
    dirs = set([output_dir, genome_dir, custom_library, tandem_repeat_library])
    bind_string = " ".join([F"-B {d}" for d in dirs])
    print("Run singularity with following bind options:")
    print(F"singularity run {bind_string} ....")



def main():
    # get arguments

    config_template="/opt/pipeline/config.yaml"
    # read config file, keep end of lines
    config_string = open(config_template, 'r').read()
    parser = argparse.ArgumentParser(
            description=
"""Repeat Annotation Pipeline using DANTE, DANTE_LTR, TideCluster adn RepeatMasker.""",
            epilog=F"""Example of config.yaml file:
            
{config_string}

custom_library and tandem_repeat_library are optional.
If custom library is provided it must be in FASTA 
format and header must followthe following format:

>unique_id#family_name/subfamily_name/..

Use following classification scheme:

Class_II/Subclass_1/TIR/EnSpm_CACTA
Class_II/Subclass_1/TIR/hAT
Class_II/Subclass_1/TIR/MITE
Class_II/Subclass_1/TIR/MITE/Stowaway
Class_II/Subclass_1/TIR/MuDR_Mutator
Class_II/Subclass_1/TIR/Tc1_Mariner
Class_II/Subclass_2/Helitron
Class_I/LINE
Class_I/pararetrovirus
rDNA_45S/18S
rDNA_45S/25S
rDNA_45S/5_8S
rDNA_45S/IGS
rDNA_45S/ITS1

            """,
            formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument('-c', '--config', required=True, help='config file')
    parser.add_argument('-t', '--threads', required=False, default=2, type=int,
                        help='Number of threads to use')
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.realpath(__file__))
    snakefile="/opt/pipeline/snakefile"
    # create output directory if it does not exist

    # get conda prefix
    CONDA_ENVS_PATH = os.environ.get('CONDA_ENVS_PATH')
    # NOTE - snake make is using --conda-prefix as path to conda envs while conda
    # CONDA_PREFIX variable points to conda installation directory!
    # run snakemake

    # for subprocess we need to set XDG_CACHE_HOME otherwise snakemake will use
    # non-writable directory
    # load yaml file from args.config
    try:
        config_object = yaml.safe_load(open(args.config))
    except FileNotFoundError:
        # the path is either wrong or path is not mounted
        print(F"Cannot open config file {args.config}")
        print(F"Check if the file exists and is accessible or if the path is mounted "
              F"using -B option in singularity run command")
        exit(1)

    output_dir = config_object['output_dir']

    # this could be relative path, so we need to get absolute path
    output_dir = os.path.abspath(output_dir)
    cache_dir = F"{output_dir}/.cache"
    # check if output_dir exists or can be created (dry-run creates directory)
    try:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    except PermissionError:
        print(F"Cannot create output directory {output_dir}")
        show_singularity_settings(config_object)
        exit(1)

    # if genome accessible
    genome_path = os.path.abspath(config_object['genome_fasta'])
    if not os.path.exists(genome_path):
        print(F"Genome fasta file {genome_path} does not exist or is not accessible")
        show_singularity_settings(config_object)
        exit(1)
    # check if custom_library is defined and files exists
    if 'custom_library' in config_object:
        custom_library = os.path.abspath(config_object['custom_library'])
        if not os.path.exists(custom_library):
            print(F"Custom library file {custom_library} does not exist or is not accessible")
            show_singularity_settings(config_object)
            exit(1)
    # check if tandem_repeat_library is defined and files exists
    if 'tandem_repeat_library' in config_object:
        tandem_repeat_library = os.path.abspath(config_object['tandem_repeat_library'])
        if not os.path.exists(tandem_repeat_library):
            print(F"Tandem repeat library file {tandem_repeat_library} does not exist or is not accessible")
            show_singularity_settings(config_object)
            exit(1)

    cmd = (F"snakemake --snakefile {script_dir}/Snakefile --configfile {args.config} "
           F"--cores {args.threads} --use-conda --conda-prefix {CONDA_ENVS_PATH} "
           F"--conda-frontend mamba --show-failed-logs")

    # append cache dir to other environment variables
    env = os.environ.copy()
    env['XDG_CACHE_HOME'] = cache_dir
    subprocess.check_call(cmd, shell=True, env=env)

if __name__ == "__main__":
    main()
