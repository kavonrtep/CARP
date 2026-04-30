#!/usr/bin/env python
""" this script take as argument config.yaml file and run snakemake pipeline """

import argparse
import os
import subprocess
import sys
import yaml

# Single source of truth for the pipeline version. Importing requires
# version.py to be next to this script, which is true in both contexts
# (repo checkout: both at repo root; container: both at /opt/pipeline/).
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from version import __version__

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

    # The config.yaml template is baked into the container at
    # /opt/pipeline/config.yaml; when running outside the container (CI,
    # local snakemake invocation) fall back to a repo-local copy, and
    # finally to a stub so --help still works.
    script_dir = os.path.dirname(os.path.realpath(__file__))
    config_string = ""
    for candidate in ("/opt/pipeline/config.yaml",
                      os.path.join(script_dir, "config.yaml")):
        if os.path.exists(candidate):
            config_string = open(candidate, "r").read()
            break
    if not config_string:
        config_string = "# (example config not found; see tests/fixtures/ for templates)"
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
    parser.add_argument('--version', action='version',
                        version=F'CARP {__version__}',
                        help='Show pipeline version and exit.')
    parser.add_argument('-c', '--config', required=True, help='config file')
    parser.add_argument('-t', '--threads', required=False, default=2, type=int,
                        help='Number of threads to use')
    parser.add_argument(
        '-S', '--snakemake_args', type=str, nargs='?', required=False, default="",
        help='Additional snakemake arguments, Usage examples: '
             '-S="--dry-run"'
             '-S="--dry-run --reason"'
             'Note: Do not use options --use-conda, --conda-prefix,'
             '--conda-frontend, --show-failed-logs, --cores, --snakefile, --configfile'
             'as they are set by the script.'
        )

    args = parser.parse_args()

    # Startup banner — first line of every run, so logs are
    # self-identifying even without the provenance JSON (which the
    # later commits in the versioning rollout add).
    print(F"CARP pipeline {__version__} starting", file=sys.stderr)

    snakefile="/opt/pipeline/snakefile"
    # create output directory if it does not exist

    # get conda prefix
    # Inside the singularity container CONDA_ENVS_PATH is set to the
    # baked-in envs dir. When running outside (CI, local snakemake), fall
    # back to the repo's .snakemake/conda so cached envs are reused.
    # NOTE - snakemake uses --conda-prefix as path to conda envs, while
    # conda's CONDA_PREFIX variable points to the conda installation dir.
    CONDA_ENVS_PATH = os.environ.get('CONDA_ENVS_PATH') or \
                      os.path.join(script_dir, ".snakemake", "conda")

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
           F"--conda-frontend conda --show-failed-logs --keep-incomplete"
           F" {args.snakemake_args}")

    # append cache dir to other environment variables
    env = os.environ.copy()
    env['XDG_CACHE_HOME'] = cache_dir
    subprocess.check_call(cmd, shell=True, env=env)

if __name__ == "__main__":
    main()
