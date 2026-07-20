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
    parser.add_argument(
        '--keep-all', action='store_true',
        help='Keep all intermediate files (overrides config cleanup_intermediates, '
             'forcing "none"). By default a successful run deletes per-tool scratch '
             'per the cleanup_intermediates config key (default "minimal").')

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

    # Phase-1 provenance: write run_provenance.json *before* snakemake
    # is invoked. If conda env solve dies, the user still has a record
    # of what was attempted. Skip on dry-runs (no real run, no need to
    # spam the user's output dir with stale provenance).
    is_dry_run = ("--dry-run" in args.snakemake_args
                  or "-n" in args.snakemake_args.split())
    finalise_fn = None
    extend_fn = None
    manifest_finalise_fn = None
    if not is_dry_run:
        try:
            sys.path.insert(0, os.path.join(script_dir, "scripts"))
            from record_provenance import (init_provenance,
                                            finalise_provenance,
                                            extend_with_envs)
            init_provenance(args.config, output_dir)
            finalise_fn = finalise_provenance
            extend_fn = extend_with_envs
        except Exception as e:
            print(F"WARNING: failed to write phase-1 provenance: {e}",
                  file=sys.stderr)
        # Output manifest (FR-3): carp_manifest.json — a SEPARATE artifact from
        # provenance (schema + file map, not a completion gate), written by the
        # same wrapper hook so manifest and provenance can't drift and both are
        # present on success and failure.
        try:
            sys.path.insert(0, os.path.join(script_dir, "scripts"))
            from manifest import init_manifest, finalise_manifest
            init_manifest(output_dir)
            manifest_finalise_fn = finalise_manifest
        except Exception as e:
            print(F"WARNING: failed to write output manifest: {e}",
                  file=sys.stderr)

    # Phase-2 provenance prep: ensure conda envs are created, then
    # amend run_provenance.json with per-env primary-package versions.
    # Splitting --conda-create-envs-only into a pre-step (rather than
    # waiting for the heavy pipeline rules) means provenance includes
    # env data even if a later rule fails. On a warm conda cache this
    # adds ~5 s of snakemake startup; on cold, it's no slower than the
    # main run would have been creating the envs lazily anyway.
    if extend_fn is not None:
        prep_cmd = (
            F"snakemake --snakefile {script_dir}/Snakefile --configfile {args.config} "
            F"--cores {args.threads} --use-conda --conda-prefix {CONDA_ENVS_PATH} "
            F"--conda-frontend conda --conda-create-envs-only --quiet")
        prep_rc = subprocess.call(prep_cmd, shell=True, env=env)
        if prep_rc == 0:
            try:
                from pathlib import Path as _Path
                extend_fn(
                    _Path(output_dir) / "run_provenance.json",
                    _Path(script_dir) / "envs",
                    _Path(CONDA_ENVS_PATH))
            except Exception as e:
                print(F"WARNING: failed to write phase-2 provenance: {e}",
                      file=sys.stderr)
        else:
            print(F"WARNING: conda env creation failed (exit {prep_rc}); "
                  F"phase-2 provenance skipped, phase-3 will record failure",
                  file=sys.stderr)

    # Run snakemake. Use subprocess.call (not check_call) so we can
    # finalise provenance with the exit status before propagating the
    # failure to our own exit code.
    rc = subprocess.call(cmd, shell=True, env=env)

    if finalise_fn is not None:
        try:
            finalise_fn(
                output_dir,
                exit_status=("completed" if rc == 0 else "failed"))
        except Exception as e:
            print(F"WARNING: failed to finalise provenance: {e}",
                  file=sys.stderr)

    # Flip the manifest's exit_status before we return (and thus before any
    # scratch->NFS copy-back and the exit-code DONE markers), so it travels
    # with the outputs. Never gates completion — a failure here only warns.
    if manifest_finalise_fn is not None:
        try:
            manifest_finalise_fn(
                output_dir,
                "completed" if rc == 0 else "failed")
        except Exception as e:
            print(F"WARNING: failed to finalise output manifest: {e}",
                  file=sys.stderr)

    # Post-run cleanup of per-tool intermediates. Runs ONLY on a clean run
    # (rc == 0, non-dry-run) — --keep-incomplete already preserves partials on
    # failure. Mode from config `cleanup_intermediates` (default "minimal"),
    # forced to "none" by --keep-all. Runs after manifest/provenance finalize
    # (both are in the keep-set) and before the exit-code / copy-back steps, so
    # it operates on the scratch copy. Best-effort: a failure only warns.
    cleanup_mode = "none" if args.keep_all else \
        str(config_object.get("cleanup_intermediates", "minimal")).strip().lower()
    if rc == 0 and not is_dry_run and cleanup_mode != "none":
        if cleanup_mode not in ("minimal", "maximal"):
            print(F"WARNING: unknown cleanup_intermediates '{cleanup_mode}'; "
                  F"skipping cleanup (use minimal | maximal | none)", file=sys.stderr)
        else:
            try:
                sys.path.insert(0, os.path.join(script_dir, "scripts"))
                from cleanup_outputs import cleanup as _cleanup
                _cleanup(output_dir, cleanup_mode,
                         log=lambda m: print(m, file=sys.stderr))
            except Exception as e:
                print(F"WARNING: post-run cleanup failed: {e}", file=sys.stderr)

    if rc != 0:
        sys.exit(rc)

if __name__ == "__main__":
    main()
