Bootstrap: docker
# PINNED base image (do not use the floating :latest). The unpinned
# continuumio/miniconda3 drifted repeatedly and broke the release build:
# first a newer conda started enforcing Anaconda's channel Terms of Service
# (CondaToSNonInteractiveError), then the base python jumped to 3.14 with an
# Anaconda helper package (anaconda-channel-guide) that hard-pinned it and
# blocked the downgrade to 3.11. 24.9.2-0 is a pre-drift conda 24.9.2 build
# (the same conda version this project's sandbox runs the pipeline on, so it is
# known-compatible). Bump this tag deliberately, and re-test the SIF build, when
# you actually want a newer base — never let it float.
From: continuumio/miniconda3:24.9.2-0

%post
    apt-get update
    # jq is a BUILD-time requirement of the bioconda bioconductor-*data post-link
    # scripts (installBiocDataPackage.sh -> yq -> jq). It must be on PATH at the
    # system level *before* the conda env creation below, because conda does not
    # guarantee jq links before genomeinfodbdata's post-link runs (there is no
    # dependency between them), so pinning jq inside the env yaml is not enough.
    # Installing it here makes it available to every post-link regardless of
    # conda's link order. Without it the SIF build fails: "yq: Error starting jq".
    apt-get install -y libxml2 libxml2-dev jq
    mkdir -p /opt/conda/config
    export CONDARC=/opt/conda/config/.condarc

    # Build ONLY from conda-forge + bioconda, never Anaconda's default channels
    # (pkgs/main, pkgs/r). This (1) keeps CARP off Anaconda's commercially-
    # licensed channels, and (2) is defence-in-depth against the channel Terms-
    # of-Service gate: newer conda refuses to install from the default channels
    # until their ToS is accepted (CondaToSNonInteractiveError), and some conda
    # builds enforce that but don't even ship the `conda tos` CLI to accept it.
    # The pinned 24.9.2-0 base predates ToS enforcement, but keeping defaults out
    # of the config means a future pin bump can't silently reintroduce the gate.
    # --override-channels covers the bootstrap installs below; stripping
    # `defaults` from BOTH the system condarc and $CONDARC covers the later
    # `snakemake --use-conda` env creation (which reads the config, not
    # --override-channels).
    conda config --system --remove channels defaults 2>/dev/null || true
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda config --remove channels defaults 2>/dev/null || true
    conda config --set channel_priority strict

    # Do NOT downgrade the base env's python. The pinned base ships python 3.12,
    # and snakemake 8.12.0 is a noarch package (build pyhdfd78af_0 — runs on
    # 3.11/3.12/3.13), so it installs into the base env as-is. Changing the base
    # python (the old `conda install python=3.11`) rewrites the base env and
    # breaks the base `conda` package, which is built for the base python; the
    # next conda operation then fails plugin discovery with
    # "PluginError: Conflicting post_solves plugins: signature-verification".
    #
    # anaconda-anon-usage (a base telemetry plugin) corrupts 'conda info --json',
    # which Snakemake 8.12+ parses, so remove it first; also drop
    # anaconda-channel-guide if present (both guarded, separate commands so a
    # missing package is a harmless no-op).
    conda remove --force -y anaconda-anon-usage    2>/dev/null || true
    conda remove --force -y anaconda-channel-guide 2>/dev/null || true

    # Install Snakemake into the base env (conda-forge + bioconda only; the base
    # python 3.12 is retained — snakemake 8.12.0 is noarch).
    conda install -y --override-channels -c conda-forge -c bioconda snakemake=8.12.0
    conda init bash

    # Source the conda.sh script to ensure conda commands are available
    echo ". /opt/conda/etc/profile.d/conda.sh" >> /etc/profile
    echo "conda activate base" >>  /etc/profile

    # Create environments using snakemake
    . /opt/conda/etc/profile.d/conda.sh
    conda activate base
    # Verify that strict channel priority is set
    conda config --show | grep channel_priority
    # (anaconda-anon-usage / anaconda-channel-guide were removed above, before
    #  the python install, so 'conda info --json' is clean for Snakemake here.)

    cd /opt/pipeline
    # make dummy data so snakemake can create the environments
    mkdir -p /opt/pipeline/data
    touch /opt/pipeline/data/CEN6_ver_220406_part.fasta
    touch /opt/pipeline/data/pisum_custom_library.fasta

    touch /opt/pipeline/data/FabTR_all_sequences_210901.db.RM_format.fasta

    snakemake --use-conda --conda-prefix /opt/conda/envs --conda-create-envs-only --cores 4 --configfile /opt/pipeline/config.yaml --conda-frontend conda
    # Clean up
    conda clean --all

    # Download Plotly.js for standalone HTML reports (no CDN needed at view time)
    wget -q -O /opt/pipeline/data/plotly.min.js \
        https://cdn.plot.ly/plotly-2.35.2.min.js

    # make root accessible for everyone
    chmod -R 777 /root

    # Stamp container labels from version.py (single source of truth)
    # so 'singularity inspect --labels <sif>' returns the version that
    # matched the source tree at build time. The %labels block above
    # only provides "unknown" fallbacks; this %post step overwrites
    # them in /.singularity.d/labels.json before image is finalised.
    PIPELINE_VERSION=$(python3 /opt/pipeline/version.py)
    BUILD_DATE=$(date -u +%Y-%m-%dT%H:%M:%SZ)
    mkdir -p /.singularity.d
    cat > /.singularity.d/labels.json <<EOF
{
    "Version": "${PIPELINE_VERSION}",
    "Build-Date": "${BUILD_DATE}"
}
EOF
    echo "Stamped container labels: Version=${PIPELINE_VERSION} Build-Date=${BUILD_DATE}"
    # remove all temp files

%files
    envs /opt/pipeline/envs
    Snakefile /opt/pipeline/Snakefile
    config.yaml /opt/pipeline/config.yaml
    classification_vocabulary.yaml /opt/pipeline/classification_vocabulary.yaml
    data/rdna_library.fasta /opt/pipeline/data/rdna_library.fasta
    run_pipeline.py /opt/pipeline/run_pipeline.py
    version.py /opt/pipeline/version.py
    scripts /opt/pipeline/scripts


%labels
    # Stamped from version.py + build environment in %post below.
    # Read at run time with: singularity inspect --labels <image.sif>
    # Hardcoded fallbacks here are only used if the %post block fails;
    # they should never appear in a successful build.
    Version unknown
    Build-Date unknown


%environment
    export PATH=/opt/pipeline/scripts:/opt/conda/bin:$PATH
    export CONDA_ENVS_PATH=/opt/conda/envs
    export CONDA_PREFIX=/opt/conda
    export CONDARC=/opt/conda/config/.condarc
    export HOME=/root


%runscript
    # Navigate to the pipeline directory
    # set cache directory

    /opt/pipeline/run_pipeline.py "$@"