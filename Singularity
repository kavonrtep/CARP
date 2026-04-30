Bootstrap: docker
From: continuumio/miniconda3

%post
    apt-get update
    apt-get install -y libxml2 libxml2-dev
    mkdir -p /opt/conda/config
    export CONDARC=/opt/conda/config/.condarc

    conda install python=3.11

    # Install Snakemake
    conda install -c bioconda -c conda-forge snakemake=8.12.0
    # mamba install h5py
    # configure strict channel priority
    conda config --set channel_priority strict
    conda init bash

    # Source the conda.sh script to ensure conda commands are available
    echo ". /opt/conda/etc/profile.d/conda.sh" >> /etc/profile
    echo "conda activate base" >>  /etc/profile

    # Create environments using snakemake
    . /opt/conda/etc/profile.d/conda.sh
    conda activate base
    # Verify that strict channel priority is set
    conda config --show | grep channel_priority

    # anaconda-anon-usage (Anaconda telemetry plugin) is incompatible with
    # conda 26.x: it prints "Error loading anaconda-anon-usage: ..." to STDOUT
    # before the JSON, corrupting 'conda info --json' output that Snakemake
    # 8.12+ parses. Remove it so conda info returns clean JSON.
    conda remove anaconda-anon-usage --force -y || true

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