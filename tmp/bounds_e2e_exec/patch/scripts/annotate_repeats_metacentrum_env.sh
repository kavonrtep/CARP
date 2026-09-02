#!/bin/bash
#PBS -N annotate_repeats
#PBS -l select=1:ncpus=25:mem=256gb:scratch_local=300gb
#PBS -q elixircz@pbs-m1.metacentrum.cz
#PBS -l walltime=172:00:00
#PBS -j oe
#PBS -m bae

## Use approximately 1 cpu per 8-10 GB of RAM
##
## Submit with -v to pass required and optional variables, e.g.:
##   qsub -v SINGULARITY_IMAGE=/path/to/image.sif,GENOME=/path/to/genome.fasta,OUTPUT_DIR=/path/to/output \
##        annotate_repeats_metacentrum_env.sh
##
## Required variables (must be passed via -v):
##   SINGULARITY_IMAGE  - path to the Singularity .sif image
##   GENOME             - path to the genome FASTA file
##
## Optional variables (passed via -v, omit to use defaults):
##   OUTPUT_DIR                - directory where output zip will be copied
##                               (default: same directory as GENOME)
##   CUSTOM_DATABASE_TAREAN    - path to custom TAREAN/tandem repeat FASTA library
##   CUSTOM_DATABASE_REPEATS   - path to custom RepeatMasker repeat FASTA library
##   REPEATMASKER_SENSITIVITY  - default | sensitive | quick  (default: default)
##   REDUCE_LIBRARY            - True | False  (default: True)


# --- Validate required variables ---
missing=()
[ -z "$SINGULARITY_IMAGE" ] && missing+=("SINGULARITY_IMAGE")
[ -z "$GENOME" ]             && missing+=("GENOME")

if [ ${#missing[@]} -gt 0 ]; then
    echo "ERROR: The following required variables are not set (pass them with qsub -v):"
    printf '  %s\n' "${missing[@]}"
    exit 1
fi

# --- Resolve OUTPUT_DIR and zip name ---
GENOME_BASENAME=$(basename "$GENOME" | sed 's/\.[^.]*$//')
OUTPUT_ZIP_NAME="${GENOME_BASENAME}_output.zip"
if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR=$(dirname "$GENOME")
fi

# --- Copy data and Singularity image to scratch ---
SINGULARITY_IMAGE_BASE_NAME=$(basename "$SINGULARITY_IMAGE")

rsync -avt "$SINGULARITY_IMAGE" "$SCRATCHDIR/"
rsync -avt "$GENOME" "${SCRATCHDIR}/genome.fasta"

if [ -n "$CUSTOM_DATABASE_TAREAN" ]; then
    rsync -avt "$CUSTOM_DATABASE_TAREAN" "${SCRATCHDIR}/custom_database_tarean.fasta"
fi

if [ -n "$CUSTOM_DATABASE_REPEATS" ]; then
    rsync -avt "$CUSTOM_DATABASE_REPEATS" "${SCRATCHDIR}/custom_database_repeats.fasta"
fi

cd "$SCRATCHDIR"

# --- Generate config.yaml ---
echo "genome_fasta: genome.fasta"  > "${SCRATCHDIR}/config.yaml"
echo "output_dir: output"         >> "${SCRATCHDIR}/config.yaml"

if [ -n "$CUSTOM_DATABASE_REPEATS" ]; then
    echo "custom_library: custom_database_repeats.fasta" >> "${SCRATCHDIR}/config.yaml"
fi

if [ -n "$CUSTOM_DATABASE_TAREAN" ]; then
    echo "tandem_repeat_library: custom_database_tarean.fasta" >> "${SCRATCHDIR}/config.yaml"
fi

if [ -n "$REPEATMASKER_SENSITIVITY" ]; then
    echo "repeatmasker_sensitivity: $REPEATMASKER_SENSITIVITY" >> "${SCRATCHDIR}/config.yaml"
fi

if [ -n "$REDUCE_LIBRARY" ]; then
    echo "reduce_library: $REDUCE_LIBRARY" >> "${SCRATCHDIR}/config.yaml"
fi

# Read pipeline version from the Singularity image's labels (stamped
# at build time by the %post block in the Singularity recipe). Falls
# back to "unknown" if the image is too old or `singularity inspect`
# fails. Cheap — runs in <1 s.
PIPELINE_VERSION=$(singularity inspect --labels "$SINGULARITY_IMAGE" 2>/dev/null \
    | awk -F': ' '/^Version:/ {print $2; exit}')
PIPELINE_VERSION="${PIPELINE_VERSION:-unknown}"

# --- Write 0INFO.txt ---
INFO_FILE="${SCRATCHDIR}/0INFO.txt"
{
    echo "===== Job configuration ====="
    echo "Date/time          : $(date --iso-8601=seconds)"
    echo "PBS job ID         : ${PBS_JOBID}"
    echo "Hostname           : $(hostname)"
    echo "PIPELINE_VERSION   : ${PIPELINE_VERSION}"
    echo ""
    echo "--- Required parameters ---"
    echo "SINGULARITY_IMAGE  : ${SINGULARITY_IMAGE}"
    echo "GENOME             : ${GENOME}"
    echo "OUTPUT_DIR         : ${OUTPUT_DIR} (zip: ${OUTPUT_ZIP_NAME})"
    echo ""
    echo "--- Optional parameters ---"
    echo "CUSTOM_DATABASE_TAREAN   : ${CUSTOM_DATABASE_TAREAN:-(not set)}"
    echo "CUSTOM_DATABASE_REPEATS  : ${CUSTOM_DATABASE_REPEATS:-(not set)}"
    echo "REPEATMASKER_SENSITIVITY : ${REPEATMASKER_SENSITIVITY:-(not set, using default)}"
    echo "REDUCE_LIBRARY           : ${REDUCE_LIBRARY:-(not set, using default)}"
    echo ""
    echo "--- PBS resource request ---"
    echo "NCPUS              : ${PBS_NCPUS}"
    echo "SCRATCHDIR         : ${SCRATCHDIR}"
    echo ""
    echo "--- Generated config.yaml ---"
    cat "${SCRATCHDIR}/config.yaml"
} > "$INFO_FILE"

# --- Set up temp directory ---
tmpdir="${SCRATCHDIR}/tmp"
mkdir -p "$tmpdir"
export TMPDIR="$tmpdir"

# --- Run the Singularity container ---
singularity run -B "$SCRATCHDIR" --env TMPDIR="$tmpdir" "$SINGULARITY_IMAGE_BASE_NAME" -c config.yaml -t "$PBS_NCPUS"

# --- Collect outputs ---
mkdir -p "${SCRATCHDIR}/output"
cp "${SCRATCHDIR}/genome.fasta"  "${SCRATCHDIR}/output/"
cp "${SCRATCHDIR}/config.yaml"   "${SCRATCHDIR}/output/"
cp "$INFO_FILE"                  "${SCRATCHDIR}/output/"
cp "$0"                          "${SCRATCHDIR}/output/pbs_script.sh"

if [ -n "$CUSTOM_DATABASE_TAREAN" ]; then
    cp "${SCRATCHDIR}/custom_database_tarean.fasta" "${SCRATCHDIR}/output/"
fi
if [ -n "$CUSTOM_DATABASE_REPEATS" ]; then
    cp "${SCRATCHDIR}/custom_database_repeats.fasta" "${SCRATCHDIR}/output/"
fi

cp -r .snakemake "${SCRATCHDIR}/output/"
env > "${SCRATCHDIR}/output/env.sh"

cp "/var/spool/pbs/spool/${PBS_JOBID}.OU" "${SCRATCHDIR}/output/pbs_log.txt"

# --- Package and transfer ---
zip -y -fz -r "${SCRATCHDIR}/${OUTPUT_ZIP_NAME}" output
mkdir -p "$OUTPUT_DIR"
cp "${SCRATCHDIR}/${OUTPUT_ZIP_NAME}" "${OUTPUT_DIR}/"