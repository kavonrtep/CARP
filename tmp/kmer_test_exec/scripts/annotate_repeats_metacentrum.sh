#!/bin/bash
#PBS -N annotate_repeats
#PBS -l select=1:ncpus=50:mem=512gb:scratch_local=300gb
#PBS -q elixircz@pbs-m1.metacentrum.cz
#PBS -l walltime=172:00:00
#PBS -j oe
#PBS -m bae

## Use approximately 1 cpu per 8-10 GB of RAM

# CONFIGURATION - use absolute paths
SINGULARITY_IMAGE="todo.sif"
GENOME="TODO.fasta"
CUSTOM_DATABASE_TAREAN=""  # optional, keep empty if not used
CUSTOM_DATABASE_REPEATS="" # optional, keep empty if not used
OUTPUT_DIR="TODO"
REPEATMASKER_SENSITIVITY="default"  # default, sensitive, quick
REDUCE_LIBRARY="True"  # True, False, if not defined, default is True
# END OF CONFIGURATION



SINGULARITY_IMAGE_BASE_NAME=$(basename $SINGULARITY_IMAGE)
# copy data and singularity image to scratch


rsync -avt $SINGULARITY_IMAGE $SCRATCHDIR/
rsync -avt $GENOME ${SCRATCHDIR}/genome.fasta


# copy is not empty
if [ -n "$CUSTOM_DATABASE_TAREAN" ]; then
     rsync -avt $CUSTOM_DATABASE_TAREAN ${SCRATCHDIR}/custom_database_tarean.fasta

fi

if [ -n "$CUSTOM_DATABASE_REPEATS" ]; then
      rsync -avt $CUSTOM_DATABASE_REPEATS ${SCRATCHDIR}/custom_database_repeats.fasta
fi

cd $SCRATCHDIR

# Generate config.yaml
# format is:
# genome_fasta: genome.fasta
# output_dir: output_dir
# custom_library: path/to/custom_library.fasta
# tandem_repeat_library: path/to/tandem_repeat_library.fasta

echo "genome_fasta: genome.fasta" > ${SCRATCHDIR}/config.yaml
echo "output_dir: output" >> ${SCRATCHDIR}/config.yaml
# if custom database is not empty
if [ -n "$CUSTOM_DATABASE_REPEATS" ]; then
    echo "custom_library: custom_database_repeats.fasta" >> ${SCRATCHDIR}/config.yaml
fi
# if custom database is not empty
if [ -n "$CUSTOM_DATABASE_TAREAN" ]; then
    echo "tandem_repeat_library: custom_database_tarean.fasta" >> ${SCRATCHDIR}/config.yaml
fi

if [ -n "$REPEATMASKER_SENSITIVITY" ]; then
    echo "repeatmasker_sensitivity: $REPEATMASKER_SENSITIVITY" >> ${SCRATCHDIR}/config.yaml
fi

if [ -n "$REDUCE_LIBRARY" ]; then
    echo "reduce_library: $REDUCE_LIBRARY" >> ${SCRATCHDIR}/config.yaml
fi

tmpdir=$SCRATCHDIR/tmp
mkdir -p $tmpdir
export TMPDIR=$tmpdir

# run the singularity container, specify temp directory variable

singularity run -B $SCRATCHDIR --env TMPDIR=$tmpdir $SINGULARITY_IMAGE_BASE_NAME  -c config.yaml -t $PBS_NCPUS


cp $SCRATCHDIR/genome.fasta $SCRATCHDIR/output
cp $SCRATCHDIR/config.yaml $SCRATCHDIR/output
if [ -n "$CUSTOM_DATABASE_TAREAN" ]; then
  cp $SCRATCHDIR/custom_database_tarean.fasta $SCRATCHDIR/output
fi
if [ -n "$CUSTOM_DATABASE_REPEATS" ]; then
  cp $SCRATCHDIR/custom_database_repeats.fasta $SCRATCHDIR/output
fi
cp $0 $SCRATCHDIR/output/pbs_script.sh
cp -r .snakemake $SCRATCHDIR/output
env > $SCRATCHDIR/output/env.sh

cp /var/spool/pbs/spool/${PBS_JOBID}.OU re_output/pbs_log.txt

zip -y -fz -r $SCRATCHDIR/output.zip output
mkdir -p $OUTPUT_DIR
cp $SCRATCHDIR/output.zip $OUTPUT_DIR/



