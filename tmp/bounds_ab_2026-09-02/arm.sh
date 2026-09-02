#!/bin/bash
set -u
REPO=/home/petr/PycharmProjects/assembly_repeat_annotation_pipeline
W=$REPO/tmp/bounds_ab_2026-09-02
E=/envs/conda/envs/582c0adf8a6de987318ca6c3ed895a14_/bin
tag=$1; D=$2; vocab=$3
out="$W/${tag}_${vocab}"
mkdir -p "$out"
export PATH="$E:$PATH"
export CARP_VOCABULARY="$W/vocab_${vocab}.yaml"
timeout 20000 "$E/python3" "$REPO/scripts/dante_line.py" \
  -g "$D/genome_cleaned.fasta" -a "$D/DANTE/DANTE_filtered.gff3" \
  -o "$out" -t 7 --max-group-size 500 --support-fraction 0.5 \
  --min-group-alignments 5 --max-extension 0 --max-element-length 8000 \
  --library-source core \
  --mask-gff3 "$D/TideCluster/default/TideCluster_tidehunter.gff3" \
  --mask-gff3 "$D/DANTE_LTR/DANTE_LTR.gff3" \
  --mask-gff3 "$D/DANTE_TIR/DANTE_TIR_final.gff3" > "$W/${tag}_${vocab}.log" 2>&1
echo $? > "$W/DONE_${tag}_${vocab}"
