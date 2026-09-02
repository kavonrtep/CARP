#!/bin/bash
set -u
W=/home/petr/PycharmProjects/assembly_repeat_annotation_pipeline/tmp/rm_ab_2026-09-02
RMBIN=/envs/conda/envs/280b2676dc55ba6e3e7da3b413edfc00_/bin
lib=$1; tag=$2
mkdir -p "$W/out_$tag"
export PATH="$RMBIN:$PATH"
timeout 20000 "$RMBIN/RepeatMasker" -lib "$W/$lib" -pa 7 -qq -nolow -no_is \
  -dir "$W/out_$tag" "$W/genome.fasta" > "$W/$tag.log" 2>&1
echo $? > "$W/DONE_$tag"
