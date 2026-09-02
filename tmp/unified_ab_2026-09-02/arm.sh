#!/bin/bash
# Run make_unified_annotation twice per genome: once with the baseline's own
# DANTE_LINE.gff3 (rc1 bounds) and once with the gated one. Every other input is
# the baseline's, unchanged. Valid because the LINE library is byte-identical in
# content between the two (library_source: core), so the RepeatMasker layer --
# and everything else feeding this step -- is unaffected.
set -u
REPO=/home/petr/PycharmProjects/assembly_repeat_annotation_pipeline
W=$REPO/tmp/unified_ab_2026-09-02
RS=/envs/conda/envs/280b2676dc55ba6e3e7da3b413edfc00_/bin/Rscript
tag=$1; D=$2; line_gff=$3; arm=$4
opt(){ [ -r "$2" ] && printf ' %s %s' "$1" "$2"; }
ARGS="--ltr $D/DANTE_LTR/DANTE_LTR.gff3"
ARGS="$ARGS$(opt --ltr_tandems $D/DANTE_LTR/DANTE_LTR_tandems.gff3)"
ARGS="$ARGS --tir $D/DANTE_TIR/DANTE_TIR_combined.gff3 --line $line_gff"
ARGS="$ARGS --dante $D/DANTE/DANTE_filtered.gff3"
ARGS="$ARGS --tc_default $D/TideCluster/default/TideCluster_clustering.gff3"
ARGS="$ARGS --tc_short $D/TideCluster/short_monomer/TideCluster_clustering.gff3"
ARGS="$ARGS --tc_rm $D/TideCluster/default/RM_on_TideCluster_Library.gff3"
ARGS="$ARGS$(opt --tc_rdna_default $D/TideCluster/default/TideCluster_rdna.tsv)"
ARGS="$ARGS$(opt --tc_rdna_short $D/TideCluster/short_monomer/TideCluster_rdna.tsv)"
ARGS="$ARGS$(opt --tc_trc_table_default $D/TideCluster/default/TideCluster_report/trc_table.tsv)"
ARGS="$ARGS$(opt --tc_trc_table_short $D/TideCluster/short_monomer/TideCluster_report/trc_table.tsv)"
ARGS="$ARGS --rm $D/RepeatMasker/RM_on_combined_library_plus_DANTE.gff3"
ARGS="$ARGS --th_default $D/TideCluster/default/TideCluster_tidehunter_short.gff3"
ARGS="$ARGS --th_short $D/TideCluster/short_monomer/TideCluster_tidehunter_short.gff3"
ARGS="$ARGS$(opt --th_raw_default $D/TideCluster/default/TideCluster_tidehunter.gff3)"
ARGS="$ARGS$(opt --th_raw_short $D/TideCluster/short_monomer/TideCluster_tidehunter.gff3)"
ARGS="$ARGS --rm_tc_tandem_gate TRUE --fai $D/genome_cleaned.fasta.fai --threads 4"
mkdir -p "$W/$tag"
export PATH="$REPO/scripts:$PATH"
timeout 20000 "$RS" "$REPO/scripts/make_unified_annotation.R" $ARGS \
  --output "$W/$tag/unified_$arm.gff3" > "$W/$tag/$arm.log" 2>&1
echo $? > "$W/DONE_${tag}_${arm}"
