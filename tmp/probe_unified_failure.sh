#!/bin/bash
# ---------------------------------------------------------------------------
# Diagnose the make_unified_annotation failure seen in the 1.7.0a1 build.
#
# Same code + same data + same batch count PASSES in the sandbox and FAILS in
# the container, so the difference is the container's R/Bioconductor stack or
# the code change — this settles which, and captures a real traceback either way.
#
# Three arms, all inside the container, all on the same inputs:
#   1. VERSIONS  — what Bioconductor generation the container actually has
#   2. OLD       — scripts/ as of f9e603a (the 1.6.2 code, before this branch)
#   3. NEW       — scripts/ from the current working tree (has traceback capture)
#
# NOTE on the first version of this script: its A/B was unsound. It varied the
# CODE between arms but both arms read the SAME inputs — inputs that the NEW code
# had produced. So "old code fails too" did NOT exonerate the change; the defect
# was in a GFF3 file, not in make_unified_annotation.R. Arm 0 below now checks the
# data first, which is where that failure actually lived.
#
# Verdict:
#   arm 0 flags a comma-valued attribute -> the upstream writer is at fault
#   OLD and NEW both clean -> environment; pin the Bioconductor stack
#   OLD passes, NEW fails  -> the code change; arm 3's traceback names the frame
#
# Usage:
#   bash tmp/probe_unified_failure.sh [OUTPUT_DIR] [THREADS]
#
# OUTPUT_DIR is a finished carp_output/ from a FAILED run (it needs everything
# except Repeat_Annotation_Unified.gff3). Default is the smaller of the two.
# ---------------------------------------------------------------------------
set -uo pipefail

REPO="${REPO:-/home/petr/PycharmProjects/assembly_repeat_annotation_pipeline}"
SIF="${SIF:-$REPO/images/assembly_repeat_annotation_pipeline.sif}"
OUT="${1:-/nfsroot/projects/darwin/runs/tmp/GCA_973357735.1/carp_output}"
THREADS="${2:-32}"
OLD_REF="${OLD_REF:-f9e603a}"
WORK="${WORK:-/tmp/carp_probe.$$}"

RUNTIME="$(command -v singularity || command -v apptainer || true)"
[ -n "$RUNTIME" ] || { echo "FATAL: neither singularity nor apptainer on PATH"; exit 1; }
[ -f "$SIF" ]     || { echo "FATAL: no image at $SIF"; exit 1; }
[ -d "$OUT" ]     || { echo "FATAL: no output dir at $OUT"; exit 1; }

mkdir -p "$WORK"/{old,new,logs}
echo "repo    : $REPO"
echo "image   : $SIF"
echo "data    : $OUT"
echo "threads : $THREADS   (batch count depends on this; the failing runs used the job's -t)"
echo "workdir : $WORK"
echo

# --- stage the two scripts/ trees -------------------------------------------
# make_unified_annotation.R sources classification.R from its OWN directory, so
# each arm needs a complete scripts/ tree, not just the one file.
echo "==> staging scripts/ trees"
git -C "$REPO" archive "$OLD_REF" scripts | tar -x -C "$WORK/old" || {
  echo "FATAL: could not export scripts/ at $OLD_REF"; exit 1; }
cp -a "$REPO/scripts/." "$WORK/new/scripts/" 2>/dev/null || {
  mkdir -p "$WORK/new/scripts" && cp -a "$REPO/scripts/." "$WORK/new/scripts/"; }
# the old tree also gets the old vocabulary, so the old parser never sees a key
# it predates (the new one adds max_consensus_length)
git -C "$REPO" show "$OLD_REF:classification_vocabulary.yaml" \
  > "$WORK/old/classification_vocabulary.yaml" 2>/dev/null || true
echo "    old: $(ls "$WORK/old/scripts" | wc -l) files from $OLD_REF"
echo "    new: $(ls "$WORK/new/scripts" | wc -l) files from the working tree"
echo

# --- arm 0: is the DATA well-formed? ----------------------------------------
# A comma is GFF3's multi-value separator, so `Key=a,b` imports as a
# CompressedCharacterList. Subsetting a DataFrame holding one dies on the
# Bioconductor 3.14 stack this container pins. Check the tier-1 inputs first —
# it is far cheaper than running the rule, and it is where the 1.7.0a1 failure was.
echo "############ ARM 0 — comma-valued GFF3 attributes in the inputs ############"
arm0_bad=0
for f in "$OUT/DANTE_LINE/DANTE_LINE.gff3" \
         "$OUT/DANTE_TIR/DANTE_TIR_combined.gff3" \
         "$OUT/DANTE_LTR/DANTE_LTR.gff3"; do
  [ -f "$f" ] || continue
  n=$(awk -F'\t' '!/^#/ && NF>=9 {print $9}' "$f" | grep -c '=[^;]*,[^;]*')
  printf "  %-34s %s feature(s) with a comma inside an attribute value\n" \
         "$(basename "$f")" "$n"
  if [ "$n" -gt 0 ]; then
    arm0_bad=1
    awk -F'\t' '!/^#/ && NF>=9 {print $9}' "$f" | grep -o '[A-Za-z_]*=[^;]*,[^;]*' \
      | sed 's/=.*//' | sort -u | sed 's/^/      offending key: /'
  fi
done
[ "$arm0_bad" -eq 0 ] && echo "  clean" || \
  echo "  ^ this alone can fail every batch on this container's S4Vectors."
echo

BINDS="-B $OUT:$OUT -B $WORK:$WORK"

# --- locate the R that the rule actually uses -------------------------------
echo "==> locating the tidecluster env inside the image"
RSCRIPT=$("$RUNTIME" exec $BINDS "$SIF" bash -lc '
  for d in /opt/conda/envs/*/bin; do
    if [ -x "$d/Rscript" ] && "$d/Rscript" -e "q(status=!requireNamespace(\"rtracklayer\",quietly=TRUE))" >/dev/null 2>&1; then
      echo "$d/Rscript"; break
    fi
  done')
[ -n "$RSCRIPT" ] || { echo "FATAL: no Rscript with rtracklayer inside the image"; exit 1; }
echo "    $RSCRIPT"
echo

# --- arm 1: what stack is in there? -----------------------------------------
echo "############ ARM 1 — Bioconductor versions in the container ############"
"$RUNTIME" exec $BINDS "$SIF" "$RSCRIPT" -e '
for (p in c("S4Vectors","IRanges","GenomicRanges","XVector","rtracklayer",
            "BiocGenerics","GenomeInfoDb")) {
  v <- tryCatch(as.character(packageVersion(p)), error=function(e) "ABSENT")
  cat(sprintf("  %-14s %s\n", p, v))
}
cat("  R              ", R.version.string, "\n")' 2>&1 | tee "$WORK/logs/versions.txt"
echo
echo "  (sandbox for comparison: S4Vectors 0.40.2  IRanges 2.36.0"
echo "                           GenomicRanges 1.54.1  XVector 0.42.0  rtracklayer 1.62.0)"
echo

# --- build the argument list ------------------------------------------------
opt_arg() { [ -f "$2" ] && printf ' %s %s' "$1" "$2"; }
ARGS="--ltr $OUT/DANTE_LTR/DANTE_LTR.gff3"
ARGS="$ARGS$(opt_arg --ltr_tandems $OUT/DANTE_LTR/DANTE_LTR_tandems.gff3)"
ARGS="$ARGS --tir $OUT/DANTE_TIR/DANTE_TIR_combined.gff3"
ARGS="$ARGS --line $OUT/DANTE_LINE/DANTE_LINE.gff3"
ARGS="$ARGS --dante $OUT/DANTE/DANTE_filtered.gff3"
ARGS="$ARGS --tc_default $OUT/TideCluster/default/TideCluster_clustering.gff3"
ARGS="$ARGS --tc_short $OUT/TideCluster/short_monomer/TideCluster_clustering.gff3"
ARGS="$ARGS --tc_rm $OUT/TideCluster/default/RM_on_TideCluster_Library.gff3"
ARGS="$ARGS$(opt_arg --tc_rdna_default $OUT/TideCluster/default/TideCluster_rdna.tsv)"
ARGS="$ARGS$(opt_arg --tc_rdna_short   $OUT/TideCluster/short_monomer/TideCluster_rdna.tsv)"
ARGS="$ARGS$(opt_arg --tc_trc_table_default $OUT/TideCluster/default/TideCluster_report/trc_table.tsv)"
ARGS="$ARGS$(opt_arg --tc_trc_table_short   $OUT/TideCluster/short_monomer/TideCluster_report/trc_table.tsv)"
ARGS="$ARGS --rm $OUT/RepeatMasker/RM_on_combined_library_plus_DANTE.gff3"
ARGS="$ARGS --th_default $OUT/TideCluster/default/TideCluster_tidehunter_short.gff3"
ARGS="$ARGS --th_short $OUT/TideCluster/short_monomer/TideCluster_tidehunter_short.gff3"
ARGS="$ARGS --th_raw_default $OUT/TideCluster/default/TideCluster_tidehunter.gff3"
ARGS="$ARGS --th_raw_short $OUT/TideCluster/short_monomer/TideCluster_tidehunter.gff3"
ARGS="$ARGS --rm_tc_tandem_gate TRUE"
ARGS="$ARGS --fai $OUT/genome_cleaned.fasta.fai"
ARGS="$ARGS --threads $THREADS"

run_arm () {                       # $1 = label, $2 = staged tree root
  local label="$1" root="$2" log="$WORK/logs/$1.log"
  echo "############ ARM — $label ############"
  local vocab_bind=""
  [ -s "$root/classification_vocabulary.yaml" ] && \
    vocab_bind="-B $root/classification_vocabulary.yaml:/opt/pipeline/classification_vocabulary.yaml"
  "$RUNTIME" exec $BINDS -B "$root/scripts:/opt/pipeline/scripts" $vocab_bind "$SIF" \
      "$RSCRIPT" /opt/pipeline/scripts/make_unified_annotation.R \
      $ARGS --output "$WORK/$label.gff3" > "$log" 2>&1
  local rc=$?
  if [ $rc -eq 0 ]; then
    echo "  RESULT: PASS  (rc=0)"
  else
    echo "  RESULT: FAIL  (rc=$rc)"
    echo "  --- first failure + stack ---"
    grep -m1 -A22 'call stack' "$log" 2>/dev/null || \
      grep -m3 -E 'failed:|Error' "$log" | head -8
  fi
  echo "  full log: $log"
  echo
  return $rc
}

run_arm old "$WORK/old"; OLD_RC=$?
run_arm new "$WORK/new"; NEW_RC=$?

# --- verdict ----------------------------------------------------------------
echo "############ VERDICT ############"
if [ $OLD_RC -ne 0 ] && [ $NEW_RC -ne 0 ]; then
  echo "  BOTH FAIL -> the container's R/Bioconductor stack, not the code change."
  echo "  envs/tidecluster.yaml pins none of rtracklayer / S4Vectors /"
  echo "  GenomicRanges / XVector; a fresh solve took a different generation."
  echo "  Fix: pin them in envs/tidecluster.yaml and rebuild."
  echo "  Compare $WORK/logs/versions.txt against the sandbox versions printed above."
elif [ $OLD_RC -eq 0 ] && [ $NEW_RC -ne 0 ]; then
  echo "  OLD PASSES, NEW FAILS -> the code change on this branch."
  echo "  The stack above names the frame; send me $WORK/logs/new.log."
elif [ $OLD_RC -ne 0 ] && [ $NEW_RC -eq 0 ]; then
  echo "  OLD FAILS, NEW PASSES -> unexpected; the branch fixed something."
  echo "  Send me both logs."
else
  echo "  BOTH PASS -> not reproducible this way. The failing rule differs from"
  echo "  this invocation: check the real rule's --threads (batch count depends"
  echo "  on it) and --max_workers / --mem_budget_gb in the run's snakemake log,"
  echo "  then re-run with the same THREADS argument."
fi
echo
echo "  logs + outputs under: $WORK"
