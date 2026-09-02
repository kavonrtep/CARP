#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# k-mer prefilter test for DANTE_LINE — exec-host runner.
#
#   ./RUN.sh              run the CALIBRATION arm (3 genomes x 2 settings)
#   ./RUN.sh --list       show what would run and what is already done
#   ./RUN.sh --score      re-score existing results without running anything
#   ./RUN.sh --validate   run the HELD-OUT genomes (only after calibration passes)
#
# WHAT IS BEING TESTED. dante_line infers how far a LINE extends past its
# protein-domain core by aligning many copies' flanks. A pair is only aligned if
# mmseqs finds a shared exact k-mer in the 30 nt next to the core. At mmseqs'
# default k that window is too short, so most genuinely related pairs are
# discarded unseen and 71% of elements get no extension at all. This measures
# whether an explicit "-k 7" fixes that WITHOUT degrading boundary quality.
#
# The decision rule was fixed in advance -- see PREREGISTERED.md. Do not change
# it after seeing results; that is the whole point of it existing.
#
# SELF-CONTAINED: uses the patched scripts/ in this directory (they carry the
# new --prefilter-kmer flag, which the shipped container does NOT have), the
# container only for its tools. Nothing outside this dir and /nfsroot is needed.
# ---------------------------------------------------------------------------
set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS=/nfsroot/projects/darwin/runs
SIF="${SIF:-$RUNS/tmp3/image/assembly_repeat_annotation_pipeline.sif}"
THREADS="${THREADS:-32}"

# genome tag | carp_output dir | role
JOBS_CAL="boechera|$RUNS/tmp3/GCA_018361405.1/carp_output
GCA_977018505|$RUNS/run-000182/carp_output
GCA_986270105|$RUNS/run-000183/carp_output"

JOBS_VAL="GCA_986264755|$RUNS/run-000181/carp_output
GCA_965641495|$RUNS/run-000185/carp_output
GCA_984574005|$RUNS/run-000178/carp_output"

MODE="${1:-run}"

# ---------------------------------------------------------------- environment
find_python() {
  # An env with python3 + parasail + yaml, and seqkit + mmseqs on the same PATH.
  local runtime
  runtime="$(command -v singularity || command -v apptainer || true)"
  if [ -n "$runtime" ] && [ -f "$SIF" ]; then
    local found
    found=$("$runtime" exec -B /nfsroot:/nfsroot "$SIF" bash -lc '
      for d in /opt/conda/envs/*/bin /usr/local/bin /opt/conda/bin; do
        [ -x "$d/python3" ] || continue
        [ -x "$d/seqkit" ] || continue
        [ -x "$d/mmseqs" ] || continue
        "$d/python3" -c "import parasail, yaml" >/dev/null 2>&1 || continue
        echo "$d"; break
      done' 2>/dev/null | tail -1)
    if [ -n "$found" ]; then echo "container:$runtime:$found"; return 0; fi
  fi
  # fallback: a host conda env with the same tools
  local d
  for d in /envs/conda/envs/*/bin /opt/conda/envs/*/bin; do
    [ -x "$d/python3" ] && [ -x "$d/seqkit" ] && [ -x "$d/mmseqs" ] || continue
    "$d/python3" -c "import parasail, yaml" >/dev/null 2>&1 || continue
    echo "host::$d"; return 0
  done
  return 1
}

ENVSPEC="$(find_python)" || {
  echo "FATAL: found no environment with python3+parasail+yaml AND seqkit AND mmseqs." >&2
  echo "       Looked inside $SIF and in /envs/conda/envs, /opt/conda/envs." >&2
  echo "       Set SIF=/path/to/image.sif, or put such an env on PATH." >&2
  exit 1
}
KIND="${ENVSPEC%%:*}"; REST="${ENVSPEC#*:}"; RUNTIME="${REST%%:*}"; BIN="${REST#*:}"

run_dante_line() {   # $1=outdir $2=carp_output $3=k
  local out="$1" D="$2" k="$3"
  local args=(
    "$BIN/python3" "$HERE/scripts/dante_line.py"
    -g "$D/genome_cleaned.fasta" -a "$D/DANTE/DANTE_filtered.gff3"
    -o "$out" -t "$THREADS"
    --max-group-size 500 --support-fraction 0.5 --min-group-alignments 5
    --max-extension 0 --max-element-length 8000 --library-source core
    --prefilter-identity 0.8 --prefilter-kmer "$k"
    --mask-gff3 "$D/TideCluster/default/TideCluster_tidehunter.gff3"
    --mask-gff3 "$D/DANTE_LTR/DANTE_LTR.gff3"
    --mask-gff3 "$D/DANTE_TIR/DANTE_TIR_final.gff3"
  )
  if [ "$KIND" = container ]; then
    "$RUNTIME" exec -B /nfsroot:/nfsroot -B "$HERE:$HERE" \
      --env "PATH=$BIN:/usr/bin:/bin" --env "CARP_VOCABULARY=$HERE/classification_vocabulary.yaml" \
      "$SIF" "${args[@]}"
  else
    PATH="$BIN:$PATH" CARP_VOCABULARY="$HERE/classification_vocabulary.yaml" "${args[@]}"
  fi
}

do_jobs() {   # $1 = job list
  local pids=()
  while IFS='|' read -r tag D; do
    [ -n "$tag" ] || continue
    for k in 0 7; do
      local out="$HERE/results/${tag}_k${k}"
      if [ -f "$HERE/results/DONE_${tag}_k${k}" ]; then
        echo "  skip   $tag k=$k (already done)"; continue
      fi
      mkdir -p "$out"
      echo "  launch $tag k=$k -> results/${tag}_k${k}.log"
      ( run_dante_line "$out" "$D" "$k" > "$HERE/results/${tag}_k${k}.log" 2>&1
        echo $? > "$HERE/results/DONE_${tag}_k${k}" ) &
      pids+=($!)
    done
  done <<< "$1"
  echo "  waiting for ${#pids[@]} job(s)..."
  wait "${pids[@]}" 2>/dev/null
}

score() {     # $1 = job list, $2 = label
  echo
  echo "==================== $2 ===================="
  while IFS='|' read -r tag D; do
    [ -n "$tag" ] || continue
    local gt="$HERE/gt/gt_${tag}.tsv"
    local a="$HERE/results/${tag}_k0/DANTE_LINE.gff3"
    local b="$HERE/results/${tag}_k7/DANTE_LINE.gff3"
    [ -f "$gt" ] && [ -f "$a" ] && [ -f "$b" ] || { echo "### $tag: incomplete, skipping"; continue; }
    "$BIN/python3" "$HERE/score_extensions.py" --truth "$gt" --label "$tag" \
        --gff "k-default=$a" "k-7=$b"
  done <<< "$1"
  echo
  echo "PRE-REGISTERED RULE (PREREGISTERED.md): adopt k=7 only if, on EVERY"
  echo "genome above, coverage rises by >= 10 percentage points AND purity stays"
  echo ">= 95% AND purity falls by no more than 1 point. Judge the WORST genome,"
  echo "never the average. Any failure -> keep the mmseqs default."
}

case "$MODE" in
  --list)
    echo "environment: $KIND ($BIN)"
    echo "threads per job: $THREADS"
    echo; echo "CALIBRATION:"; echo "$JOBS_CAL" | sed 's/|/  ->  /'
    echo; echo "VALIDATION (held out):"; echo "$JOBS_VAL" | sed 's/|/  ->  /'
    echo; echo "already done:"; ls "$HERE/results" 2>/dev/null | grep '^DONE_' || echo "  (none)"
    ;;
  --score)   score "$JOBS_CAL" "CALIBRATION" ;;
  --validate)
    echo "environment: $KIND ($BIN); threads=$THREADS"
    do_jobs "$JOBS_VAL"; score "$JOBS_VAL" "VALIDATION (held out)" ;;
  *)
    echo "environment: $KIND ($BIN); threads=$THREADS"
    echo "running CALIBRATION (3 genomes x 2 settings)"
    do_jobs "$JOBS_CAL"; score "$JOBS_CAL" "CALIBRATION" ;;
esac
