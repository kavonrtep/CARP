#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# k-mer prefilter test for DANTE_LINE — exec-host runner.
#
#   ./RUN.sh --check      verify environment + inputs, run nothing   <-- DO THIS FIRST
#   ./RUN.sh --setup      build the conda env if none was found
#   ./RUN.sh              run the CALIBRATION arm and score it
#   ./RUN.sh --list       what would run, what is already done
#   ./RUN.sh --score      re-score existing results, run nothing
#   ./RUN.sh --validate   the HELD-OUT genomes -- only after calibration passes
#
# WHAT IS BEING TESTED. dante_line infers how far a LINE extends past its
# protein-domain core by aligning many copies' flanks. A pair is only aligned if
# mmseqs finds a shared exact k-mer in the 30 nt next to the core. At mmseqs'
# default k that window is too short, so most genuinely related pairs are
# discarded unseen and 71% of elements get no extension. This measures whether
# an explicit "-k 7" fixes that WITHOUT degrading boundary quality.
#
# The decision rule was fixed in advance -- see PREREGISTERED.md. Do not change
# it after seeing results; that is the whole point of it existing.
#
# SELF-CONTAINED: the shipped container does NOT have --prefilter-kmer, so this
# runs the patched scripts/ in this directory and uses an environment only for
# its tools (python3+parasail+yaml, seqkit, mmseqs).
# ---------------------------------------------------------------------------
set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS=/nfsroot/projects/darwin/runs
SIF="${SIF:-$RUNS/tmp3/image/assembly_repeat_annotation_pipeline.sif}"
THREADS="${THREADS:-32}"
ENV_PREFIX="${ENV_PREFIX:-$HERE/conda_env}"

JOBS_CAL="boechera|$RUNS/tmp3/GCA_018361405.1/carp_output
GCA_977018505|$RUNS/run-000182/carp_output
GCA_986270105|$RUNS/run-000183/carp_output"

JOBS_VAL="GCA_986264755|$RUNS/run-000181/carp_output
GCA_965641495|$RUNS/run-000185/carp_output
GCA_984574005|$RUNS/run-000178/carp_output"

MODE="${1:-run}"

# ------------------------------------------------------- scoring interpreter
# score_extensions.py is pure stdlib, so ANY python3 runs it. It must NOT use
# the tool env's python, which may live inside a container and not exist on the
# host -- that bug made an earlier version die at the scoring step.
SCORE_PY="$(command -v python3 || true)"

# --------------------------------------------------------------- environment
# Returns "KIND:RUNTIME:BIN". KIND is container|host. For host, RUNTIME is empty.
probe_bin() {   # $1 = bin dir on THIS host
  [ -x "$1/python3" ] && [ -x "$1/seqkit" ] && [ -x "$1/mmseqs" ] || return 1
  "$1/python3" -c "import parasail, yaml" >/dev/null 2>&1 || return 1
  return 0
}

find_env() {
  # 1. an env this script built earlier
  if probe_bin "$ENV_PREFIX/bin"; then echo "host::$ENV_PREFIX/bin"; return 0; fi
  # 2. inside the container
  local runtime; runtime="$(command -v singularity || command -v apptainer || true)"
  if [ -n "$runtime" ] && [ -f "$SIF" ]; then
    local found
    found=$("$runtime" exec -B /nfsroot:/nfsroot "$SIF" bash -lc '
      for d in /opt/conda/envs/*/bin /opt/conda/bin /usr/local/bin; do
        [ -x "$d/python3" ] && [ -x "$d/seqkit" ] && [ -x "$d/mmseqs" ] || continue
        "$d/python3" -c "import parasail, yaml" >/dev/null 2>&1 || continue
        echo "$d"; break
      done' 2>/dev/null | tail -1)
    if [ -n "$found" ]; then echo "container:$runtime:$found"; return 0; fi
  fi
  # 3. any env already on this host
  local d
  for d in /envs/conda/envs/*/bin /opt/conda/envs/*/bin "$HOME"/miniconda3/envs/*/bin \
           "$HOME"/mambaforge/envs/*/bin "$HOME"/miniforge3/envs/*/bin; do
    probe_bin "$d" && { echo "host::$d"; return 0; }
  done
  return 1
}

setup_env() {
  local conda; conda="$(command -v mamba || command -v micromamba || command -v conda || true)"
  [ -n "$conda" ] || { echo "FATAL: no conda/mamba/micromamba on PATH; cannot build the env." >&2; return 1; }
  echo "building env at $ENV_PREFIX using $conda (a few minutes)..."
  "$conda" env create --yes -p "$ENV_PREFIX" -f "$HERE/dante_line_env.yaml" 2>&1 | tail -20 \
    || "$conda" create --yes -p "$ENV_PREFIX" -c conda-forge -c bioconda \
         seqkit parasail-python mmseqs2 pyyaml 2>&1 | tail -20
  probe_bin "$ENV_PREFIX/bin" || { echo "FATAL: env built but still missing tools." >&2; return 1; }
  echo "env ready: $ENV_PREFIX/bin"
}

resolve_env() {
  ENVSPEC="$(find_env)" || return 1
  KIND="${ENVSPEC%%:*}"; local rest="${ENVSPEC#*:}"
  RUNTIME="${rest%%:*}"; BIN="${rest#*:}"
  return 0
}


selftest_scripts() {
  # The environment can be perfect and the run still die instantly because a
  # module the scripts import is missing. Prove dante_line.py LOADS, in the
  # environment that will run it, before committing hours to the jobs.
  local out
  if [ "$KIND" = container ]; then
    out=$("$RUNTIME" exec -B /nfsroot:/nfsroot -B "$HERE:$HERE" \
          --env "PATH=$BIN:/usr/bin:/bin" \
          --env "CARP_VOCABULARY=$HERE/classification_vocabulary.yaml" \
          "$SIF" "$BIN/python3" "$HERE/scripts/dante_line.py" --help 2>&1)
  else
    out=$(PATH="$BIN:$PATH" CARP_VOCABULARY="$HERE/classification_vocabulary.yaml" \
          "$BIN/python3" "$HERE/scripts/dante_line.py" --help 2>&1)
  fi
  if [ $? -ne 0 ]; then
    echo "  SCRIPTS FAILED TO LOAD:"; echo "$out" | tail -6 | sed 's/^/    /'; return 1
  fi
  case "$out" in *--prefilter-kmer*) echo "  scripts load, and carry --prefilter-kmer";;
                 *) echo "  SCRIPTS LOAD BUT LACK --prefilter-kmer (stale copy?)"; return 1;; esac
  return 0
}

# ---------------------------------------------------------------- the runner
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
      --env "PATH=$BIN:/usr/bin:/bin" \
      --env "CARP_VOCABULARY=$HERE/classification_vocabulary.yaml" \
      "$SIF" "${args[@]}"
  else
    PATH="$BIN:$PATH" CARP_VOCABULARY="$HERE/classification_vocabulary.yaml" "${args[@]}"
  fi
}

check_inputs() {
  local bad=0
  while IFS='|' read -r tag D; do
    [ -n "$tag" ] || continue
    for f in genome_cleaned.fasta genome_cleaned.fasta.fai DANTE/DANTE_filtered.gff3 \
             TideCluster/default/TideCluster_tidehunter.gff3 DANTE_LTR/DANTE_LTR.gff3 \
             DANTE_TIR/DANTE_TIR_final.gff3; do
      [ -r "$D/$f" ] || { echo "  MISSING  $tag: $D/$f"; bad=1; }
    done
    [ -r "$HERE/gt/gt_${tag}.tsv" ] || { echo "  MISSING  $tag: ground truth gt_${tag}.tsv"; bad=1; }
  done <<< "$1"
  return $bad
}

do_jobs() {
  local pids=() tags=()
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
      pids+=($!); tags+=("${tag}_k${k}")
    done
  done <<< "$1"
  [ ${#pids[@]} -eq 0 ] && { echo "  nothing to run"; return 0; }
  echo "  waiting for ${#pids[@]} job(s)..."
  wait "${pids[@]}" 2>/dev/null
  local i
  for i in "${!tags[@]}"; do
    local rc; rc="$(cat "$HERE/results/DONE_${tags[$i]}" 2>/dev/null || echo ?)"
    if [ "$rc" != 0 ]; then
      echo; echo "  !! ${tags[$i]} FAILED (rc=$rc). Last lines of its log:"
      tail -12 "$HERE/results/${tags[$i]}.log" 2>/dev/null | sed 's/^/     /'
      rm -f "$HERE/results/DONE_${tags[$i]}"   # so a re-run retries it
    fi
  done
}

score() {
  echo; echo "==================== $2 ===================="
  [ -n "$SCORE_PY" ] || { echo "no python3 on PATH for scoring" >&2; return 1; }
  while IFS='|' read -r tag D; do
    [ -n "$tag" ] || continue
    local gt="$HERE/gt/gt_${tag}.tsv"
    local a="$HERE/results/${tag}_k0/DANTE_LINE.gff3"
    local b="$HERE/results/${tag}_k7/DANTE_LINE.gff3"
    [ -f "$gt" ] && [ -f "$a" ] && [ -f "$b" ] || { echo "### $tag: incomplete, skipping"; continue; }
    "$SCORE_PY" "$HERE/score_extensions.py" --truth "$gt" --label "$tag" \
        --gff "k-default=$a" "k-7=$b"
  done <<< "$1"
  cat <<'TXT'

PRE-REGISTERED RULE (PREREGISTERED.md): adopt k=7 only if, on EVERY genome
above, coverage rises by >= 10 percentage points AND purity stays >= 95% AND
purity falls by no more than 1 point. Judge the WORST genome, never the
average. Any failure -> keep the mmseqs default.
TXT
}

# --------------------------------------------------------------------- modes
case "$MODE" in
  --setup)
    if resolve_env; then echo "an environment already works: $KIND $BIN"; else setup_env; fi ;;
  --check)
    echo "== environment =="
    if resolve_env; then echo "  OK  $KIND  $BIN"
    else echo "  NOT FOUND — run:  ./RUN.sh --setup"; fi
    echo "  scoring python3: ${SCORE_PY:-MISSING}"
    [ -n "${BIN:-}" ] && selftest_scripts
    echo "  container: $SIF $([ -f "$SIF" ] && echo '(present)' || echo '(ABSENT)')"
    echo "== inputs =="
    ok=0
    check_inputs "$JOBS_CAL" || ok=1
    check_inputs "$JOBS_VAL" || ok=1
    [ $ok -eq 0 ] && echo "  all inputs readable"
    echo "== already done =="
    ls "$HERE/results" 2>/dev/null | grep '^DONE_' | sed 's/^/  /' || echo "  (none)" ;;
  --list)
    resolve_env && echo "environment: $KIND ($BIN)" || echo "environment: NOT FOUND (./RUN.sh --setup)"
    echo "threads per job: $THREADS"
    echo; echo "CALIBRATION:"; echo "$JOBS_CAL" | sed 's/|/  ->  /'
    echo; echo "VALIDATION (held out):"; echo "$JOBS_VAL" | sed 's/|/  ->  /'
    echo; echo "already done:"; ls "$HERE/results" 2>/dev/null | grep '^DONE_' | sed 's/^/  /' || echo "  (none)" ;;
  --score)   score "$JOBS_CAL" "CALIBRATION" ;;
  --validate)
    resolve_env || { setup_env && resolve_env; } || exit 1
    echo "environment: $KIND ($BIN); threads=$THREADS"
    selftest_scripts || exit 1
    do_jobs "$JOBS_VAL"; score "$JOBS_VAL" "VALIDATION (held out)" ;;
  *)
    resolve_env || { echo "no usable environment; building one..."; setup_env && resolve_env; } || exit 1
    echo "environment: $KIND ($BIN); threads=$THREADS"
    check_inputs "$JOBS_CAL" || { echo "FATAL: inputs missing (see above)"; exit 1; }
    selftest_scripts || { echo "FATAL: the scripts do not load in this environment"; exit 1; }
    echo "running CALIBRATION (3 genomes x 2 settings)"
    do_jobs "$JOBS_CAL"; score "$JOBS_CAL" "CALIBRATION" ;;
esac
