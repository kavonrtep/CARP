#!/usr/bin/env bash
# =============================================================================
# Measure make_unified_annotation on run-000156's real inputs.
#
# WHY: run-000156 (94.26 Gbp) died in this rule with 15 of 55 mclapply workers
# killed, and left NO memory data behind -- the rule's benchmark file is only
# written on success. The rule now logs /proc-based [mem] lines (parent RSS at
# each phase, per-worker VmHWM, host MemAvailable at the fork point). This script
# re-runs ONLY that rule against the already-computed carp_output/ so we get real
# numbers to size the batching + concurrency fixes from, without re-running the
# 10 days of upstream work.
#
# Run it at a LOW thread count first (default 4). Two questions get answered at
# once:
#   * does one batch finish quickly with no memory pressure?  -> confirms the
#     17,045 s batches were host thrashing, not per-batch compute
#   * what does one worker actually cost (VmHWM)?             -> sizes the pool
#
# USAGE
#   ./measure_unified_run156.sh              # threads=4  (recommended first run)
#   ./measure_unified_run156.sh 16           # threads=16 (only after 4 succeeds)
#
# Detach it -- loading alone reads ~40 GB of GFF3 and took 64 min in the original
# run, before any batch starts:
#   nohup ./measure_unified_run156.sh 4 > /dev/null 2>&1 &
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------- EDIT THESE --
REPO="${CARP_REPO:-$HOME/assembly_repeat_annotation_pipeline}"   # git checkout on this host
RUN="${CARP_RUN:-/nfsroot/projects/darwin/runs/run-000156/carp_output}"
WORK="${CARP_WORK:-/nfsroot/projects/darwin/runs/run-000156/unified_measure}"
ENV_NAME="${CARP_ENV:-carp-unified}"       # conda env created from envs/tidecluster.yaml
# ------------------------------------------------------------------------------

THREADS="${1:-4}"
STAMP="$(date +%Y%m%d-%H%M%S)"
LOG="$WORK/measure-t${THREADS}-${STAMP}.err"     # the [mem] lines land here (stderr)
OUT="$WORK/Repeat_Annotation_Unified.t${THREADS}.gff3"
TIMING="$WORK/measure-t${THREADS}-${STAMP}.time"

say() { printf '\n=== %s ===\n' "$*"; }

# --- 1. update the checkout ----------------------------------------------------
say "1/5  updating $REPO"
[ -d "$REPO/.git" ] || { echo "ERROR: no git checkout at $REPO (set CARP_REPO)"; exit 1; }
cd "$REPO"
git pull --ff-only
echo "at commit: $(git rev-parse --short HEAD)  $(git log -1 --format=%s)"
# The instrumentation must be present, or the run tells us nothing.
grep -q 'log_mem <- function' scripts/make_unified_annotation.R \
  || { echo "ERROR: scripts/make_unified_annotation.R has no [mem] logging — pull did not bring the fix"; exit 1; }

# --- 2. preflight: inputs + writable output ------------------------------------
say "2/5  checking inputs under $RUN"
inputs=(
  DANTE_LTR/DANTE_LTR.gff3
  DANTE_LTR/DANTE_LTR_tandems.gff3
  DANTE_TIR/DANTE_TIR_combined.gff3
  DANTE_LINE/DANTE_LINE.gff3
  DANTE/DANTE_filtered.gff3
  TideCluster/default/TideCluster_clustering.gff3
  TideCluster/short_monomer/TideCluster_clustering.gff3
  TideCluster/default/RM_on_TideCluster_Library.gff3
  TideCluster/default/TideCluster_rdna.tsv
  TideCluster/short_monomer/TideCluster_rdna.tsv
  TideCluster/default/TideCluster_report/data/trc_table.tsv
  TideCluster/short_monomer/TideCluster_report/data/trc_table.tsv
  RepeatMasker/RM_on_combined_library_plus_DANTE.gff3
  TideCluster/default/TideCluster_tidehunter_short.gff3
  TideCluster/short_monomer/TideCluster_tidehunter_short.gff3
  TideCluster/default/TideCluster_tidehunter.gff3
  TideCluster/short_monomer/TideCluster_tidehunter.gff3
  genome_cleaned.fasta.fai
)
missing=0
for f in "${inputs[@]}"; do
  if [ -s "$RUN/$f" ]; then printf '  ok       %s\n' "$f"
  else printf '  MISSING  %s\n' "$f"; missing=1; fi
done
[ "$missing" -eq 0 ] || { echo "ERROR: inputs missing under $RUN"; exit 1; }

# carp_output/ is read-only in place -- everything we write goes to $WORK.
mkdir -p "$WORK"
[ -w "$WORK" ] || { echo "ERROR: $WORK is not writable (set CARP_WORK)"; exit 1; }
echo "output dir: $WORK  ($(df -h "$WORK" | awk 'NR==2{print $4}') free)"

# --- 3. conda env --------------------------------------------------------------
# Needs: optparse, rtracklayer, GenomicRanges, parallel, yaml -- i.e. exactly the
# rule's env. Build it from the repo's own envs/tidecluster.yaml so the R stack
# matches what the pipeline uses. ~234 packages; one-off, then reused.
say "3/5  resolving conda env '$ENV_NAME'"
command -v conda >/dev/null || { echo "ERROR: conda not on PATH"; exit 1; }
CONDA_BASE="$(conda info --base)"
# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"
SOLVER=mamba; command -v mamba >/dev/null || SOLVER=conda

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "env '$ENV_NAME' already exists — reusing (delete it to rebuild:"
  echo "  conda env remove -n $ENV_NAME )"
else
  echo "creating '$ENV_NAME' from envs/tidecluster.yaml with $SOLVER (this takes a while)..."
  $SOLVER env create -n "$ENV_NAME" -f envs/tidecluster.yaml
fi
conda activate "$ENV_NAME"

echo "checking R packages..."
Rscript -e 'for (p in c("optparse","rtracklayer","GenomicRanges","parallel","yaml")) {
              ok <- requireNamespace(p, quietly=TRUE)
              cat(sprintf("  %-16s %s\n", p, if (ok) "ok" else "MISSING"))
              if (!ok) quit(status=1) }' \
  || { echo "ERROR: R dependencies missing in '$ENV_NAME'"; exit 1; }

# --- 4. run --------------------------------------------------------------------
# /usr/bin/time -v gives the whole-tree peak as an independent cross-check on the
# script's own [mem] numbers.
say "4/5  running make_unified_annotation.R (threads=$THREADS)"
echo "log:    $LOG"
echo "output: $OUT"
echo "started $(date)"
TIMEBIN=""; [ -x /usr/bin/time ] && TIMEBIN="/usr/bin/time -v -o $TIMING"

set +e
$TIMEBIN "$REPO/scripts/make_unified_annotation.R" \
  --ltr                  "$RUN/DANTE_LTR/DANTE_LTR.gff3" \
  --ltr_tandems          "$RUN/DANTE_LTR/DANTE_LTR_tandems.gff3" \
  --tir                  "$RUN/DANTE_TIR/DANTE_TIR_combined.gff3" \
  --line                 "$RUN/DANTE_LINE/DANTE_LINE.gff3" \
  --dante                "$RUN/DANTE/DANTE_filtered.gff3" \
  --tc_default           "$RUN/TideCluster/default/TideCluster_clustering.gff3" \
  --tc_short             "$RUN/TideCluster/short_monomer/TideCluster_clustering.gff3" \
  --tc_rm                "$RUN/TideCluster/default/RM_on_TideCluster_Library.gff3" \
  --tc_rdna_default      "$RUN/TideCluster/default/TideCluster_rdna.tsv" \
  --tc_rdna_short        "$RUN/TideCluster/short_monomer/TideCluster_rdna.tsv" \
  --tc_trc_table_default "$RUN/TideCluster/default/TideCluster_report/data/trc_table.tsv" \
  --tc_trc_table_short   "$RUN/TideCluster/short_monomer/TideCluster_report/data/trc_table.tsv" \
  --rm                   "$RUN/RepeatMasker/RM_on_combined_library_plus_DANTE.gff3" \
  --th_default           "$RUN/TideCluster/default/TideCluster_tidehunter_short.gff3" \
  --th_short             "$RUN/TideCluster/short_monomer/TideCluster_tidehunter_short.gff3" \
  --th_raw_default       "$RUN/TideCluster/default/TideCluster_tidehunter.gff3" \
  --th_raw_short         "$RUN/TideCluster/short_monomer/TideCluster_tidehunter.gff3" \
  --rm_tc_tandem_gate    True \
  --fai                  "$RUN/genome_cleaned.fasta.fai" \
  --output               "$OUT" \
  --threads              "$THREADS" \
  2> "$LOG"
RC=$?
set -e
echo "finished $(date)  exit=$RC"

# --- 5. summary ----------------------------------------------------------------
say "5/5  summary"
echo "--- memory ([mem] lines) ---"
grep -F '[mem]' "$LOG" || echo "  (none — is this the instrumented version?)"

echo
echo "--- batching ---"
grep -E 'Batching:|batch sizes:' "$LOG" || true

echo
# The completion line is "◀ batch[SEQ]  (N seqs)  (T s) — final: ...", so take the
# LAST parenthesised number (greedy .*), not the first — the first is "(N seqs)".
echo "--- per-batch wall time in seconds (slowest 5) ---"
grep -F '◀ batch[' "$LOG" 2>/dev/null \
  | sed -nE 's/.*◀ (batch\[[^]]*\]).*\(([0-9.]+)s\).*/\2 \1/p' | sort -n | tail -5 || true
printf 'batches started : %s\n' "$(grep -c '▶ batch\[' "$LOG" || echo 0)"
printf 'batches finished: %s\n' "$(grep -c '◀ batch\[' "$LOG" || echo 0)"

echo
echo "--- worker loss (must be silent) ---"
grep -F 'did not deliver results' "$LOG" && echo "  ^^ WORKERS WERE KILLED (very likely OOM)" || echo "  none"

if [ -n "$TIMEBIN" ] && [ -f "$TIMING" ]; then
  echo
  echo "--- /usr/bin/time -v cross-check ---"
  grep -E 'Maximum resident set size|Elapsed \(wall clock\)|Major .page faults' "$TIMING" || true
fi

echo
if [ "$RC" -eq 0 ]; then
  echo "RESULT: completed. features=$(grep -vc '^#' "$OUT" 2>/dev/null || echo '?')  -> $OUT"
  echo "Report back the [mem] block above; it sizes the batching + concurrency fix."
else
  echo "RESULT: FAILED (exit $RC). Last 30 log lines:"
  tail -30 "$LOG"
fi
exit "$RC"
