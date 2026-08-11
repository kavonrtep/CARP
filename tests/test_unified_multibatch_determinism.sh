#!/usr/bin/env bash
# Multi-batch determinism gate for make_unified_annotation.R.
#
# make_unified splits the genome into per-thread batches whose COMPOSITION depends
# on the thread count (threads=1 -> one batch; threads>1 -> make_batches). Several
# per-batch operations were batch-composition-dependent -- most severely the
# TE-derived-satellite (te_sat) pre-pass, which trimmed still-self-overlapping t1
# against the satellite regions and thereby DISJOINED t1's internal overlaps, but
# only in batches that contained a te_sat array (see
# docs/archive/tier1_resolution_determinism_audit.md). On a large genome the same
# sequence resolved differently at threads=1 vs threads=N (816/2.82M features).
#
# The full-pipeline determinism gate (pipeline.yml) can't catch this: its small /
# medium fixtures are below the batch-target floor, so make_unified always takes
# the single-batch path regardless of thread count. This test drives make_unified
# DIRECTLY on a tiny purpose-built fixture whose .fai declares 3x 2 Mb sequences
# (make_unified reads only lengths from the .fai, never the FASTA) so
# --batch_size 1000000 forces 3 batches at threads>1 while threads=1 stays single-
# batch. seq_te carries a TE-derived-TRC array (the te_sat trigger); seq_ovl
# carries overlapping tier-1 elements in a DIFFERENT batch (the bug's victim).
#
# PASS = the unified GFF3 is byte-identical single-batch vs multi-batch.
# Verified to FAIL against the pre-fix code (commit 67acf4c): 9 differing lines.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
FX="$REPO/tests/fixtures/unified_multibatch"
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT

# Locate an Rscript with the bioconductor stack. Prefer $RSCRIPT; else scan the
# snakemake-managed conda envs by filesystem (no R startup) for the needed libs.
RSCRIPT="${RSCRIPT:-}"
if [[ -z "$RSCRIPT" ]]; then
  for d in "$REPO"/.snakemake/conda/*; do
    lib="$d/lib/R/library"
    if [[ -d "$lib/rtracklayer" && -d "$lib/GenomicRanges" && -d "$lib/optparse" ]]; then
      RSCRIPT="$d/bin/Rscript"; break
    fi
  done
fi
if [[ -z "$RSCRIPT" || ! -x "$RSCRIPT" ]]; then
  echo "SKIP: no Rscript with rtracklayer/GenomicRanges/optparse found" >&2
  echo "      (set RSCRIPT=... to run this test locally)" >&2
  exit 0
fi

export PATH="$REPO/scripts:$PATH"
export CARP_VOCABULARY="$REPO/classification_vocabulary.yaml"

run() {  # tag threads batch_size [extra args...]
  local tag="$1" threads="$2" bsize="$3"; shift 3
  "$RSCRIPT" "$REPO/scripts/make_unified_annotation.R" \
    --ltr        "$FX/DANTE_LTR.gff3" \
    --ltr_tandems "$FX/DANTE_LTR_tandems.gff3" \
    --tir        "$FX/DANTE_TIR_combined.gff3" \
    --line       "$FX/DANTE_LINE.gff3" \
    --dante      "$FX/DANTE_filtered.gff3" \
    --tc_default "$FX/TideCluster_clustering.gff3" \
    --tc_short   "$FX/TideCluster_clustering_short.gff3" \
    --tc_rm      "$FX/RM_on_TideCluster_Library.gff3" \
    --tc_rdna_default "$FX/rdna_default.tsv" \
    --tc_rdna_short   "$FX/rdna_short.tsv" \
    --tc_trc_table_default "$FX/trc_table_default.tsv" \
    --tc_trc_table_short   "$FX/trc_table_short.tsv" \
    --rm         "$FX/RM_on_combined_library_plus_DANTE.gff3" \
    --th_default "$FX/TideHunter_short_default.gff3" \
    --th_short   "$FX/TideHunter_short_short.gff3" \
    --th_raw_default "$FX/TideHunter_raw_default.gff3" \
    --th_raw_short   "$FX/TideHunter_raw_short.gff3" \
    --rm_tc_tandem_gate TRUE \
    --fai        "$FX/genome.fai" \
    --output     "$OUT/unified_${tag}.gff3" \
    --threads    "$threads" \
    --batch_size "$bsize" \
    "$@" \
    > "$OUT/log_${tag}.txt" 2>&1
}

run single 1 200000000
run multi  3 1000000
# Worker-pool cap must not change the result. --threads (and therefore batch
# COMPOSITION) is held at 3 while only the CONCURRENCY is squeezed to 1, so this
# isolates the mclapply pool size introduced by make_unified_max_workers /
# --mem_budget_gb. Guards the property those knobs are only safe because of.
run serialw 3 1000000 --max_workers 1

fail() { echo "FAIL: $1" >&2; echo "--- single log tail ---" >&2; tail -20 "$OUT/log_single.txt" >&2
         echo "--- multi log tail ---" >&2; tail -20 "$OUT/log_multi.txt" >&2; exit 1; }

# The fixture must actually exercise the paths it claims to guard, else a green
# result is meaningless: te_sat non-empty (single batch), and a genuine >1-batch
# split with a te_sat batch AND a te_sat-free batch (multi).
grep -q "TE-derived TRC(s) → te_sat=1" "$OUT/log_single.txt" \
  || fail "te_sat did not fire in the single-batch run (fixture no longer triggers it)"
grep -qE "Batching: [0-9]+ batch" "$OUT/log_multi.txt" \
  || fail "multi run did not split into >1 batch (fixture no longer forces batching)"
grep -qE "batch\[seq_te\].*te_sat=1" "$OUT/log_multi.txt" \
  || fail "no te_sat batch in the multi-batch run"
grep -qE "batch\[seq_ovl\].*te_sat=0" "$OUT/log_multi.txt" \
  || fail "seq_ovl did not land in a te_sat-free batch (the split that reproduces the bug)"

# The worker cap must actually have bitten, or the comparison below proves nothing.
grep -qE "capping workers [0-9]+ -> 1 \(--max_workers\)" "$OUT/log_serialw.txt" \
  || fail "--max_workers 1 did not cap the pool (gating not reached — comparison would be vacuous)"
grep -qE "Batching: 3 batch" "$OUT/log_serialw.txt" \
  || fail "capped run did not still produce 3 batches (concurrency cap must not change batch composition)"

rc=0
if diff <(grep -v '^#' "$OUT/unified_single.gff3") \
        <(grep -v '^#' "$OUT/unified_multi.gff3") > "$OUT/diff.txt"; then
  echo "test_unified_multibatch_determinism: PASSED (byte-identical single-batch vs 3-batch)"
else
  n=$(grep -c '^[<>]' "$OUT/diff.txt" || true)
  echo "FAIL: unified GFF3 differs across batchings ($n lines):" >&2
  head -40 "$OUT/diff.txt" >&2
  rc=1
fi

if diff <(grep -v '^#' "$OUT/unified_multi.gff3") \
        <(grep -v '^#' "$OUT/unified_serialw.gff3") > "$OUT/diff_workers.txt"; then
  echo "test_unified_multibatch_determinism: PASSED (byte-identical 3 workers vs 1 worker, same batching)"
else
  n=$(grep -c '^[<>]' "$OUT/diff_workers.txt" || true)
  echo "FAIL: unified GFF3 differs by WORKER COUNT at fixed batching ($n lines)." >&2
  echo "      make_unified_max_workers / --mem_budget_gb are only safe if this stays identical." >&2
  head -40 "$OUT/diff_workers.txt" >&2
  rc=1
fi
exit $rc
