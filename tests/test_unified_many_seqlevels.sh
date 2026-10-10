#!/usr/bin/env bash
# Regression test for GitHub issue #4: make_unified_annotation.R died with
#   Error: C stack usage ... is too close to the limit
# on a 305,167-sequence assembly, in identify_te_derived_trcs().
#
# Cause: every overlap / c() between two GRanges merges their Seqinfo. When EACH
# side has seqlevels the other lacks, GenomeInfoDb 1.30.x warns with a message
# that pastes every non-shared seqlevel name, and R's warning() overflows the C
# stack once that message is a few MB (suppressWarnings() does not help). The
# trigger needs (a) a few hundred thousand seqlevels and (b) realistic INSDC-style
# names — generic short names (s000001...) keep the message below the limit.
#
# This test builds that shape on the fly from the unified_multibatch fixture:
#   * 210,000 extra 10 kb sequences named JBXXXX01<7 digits>.1 in the .fai,
#   * one structural TE (DANTE_LTR) on each of them  -> t1 carries ~210k levels,
#   * one extra TideCluster array on a sequence with NO structural TE
#                                                    -> t3 has a level t1 lacks.
# PASS = make_unified completes at the default 8 MiB stack and its outputs on the
# original three sequences match a run on the unmodified fixture.
# Verified to FAIL against the pre-fix code (1.9.3): C stack overflow.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
FX="$REPO/tests/fixtures/unified_multibatch"
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT

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
SCRIPT="${MAKE_UNIFIED:-$REPO/scripts/make_unified_annotation.R}"

export PATH="$REPO/scripts:$PATH"
export CARP_VOCABULARY="$REPO/classification_vocabulary.yaml"

# The bug is a C stack overflow, so the stack limit is part of the test.
ulimit -s 8192

N=210000
BIG="$OUT/big"
cp -r "$FX" "$BIG"
awk -v n="$N" 'BEGIN { off = 6000003
  for (i = 1; i <= n + 1; i++) { printf "JBXXXX01%07d.1\t10000\t%d\t80\t81\n", i, off; off += 10126 } }' \
  >> "$BIG/genome.fai"
awk -v n="$N" 'BEGIN { for (i = 1; i <= n; i++)
  printf "JBXXXX01%07d.1\tDANTE_LTR\ttransposable_element\t1001\t8000\t.\t+\t.\tID=TE_bulk_%d;Final_Classification=Class_I|LTR|Ty3/gypsy|chromovirus|CRM;Name=Class_I|LTR|Ty3/gypsy|chromovirus|CRM;Rank=DL\n", i, i }' \
  >> "$BIG/DANTE_LTR.gff3"
# Array on the one extra sequence that carries no structural TE.
printf 'JBXXXX01%07d.1\tTideCluster\trepeat_region\t2001\t6000\t.\t+\t.\tName=TRC_2\n' $((N + 1)) \
  >> "$BIG/TideCluster_clustering.gff3"

run() {  # fixture_dir out_prefix
  "$RSCRIPT" "$SCRIPT" \
    --ltr        "$1/DANTE_LTR.gff3" \
    --ltr_tandems "$1/DANTE_LTR_tandems.gff3" \
    --tir        "$1/DANTE_TIR_combined.gff3" \
    --line       "$1/DANTE_LINE.gff3" \
    --dante      "$1/DANTE_filtered.gff3" \
    --tc_default "$1/TideCluster_clustering.gff3" \
    --tc_short   "$1/TideCluster_clustering_short.gff3" \
    --tc_rm      "$1/RM_on_TideCluster_Library.gff3" \
    --tc_rdna_default "$1/rdna_default.tsv" \
    --tc_rdna_short   "$1/rdna_short.tsv" \
    --tc_trc_table_default "$1/trc_table_default.tsv" \
    --tc_trc_table_short   "$1/trc_table_short.tsv" \
    --rm         "$1/RM_on_combined_library_plus_DANTE.gff3" \
    --th_default "$1/TideHunter_short_default.gff3" \
    --th_short   "$1/TideHunter_short_short.gff3" \
    --th_raw_default "$1/TideHunter_raw_default.gff3" \
    --th_raw_short   "$1/TideHunter_raw_short.gff3" \
    --rm_tc_tandem_gate TRUE \
    --fai        "$1/genome.fai" \
    --output     "$OUT/$2.gff3" \
    --threads    2 \
    > "$OUT/$2.log" 2>&1
}

fail() { echo "FAIL: $1" >&2; tail -20 "$OUT/big.log" >&2; exit 1; }

run "$FX"  small || { tail -20 "$OUT/small.log" >&2; echo "FAIL: baseline run failed" >&2; exit 1; }
run "$BIG" big   || fail "make_unified failed on the ${N}-sequence fixture"

grep -q "TE-derived TRCs (global): 1 default" "$OUT/big.log" \
  || fail "TE-derived TRC pass did not run / changed its decision"

# The bulk sequences must not change anything on the original three sequences.
# Feature IDs are a genome-wide running number (the bulk features sort first),
# so compare with ID/Parent values masked.
orig='^(seq_te|seq_ovl|seq_pad)\t'
pick() { grep -P "$orig" "$1" | sed -E 's/(ID|Parent)=[^;]*/\1=_/g'; }
diff <(pick "$OUT/small.gff3") <(pick "$OUT/big.gff3") > "$OUT/diff.txt" \
  || { head -20 "$OUT/diff.txt" >&2; fail "features on the original sequences differ"; }
cmp -s "$OUT/small.te_derived_trc.csv" "$OUT/big.te_derived_trc.csv" \
  || fail "te_derived_trc.csv differs"

echo "test_unified_many_seqlevels: PASSED ($N extra seqlevels, default 8 MiB stack)"
