#!/usr/bin/env bash
# Assert that a fixture run produced the expected minimum outputs.
#
# Usage:  scripts/assert_fixture_outputs.sh <fixture-name>
#   fixture-name ∈ {micro, tiny, small, medium}
#
# Thresholds are tight enough to catch regressions but loose enough to tolerate
# harmless drift in upstream tool versions. Update when fixture content changes.
set -euo pipefail

fixture="${1:?usage: $0 <fixture>}"
out="tests/fixtures/output_${fixture}"

if [[ ! -d "$out" ]]; then
  echo "FAIL: output directory missing: $out" >&2
  exit 1
fi

# Every fixture must produce these top-level artefacts
require_nonempty() {
  [[ -s "$1" ]] || { echo "FAIL: $1 missing or empty" >&2; exit 1; }
}
require_exists() {
  [[ -e "$1" ]] || { echo "FAIL: $1 missing" >&2; exit 1; }
}
require_nonempty "$out/summary_statistics.csv"
require_nonempty "$out/Repeat_Annotation_Unified.gff3"
# Marker file intentionally 0 bytes — only its existence matters.
require_exists   "$out/.classifications_validated"

count_fasta() {
  [[ -f "$1" ]] || { echo 0; return; }
  # grep -c exits non-zero when count is 0; normalise to just the number.
  grep -c '^>' "$1" || true
}
count_gff_type() {
  [[ -f "$1" ]] || { echo 0; return; }
  awk -F'\t' -v t="$2" '$3==t' "$1" | wc -l
}

primary_tir=$(count_fasta "$out/DANTE_TIR/DANTE_TIR_final.fasta")
combined_tir=$(count_gff_type "$out/DANTE_TIR/DANTE_TIR_combined.gff3" sequence_feature)
line_n=$(count_gff_type "$out/DANTE_LINE/DANTE_LINE.gff3" LINE_element)
ltr_lib=$(count_fasta "$out/DANTE_LTR/LTR_RTs_library.fasta")
sumstats_rows=$(($(wc -l < "$out/summary_statistics.csv") - 1))

case "$fixture" in
  micro|tiny)
    # Smoke fixtures — only verify completion, no positive-count checks.
    ;;
  small)
    # Small exercises fallback TIR + RepeatMasker paths end-to-end.
    [[ "$combined_tir" -ge 1 ]]  || { echo "FAIL: $fixture expected ≥1 combined TIR, got $combined_tir" >&2; exit 1; }
    [[ "$ltr_lib"      -ge 5 ]]  || { echo "FAIL: $fixture expected ≥5 LTR lib seqs, got $ltr_lib" >&2; exit 1; }
    [[ "$sumstats_rows" -ge 10 ]] || { echo "FAIL: $fixture expected ≥10 summary classes, got $sumstats_rows" >&2; exit 1; }
    ;;
  medium)
    # Medium is the release gate — primary DANTE_TIR + DANTE_LINE must fire.
    [[ "$primary_tir"  -ge 3  ]] || { echo "FAIL: $fixture expected ≥3 primary TIR, got $primary_tir" >&2; exit 1; }
    [[ "$combined_tir" -ge 15 ]] || { echo "FAIL: $fixture expected ≥15 combined TIR, got $combined_tir" >&2; exit 1; }
    [[ "$line_n"       -ge 5  ]] || { echo "FAIL: $fixture expected ≥5 LINE elements, got $line_n" >&2; exit 1; }
    [[ "$ltr_lib"      -ge 5  ]] || { echo "FAIL: $fixture expected ≥5 LTR lib seqs, got $ltr_lib" >&2; exit 1; }
    [[ "$sumstats_rows" -ge 12 ]] || { echo "FAIL: $fixture expected ≥12 summary classes, got $sumstats_rows" >&2; exit 1; }
    ;;
  *)
    echo "FAIL: unknown fixture '$fixture'" >&2
    exit 2
    ;;
esac

cat <<EOT
OK: $fixture outputs meet expectations
  primary_TIR  = $primary_tir
  combined_TIR = $combined_tir
  LINE         = $line_n
  LTR_lib      = $ltr_lib
  summary rows = $sumstats_rows
EOT
