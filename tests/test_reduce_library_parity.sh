#!/bin/bash
# Parity check: scripts/reduce_library_size.py vs scripts/reduce_library_size.R.
#
# Runs both implementations on the same input combined library FASTA and
# verifies the output FASTAs are byte-identical. The Python implementation
# is supposed to be a drop-in CLI-compatible replacement; this test is the
# proof. Until it runs green on a representative input, the Snakefile rule
# `reduce_library` should keep calling the R version.
#
# Usage:
#   ./tests/test_reduce_library_parity.sh <input_combined_library.fasta>
#
# Optional environment variables:
#   THREADS  — how many workers to give each implementation (default 4)
#   RSCRIPT  — Rscript binary to use (default: from PATH)
#   PYTHON3  — Python 3 binary to use (default: from PATH)
#   KEEP_WORKDIR — if set, do not delete the workdir on success
#
# Exit codes:
#   0  — outputs are byte-identical
#   1  — outputs differ (workdir preserved for inspection)
#   2  — either implementation failed to produce an output FASTA
#   3  — usage error
#
# Recommended invocations during development:
#   ./tests/test_reduce_library_parity.sh tests/fixtures/output_micro/Libraries/combined_library.fasta
#   ./tests/test_reduce_library_parity.sh tests/fixtures/output_small/Libraries/combined_library.fasta
#   ./tests/test_reduce_library_parity.sh tests/fixtures/output_medium/Libraries/combined_library.fasta
#
# Each fixture exercises a different mix of CAP3 / mmseqs2 / BLAST paths.
# The medium fixture is the one CI's release.yml test-in-container step
# uses; passing on it is the green-light signal for the Snakefile cutover.

set -euo pipefail

if [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ] || [ -z "${1:-}" ]; then
    sed -n '2,40p' "$0" | sed 's/^# \{0,1\}//'
    exit 3
fi

INPUT="$1"
THREADS="${THREADS:-4}"
RSCRIPT="${RSCRIPT:-Rscript}"
PYTHON3="${PYTHON3:-python3}"

if [ ! -f "$INPUT" ]; then
    echo "ERROR: input not found: $INPUT" >&2
    exit 3
fi

# Resolve to absolute paths so the implementations don't get confused by
# their own per-run working dirs.
INPUT_ABS="$(readlink -f "$INPUT")"
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
R_SCRIPT="$REPO_ROOT/scripts/reduce_library_size.R"
PY_SCRIPT="$REPO_ROOT/scripts/reduce_library_size.py"

for s in "$R_SCRIPT" "$PY_SCRIPT"; do
    [ -f "$s" ] || { echo "ERROR: missing $s" >&2; exit 3; }
done

# Sanity-check the binaries are reachable. cap3, mmseqs, makeblastdb, and
# blastn are needed by both implementations; if any is missing both will
# fail identically (which would be misreported as parity).
for tool in cap3 mmseqs makeblastdb blastn; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: '$tool' not on PATH; activate the conda env that ships it" >&2
        echo "       (envs/tidecluster.yaml provides all four)" >&2
        exit 3
    fi
done

WORKDIR="$(mktemp -d -t reduce_lib_parity.XXXXXX)"
echo "parity test: input=$INPUT_ABS  threads=$THREADS  workdir=$WORKDIR"
echo "  R script:  $R_SCRIPT"
echo "  Py script: $PY_SCRIPT"
echo

# ── R reference run ──────────────────────────────────────────────────
echo "[R] running…"
if ! "$RSCRIPT" "$R_SCRIPT" \
        -i "$INPUT_ABS" \
        -o "$WORKDIR/r_out.fasta" \
        -t "$THREADS" \
        -d "$WORKDIR/r_wd" \
        > "$WORKDIR/r.stdout" 2> "$WORKDIR/r.stderr"; then
    echo "ERROR: R implementation failed; see $WORKDIR/r.stderr" >&2
    tail -20 "$WORKDIR/r.stderr" >&2
    exit 2
fi
if [ ! -s "$WORKDIR/r_out.fasta" ]; then
    echo "ERROR: R produced no output (empty file at $WORKDIR/r_out.fasta)" >&2
    exit 2
fi

# ── Python implementation under test ─────────────────────────────────
echo "[Py] running…"
if ! "$PYTHON3" "$PY_SCRIPT" \
        -i "$INPUT_ABS" \
        -o "$WORKDIR/py_out.fasta" \
        -t "$THREADS" \
        -d "$WORKDIR/py_wd" \
        > "$WORKDIR/py.stdout" 2> "$WORKDIR/py.stderr"; then
    echo "ERROR: Python implementation failed; see $WORKDIR/py.stderr" >&2
    tail -20 "$WORKDIR/py.stderr" >&2
    exit 2
fi
if [ ! -s "$WORKDIR/py_out.fasta" ]; then
    echo "ERROR: Python produced no output (empty file at $WORKDIR/py_out.fasta)" >&2
    exit 2
fi

# ── Byte-identity verdict ────────────────────────────────────────────
R_HASH="$(md5sum "$WORKDIR/r_out.fasta" | awk '{print $1}')"
PY_HASH="$(md5sum "$WORKDIR/py_out.fasta" | awk '{print $1}')"
R_BYTES="$(wc -c < "$WORKDIR/r_out.fasta")"
PY_BYTES="$(wc -c < "$WORKDIR/py_out.fasta")"

echo
echo "  R  output: $R_HASH  ($R_BYTES bytes)"
echo "  Py output: $PY_HASH  ($PY_BYTES bytes)"

if cmp -s "$WORKDIR/r_out.fasta" "$WORKDIR/py_out.fasta"; then
    echo
    echo "PARITY: byte-identical output ✓"
    if [ -z "${KEEP_WORKDIR:-}" ]; then
        rm -rf "$WORKDIR"
    else
        echo "  workdir preserved at $WORKDIR (KEEP_WORKDIR set)"
    fi
    exit 0
fi

# ── Diff diagnostic — outputs differ ─────────────────────────────────
echo
echo "PARITY FAILED: outputs differ" >&2

# Structural quick-look: same record count? same name set?
echo "  R record count:  $(grep -c '^>' "$WORKDIR/r_out.fasta")" >&2
echo "  Py record count: $(grep -c '^>' "$WORKDIR/py_out.fasta")" >&2
echo

# Diff sample (first 40 lines that differ).
echo "  diff sample (first 40 lines):" >&2
diff "$WORKDIR/r_out.fasta" "$WORKDIR/py_out.fasta" 2>/dev/null \
    | head -40 \
    | sed 's/^/    /' >&2

echo >&2
echo "  workdir preserved for inspection: $WORKDIR" >&2
echo "    R output:  $WORKDIR/r_out.fasta" >&2
echo "    Py output: $WORKDIR/py_out.fasta" >&2
echo "    R stderr:  $WORKDIR/r.stderr" >&2
echo "    Py stderr: $WORKDIR/py.stderr" >&2
echo "    R workdir: $WORKDIR/r_wd/" >&2
echo "    Py workdir: $WORKDIR/py_wd/" >&2

exit 1
