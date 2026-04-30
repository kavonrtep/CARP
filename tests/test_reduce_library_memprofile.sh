#!/bin/bash
# Memory-profile comparison: scripts/reduce_library_size.py vs .R.
#
# Wraps both implementations under /usr/bin/time -v (the GNU time, not
# the bash builtin) which reports honest peak resident-set size for the
# whole process tree — no COW double-counting like snakemake's
# benchmark TSV. Output is a small two-row table comparing peak RSS,
# wall time, and CPU time.
#
# Usage:
#   ./tests/test_reduce_library_memprofile.sh <input_combined_library.fasta>
#
# Optional environment variables:
#   THREADS  — workers per implementation (default 4)
#   RSCRIPT  — Rscript binary (default: from PATH)
#   PYTHON3  — Python 3 binary (default: from PATH)
#
# Requires: /usr/bin/time (Debian/Ubuntu: apt install time).
# Note: the bash builtin `time` does not support -v; this script
# explicitly invokes /usr/bin/time.
#
# Exit codes:
#   0 — both runs completed; report printed
#   2 — either implementation failed
#   3 — usage / dependency error

set -euo pipefail

if [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ] || [ -z "${1:-}" ]; then
    sed -n '2,25p' "$0" | sed 's/^# \{0,1\}//'
    exit 3
fi

INPUT="$1"
THREADS="${THREADS:-4}"
RSCRIPT="${RSCRIPT:-Rscript}"
PYTHON3="${PYTHON3:-python3}"
TIME_BIN=/usr/bin/time

if [ ! -f "$INPUT" ]; then
    echo "ERROR: input not found: $INPUT" >&2
    exit 3
fi
if [ ! -x "$TIME_BIN" ]; then
    echo "ERROR: $TIME_BIN not available — install GNU time:" >&2
    echo "         sudo apt install time   (Debian/Ubuntu)" >&2
    echo "         brew install gnu-time   (macOS)" >&2
    exit 3
fi

INPUT_ABS="$(readlink -f "$INPUT")"
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
R_SCRIPT="$REPO_ROOT/scripts/reduce_library_size.R"
PY_SCRIPT="$REPO_ROOT/scripts/reduce_library_size.py"

WORKDIR="$(mktemp -d -t reduce_lib_memprofile.XXXXXX)"
trap 'rm -rf "$WORKDIR"' EXIT

# Custom format — one line per probe, no preamble:
#   peak_rss_kb wall_seconds user_seconds system_seconds
TIME_FMT='%M %e %U %S'

echo "memprofile: input=$INPUT_ABS  threads=$THREADS"
echo "          : workdir=$WORKDIR (cleaned up on exit)"
echo

run_probe() {
    local label="$1" timing_file="$2"
    shift 2
    if ! "$TIME_BIN" -f "$TIME_FMT" -o "$timing_file" "$@" \
            > "$WORKDIR/$label.stdout" 2> "$WORKDIR/$label.stderr"; then
        echo "ERROR: $label run failed; see $WORKDIR/$label.stderr" >&2
        tail -20 "$WORKDIR/$label.stderr" >&2
        exit 2
    fi
}

echo "[R]  probing…"
run_probe r "$WORKDIR/r.time" \
    "$RSCRIPT" "$R_SCRIPT" \
    -i "$INPUT_ABS" -o "$WORKDIR/r.fasta" \
    -t "$THREADS" -d "$WORKDIR/r_wd"

echo "[Py] probing…"
run_probe py "$WORKDIR/py.time" \
    "$PYTHON3" "$PY_SCRIPT" \
    -i "$INPUT_ABS" -o "$WORKDIR/py.fasta" \
    -t "$THREADS" -d "$WORKDIR/py_wd"

read -r r_rss_kb r_wall r_user r_sys < "$WORKDIR/r.time"
read -r py_rss_kb py_wall py_user py_sys < "$WORKDIR/py.time"

# Convert peak RSS from KB to MB for readability
r_rss_mb=$(awk -v k="$r_rss_kb" 'BEGIN{printf "%.1f", k/1024}')
py_rss_mb=$(awk -v k="$py_rss_kb" 'BEGIN{printf "%.1f", k/1024}')

# Reduction ratio (R peak / Py peak), useful headline number.
ratio=$(awk -v r="$r_rss_kb" -v p="$py_rss_kb" \
    'BEGIN { if (p > 0) printf "%.1f", r/p; else print "NA" }')

# Output table — fixed-width so it's easy to copy into a commit message
# or release note.
echo
printf '%-15s | %12s | %12s | %12s | %12s\n' \
    "implementation" "peak RSS (MB)" "wall (s)" "user (s)" "sys (s)"
printf '%s\n' "------------------------------------------------------------------------------"
printf '%-15s | %12s | %12s | %12s | %12s\n' "R"      "$r_rss_mb"  "$r_wall"  "$r_user"  "$r_sys"
printf '%-15s | %12s | %12s | %12s | %12s\n' "Python" "$py_rss_mb" "$py_wall" "$py_user" "$py_sys"
echo
echo "R peak / Python peak = ${ratio}x"

# Sanity: both produced output. We don't compare bytes here — that's
# what test_reduce_library_parity.sh is for. But report sizes so any
# accidental empty output is obvious.
r_sz=$(wc -c < "$WORKDIR/r.fasta" 2>/dev/null || echo 0)
py_sz=$(wc -c < "$WORKDIR/py.fasta" 2>/dev/null || echo 0)
echo "Output FASTA sizes: R=$r_sz B  Py=$py_sz B  (parity check: run test_reduce_library_parity.sh separately)"
