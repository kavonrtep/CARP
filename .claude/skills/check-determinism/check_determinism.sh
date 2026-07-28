#!/usr/bin/env bash
# Manually-triggered FULL-SIZE determinism check. The every-PR CI gate
# (pipeline.yml `determinism` job) runs the small fixture twice, but several
# non-determinism bugs only surfaced on a real, full-size genome — so run this on
# such a genome before trusting a release / a big change.
#
# It runs the whole pipeline TWICE on the SAME config — the two passes under a
# different PYTHONHASHSEED (exposes Python set / hash-order leaks) and a different
# thread count (exposes thread-order non-determinism) — then diffs the manifest
# data outputs with scripts/assert_run_determinism.py. Byte-identical => pass.
#
# Usage (run where run_pipeline.py works — host/HPC, NOT the sandbox where the
# conda wrapper hangs):
#   .claude/skills/check-determinism/check_determinism.sh -c my_genome.yaml [-a T1] [-b T2]
# -c  CARP config (its output_dir is used; two sibling dirs _det_a / _det_b are made)
# -a  threads for run 1 (default: all cores)
# -b  threads for run 2 (default: half, min 1 — deliberately different from -a)
#
# NOTE: a full-size genome takes ~30-60 min per pass, so this is ~1-2 h. It is a
# manual gate, never wired into per-PR CI.
set -euo pipefail

CONFIG="" ; T1="" ; T2=""
while getopts "c:a:b:h" opt; do
  case "$opt" in
    c) CONFIG="$OPTARG" ;;
    a) T1="$OPTARG" ;;
    b) T2="$OPTARG" ;;
    h) grep '^#' "$0" | sed 's/^# \?//'; exit 0 ;;
    *) echo "see -h" >&2; exit 2 ;;
  esac
done
[ -n "$CONFIG" ] || { echo "usage: check_determinism.sh -c <config.yaml> [-a t1] [-b t2]" >&2; exit 2; }
[ -f "$CONFIG" ] || { echo "config not found: $CONFIG" >&2; exit 2; }

ncores=$(nproc 2>/dev/null || echo 4)
T1="${T1:-$ncores}"
T2="${T2:-$(( ncores / 2 > 0 ? ncores / 2 : 1 ))}"

repo_root=$(cd "$(dirname "$0")/../../.." && pwd)   # .claude/skills/check-determinism -> repo root
outdir=$(python3 -c "import yaml,sys; print(yaml.safe_load(open(sys.argv[1]))['output_dir'])" "$CONFIG")
a="${outdir}_det_a" ; b="${outdir}_det_b"

echo "== determinism check =="
echo "  config:  $CONFIG   (output_dir: $outdir)"
echo "  run A:   PYTHONHASHSEED=1  -t $T1  -> $a"
echo "  run B:   PYTHONHASHSEED=2  -t $T2  -> $b"

rm -rf "$outdir" "$a" "$b"
PYTHONHASHSEED=1 python3 run_pipeline.py -c "$CONFIG" -t "$T1"
mv "$outdir" "$a"
PYTHONHASHSEED=2 python3 run_pipeline.py -c "$CONFIG" -t "$T2"
mv "$outdir" "$b"

python3 "$repo_root/scripts/assert_run_determinism.py" "$a" "$b"
