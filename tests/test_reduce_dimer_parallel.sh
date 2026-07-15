#!/bin/bash
# Regression test for the reduce_dimer_library.py per-TRC parallelization.
# The reduction now runs one TRC group per worker, each mmseqs single-threaded;
# it is output-identical to the old serial multi-threaded loop ONLY IF the result
# is invariant to (a) the mmseqs thread count and (b) the worker count / order.
# This asserts both by comparing -t 1 vs -t 4 output on a self-generated dimer
# library that exercises the reduce path AND the short-monomer skip path.
# Requires mmseqs on PATH (skips cleanly if absent).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
command -v mmseqs >/dev/null 2>&1 || { echo "mmseqs not on PATH; skipping"; exit 0; }
tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT
python3 - "$tmp/lib.fasta" <<'PY'
import random, sys
random.seed(7)
mono = lambda n: "".join(random.choice("ACGT") for _ in range(n))
recs = []
m = mono(60)
for i in range(4):                    # TRC_1: rotational variants -> reducible
    rot = m[i*7:] + m[:i*7]; recs.append((f"s1_{i}#TRC_1", rot + rot))
for i in range(3):                    # TRC_3: distinct dimers
    m2 = mono(80); recs.append((f"s3_{i}#TRC_3", m2 + m2))
for i in range(3):                    # TRC_2: 10 bp monomer -> short-skip path
    ms = mono(10); recs.append((f"s2_{i}#TRC_2", ms + ms))
with open(sys.argv[1], "w") as f:
    for h, s in recs:
        f.write(f">{h}\n")
        for j in range(0, len(s), 80): f.write(s[j:j+80] + "\n")
PY
python3 "$ROOT/scripts/reduce_dimer_library.py" -i "$tmp/lib.fasta" -o "$tmp/t1.fasta" -t 1 2>/dev/null
python3 "$ROOT/scripts/reduce_dimer_library.py" -i "$tmp/lib.fasta" -o "$tmp/t4.fasta" -t 4 2>/dev/null
if diff -q "$tmp/t1.fasta" "$tmp/t4.fasta" >/dev/null; then
  echo "test_reduce_dimer_parallel: PASSED (-t1 == -t4, $(grep -c '>' "$tmp/t1.fasta") reps)"
else
  echo "test_reduce_dimer_parallel: FAILED (-t1 != -t4)"; diff "$tmp/t1.fasta" "$tmp/t4.fasta" | head; exit 1
fi
