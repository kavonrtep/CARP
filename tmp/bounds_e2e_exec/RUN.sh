#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# End-to-end validation of the new LINE boundary bounds.
#
#   ./RUN.sh --check     verify image, patch, genomes, baselines. RUNS NOTHING.
#   ./RUN.sh             run both genomes, then compare against their baselines
#   ./RUN.sh --compare   re-print the comparison without running
#   ./RUN.sh <tag>       run one genome only
#
# WHAT IS BEING TESTED. The shipped 1.7.0rc1 container caps a LINE's flank
# extension at 2000 bp (5') and 800 bp (3'). Measured on 218 LINE loci across 7
# genomes where BOTH ends were observed, those caps truncate real element: the
# 3' extension is not a stable distance at all (it varies 13.8x) because it is
# only what is left over after the domain annotation stops. What IS stable is
# the span from the START of the domain core to the element's 3' end (1.16x).
# So main replaces the flat 3' cap with a 4500 bp bound on that span, and raises
# the 5' cap to 3400 bp.
#
# This is the only behaviour change on main that has never been through a full
# pipeline run. A local A/B on one small genome moved 4 of 343 elements
# (+0.25% extension), so the expected effect is SMALL. A large move here means
# something is wrong, not that the change is working.
#
# HOW. The container is 1.7.0rc1 and does not carry the change. Rather than
# rebuild it, the patched scripts/ and classification_vocabulary.yaml are
# bind-mounted over /opt/pipeline/. The Snakefile is IDENTICAL between the tag
# and main, so the pipeline structure is untouched -- only these files differ.
# The run goes through the container runscript (run_pipeline.py), not bare
# snakemake, because provenance-gated code is skipped otherwise.
#
# Settings and invocation mirror the baseline runs exactly (same four
# carp_params), so the only difference is the bounds.
# ---------------------------------------------------------------------------
set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS=/nfsroot/projects/darwin/runs
SIF="${SIF:-$RUNS/tmp3/image/assembly_repeat_annotation_pipeline.sif}"
THREADS="${THREADS:-32}"
SCRATCH="${SCRATCH:-$HERE/scratch}"
RUNTIME="$(command -v singularity || command -v apptainer || true)"

list_jobs() { grep -v '^#' "$HERE/genomes.tsv" | grep -v '^[[:space:]]*$'; }

check() {
  local bad=0
  echo "== container =="
  [ -n "$RUNTIME" ] && echo "  runtime: $RUNTIME" || { echo "  NO singularity/apptainer"; bad=1; }
  [ -f "$SIF" ] && echo "  image:   $SIF" || { echo "  MISSING image: $SIF"; bad=1; }
  echo "== patch payload =="
  for f in patch/classification_vocabulary.yaml patch/scripts/dante_line.py \
           patch/scripts/classification.py patch/scripts/global_local_aln.py; do
    [ -r "$HERE/$f" ] && echo "  ok      $f" || { echo "  MISSING $f"; bad=1; }
  done
  if grep -q "max_core_to_3prime_end" "$HERE/patch/classification_vocabulary.yaml" 2>/dev/null; then
    echo "  ok      vocabulary carries max_core_to_3prime_end (the change under test)"
  else
    echo "  PATCH IS STALE: no max_core_to_3prime_end in the vocabulary"; bad=1
  fi
  echo "== genomes and baselines =="
  while IFS=$'\t' read -r tag asm base why; do
    [ -n "${tag:-}" ] || continue
    [ -r "$asm" ] && echo "  ok      $tag assembly" || { echo "  MISSING $tag assembly: $asm"; bad=1; }
    local b="$RUNS/$base/carp_output/summary_statistics.csv"
    [ -r "$b" ] && echo "  ok      $tag baseline $base" || { echo "  MISSING $tag baseline: $b"; bad=1; }
  done < <(list_jobs)
  echo "== scripts load inside the container =="
  if [ -n "$RUNTIME" ] && [ -f "$SIF" ]; then
    if "$RUNTIME" exec -B /nfsroot:/nfsroot -B "$HERE:$HERE" \
         -B "$HERE/patch/scripts:/opt/pipeline/scripts" \
         -B "$HERE/patch/classification_vocabulary.yaml:/opt/pipeline/classification_vocabulary.yaml" \
         "$SIF" bash -lc 'dante_line.py --help' 2>&1 | grep -q -- --prefilter-kmer; then
      echo "  ok      patched dante_line.py loads and is the new version"
    else
      echo "  FAILED  patched scripts do not load in the container"; bad=1
    fi
  fi
  [ $bad -eq 0 ] && echo "ALL CHECKS PASSED" || echo "CHECKS FAILED — fix the above before running"
  return $bad
}

run_one() {   # $1=tag $2=assembly
  local tag="$1" asm="$2"
  local out="$HERE/results/$tag"
  if [ -f "$HERE/results/DONE_$tag" ]; then echo "  skip $tag (done)"; return 0; fi
  mkdir -p "$out" "$SCRATCH/$tag/work" "$SCRATCH/$tag/tmp"
  # config identical to the baseline runs' four carp_params
  cat > "$SCRATCH/$tag/carp_config.yaml" <<CFG
genome_fasta: $asm
output_dir: $out
repeatmasker_sensitivity: rush
repeatmasker_culling_limit: 2
tidecluster_reannotate_culling_limit: 2
reduce_library: True
cleanup_intermediates: maximal
CFG
  echo "  run  $tag  ($THREADS threads) -> results/$tag.log"
  "$RUNTIME" run \
    -B "$HERE" -B "$(dirname "$asm")" -B "$SCRATCH" \
    -B "$HERE/patch/scripts:/opt/pipeline/scripts" \
    -B "$HERE/patch/classification_vocabulary.yaml:/opt/pipeline/classification_vocabulary.yaml" \
    --pwd "$SCRATCH/$tag/work" --env TMPDIR="$SCRATCH/$tag/tmp" \
    "$SIF" -c "$SCRATCH/$tag/carp_config.yaml" -t "$THREADS" \
    > "$HERE/results/$tag.log" 2>&1
  local rc=$?
  echo $rc > "$HERE/results/DONE_$tag"
  if [ $rc -ne 0 ]; then
    echo "  !! $tag FAILED (rc=$rc); last lines:"
    tail -15 "$HERE/results/$tag.log" | sed 's/^/     /'
    rm -f "$HERE/results/DONE_$tag"
  fi
}

compare() {
  echo
  echo "================= NEW BOUNDS vs 1.7.0rc1 BASELINE ================="
  while IFS=$'\t' read -r tag asm base why; do
    [ -n "${tag:-}" ] || continue
    local new="$HERE/results/$tag" old="$RUNS/$base/carp_output"
    [ -f "$new/summary_statistics.csv" ] || { echo "### $tag: not finished"; continue; }
    echo
    echo "### $tag   (baseline $base)"
    echo "    $why"
    python3 - "$old" "$new" <<'PY'
import csv, sys, os
old, new = sys.argv[1], sys.argv[2]
def load(p):
    d = {}
    f = os.path.join(p, "summary_statistics.csv")
    if not os.path.exists(f): return d
    with open(f) as fh:
        r = csv.reader(fh, delimiter="\t"); next(r, None)
        for row in r:
            if len(row) >= 3: d[row[0]] = float(row[2])
    return d
a, b = load(old), load(new)
line_a, line_b = a.get("Class_I/LINE", 0.0), b.get("Class_I/LINE", 0.0)
tot_a, tot_b = sum(a.values()), sum(b.values())
print(f"    Class_I/LINE     {line_a:6.2f}% -> {line_b:6.2f}%   ({line_b-line_a:+.2f} pts)")
print(f"    total annotated  {tot_a:6.2f}% -> {tot_b:6.2f}%   ({tot_b-tot_a:+.2f} pts)")
moved = sorted(((abs(b.get(k,0)-a.get(k,0)), k) for k in set(a)|set(b)), reverse=True)[:5]
print("    largest class moves:")
for d_, k in moved:
    if d_ < 0.01: continue
    print(f"      {k[:48]:<48} {a.get(k,0):6.2f} -> {b.get(k,0):6.2f}  ({b.get(k,0)-a.get(k,0):+.2f})")
def health(p, key):
    f = os.path.join(p, "Libraries", "library_health.tsv")
    if not os.path.exists(f): return "n/a"
    for l in open(f):
        c = l.rstrip("\n").split("\t")
        if len(c) >= 4 and c[0] == "boundary" and c[1] == "DANTE_LINE" and c[2] == key:
            return c[3]
    return "n/a"
for k in ("n_elements", "n_capped", "extension_bp", "extension_fraction", "n_at_flank_ceiling"):
    print(f"    {k:<20} {health(old,k):>12} -> {health(new,k):>12}")
PY
  done < <(list_jobs)
  cat <<'TXT'

HOW TO READ THIS. The change raises the 5' cap and replaces the flat 3' cap with
a span bound, so LINE elements get LONGER. Expect Class_I/LINE to rise SLIGHTLY.
A local A/B moved 4 of 343 elements (+0.25% extension), so single-digit-percent
moves are the expectation.

  * a small LINE rise, total annotated roughly flat   -> as designed
  * a large LINE rise, or LTR/Tandem falling to pay   -> INVESTIGATE. That is
    the displacement pattern the whole boundary investigation started from.
  * n_capped falling sharply is expected (the caps bind less often now).
  * n_at_flank_ceiling rising above ~0 is a red flag: extensions running to the
    10 kb search limit is the runaway signature.
TXT
}

case "${1:---run}" in
  --check)   check ;;
  --compare) compare ;;
  --run|"")  check || exit 1
             while IFS=$'\t' read -r tag asm base why; do
               [ -n "${tag:-}" ] || continue
               run_one "$tag" "$asm"
             done < <(list_jobs)
             compare ;;
  *)         check || exit 1
             while IFS=$'\t' read -r tag asm base why; do
               [ "$tag" = "$1" ] || continue
               run_one "$tag" "$asm"
             done < <(list_jobs)
             compare ;;
esac
