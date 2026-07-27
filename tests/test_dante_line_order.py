#!/usr/bin/env python3
"""Item 8 guard: LINE pattern/FASTA ordering must not depend on PYTHONHASHSEED.

mmseqs clustering/grouping (used by both the dante_line and dante_tir_fallback
rules via the shared engine) is order-sensitive — the same sequences in a
different order yield different clusters. The order of the FASTA handed to mmseqs
comes from dante_line's GFF3 -> feature -> pattern path. That path is deterministic
by construction today (position-sorted parse, insertion-ordered dict, a
membership-only `set`), but nothing asserted it. If a `set` ever leaked into the
order, clusters would silently differ between identical runs (Python randomises
string hashing per process), which is exactly the kind of bug nobody thinks to
check for.

This runs the parse -> get_line_features -> find_domain_patterns path in separate
interpreters under different PYTHONHASHSEED values and asserts byte-identical
ordered output. Data-free, no parasail/mmseqs — runs in the cheap CI env.
Mirrors DANTE_TIR's tests/test_aa_fasta_order.py.
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"

# chr1 (+): ENDO,RT -> a 2-feature pattern; chr2 (+): ENDO,RT,RH -> a 3-feature
# pattern. Two groups across two seqnames so the ORDER of the emitted patterns is
# actually exercised (not a single trivial group).
GFF3 = (
    "chr1\tdante\tprotein_domain\t1000\t1100\t.\t+\t.\tName=ENDO;Final_Classification=Class_I|LINE\n"
    "chr1\tdante\tprotein_domain\t1150\t1250\t.\t+\t.\tName=RT;Final_Classification=Class_I|LINE\n"
    "chr2\tdante\tprotein_domain\t100\t200\t.\t+\t.\tName=ENDO;Final_Classification=Class_I|LINE\n"
    "chr2\tdante\tprotein_domain\t250\t350\t.\t+\t.\tName=RT;Final_Classification=Class_I|LINE\n"
    "chr2\tdante\tprotein_domain\t400\t500\t.\t+\t.\tName=RH;Final_Classification=Class_I|LINE\n"
)

SNIPPET = r"""
import sys, importlib.util
scripts, gff = sys.argv[1], sys.argv[2]
spec = importlib.util.spec_from_file_location("dante_line", scripts + "/dante_line.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
patterns = m.find_domain_patterns(m.get_line_features(m.parse_gff3_features(gff)), 10000)
for p in patterns:
    lo, hi = p.get_region_bounds()
    print(p.group_id, p.seqname, p.strand, p.pattern_type, lo, hi)
"""


def _run(gff, seed):
    env = dict(os.environ, PYTHONHASHSEED=str(seed))
    r = subprocess.run([sys.executable, "-c", SNIPPET, str(SCRIPTS), gff],
                       env=env, capture_output=True, text=True)
    assert r.returncode == 0, f"seed={seed} failed:\n{r.stderr}"
    return r.stdout


def test_pattern_order_is_hashseed_invariant():
    with tempfile.TemporaryDirectory() as d:
        gff = os.path.join(d, "line.gff3")
        Path(gff).write_text(GFF3)
        outputs = [_run(gff, s) for s in (0, 1, 42, 12345)]
        # the input must actually produce patterns, else the guard is vacuous
        assert outputs[0].strip(), "test input produced no patterns — not exercising ordering"
        assert len(outputs[0].splitlines()) == 2, f"expected 2 patterns, got:\n{outputs[0]}"
        for seed, o in zip((1, 42, 12345), outputs[1:]):
            assert o == outputs[0], (
                f"pattern order varies with PYTHONHASHSEED (seed 0 vs {seed}):\n"
                f"{outputs[0]!r}\n{o!r}"
            )
    print("  test_pattern_order_is_hashseed_invariant: PASS")


if __name__ == "__main__":
    test_pattern_order_is_hashseed_invariant()
    print("test_dante_line_order: ALL PASS")
