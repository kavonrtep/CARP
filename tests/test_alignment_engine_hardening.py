#!/usr/bin/env python3
"""Unit tests for the alignment-engine hardening (PR-A: items 3, 6).

Covers only the pure-Python writer/analysis helpers — no parasail/mmseqs — so it
runs in the cheap CI unit env:
  * global_local_aln.write_results_table — atomic write, no .tmp remnant, empty=no file
  * dante_line.analyze_alignment_lengths — correctness, determinism, atomic write

The end-to-end determinism of the alignment TSV row order (item 3's sort in
run_all_vs_all_alignment) needs parasail and is validated by the fixture run.
"""
import importlib.util
import os
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS))

import global_local_aln as gla  # noqa: E402


def _load_dante_line():
    spec = importlib.util.spec_from_file_location("dante_line", SCRIPTS / "dante_line.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


dl = _load_dante_line()


def test_write_results_table_atomic_and_ordered():
    results = [
        {"query_id": "A", "ref_id": "B", "max_score": 42},
        {"query_id": "A", "ref_id": "C", "max_score": 7},
    ]
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "aln.tsv")
        gla.write_results_table(results, out)
        assert os.path.exists(out), "output not written"
        assert not os.path.exists(out + ".tmp"), "atomic write left a .tmp remnant"
        lines = open(out).read().splitlines()
        assert lines[0].split("\t")[0] == "query_id", lines[0]
        # rows written in the given order (caller sorts upstream)
        assert lines[1].startswith("A\tB\t"), lines[1]
        assert lines[2].startswith("A\tC\t"), lines[2]

    # empty results -> no file at all (unchanged semantics; nothing to be partial)
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "empty.tsv")
        gla.write_results_table([], out)
        assert not os.path.exists(out), "empty results should not create a file"
    print("  test_write_results_table_atomic_and_ordered: PASS")


def test_analyze_alignment_lengths_correct_deterministic_atomic():
    # Group A accrues lengths [10,8,6,4]; with N=3 the Nth-largest is 6 and the
    # number of alignments <= 6 (excluding the selected one) is 1. Groups B..E get
    # a single length each, below N, so they are not reported.
    tsv = (
        "query_id\tref_id\tdegapped_query_len\tdegapped_ref_len\n"
        "A\tB\t10\t1\n"
        "A\tC\t8\t1\n"
        "A\tD\t6\t1\n"
        "A\tE\t4\t1\n"
    )
    with tempfile.TemporaryDirectory() as d:
        aln = Path(d) / "aln.tsv"
        out = Path(d) / "len.tsv"
        aln.write_text(tsv)
        dl.analyze_alignment_lengths(aln, out, 3)

        assert not out.with_name(out.name + ".tmp").exists(), "atomic write left a .tmp remnant"
        rows = out.read_text().splitlines()
        assert rows[0] == "Group_ID\tSelected_Length\tNum_Shorter", rows[0]
        data = [r for r in rows[1:] if r]
        assert data == ["A\t6\t1"], data

        # determinism: a second run on the same input is byte-identical
        out2 = Path(d) / "len2.tsv"
        dl.analyze_alignment_lengths(aln, out2, 3)
        assert out.read_text() == out2.read_text(), "analyze output not deterministic"
    print("  test_analyze_alignment_lengths_correct_deterministic_atomic: PASS")


def test_analyze_alignment_lengths_heap_equivalence():
    """PR-B item 4: the bounded-heap rewrite must produce byte-identical output to
    the old 'store every length, sort, index' reference, across random inputs."""
    import random
    from collections import defaultdict

    def reference(rows, N):
        gl = defaultdict(list)
        for q, r, dq, dr in rows:
            gl[q].append(dq)
            gl[r].append(dr)
        out = []
        for g in sorted(gl):
            lengths = sorted(gl[g], reverse=True)
            if len(lengths) >= N:
                sel = lengths[N - 1]
                ns = sum(1 for x in lengths if x <= sel) - 1
                out.append((g, sel, ns))
        return out

    rng = random.Random(1234)
    for trial in range(300):
        N = rng.randint(1, 4)
        groups = [f"G{i}" for i in range(rng.randint(1, 6))]
        rows = [
            (rng.choice(groups), rng.choice(groups), rng.randint(1, 20), rng.randint(1, 20))
            for _ in range(rng.randint(0, 40))
        ]
        with tempfile.TemporaryDirectory() as d:
            aln = Path(d) / "a.tsv"
            out = Path(d) / "o.tsv"
            with open(aln, "w") as f:
                f.write("query_id\tref_id\tdegapped_query_len\tdegapped_ref_len\n")
                for q, r, dq, dr in rows:
                    f.write(f"{q}\t{r}\t{dq}\t{dr}\n")
            dl.analyze_alignment_lengths(aln, out, N)
            got = []
            for line in out.read_text().splitlines()[1:]:
                if not line:
                    continue
                g, sel, ns = line.split("\t")
                got.append((g, int(sel), int(ns)))
        assert got == reference(rows, N), (
            f"trial {trial} N={N}: heap != reference\ngot={got}\nexp={reference(rows, N)}\nrows={rows}"
        )
    print("  test_analyze_alignment_lengths_heap_equivalence: PASS")


def test_dante_line_exit_code_split():
    """Item 7: dante_line.py must exit 3 for the benign 'no LINE content' case and
    a non-3 code for a real error, so the Snakefile can tell them apart (and stop
    swallowing crashes into empty LINE outputs). Reaches the exit before any
    alignment, so it needs no parasail/mmseqs."""
    dl = SCRIPTS / "dante_line.py"
    with tempfile.TemporaryDirectory() as d:
        genome = Path(d) / "g.fa"
        genome.write_text(">chr1\nACGTACGT\n")
        # a valid GFF3 with features but NO Class_I|LINE -> benign "no LINE content"
        gff = Path(d) / "noline.gff3"
        gff.write_text(
            "chr1\tdante\tprotein_domain\t10\t20\t.\t+\t.\tName=GAG;Final_Classification=Class_I|LTR\n"
        )
        out = Path(d) / "out"
        r = subprocess.run(
            [sys.executable, str(dl), "-g", str(genome), "-a", str(gff), "-o", str(out)],
            capture_output=True, text=True,
        )
        assert r.returncode == 3, (
            f"expected exit 3 for no-LINE input, got {r.returncode}\n{r.stdout}\n{r.stderr}"
        )

        # a real error (missing genome) must NOT masquerade as the benign code
        r2 = subprocess.run(
            [sys.executable, str(dl), "-g", str(Path(d) / "missing.fa"),
             "-a", str(gff), "-o", str(out)],
            capture_output=True, text=True,
        )
        assert r2.returncode not in (0, 3), (
            f"missing genome should be a real error (not 0/3), got {r2.returncode}"
        )
    print("  test_dante_line_exit_code_split: PASS")


if __name__ == "__main__":
    test_write_results_table_atomic_and_ordered()
    test_analyze_alignment_lengths_correct_deterministic_atomic()
    test_analyze_alignment_lengths_heap_equivalence()
    test_dante_line_exit_code_split()
    print("test_alignment_engine_hardening: ALL PASS")
