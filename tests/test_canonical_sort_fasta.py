#!/usr/bin/env python3
"""Determinism test for canonical_sort_fasta.py.

The greedy clustering steps CARP drives (mmseqs easy-cluster, CAP3) are
order-sensitive: the same sequences in a different order give different
representatives and a different cluster count, which made the LINE/TIR/LTR
libraries — and the RepeatMasker annotation — vary run-to-run. The fix sorts
every clustering input into one canonical order first. This test proves the
sorter delivers that guarantee: the output is invariant to input order,
loss-less, and totally ordered (so it is a deterministic function of the input
set).

This is the CARP-side half of the reproducibility fix. The end-to-end
"shuffle -> identical mmseqs representatives" behaviour needs mmseqs and is
covered by the fixture pipeline runs; here we test the ordering primitive.

Run: python3 tests/test_canonical_sort_fasta.py
"""
import os
import random
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from canonical_sort_fasta import sort_fasta_by_sequence  # noqa: E402


def _write_fasta(path, records):
    """records: list of (header_without_gt, sequence)."""
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(">" + header + "\n")
            # Vary the wrapping deliberately so the test doesn't assume a layout.
            width = random.choice((0, 60, 70, 80))
            if width and len(seq) > width:
                for i in range(0, len(seq), width):
                    fh.write(seq[i:i + width] + "\n")
            else:
                fh.write(seq + "\n")


def _read_fasta(path):
    """Return list of (header_without_gt, sequence)."""
    out = []
    header = None
    parts = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    out.append((header, "".join(parts)))
                header = line[1:].rstrip("\n")
                parts = []
            else:
                parts.append(line.strip())
    if header is not None:
        out.append((header, "".join(parts)))
    return out


def _random_records(n, rng):
    recs = []
    for i in range(n):
        length = rng.randint(1, 400)
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        recs.append((f"elem_{i}#Class_I/LTR/Ty1_copia/Angela", seq))
    return recs


def _sort(records, tmpdir, tag):
    src = os.path.join(tmpdir, f"in_{tag}.fasta")
    dst = os.path.join(tmpdir, f"out_{tag}.fasta")
    _write_fasta(src, records)
    sort_fasta_by_sequence(src, dst, threads=1, tmpdir=tmpdir)
    return dst


def test_order_invariance():
    """Two different shuffles of the same records -> byte-identical output."""
    rng = random.Random(1234)
    records = _random_records(500, rng)
    with tempfile.TemporaryDirectory() as td:
        a = records[:]
        b = records[:]
        rng.shuffle(a)
        rng.shuffle(b)
        out_a = _sort(a, td, "a")
        out_b = _sort(b, td, "b")
        with open(out_a) as fa, open(out_b) as fb:
            assert fa.read() == fb.read(), "sort is not order-invariant"


def test_lossless_and_total_order():
    """No record lost/duplicated; output is sorted by (sequence, header)."""
    rng = random.Random(99)
    records = _random_records(300, rng)
    # Add exact-duplicate sequences with different headers to exercise the
    # header tie-break (must still be a deterministic total order).
    dup_seq = "ACGTACGTACGTACGT"
    records += [(f"dup_{k}#Class_I/LINE", dup_seq) for k in range(5)]
    with tempfile.TemporaryDirectory() as td:
        shuffled = records[:]
        rng.shuffle(shuffled)
        out = _read_fasta(_sort(shuffled, td, "x"))
        assert sorted(out) == sorted(records), "records lost/changed by sort"
        keys = [(seq, header) for header, seq in out]
        assert keys == sorted(keys), "output not in canonical (sequence, header) order"


def test_empty_input():
    with tempfile.TemporaryDirectory() as td:
        src = os.path.join(td, "empty.fasta")
        dst = os.path.join(td, "empty_out.fasta")
        open(src, "w").close()
        sort_fasta_by_sequence(src, dst, tmpdir=td)
        assert os.path.getsize(dst) == 0, "empty input must give empty output"


def main():
    test_order_invariance()
    test_lossless_and_total_order()
    test_empty_input()
    print("OK  test_canonical_sort_fasta: order-invariant, lossless, total-order, empty-safe")


if __name__ == "__main__":
    main()
