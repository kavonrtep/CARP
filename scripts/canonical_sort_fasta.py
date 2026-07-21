#!/usr/bin/env python3
"""Canonically sort a FASTA by sequence content — deterministic, out-of-core.

Why
---
The greedy clustering steps CARP drives directly (``mmseqs easy-cluster`` and
CAP3) are **order-sensitive**: the same set of sequences presented in a
different order yields a different set of representative consensi *and* a
different cluster count (measured: ~19% of mmseqs reps and the cluster count
change under a pure shuffle; ~3.7% for CAP3). Upstream tools can hand CARP the
same sequences in a run-varying order (e.g. DANTE's ~0.1% boundary jitter, or
environment-sensitive chunk grouping), and the order-sensitive clustering then
amplifies that into large run-to-run library churn.

Sorting every clustering input into one canonical order — by sequence content,
so it is invariant to upstream coordinate/ID reordering — makes the clustering a
deterministic function of the input *set*. Verified: sorting a shuffled input
before ``mmseqs easy-cluster`` reproduces byte-identical representatives.

Performance
-----------
The sort is delegated to GNU ``sort`` with ``LC_ALL=C`` (same approach TideCluster
uses for its deterministic comparative analysis). GNU sort is disk-backed: the
in-RAM footprint is capped by ``--buffer`` and overflow spills to ``--tmpdir``,
so this stays safe on the multi-GB intermediates produced by 30-90 Gbp genomes —
it never loads the library into memory. The work is O(n log n) over bytes the
clustering step reads anyway, so it is a small fraction of clustering wall-time.

Usage
-----
    canonical_sort_fasta.py IN.fasta OUT.fasta [--threads N] [--buffer 2G] [--tmpdir DIR]
    canonical_sort_fasta.py -  -                       # stdin -> stdout

The ordering key is (sequence, header), giving a total, reproducible order.
Empty input yields empty output. Sequences are re-wrapped at 60 columns.
"""
import argparse
import os
import subprocess
import sys
import tempfile


def _iter_records(fh):
    """Yield (header_without_gt, sequence) for each FASTA record, streaming."""
    header = None
    seq_parts = []
    for line in fh:
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq_parts)
            header = line[1:].rstrip("\n")
            seq_parts = []
        else:
            seq_parts.append(line.strip())
    if header is not None:
        yield header, "".join(seq_parts)


def sort_fasta_by_sequence(input_path, output_path, threads=1, buffer="1G",
                           tmpdir=None):
    """Write ``input_path`` to ``output_path`` sorted canonically by sequence.

    Uses a linearize -> GNU sort (out-of-core) -> delinearize pipeline. ``tmpdir``
    defaults to $TMPDIR then the output directory so large spills land on scratch.
    """
    in_fh = sys.stdin if input_path == "-" else open(input_path)
    if tmpdir is None:
        tmpdir = os.environ.get("TMPDIR")
    if tmpdir is None:
        tmpdir = (os.path.dirname(os.path.abspath(output_path))
                  if output_path != "-" else ".")
    os.makedirs(tmpdir, exist_ok=True)

    # Pass 1: linearize to "SEQ\tHEADER" (streaming, bounded memory).
    lin = tempfile.NamedTemporaryFile("w", dir=tmpdir, delete=False,
                                      suffix=".canonsort.lin")
    n = 0
    try:
        with in_fh:
            for header, seq in _iter_records(in_fh):
                # FASTA headers/sequences contain no tabs; guard defensively so a
                # stray tab can never corrupt the field split.
                if "\t" in header:
                    header = header.replace("\t", " ")
                lin.write(seq)
                lin.write("\t")
                lin.write(header)
                lin.write("\n")
                n += 1
    finally:
        lin.close()

    out_fh = sys.stdout if output_path == "-" else open(output_path, "w")
    try:
        if n == 0:
            # Empty input -> empty output; nothing to sort.
            return
        # Pass 2: sort by sequence (k1) then header (k2) -> total order.
        sorted_path = lin.name + ".sorted"
        env = dict(os.environ, LC_ALL="C")
        with open(sorted_path, "w") as sf:
            subprocess.check_call(
                ["sort", "-t", "\t", "-k1,1", "-k2,2",
                 "-S", str(buffer), "--parallel", str(max(1, int(threads))),
                 "-T", tmpdir, lin.name],
                stdout=sf, env=env,
            )
        # Pass 3: delinearize -> wrapped FASTA (streaming).
        with open(sorted_path) as sf:
            for line in sf:
                seq, _, header = line.rstrip("\n").partition("\t")
                out_fh.write(">")
                out_fh.write(header)
                out_fh.write("\n")
                for i in range(0, len(seq), 60):
                    out_fh.write(seq[i:i + 60])
                    out_fh.write("\n")
        os.unlink(sorted_path)
    finally:
        if out_fh is not sys.stdout:
            out_fh.close()
        try:
            os.unlink(lin.name)
        except OSError:
            pass


def main():
    ap = argparse.ArgumentParser(
        description="Canonically sort a FASTA by sequence content (out-of-core).")
    ap.add_argument("input", help="input FASTA ('-' for stdin)")
    ap.add_argument("output", help="output FASTA ('-' for stdout)")
    ap.add_argument("--threads", type=int, default=1,
                    help="threads for GNU sort --parallel (default: 1)")
    ap.add_argument("--buffer", default="1G",
                    help="GNU sort -S in-RAM buffer before spilling (default: 1G)")
    ap.add_argument("--tmpdir", default=None,
                    help="spill directory for GNU sort (default: $TMPDIR or output dir)")
    args = ap.parse_args()
    sort_fasta_by_sequence(args.input, args.output, threads=args.threads,
                           buffer=args.buffer, tmpdir=args.tmpdir)


if __name__ == "__main__":
    main()
