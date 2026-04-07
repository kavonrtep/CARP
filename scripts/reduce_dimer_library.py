#!/usr/bin/env python3
"""
Reduce a TideCluster dimer library by clustering sequences within each TRC group
independently using mmseqs2 easy-linclust.

Clustering is performed per-TRC to guarantee that every TRC cluster retains at
least one representative in the output, regardless of cross-TRC similarity.

Sequence headers follow the RepeatMasker naming convention:
    >name#classification
e.g. >TRC_1#TRC_1

Duplicate headers (same TRC, different phase-rotations of the dimer) are handled
by temporarily prepending a numeric counter before clustering, then restoring the
original names in the output.
"""
import argparse
import os
import subprocess
import sys
import tempfile
from collections import defaultdict


def parse_fasta(path):
    """Yield (header_without_gt, sequence) pairs from a FASTA file."""
    header = None
    parts = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:]
                parts = []
            else:
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def write_fasta(records, path):
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def trc_id_of(header):
    """Return the classification part (after '#') of a FASTA header."""
    idx = header.find("#")
    return header[idx + 1:] if idx != -1 else header


def run_mmseqs(input_fasta, output_prefix, tmp_dir, min_seq_id, coverage, threads):
    cmd = [
        "mmseqs", "easy-linclust",
        input_fasta,
        output_prefix,
        tmp_dir,
        "--min-seq-id", str(min_seq_id),
        "--cov-mode", "1",
        "-c", str(coverage),
        "--threads", str(threads),
        "-v", "1",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.stderr.write(result.stderr)
        raise RuntimeError(f"mmseqs2 failed for {input_fasta}")
    sys.stderr.write(result.stderr)


def reduce_trc_group(trc_id, records, work_dir, min_seq_id, coverage, threads):
    """
    Cluster one TRC group.  Returns list of (original_header, seq) representatives.
    If the group has only one sequence, return it unchanged.
    """
    if len(records) <= 1:
        return records

    # Uniquify headers: prepend "N__" so mmseqs2 sees distinct IDs
    unique_records = [(f"{i + 1}__{hdr}", seq) for i, (hdr, seq) in enumerate(records)]

    safe = trc_id.replace("/", "_")
    input_fasta  = os.path.join(work_dir, f"{safe}_input.fasta")
    output_prefix = os.path.join(work_dir, f"{safe}_clust")
    mmseqs_tmp   = os.path.join(work_dir, f"{safe}_mmseqs_tmp")
    rep_fasta    = output_prefix + "_rep_seq.fasta"

    write_fasta(unique_records, input_fasta)
    run_mmseqs(input_fasta, output_prefix, mmseqs_tmp, min_seq_id, coverage, threads)

    # Read representatives; strip the "N__" prefix to restore original header
    representatives = []
    for uniq_hdr, seq in parse_fasta(rep_fasta):
        original_hdr = uniq_hdr.split("__", 1)[1]
        representatives.append((original_hdr, seq))
    return representatives


def main():
    parser = argparse.ArgumentParser(
        description="Reduce TideCluster dimer library via per-TRC mmseqs2 clustering."
    )
    parser.add_argument("-i", "--input",     required=True,  help="Input dimer library FASTA")
    parser.add_argument("-o", "--output",    required=True,  help="Output reduced library FASTA")
    parser.add_argument("-t", "--threads",   type=int, default=4, help="Threads for mmseqs2")
    parser.add_argument("--min_seq_id",      type=float, default=0.8,
                        help="Minimum sequence identity for clustering (default: 0.8)")
    parser.add_argument("--coverage",        type=float, default=0.8,
                        help="Bidirectional coverage threshold (default: 0.8)")
    args = parser.parse_args()

    # Group sequences by TRC classification
    groups = defaultdict(list)
    for header, seq in parse_fasta(args.input):
        groups[trc_id_of(header)].append((header, seq))

    n_input = sum(len(v) for v in groups.values())
    sys.stderr.write(
        f"Reducing dimer library: {n_input} sequences in {len(groups)} TRC groups\n"
    )

    with tempfile.TemporaryDirectory() as work_dir:
        output_records = []
        for trc_id in sorted(groups):
            records = groups[trc_id]
            reps = reduce_trc_group(
                trc_id, records, work_dir, args.min_seq_id, args.coverage, args.threads
            )
            sys.stderr.write(
                f"  {trc_id}: {len(records)} → {len(reps)} sequences\n"
            )
            output_records.extend(reps)

    write_fasta(output_records, args.output)
    sys.stderr.write(
        f"Dimer library: {n_input} → {len(output_records)} sequences after per-TRC reduction\n"
    )


if __name__ == "__main__":
    main()
