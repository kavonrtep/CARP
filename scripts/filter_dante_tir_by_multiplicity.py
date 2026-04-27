#!/usr/bin/env python3
"""
Drop primary DANTE_TIR FASTA records whose parent ``sequence_feature`` row
in DANTE_TIR_final.gff3 has ``Multiplicity`` below a threshold.

Used by ``rule make_tir_combined_library`` when the user raises
``dante_tir_min_multiplicity`` above 1. Affects only the library FASTA;
the GFF and the unified annotation still carry every primary element.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


def load_keep_ids(gff_path: Path, threshold: int) -> set:
    keep = set()
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "sequence_feature":
                continue
            attrs = fields[8]
            m_id = re.search(r"(?:^|;)ID=([^;]*)", attrs)
            m_mul = re.search(r"(?:^|;)Multiplicity=([^;]*)", attrs)
            if not m_id:
                continue
            try:
                mul = int(m_mul.group(1)) if m_mul else 1
            except ValueError:
                mul = 1
            if mul >= threshold:
                keep.add(m_id.group(1))
    return keep


def filter_fasta(fa_in: Path, fa_out: Path, keep: set) -> tuple[int, int]:
    n_in = n_out = 0
    writing = False
    with open(fa_in) as fin, open(fa_out, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                n_in += 1
                seq_id = line[1:].split("#", 1)[0].split()[0].replace("_revcomp", "")
                writing = seq_id in keep
                if writing:
                    fout.write(line)
                    n_out += 1
            elif writing:
                fout.write(line)
    return n_in, n_out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gff", required=True, type=Path,
                        help="DANTE_TIR_final.gff3 (provides Multiplicity per ID)")
    parser.add_argument("--fasta-in", required=True, type=Path,
                        help="Input FASTA, one record per primary element")
    parser.add_argument("--fasta-out", required=True, type=Path,
                        help="Output FASTA, filtered by Multiplicity threshold")
    parser.add_argument("--min-multiplicity", required=True, type=int,
                        help="Keep records whose parent has Multiplicity >= this value")
    args = parser.parse_args()

    if args.min_multiplicity < 1:
        parser.error("--min-multiplicity must be >= 1")

    keep = load_keep_ids(args.gff, args.min_multiplicity)
    n_in, n_out = filter_fasta(args.fasta_in, args.fasta_out, keep)
    sys.stderr.write(
        f"multiplicity filter: kept {n_out}/{n_in} primary records "
        f"(threshold={args.min_multiplicity})\n"
    )


if __name__ == "__main__":
    main()
