#!/usr/bin/env python3
"""
Drop primary DANTE_TIR FASTA records whose parent ``sequence_feature`` row
in DANTE_TIR_final.gff3 has ``Multiplicity`` below a threshold.

Used by ``rule make_tir_combined_library`` when the user raises
``dante_tir_min_multiplicity`` above 1. Affects only the library FASTA;
the GFF and the unified annotation still carry every primary element.

Join key
--------
DANTE_TIR spells the same element two different ways -- the FASTA name is
the GFF3 ``ID`` with the ``Class_II_Subclass_1_TIR_`` prefix stripped
(``detect_tirs.R``, ``names(tir_seqs)``)::

    GFF3    ID=Class_II_Subclass_1_TIR_EnSpm_CACTA_1611
    FASTA   >EnSpm_CACTA_1611#Class_II/Subclass_1/TIR/EnSpm_CACTA

so the FASTA name is never equal to the GFF3 ID. DANTE_TIR >= 0.3.1
publishes the FASTA name directly as ``Name=``; key on that when present
and fall back to stripping the prefix off ``ID=`` so older DANTE_TIR
outputs keep working.

Getting this wrong is silent -- the join simply matches nothing and the
library comes out empty -- so a run where no record joins is treated as a
bug and exits non-zero rather than emitting an empty FASTA.

Multiplicity is deliberately NOT defaulted to 1 when the attribute is
absent. DANTE_TIR < 0.3.1 ran ``cluster_tir_sequences()`` before
``round4()``, so round-4 elements carry no ``Multiplicity`` at all (25% of
elements on a real run); defaulting to 1 dropped every one of them at any
threshold >= 2, for a reason unrelated to their copy number. They are
still dropped by default, but counted and reported, and
``--keep-missing-multiplicity`` keeps them instead.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ID_PREFIX = "Class_II_Subclass_1_TIR_"


def fasta_name_from_attrs(attrs: str) -> str | None:
    """FASTA name for a GFF3 attribute string, or None if it has no ID.

    Prefers ``Name=`` (DANTE_TIR >= 0.3.1, which publishes the FASTA name
    verbatim); falls back to stripping the classification prefix off ``ID=``.
    """
    m_name = re.search(r"(?:^|;)Name=([^;]*)", attrs)
    if m_name and m_name.group(1):
        return m_name.group(1)
    m_id = re.search(r"(?:^|;)ID=([^;]*)", attrs)
    if not m_id:
        return None
    return re.sub(rf"^{re.escape(ID_PREFIX)}", "", m_id.group(1))


def load_multiplicity(gff_path: Path) -> dict:
    """Map FASTA name -> Multiplicity, or None where the attribute is absent."""
    mult = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "sequence_feature":
                continue
            attrs = fields[8]
            name = fasta_name_from_attrs(attrs)
            if name is None:
                continue
            m_mul = re.search(r"(?:^|;)Multiplicity=([^;]*)", attrs)
            try:
                value = int(m_mul.group(1)) if m_mul else None
            except ValueError:
                value = None
            mult[name] = value
    return mult


def filter_fasta(fa_in: Path, fa_out: Path, mult: dict, threshold: int,
                 keep_missing: bool) -> tuple[int, int, int, int]:
    """Returns (n_in, n_out, n_joined, n_missing_multiplicity)."""
    n_in = n_out = n_joined = n_missing = 0
    writing = False
    with open(fa_in) as fin, open(fa_out, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                n_in += 1
                seq_id = line[1:].split("#", 1)[0].split()[0].replace("_revcomp", "")
                if seq_id not in mult:
                    writing = False
                    continue
                n_joined += 1
                value = mult[seq_id]
                if value is None:
                    n_missing += 1
                    writing = keep_missing
                else:
                    writing = value >= threshold
                if writing:
                    fout.write(line)
                    n_out += 1
            elif writing:
                fout.write(line)
    return n_in, n_out, n_joined, n_missing


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
    parser.add_argument("--keep-missing-multiplicity", action="store_true",
                        help="Keep elements whose parent row carries no Multiplicity "
                             "attribute (DANTE_TIR < 0.3.1 round-4 elements). "
                             "Default: drop them, but report the count.")
    args = parser.parse_args()

    if args.min_multiplicity < 1:
        parser.error("--min-multiplicity must be >= 1")

    mult = load_multiplicity(args.gff)
    n_in, n_out, n_joined, n_missing = filter_fasta(
        args.fasta_in, args.fasta_out, mult, args.min_multiplicity,
        args.keep_missing_multiplicity)

    sys.stderr.write(
        f"multiplicity filter: kept {n_out}/{n_in} primary records "
        f"(threshold={args.min_multiplicity})\n"
        f"  GFF3 sequence_feature rows : {len(mult)}\n"
        f"  joined to a GFF3 parent    : {n_joined}\n"
        f"  no GFF3 parent found       : {n_in - n_joined}\n"
        f"  no Multiplicity attribute  : {n_missing} "
        f"({'kept' if args.keep_missing_multiplicity else 'dropped'})\n"
    )

    if n_missing and not args.keep_missing_multiplicity:
        sys.stderr.write(
            f"WARNING: {n_missing} element(s) carry no Multiplicity attribute and "
            f"were dropped regardless of copy number. DANTE_TIR >= 0.3.1 populates "
            f"it for every element; on older versions round-4 elements never get "
            f"one. Use --keep-missing-multiplicity to keep them.\n"
        )

    # A join that matches nothing is a bug, not an empty result. Fail loudly so
    # the calling rule cannot touch an empty library and exit 0.
    if n_in and not n_joined:
        sys.exit(
            "ERROR: no FASTA record joined to a GFF3 parent. Expected the FASTA "
            "name to equal Name= (DANTE_TIR >= 0.3.1) or ID= with the "
            f"'{ID_PREFIX}' prefix stripped. Check {args.gff} against {args.fasta_in}."
        )


if __name__ == "__main__":
    main()
