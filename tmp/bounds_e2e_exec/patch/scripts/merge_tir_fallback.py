#!/usr/bin/env python3
"""
Merge DANTE_TIR and DANTE_TIR_FALLBACK annotations.

This script:
1. Loads primary DANTE_TIR GFF3 and fallback GFF3.
2. Removes fallback elements that have ANY coordinate overlap with a
   primary DANTE_TIR element on the same sequence.
3. Labels surviving fallback elements as partial (Status=partial,
   source=DANTE_TIR_fallback).
4. Writes a combined GFF3 (primary + non-overlapping fallback).
5. Writes a FASTA of non-overlapping fallback element sequences
   (extracted from the fallback extended FASTA) for downstream
   re-clustering with primary element sequences.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, NamedTuple, Set


class Interval(NamedTuple):
    seqname: str
    start: int
    end: int


def parse_gff3_intervals(gff3_file: str, feature_type: str = "sequence_feature") -> List[Interval]:
    """Extract intervals for parent-level features from a GFF3 file."""
    intervals = []
    with open(gff3_file) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            if fields[2] != feature_type:
                continue
            intervals.append(Interval(fields[0], int(fields[3]), int(fields[4])))
    return intervals


def parse_fallback_parent_ids(gff3_file: str, feature_type: str = "sequence_feature") -> Dict[str, Interval]:
    """Return mapping of fallback parent ID -> Interval."""
    id_to_interval = {}
    with open(gff3_file) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            if fields[2] != feature_type:
                continue
            attrs = dict(
                kv.split("=", 1) for kv in fields[8].split(";") if "=" in kv
            )
            parent_id = attrs.get("ID", "")
            if parent_id:
                id_to_interval[parent_id] = Interval(fields[0], int(fields[3]), int(fields[4]))
    return id_to_interval


def intervals_overlap(a: Interval, b: Interval) -> bool:
    """Check if two intervals on the same sequence overlap."""
    if a.seqname != b.seqname:
        return False
    return a.start <= b.end and b.start <= a.end


def find_overlapping_ids(
    primary_intervals: List[Interval],
    fallback_id_intervals: Dict[str, Interval],
) -> Set[str]:
    """Find fallback parent IDs that overlap any primary interval."""
    # Build per-sequence index for primary intervals
    primary_by_seq: Dict[str, List[Interval]] = {}
    for iv in primary_intervals:
        primary_by_seq.setdefault(iv.seqname, []).append(iv)

    overlapping = set()
    for parent_id, fb_iv in fallback_id_intervals.items():
        for prim_iv in primary_by_seq.get(fb_iv.seqname, []):
            if intervals_overlap(fb_iv, prim_iv):
                overlapping.add(parent_id)
                break

    return overlapping


def write_combined_gff(
    primary_gff: str,
    fallback_gff: str,
    output_gff: str,
    overlapping_ids: Set[str],
) -> int:
    """Write combined GFF3: all primary lines + non-overlapping fallback lines.

    Fallback lines are relabeled: source -> DANTE_TIR_fallback, Status=partial added.
    Returns number of non-overlapping fallback parent features written.
    """
    n_fallback = 0

    with open(output_gff, "w") as out:
        # Write header
        out.write("##gff-version 3\n")

        # Copy primary GFF body
        with open(primary_gff) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                out.write(line)

        # Copy non-overlapping fallback lines with relabeling
        with open(fallback_gff) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 9:
                    continue

                attrs = dict(
                    kv.split("=", 1) for kv in fields[8].split(";") if "=" in kv
                )

                # Determine the parent ID for this line
                parent_id = attrs.get("ID") if fields[2] == "sequence_feature" else attrs.get("Parent")
                if not parent_id:
                    # Try Group_ID as fallback
                    parent_id = attrs.get("Group_ID", "")

                if parent_id in overlapping_ids:
                    continue

                # Relabel source column
                fields[1] = "DANTE_TIR_fallback"

                # IDs in the fallback GFF3 already carry the "_partial" suffix
                # (embedded upstream by dante_tir_fallback.py), so no rewriting
                # of ID / Parent / Group_ID is needed here. Just add a
                # Status=partial marker to parent features.
                if fields[2] == "sequence_feature":
                    if "Status" not in attrs:
                        fields[8] = fields[8] + ";Status=partial"
                    n_fallback += 1

                out.write("\t".join(fields) + "\n")

    return n_fallback


def filter_fallback_fasta(
    fallback_fasta: str,
    output_fasta: str,
    overlapping_ids: Set[str],
) -> int:
    """Write fallback FASTA sequences excluding those with overlapping parent IDs.

    The fallback FASTA is expected to use the canonical
    ``>unique_id#classification`` naming produced upstream by
    ``dante_tir_fallback.py`` (where ``unique_id`` already ends in
    ``_partial``). This function just filters records whose parent ID
    appears in ``overlapping_ids`` and copies the rest through verbatim.

    Returns number of sequences written.
    """
    n_written = 0
    writing = False

    with open(fallback_fasta) as fin, open(output_fasta, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                header = line[1:].strip()
                seq_id = header.split("#", 1)[0].split()[0]
                # Strip optional strand-orientation suffix so the ID matches
                # the parent-ID keys loaded from the GFF3 (which do NOT carry
                # _revcomp but DO already carry _partial).
                base_id = seq_id.replace("_revcomp", "")

                if base_id in overlapping_ids:
                    writing = False
                else:
                    writing = True
                    fout.write(line)
                    n_written += 1
            elif writing:
                fout.write(line)

    return n_written


def main():
    parser = argparse.ArgumentParser(
        description="Merge DANTE_TIR and DANTE_TIR_FALLBACK, removing overlapping fallback elements"
    )
    parser.add_argument(
        "--primary-gff", required=True,
        help="Primary DANTE_TIR GFF3 (DANTE_TIR_final.gff3)",
    )
    parser.add_argument(
        "--fallback-gff", required=True,
        help="Fallback DANTE_TIR_FALLBACK GFF3",
    )
    parser.add_argument(
        "--fallback-fasta", required=True,
        help="Fallback extended element FASTA (for filtering)",
    )
    parser.add_argument(
        "--output-gff", required=True,
        help="Output combined GFF3",
    )
    parser.add_argument(
        "--output-fasta", required=True,
        help="Output filtered fallback FASTA (non-overlapping elements only)",
    )

    args = parser.parse_args()

    for path, label in [
        (args.primary_gff, "Primary GFF"),
        (args.fallback_gff, "Fallback GFF"),
        (args.fallback_fasta, "Fallback FASTA"),
    ]:
        if not Path(path).exists():
            print(f"Error: {label} file not found: {path}", file=sys.stderr)
            sys.exit(1)

    # Parse primary DANTE_TIR intervals
    primary_intervals = parse_gff3_intervals(args.primary_gff)
    print(f"Primary DANTE_TIR elements: {len(primary_intervals)}")

    # Parse fallback parent IDs and intervals
    fallback_id_intervals = parse_fallback_parent_ids(args.fallback_gff)
    print(f"Fallback elements (total): {len(fallback_id_intervals)}")

    # Find overlapping fallback IDs
    overlapping_ids = find_overlapping_ids(primary_intervals, fallback_id_intervals)
    print(f"Fallback elements overlapping primary: {len(overlapping_ids)}")
    print(f"Fallback elements retained: {len(fallback_id_intervals) - len(overlapping_ids)}")

    # Write combined GFF3
    n_fallback_written = write_combined_gff(
        args.primary_gff, args.fallback_gff, args.output_gff, overlapping_ids,
    )
    print(f"Combined GFF3 written: {args.output_gff} ({len(primary_intervals)} primary + {n_fallback_written} fallback)")

    # Filter fallback FASTA
    n_fasta_written = filter_fallback_fasta(
        args.fallback_fasta, args.output_fasta, overlapping_ids,
    )
    print(f"Filtered fallback FASTA: {args.output_fasta} ({n_fasta_written} sequences)")


if __name__ == "__main__":
    main()