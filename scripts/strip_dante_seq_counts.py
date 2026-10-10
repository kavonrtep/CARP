#!/usr/bin/env python3
"""Remove DANTE's per-sequence SEQ/DOMAIN/COUNTS table from GFF3 headers.

dante_gff_output_filtering.py writes, between a ``##SEQ<TAB>DOMAIN<TAB>COUNTS``
header and the next ``##-----`` separator, one ``##<seqid><TAB><domain><TAB><n>``
comment row per sequence x domain. On a fragmented assembly that is ~1 M rows
(issue #5), and rtracklayer's import.gff3 is quadratic in comment lines (6.2 h
on that genome). Nothing in CARP reads the table. Every other line — the other
header lines and all features — is copied unchanged.

Only a table in the HEADER costs time: rtracklayer reads pragmas up to the
first feature line, so the same rows after a feature are free (0.13 s vs 127 s
for 160k rows) — which is why concatenating readers such as
merge_repeat_annotations.R (RepeatMasker GFF3 first) were never affected.

Usage: strip_dante_seq_counts.py IN.gff3 [IN2.gff3 ...] -o OUT.gff3
Several inputs are concatenated in order.
Also importable: ``strip_lines(iterable_of_lines)`` yields the kept lines.
"""
import argparse
import os
import re

TABLE_HEADER = "##SEQ\tDOMAIN\tCOUNTS"
SEPARATOR = re.compile(r"^##-+$")


def strip_lines(lines):
    """Yield every line except the SEQ/DOMAIN/COUNTS block (header, rows and
    its closing separator). Only comment lines can be dropped: the block is
    also closed by the first non-comment line, so a feature is never lost."""
    in_table = False
    for line in lines:
        bare = line.rstrip("\r\n")
        if in_table:
            if SEPARATOR.match(bare):
                in_table = False
                continue
            if bare.startswith("#"):
                continue
            in_table = False          # unterminated table: keep the feature
        elif bare == TABLE_HEADER:
            in_table = True
            continue
        yield line


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("inputs", nargs="+", help="input GFF3 file(s)")
    ap.add_argument("-o", "--output", required=True, help="output GFF3")
    args = ap.parse_args()
    tmp = args.output + ".tmp"
    with open(tmp, "w") as fout:
        for path in args.inputs:
            with open(path) as fin:
                fout.writelines(strip_lines(fin))
    os.replace(tmp, args.output)


if __name__ == "__main__":
    main()
