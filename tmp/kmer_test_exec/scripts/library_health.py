#!/usr/bin/env python3
"""
Summarise the health of the repeat library and the inferred element boundaries.

Why
---
Two defects shipped in CARP for several releases and were invisible in every
output the pipeline produced:

  * the DANTE_TIR library was **empty** on every run (a GFF3<->FASTA join bug),
    so the RepeatMasker library held no Class_II sequences at all;
  * DANTE_LINE built 16-22 kb "LINE" consensi whose flanks were a different,
    far more abundant repeat, and RepeatMasker then masked those flanks
    genome-wide under the wrong label.

Neither was hard to see once someone looked at the library; nothing in the run
invited anyone to look. This writes the numbers that would have made both
obvious, as a TSV that ships with the run and a set of warnings on stderr.

It never fails a run. Everything here is advisory: a warning is a prompt to look
at `Libraries/cross_class_screen.tsv` or at the element GFF3, not a verdict.

Output
------
`Libraries/library_health.tsv`, long format, four columns::

    section   item              metric               value

    library   all               n_consensi           11209
    class     Class_I/LINE      n_consensi           2141
    class     Class_I/LINE      max_len              7998
    boundary  DANTE_LINE        n_at_flank_ceiling   37
    boundary  DANTE_LINE        extension_fraction   0.2984
    screen    cross_class       n_trimmed            156

Rows are sorted, so the file is a deterministic function of its inputs.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

HEADER = ["section", "item", "metric", "value"]

# An inferred extension this close to dante_line's --flank ceiling means the
# boundary search never found an end. Calibrated across 87 assemblies: genomes
# with no such elements take a median 16% of their LINE base pairs from
# extensions, genomes with twenty or more take 80%.
FLANK_CEILING = 10000
CEILING_MARGIN = 100

# Warn when this much of a builder's output is inferred flank rather than the
# domain core it is anchored on.
EXTENSION_FRACTION_WARN = 0.50

# Below this many consensi the library is a test fixture, not a genome: "no
# Class_II sequences" is then an unremarkable consequence of the multiplicity
# floor having nothing to keep, and warning about it only trains people to
# ignore the warning.
MIN_CONSENSI_FOR_CLASS_WARN = 50


def read_fasta_classes(path: Path):
    """[(classification, length)] for every record."""
    out = []
    cls = None
    length = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cls is not None:
                    out.append((cls, length))
                header = line[1:].strip()
                cls = header.partition("#")[2] or "(unclassified)"
                length = 0
            else:
                length += len(line.strip())
    if cls is not None:
        out.append((cls, length))
    return out


def read_elements(path: Path, feature_types):
    """[(span, ext5, ext3, inferred5, inferred3)] for each element feature."""
    id_re = re.compile(r"(?:^|;)ID=")
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] not in feature_types:
                continue
            attrs = cols[8]
            if not id_re.search(attrs):
                continue

            def attr(name, default=0):
                m = re.search(rf"(?:^|;){name}=(\d+)", attrs)
                return int(m.group(1)) if m else default

            e5, e3 = attr("Extension_5prime"), attr("Extension_3prime")
            out.append((int(cols[4]) - int(cols[3]) + 1, e5, e3,
                        attr("Extension_5prime_inferred", e5),
                        attr("Extension_3prime_inferred", e3)))
    return out


def load_class_bounds(vocabulary=None):
    try:
        from classification import load_vocabulary
        vocab = load_vocabulary(vocabulary)
        return dict(getattr(vocab, "max_consensus_length", {}) or {})
    except Exception as exc:
        sys.stderr.write(f"NOTE: class length bounds unavailable ({exc}).\n")
        return {}


def bound_for(cls, table):
    best, best_depth = 0, -1
    for prefix, value in table.items():
        if cls == prefix or cls.startswith(prefix + "/"):
            depth = len(prefix.split("/"))
            if depth > best_depth:
                best, best_depth = value, depth
    return best


def median(values):
    if not values:
        return 0
    ordered = sorted(values)
    return ordered[len(ordered) // 2]


def collect(args):
    rows = []
    warnings = []

    def add(section, item, metric, value):
        rows.append((section, item, metric, value))

    # ── the library RepeatMasker will actually search ────────────────
    bounds = load_class_bounds(args.vocabulary)
    if args.library and Path(args.library).is_file():
        records = read_fasta_classes(Path(args.library))
        add("library", "all", "n_consensi", len(records))
        add("library", "all", "total_bp", sum(l for _, l in records))

        per_class = {}
        for cls, length in records:
            per_class.setdefault(cls, []).append(length)
        for cls in sorted(per_class):
            lengths = per_class[cls]
            limit = bound_for(cls, bounds)
            over = sum(1 for l in lengths if limit and l > limit)
            add("class", cls, "n_consensi", len(lengths))
            add("class", cls, "total_bp", sum(lengths))
            add("class", cls, "median_len", median(lengths))
            add("class", cls, "max_len", max(lengths))
            add("class", cls, "max_allowed_len", limit)
            add("class", cls, "n_over_bound", over)
            if over:
                warnings.append(
                    f"{cls}: {over} consensus/consensi longer than the "
                    f"{limit} bp bound for this class (longest {max(lengths)} bp)")

        # An empty class is the shape of the pre-1.4.0 DANTE_TIR join bug.
        if (len(records) >= MIN_CONSENSI_FOR_CLASS_WARN
                and not any(c.startswith("Class_II") for c in per_class)):
            warnings.append(
                "the library contains NO Class_II sequences — DNA transposons "
                "will only be annotated where a structural element was found. "
                "Check DANTE_TIR/all_representative_elements_combined.fasta is "
                "non-empty.")
    else:
        add("library", "all", "n_consensi", 0)

    # ── inferred element boundaries, per builder ─────────────────────
    for label, path, types in (
            ("DANTE_LINE", args.dante_line_gff3, {"LINE_element"}),
            ("DANTE_TIR_FALLBACK", args.fallback_gff3,
             {"sequence_feature", "TIR_fallback_element"})):
        if not path or not Path(path).is_file():
            continue
        elements = read_elements(Path(path), types)
        if not elements:
            add("boundary", label, "n_elements", 0)
            continue
        total = sum(e[0] for e in elements)
        ext = sum(e[1] + e[2] for e in elements)
        ceiling = sum(1 for e in elements
                      if max(e[3], e[4]) >= FLANK_CEILING - CEILING_MARGIN)
        capped = sum(1 for e in elements if (e[1], e[2]) != (e[3], e[4]))
        fraction = ext / total if total else 0.0
        add("boundary", label, "n_elements", len(elements))
        add("boundary", label, "total_bp", total)
        add("boundary", label, "extension_bp", ext)
        add("boundary", label, "extension_fraction", f"{fraction:.4f}")
        add("boundary", label, "n_at_flank_ceiling", ceiling)
        add("boundary", label, "n_capped", capped)
        if ceiling:
            warnings.append(
                f"{label}: {ceiling} element(s) had a flank alignment reaching the "
                f"{FLANK_CEILING} bp ceiling — the boundary search found no end for "
                f"them. Across 87 assemblies this count tracks how much of the "
                f"annotation comes from inferred flank rather than the element.")
        if fraction > EXTENSION_FRACTION_WARN:
            warnings.append(
                f"{label}: {100 * fraction:.0f}% of its base pairs are inferred "
                f"flank rather than the anchoring domain(s).")

    # ── what the cross-class screen did ──────────────────────────────
    if args.screen_audit and Path(args.screen_audit).is_file():
        import csv as _csv
        trimmed = dropped = removed = 0
        with open(args.screen_audit) as fh:
            for row in _csv.DictReader(fh, delimiter="\t"):
                action = row.get("action", "")
                try:
                    delta = int(row["original_length"]) - int(row["retained_length"])
                except (KeyError, ValueError):
                    delta = 0
                if action == "trimmed":
                    trimmed += 1
                    removed += delta
                elif action == "dropped":
                    dropped += 1
                    removed += delta
        add("screen", "cross_class", "n_trimmed", trimmed)
        add("screen", "cross_class", "n_dropped", dropped)
        add("screen", "cross_class", "bp_removed", removed)

    return rows, warnings


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Summarise repeat-library and element-boundary health.")
    p.add_argument("-o", "--output", required=True, help="TSV to write")
    p.add_argument("--library", help="Final RepeatMasker library FASTA")
    p.add_argument("--dante-line-gff3", help="DANTE_LINE/DANTE_LINE.gff3")
    p.add_argument("--fallback-gff3", help="DANTE_TIR_FALLBACK/DANTE_TIR_FALLBACK.gff3")
    p.add_argument("--screen-audit", help="Libraries/cross_class_screen.tsv")
    p.add_argument("--vocabulary", help="classification_vocabulary.yaml")
    args = p.parse_args(argv)

    rows, warnings = collect(args)

    out = Path(args.output)
    tmp = out.with_suffix(out.suffix + ".tmp")
    with open(tmp, "w") as fh:
        fh.write("\t".join(HEADER) + "\n")
        for row in sorted(rows, key=lambda r: (r[0], r[1], r[2])):
            fh.write("\t".join(str(x) for x in row) + "\n")
    os.replace(tmp, out)

    print(f"library health summary -> {out}")
    for warning in warnings:
        print(f"WARNING: {warning}", file=sys.stderr)
    if not warnings:
        print("no library-health warnings")
    return 0


if __name__ == "__main__":
    sys.exit(main())
