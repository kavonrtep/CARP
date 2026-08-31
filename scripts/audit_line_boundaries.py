#!/usr/bin/env python3
"""
Audit DANTE_LINE consensus boundaries in one or more finished CARP output trees.

Why
---
``dante_line.py`` fixes each element's end from an all-vs-all alignment of its
flanks and keeps the ``min_num_alignments``-th largest such alignment, bounded
only by a neighbouring DANTE domain, a TideHunter array, or ``--flank``
(default 10000). In a repeat-dense genome the flank of a LINE is frequently
another, far more abundant repeat, so the extension runs to the ceiling and the
"LINE consensus" becomes mostly foreign sequence. RepeatMasker then masks that
foreign part genome-wide and labels it ``Class_I/LINE``.

Nothing downstream catches it: ``reduce_library_size.py`` and
``containment_reduce_library.py`` compare only within a class, and
``filter_ltr_rt_library`` screens only the LTR library.

What this measures
------------------
Every RepeatMasker hit records which interval of the consensus it matched, and
every LINE consensus has a recoverable core::

    core = [Extension_5prime + 1, rep_length - Extension_3prime]

(``*_revcomp`` library names are already in element orientation, so no
coordinate flip is needed). Joining the two attributes each LINE base pair
either to the element's own ENDO/RT core or to the appended extension.

The RepeatMasker ``.out`` contains overlapping fragment hits, so its raw base
pairs exceed the non-overlapping total in ``summary_statistics.csv``. Where a
``summary_statistics.csv`` is present the raw split is rescaled by
``unified_LINE_bp / raw_LINE_bp`` to give a comparable ``*_scaled`` figure.

Also flagged, because it was silent in the same runs: an empty DANTE_TIR
library (``all_representative_elements_combined.fasta`` of zero bytes), the
join-key bug fixed in CARP 1.4.0, which leaves the RepeatMasker library with no
Class_II sequences at all.

This script is read-only and never touches a run directory.

Usage
-----
    audit_line_boundaries.py -o line_audit.tsv RUNDIR [RUNDIR ...]
    audit_line_boundaries.py -o line_audit.tsv --detail-dir det/ RUNDIR
    audit_line_boundaries.py --jobs 8 -o line_audit.tsv /path/to/runs/*/carp_output

RUNDIR is a CARP output directory (the one holding ``DANTE_LINE/``,
``RepeatMasker/`` and ``Libraries/``), or a run directory containing
``carp_output/``.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

# A hit counts as core-driven when this much of the consensus interval it
# matched lies inside the element's own ENDO/RT core.
CORE_FRACTION_STRICT = 0.80
CORE_FRACTION_LOOSE = 0.50

# A consensus whose extension reached this close to dante_line's --flank
# ceiling never found a boundary at all.
FLANK_CEILING = 10000
CEILING_MARGIN = 100

# No known plant LINE is this long; the ENDO+RT core alone is ~2.1 kb.
IMPLAUSIBLE_REP_LEN = 8000

LIBRARY_CANDIDATES = (
    "RepeatMasker/combined_library_reduced_containment.fasta",
    "Libraries/combined_library_reduced_containment.fasta",
    "Libraries/combined_library_reduced.fasta",
    "Libraries/combined_library.fasta",
)

FAI_CANDIDATES = (
    "genome_cleaned.fasta.fai",
    "genome_cleaned.fasta.seqkit.fai",
)

FIELDS = [
    "run",
    "genome_id",
    "carp_version",
    "genome_bp",
    "line_unified_bp",
    "line_unified_pct",
    "line_rm_raw_bp",
    "core_strict_bp",
    "core_loose_bp",
    "straddle_bp",
    "extension_only_bp",
    "extension_only_pct_of_line",
    "line_corrected_pct_low",
    "line_corrected_pct_high",
    "n_line_reps",
    "median_rep_len",
    "max_rep_len",
    "n_reps_over_8kb",
    "pct_line_bp_from_reps_over_8kb",
    "n_reps_at_flank_ceiling",
    "pct_line_bp_from_reps_at_ceiling",
    "worst_rep",
    "worst_rep_len",
    "worst_rep_pct_genome",
    "class_ii_library_empty",
    "status",
]


def resolve_output_dir(path: Path) -> Path:
    """Accept either a CARP output dir or a run dir holding carp_output/."""
    if (path / "DANTE_LINE").is_dir():
        return path
    nested = path / "carp_output"
    if (nested / "DANTE_LINE").is_dir():
        return nested
    return path


def read_rep_lengths(out_dir: Path):
    """(rep name -> length, rep name -> classification) from the library RepeatMasker used."""
    for rel in LIBRARY_CANDIDATES:
        fasta = out_dir / rel
        if fasta.is_file() and fasta.stat().st_size > 0:
            lengths, classes = {}, {}
            name = None
            total = 0
            with open(fasta) as fh:
                for line in fh:
                    if line.startswith(">"):
                        if name is not None:
                            lengths[name] = total
                        header = line[1:].strip()
                        name = header.split("#", 1)[0]
                        classes[name] = header.split("#", 1)[1] if "#" in header else ""
                        total = 0
                    else:
                        total += len(line.strip())
            if name is not None:
                lengths[name] = total
            return lengths, classes, rel
    return {}, {}, ""


def read_extensions(gff3: Path):
    """LINE_element ID -> (Extension_5prime, Extension_3prime)."""
    extensions = {}
    id_re = re.compile(r"(?:^|;)ID=([^;\n]+)")
    e5_re = re.compile(r"(?:^|;)Extension_5prime=(\d+)")
    e3_re = re.compile(r"(?:^|;)Extension_3prime=(\d+)")
    with open(gff3) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9 or cols[2] != "LINE_element":
                continue
            m = id_re.search(cols[8])
            if not m:
                continue
            m5 = e5_re.search(cols[8])
            m3 = e3_re.search(cols[8])
            extensions[m.group(1)] = (
                int(m5.group(1)) if m5 else 0,
                int(m3.group(1)) if m3 else 0,
            )
    return extensions


def build_cores(rep_lengths, rep_classes, extensions):
    """rep name -> (core_start, core_end, ext5, ext3) in consensus coordinates."""
    cores = {}
    for name, length in rep_lengths.items():
        if rep_classes.get(name) != "Class_I/LINE":
            continue
        group = name[:-8] if name.endswith("_revcomp") else name
        ext5, ext3 = extensions.get(group, (0, 0))
        start, end = ext5 + 1, length - ext3
        if end < start:  # extensions inconsistent with the library length
            start, end = 1, length
        cores[name] = (start, end, ext5, ext3)
    return cores


def matched_repeat_interval(cols):
    """(begin, end) of the consensus interval a RepeatMasker .out row matched.

    Column 9 is '+' or 'C'. On '+' the trailing triplet is
    ``begin end (left)``; on 'C' it is ``(left) end begin`` -- so the first of
    the three is the parenthesised remainder, not a coordinate, and reading it
    as one silently mixes a length into the interval.
    """
    end = int(cols[12])
    begin = int(cols[13].strip("()")) if cols[8] == "C" else int(cols[11])
    return (begin, end) if begin <= end else (end, begin)


def scan_repeatmasker(out_file: Path, cores):
    """Split LINE base pairs into core / straddle / extension, and per consensus."""
    buckets = {"strict": 0, "loose": 0, "straddle": 0, "extension": 0}
    per_rep = defaultdict(float)
    hits = 0
    with open(out_file) as fh:
        for line in fh:
            cols = line.split()
            if len(cols) < 14 or not cols[0].isdigit():
                continue
            if cols[10] != "Class_I/LINE":
                continue
            core = cores.get(cols[9])
            if core is None:
                continue
            span = int(cols[6]) - int(cols[5]) + 1
            begin, end = matched_repeat_interval(cols)
            core_start, core_end = core[0], core[1]
            overlap = min(end, core_end) - max(begin, core_start) + 1
            matched = end - begin + 1
            fraction = (overlap / matched) if (matched > 0 and overlap > 0) else 0.0
            if fraction >= CORE_FRACTION_STRICT:
                buckets["strict"] += span
            elif fraction >= CORE_FRACTION_LOOSE:
                buckets["loose"] += span
            elif fraction > 0:
                buckets["straddle"] += span
            else:
                buckets["extension"] += span
            per_rep[cols[9]] += span
            hits += 1
    return buckets, per_rep, hits


def read_genome_bp(out_dir: Path):
    for rel in FAI_CANDIDATES:
        fai = out_dir / rel
        if fai.is_file():
            total = 0
            with open(fai) as fh:
                for line in fh:
                    parts = line.split("\t")
                    if len(parts) >= 2:
                        total += int(parts[1])
            if total:
                return total
    return 0


def read_unified_line_bp(out_dir: Path):
    stats = out_dir / "summary_statistics.csv"
    if not stats.is_file():
        return 0, 0.0
    with open(stats) as fh:
        for line in fh:
            parts = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
            if len(parts) >= 3 and parts[0] == "Class_I/LINE":
                try:
                    return int(parts[1]), float(parts[2])
                except ValueError:
                    return 0, 0.0
    return 0, 0.0


def read_run_metadata(out_dir: Path):
    import json

    genome_id = version = ""
    provenance = out_dir / "run_provenance.json"
    if provenance.is_file():
        try:
            data = json.loads(provenance.read_text())
            version = str(data.get("pipeline_version", ""))
            fasta = data.get("config", {}).get("genome_fasta", "")
            genome_id = Path(fasta).parent.name if fasta else ""
        except (ValueError, OSError):
            pass
    params = out_dir.parent / "params"
    if params.is_file():
        try:
            genome_id = json.loads(params.read_text()).get("genome_id", genome_id)
        except (ValueError, OSError):
            pass
    return genome_id, version


def class_ii_library_empty(out_dir: Path):
    lib = out_dir / "DANTE_TIR" / "all_representative_elements_combined.fasta"
    if not lib.is_file():
        return ""
    return "YES" if lib.stat().st_size == 0 else "no"


def percentile(values, q):
    if not values:
        return 0
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, int(len(ordered) * q))]


def audit_one(run_path: str, detail_dir: str | None = None):
    path = Path(run_path)
    out_dir = resolve_output_dir(path)
    label = path.name if path.name != "carp_output" else path.parent.name
    row = {field: "" for field in FIELDS}
    row["run"] = label

    gff3 = out_dir / "DANTE_LINE" / "DANTE_LINE.gff3"
    rm_out = out_dir / "RepeatMasker" / "RM_on_combined_library.out"
    if not gff3.is_file() or not rm_out.is_file():
        row["status"] = "skipped: no DANTE_LINE.gff3 or RM_on_combined_library.out"
        return row

    genome_id, version = read_run_metadata(out_dir)
    row["genome_id"] = genome_id
    row["carp_version"] = version
    row["class_ii_library_empty"] = class_ii_library_empty(out_dir)

    rep_lengths, rep_classes, _ = read_rep_lengths(out_dir)
    if not rep_lengths:
        row["status"] = "skipped: no RepeatMasker library found"
        return row

    cores = build_cores(rep_lengths, rep_classes, read_extensions(gff3))
    if not cores:
        row["status"] = "no Class_I/LINE consensi in the library"
        return row

    buckets, per_rep, hits = scan_repeatmasker(rm_out, cores)
    raw_total = sum(buckets.values())
    if raw_total == 0:
        row["status"] = "no Class_I/LINE hits in RepeatMasker output"
        return row

    genome_bp = read_genome_bp(out_dir)
    unified_bp, unified_pct = read_unified_line_bp(out_dir)
    scale = (unified_bp / raw_total) if unified_bp else 1.0

    def as_pct_of_genome(raw_bp):
        if not genome_bp:
            return ""
        return round(100.0 * raw_bp * scale / genome_bp, 4)

    lengths = [rep_lengths[name] for name in cores]
    long_reps = [n for n in cores if rep_lengths[n] > IMPLAUSIBLE_REP_LEN]
    ceiling_reps = [
        n for n, c in cores.items()
        if max(c[2], c[3]) >= FLANK_CEILING - CEILING_MARGIN
    ]
    worst = max(per_rep, key=per_rep.get) if per_rep else ""

    row.update({
        "genome_bp": genome_bp or "",
        "line_unified_bp": unified_bp or "",
        "line_unified_pct": round(unified_pct, 4) if unified_pct else "",
        "line_rm_raw_bp": raw_total,
        "core_strict_bp": buckets["strict"],
        "core_loose_bp": buckets["loose"],
        "straddle_bp": buckets["straddle"],
        "extension_only_bp": buckets["extension"],
        "extension_only_pct_of_line": round(100.0 * buckets["extension"] / raw_total, 2),
        "line_corrected_pct_low": as_pct_of_genome(buckets["strict"]),
        "line_corrected_pct_high": as_pct_of_genome(
            buckets["strict"] + buckets["loose"] + buckets["straddle"]
        ),
        "n_line_reps": len(cores),
        "median_rep_len": percentile(lengths, 0.5),
        "max_rep_len": max(lengths),
        "n_reps_over_8kb": len(long_reps),
        "pct_line_bp_from_reps_over_8kb": round(
            100.0 * sum(per_rep.get(n, 0) for n in long_reps) / raw_total, 2),
        "n_reps_at_flank_ceiling": len(ceiling_reps),
        "pct_line_bp_from_reps_at_ceiling": round(
            100.0 * sum(per_rep.get(n, 0) for n in ceiling_reps) / raw_total, 2),
        "worst_rep": worst,
        "worst_rep_len": rep_lengths.get(worst, ""),
        "worst_rep_pct_genome": as_pct_of_genome(per_rep.get(worst, 0)),
        "status": "ok",
    })

    if detail_dir:
        write_detail(Path(detail_dir) / f"{label}.line_consensi.tsv",
                     cores, rep_lengths, per_rep, genome_bp, scale)
    return row


def write_detail(dest: Path, cores, rep_lengths, per_rep, genome_bp, scale):
    """Per-consensus detail, ordered by the genomic base pairs each one masks."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    rows = sorted(
        cores.items(),
        key=lambda kv: (-per_rep.get(kv[0], 0), kv[0]),
    )
    with open(tmp, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "consensus", "rep_len", "ext_5prime", "ext_3prime", "core_len",
            "at_flank_ceiling", "genomic_bp_raw", "pct_of_genome",
        ])
        for name, (start, end, ext5, ext3) in rows:
            raw = per_rep.get(name, 0)
            writer.writerow([
                name,
                rep_lengths[name],
                ext5,
                ext3,
                end - start + 1,
                "YES" if max(ext5, ext3) >= FLANK_CEILING - CEILING_MARGIN else "",
                int(raw),
                round(100.0 * raw * scale / genome_bp, 5) if genome_bp else "",
            ])
    os.replace(tmp, dest)


def _worker(args):
    run_path, detail_dir = args
    try:
        return audit_one(run_path, detail_dir)
    except Exception as exc:  # one bad run must not lose the sweep
        row = {field: "" for field in FIELDS}
        row["run"] = Path(run_path).name
        row["status"] = f"error: {type(exc).__name__}: {exc}"
        return row


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Audit DANTE_LINE consensus boundaries in finished CARP runs.")
    parser.add_argument("rundirs", nargs="+",
                        help="CARP output directories (or run dirs containing carp_output/)")
    parser.add_argument("-o", "--output", required=True,
                        help="TSV summary, one row per run")
    parser.add_argument("--detail-dir",
                        help="also write a per-consensus TSV per run into this directory")
    parser.add_argument("-j", "--jobs", type=int, default=1,
                        help="runs to scan concurrently (default 1)")
    args = parser.parse_args(argv)

    work = [(d, args.detail_dir) for d in args.rundirs]
    if args.jobs > 1 and len(work) > 1:
        with Pool(min(args.jobs, len(work))) as pool:
            rows = pool.map(_worker, work)
    else:
        rows = [_worker(item) for item in work]

    rows.sort(key=lambda r: r["run"])
    out = Path(args.output)
    tmp = out.with_suffix(out.suffix + ".tmp")
    with open(tmp, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, out)

    audited = [r for r in rows if r["status"] == "ok"]
    print(f"audited {len(audited)} of {len(rows)} run(s) -> {out}", file=sys.stderr)
    if audited:
        inflated = [r for r in audited if float(r["extension_only_pct_of_line"] or 0) >= 50]
        print(f"  {len(inflated)} run(s) with >=50% of LINE bp from consensus extensions",
              file=sys.stderr)
        empty = [r for r in audited if r["class_ii_library_empty"] == "YES"]
        print(f"  {len(empty)} run(s) with an empty DANTE_TIR library "
              f"(pre-1.4.0 join bug)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
