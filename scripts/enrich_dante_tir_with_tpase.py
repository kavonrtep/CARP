#!/usr/bin/env python3
"""
Inject TPase ``protein_domain`` children into the parent ``sequence_feature``
rows of ``DANTE_TIR_final.gff3``.

Background
----------
The current ``dante_tir.py`` upstream tool emits one ``sequence_feature``
parent per identified TIR element but does not include the DANTE TPase
domain that triggered the structural call as a GFF3 child. This script
fills that gap by reading the raw DANTE GFF, applying the standard DANTE
quality filter, and writing the matching TPase domain row underneath each
parent.

Idempotency
-----------
A future ``dante_tir.py`` release may add TPase children itself. To stay
safe, the script runs a per-parent guard: any parent that already has at
least one ``protein_domain`` row carrying ``Parent=<this parent ID>`` is
left untouched (parent + its existing children pass through unchanged).

Defensive handling
------------------
TIR elements are expected to contain exactly one TPase domain that matches
the parent's ``Classification``. None of the following should happen in
practice; if they do, we log a ``WARN`` line and emit the parent without a
new child rather than fabricating data:

* zero matching TPase candidates,
* two or more matching TPase candidates,
* a candidate whose ``Final_Classification`` (canonicalised) does not equal
  the parent's ``Classification`` (canonicalised).

CLI
---
    enrich_dante_tir_with_tpase.py \\
        --dante-gff DANTE/DANTE.gff3 \\
        --dante-tir-gff DANTE_TIR/DANTE_TIR_final.gff3

The input ``DANTE_TIR_final.gff3`` is rewritten in place via an atomic
``os.replace`` from a sibling temp file. If the file is empty (no primary
elements) the script is a no-op.
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import classification  # noqa: E402

# DANTE per-domain quality thresholds (mirrors DANTE_filtered.gff3 defaults
# and dante_tir_fallback.py constants — kept inline because the upstream
# rule consumes raw DANTE.gff3 by design and we don't want to add a
# Snakemake edge to filter_dante just for this step).
QUALITY_MIN_IDENTITY = 0.35
QUALITY_MIN_SIMILARITY = 0.45
QUALITY_MAX_RELAT_INTERRUPTIONS = 3.0
QUALITY_RELAT_LENGTH_RANGE = (0.8, 1.2)


@dataclass
class GFFRow:
    seqname: str
    source: str
    feature: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: str
    raw: str = ""

    @classmethod
    def from_line(cls, line: str) -> Optional["GFFRow"]:
        if not line or line.startswith("#"):
            return None
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            return None
        try:
            start = int(fields[3])
            end = int(fields[4])
        except ValueError:
            return None
        return cls(
            seqname=fields[0],
            source=fields[1],
            feature=fields[2],
            start=start,
            end=end,
            score=fields[5],
            strand=fields[6],
            phase=fields[7],
            attributes=fields[8],
            raw=line.rstrip("\n"),
        )

    def attr(self, key: str) -> Optional[str]:
        m = re.search(rf"(?:^|;){re.escape(key)}=([^;]*)", self.attributes)
        return m.group(1) if m else None


@dataclass
class Parent:
    row: GFFRow
    parent_id: str
    classification: str
    children: List[GFFRow] = field(default_factory=list)


def parse_attributes_id(attr: str) -> Optional[str]:
    m = re.search(r"(?:^|;)ID=([^;]*)", attr)
    return m.group(1) if m else None


def parse_attributes_parent(attr: str) -> Optional[str]:
    m = re.search(r"(?:^|;)Parent=([^;]*)", attr)
    return m.group(1) if m else None


def gff3_quality_ok(row: GFFRow) -> bool:
    """Per-domain DANTE quality filter (same thresholds as DANTE_filtered.gff3)."""
    try:
        identity = float(row.attr("Identity") or "nan")
        similarity = float(row.attr("Similarity") or "nan")
        relat_interruptions = float(row.attr("Relat_Interruptions") or "nan")
        relat_length = float(row.attr("Relat_Length") or "nan")
    except ValueError:
        return False
    if not (identity >= QUALITY_MIN_IDENTITY):
        return False
    if not (similarity >= QUALITY_MIN_SIMILARITY):
        return False
    if not (relat_interruptions <= QUALITY_MAX_RELAT_INTERRUPTIONS):
        return False
    lo, hi = QUALITY_RELAT_LENGTH_RANGE
    if not (lo <= relat_length <= hi):
        return False
    return True


def load_quality_filtered_tpase(dante_gff: Path) -> List[GFFRow]:
    """Return DANTE TPase rows (Class_II|Subclass_1|TIR|*) that clear the
    standard quality filter, indexed externally by caller.
    """
    out: List[GFFRow] = []
    with open(dante_gff) as fh:
        for line in fh:
            row = GFFRow.from_line(line)
            if row is None:
                continue
            if row.feature != "protein_domain":
                continue
            if row.attr("Name") != "TPase":
                continue
            fc = row.attr("Final_Classification") or ""
            if "Subclass_1" not in fc:
                continue
            if not gff3_quality_ok(row):
                continue
            out.append(row)
    return out


def index_tpase_by_seqname(rows: List[GFFRow]) -> Dict[str, List[GFFRow]]:
    """Group TPase rows by chromosome / contig for fast per-parent lookup."""
    idx: Dict[str, List[GFFRow]] = defaultdict(list)
    for r in rows:
        idx[r.seqname].append(r)
    return dict(idx)


def find_tpase_for_parent(
    parent: Parent, tpase_by_seq: Dict[str, List[GFFRow]]
) -> List[GFFRow]:
    """Return TPase candidates that lie inside the parent's interval on the
    same strand. Strand match is required because TIR elements have a
    well-defined orientation derived from the TIR pair; a TPase on the
    opposite strand cannot be the correct anchor.
    """
    candidates: List[GFFRow] = []
    for r in tpase_by_seq.get(parent.row.seqname, ()):
        if r.start < parent.row.start or r.end > parent.row.end:
            continue
        if parent.row.strand in ("+", "-") and r.strand != parent.row.strand:
            continue
        candidates.append(r)
    return candidates


def parse_dante_tir_gff(path: Path) -> Tuple[List[str], List[Parent]]:
    """Parse DANTE_TIR_final.gff3 into header lines + parent records that
    carry their already-existing children. Non-parent rows without a
    ``Parent`` attribute are also stored on the previous parent's children
    list to keep ordering stable on rewrite.
    """
    headers: List[str] = []
    parents: List[Parent] = []
    parent_by_id: Dict[str, Parent] = {}

    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                headers.append(line.rstrip("\n"))
                continue
            row = GFFRow.from_line(line)
            if row is None:
                continue
            if row.feature == "sequence_feature":
                pid = parse_attributes_id(row.attributes) or ""
                cls = row.attr("Classification") or ""
                p = Parent(row=row, parent_id=pid, classification=cls)
                parents.append(p)
                if pid:
                    parent_by_id[pid] = p
            else:
                parent_id = parse_attributes_parent(row.attributes)
                if parent_id and parent_id in parent_by_id:
                    parent_by_id[parent_id].children.append(row)
                elif parents:
                    # Orphan row — stash on the most recent parent so order
                    # survives rewrite. Should not happen in practice.
                    parents[-1].children.append(row)
    return headers, parents


def child_row_for_tpase(parent: Parent, tpase: GFFRow) -> str:
    """Format the DANTE TPase row as a child of the DANTE_TIR parent.
    The original DANTE attributes (Name, Final_Classification, Identity,
    Similarity, Relat_Length, etc.) are preserved verbatim. ``Parent``
    pointer is appended.
    """
    attrs = tpase.attributes
    if "Parent=" in attrs:
        attrs = re.sub(r"(?:^|;)Parent=[^;]*", "", attrs).lstrip(";")
    if attrs and not attrs.endswith(";"):
        attrs = attrs + ";"
    attrs = attrs + f"Parent={parent.parent_id};Source=DANTE_TIR_enriched"
    return "\t".join([
        tpase.seqname, tpase.source, tpase.feature,
        str(tpase.start), str(tpase.end),
        tpase.score, tpase.strand, tpase.phase, attrs,
    ])


def enrich(parents: List[Parent], tpase_by_seq: Dict[str, List[GFFRow]],
           log: logging.Logger) -> Tuple[int, int, int, int, int]:
    """Inject TPase children into each parent that has none. Returns counts:
    (parents_total, parents_already_enriched, parents_injected,
     parents_skipped_zero, parents_skipped_multi_or_mismatch).
    """
    n_total = len(parents)
    n_already = 0
    n_injected = 0
    n_zero = 0
    n_skip_anomaly = 0

    for p in parents:
        existing_protein_domains = [c for c in p.children if c.feature == "protein_domain"]
        if existing_protein_domains:
            n_already += 1
            continue

        candidates = find_tpase_for_parent(p, tpase_by_seq)
        if len(candidates) == 0:
            log.warning(
                "no TPase candidate for parent %s at %s:%d-%d (%s); "
                "leaving parent without protein_domain child",
                p.parent_id, p.row.seqname, p.row.start, p.row.end,
                p.row.strand,
            )
            n_zero += 1
            continue

        if len(candidates) > 1:
            log.warning(
                "%d TPase candidates inside parent %s at %s:%d-%d — "
                "expected exactly 1; skipping injection",
                len(candidates), p.parent_id,
                p.row.seqname, p.row.start, p.row.end,
            )
            n_skip_anomaly += 1
            continue

        cand = candidates[0]

        # Classification cross-check. Both sides go through canonicalise()
        # to absorb the DANTE pipe form vs DANTE_TIR slash form.
        try:
            parent_canon = classification.canonicalise(p.classification, source="DANTE_TIR")
        except Exception as e:  # noqa: BLE001
            log.warning(
                "parent %s has unparseable Classification=%r (%s); skipping",
                p.parent_id, p.classification, e,
            )
            n_skip_anomaly += 1
            continue
        try:
            cand_canon = classification.canonicalise(
                cand.attr("Final_Classification") or "", source="DANTE",
            )
        except Exception as e:  # noqa: BLE001
            log.warning(
                "candidate TPase for parent %s has unparseable "
                "Final_Classification=%r (%s); skipping",
                p.parent_id, cand.attr("Final_Classification"), e,
            )
            n_skip_anomaly += 1
            continue

        if parent_canon != cand_canon:
            log.warning(
                "classification mismatch for parent %s: parent=%s vs TPase=%s; "
                "skipping injection",
                p.parent_id, parent_canon, cand_canon,
            )
            n_skip_anomaly += 1
            continue

        # All checks passed — inject.
        new_child = GFFRow.from_line(child_row_for_tpase(p, cand))
        if new_child is not None:
            p.children.append(new_child)
            n_injected += 1

    return n_total, n_already, n_injected, n_zero, n_skip_anomaly


def write_output(headers: List[str], parents: List[Parent], dest: Path) -> None:
    """Atomic in-place rewrite: write to a sibling temp file, then replace."""
    tmp_dir = dest.parent
    tmp_dir.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(suffix=".gff3", dir=str(tmp_dir))
    os.close(fd)
    try:
        with open(tmp, "w") as out:
            for h in headers:
                out.write(h + "\n")
            for p in parents:
                out.write(p.row.raw + "\n")
                for c in p.children:
                    out.write(c.raw + "\n")
        os.replace(tmp, dest)
    except Exception:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Inject DANTE TPase protein_domain children into "
                    "DANTE_TIR_final.gff3 (idempotent).",
    )
    parser.add_argument("--dante-gff", required=True, type=Path,
                       help="Path to DANTE.gff3 (raw, unfiltered)")
    parser.add_argument("--dante-tir-gff", required=True, type=Path,
                       help="Path to DANTE_TIR_final.gff3 (rewritten in place)")
    parser.add_argument("--verbose", action="store_true",
                       help="Log per-parent decisions at DEBUG level")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        stream=sys.stderr,
    )
    log = logging.getLogger("enrich_dante_tir")

    if not args.dante_tir_gff.exists() or args.dante_tir_gff.stat().st_size == 0:
        log.info("%s is missing or empty — nothing to enrich, exiting", args.dante_tir_gff)
        return

    if not args.dante_gff.exists():
        log.error("DANTE GFF not found: %s", args.dante_gff)
        sys.exit(2)

    headers, parents = parse_dante_tir_gff(args.dante_tir_gff)
    if not parents:
        log.info("no sequence_feature parents in %s — nothing to enrich",
                 args.dante_tir_gff)
        return

    log.info("loaded %d primary DANTE_TIR parents from %s",
             len(parents), args.dante_tir_gff)

    tpase_rows = load_quality_filtered_tpase(args.dante_gff)
    log.info("loaded %d quality-passing Subclass_1 TPase rows from %s",
             len(tpase_rows), args.dante_gff)
    tpase_by_seq = index_tpase_by_seqname(tpase_rows)

    n_total, n_already, n_injected, n_zero, n_skip = enrich(parents, tpase_by_seq, log)
    log.info(
        "enrichment summary: total=%d already_had_children=%d "
        "injected=%d no_candidate=%d skipped_anomaly=%d",
        n_total, n_already, n_injected, n_zero, n_skip,
    )

    write_output(headers, parents, args.dante_tir_gff)
    log.info("wrote %s", args.dante_tir_gff)


if __name__ == "__main__":
    main()
