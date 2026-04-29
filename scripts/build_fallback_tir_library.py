#!/usr/bin/env python3
"""Build the optional DANTE_TIR_FALLBACK-derived library.

When the pipeline-level flag ``include_dante_tir_fallback_in_library``
is on, this script promotes a subset of fallback elements into the
RepeatMasker library. Default behaviour (flag off) is to emit an empty
library so concatenate_libraries can append unconditionally.

Stages (when enabled):

  1. Re-cluster the post-overlap fallback survivors with mmseqs2. The
     re-cluster is essential — the cluster sizes recorded earlier in
     dante_tir_fallback.py were computed before primary-overlap
     filtering, so they include members that have since been dropped.
  2. Apply a Multiplicity floor on cluster size (--min-multiplicity).
  3. Canonicalise the default-library headers (LTR / DANTE_TIR primary
     / LINE / optional custom) via classification.canonicalise_fasta_headers,
     concatenate them into a single BLAST DB, and remember each subject
     sequence's classification.
  4. blastn the fallback rep candidates against that DB
     (-evalue 1e-19, -max_target_seqs 10 — matching filter_ltr_rt_library).
  5. Strict class-aware filter: drop any candidate that has at least one
     hit to a subject whose classification is not on the candidate's
     canonical lineage chain (same path or one is an ancestor of the
     other). A CACTA fallback hitting a hAT primary is dropped — same
     parent (TIR) but different superfamily.

The audit TSV records one row per kept rep and one row per
(dropped rep × conflicting subject) so a reviewer can re-trace why a
rep was rejected. Dropping is "better safe than sorry": a single
misclassified library entry can cascade into a genome-wide annotation
shift, so we err toward under-inclusion.
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path

# Sibling import — classification.py lives next to this script and is
# already on PATH inside the conda env (rules prepend scripts/).
sys.path.insert(0, str(Path(__file__).resolve().parent))
from classification import canonicalise_fasta_headers  # noqa: E402


def parse_class_from_header(header: str) -> str:
    """``>name#Class_X/Y/Z extra`` → ``Class_X/Y/Z``. Empty / no-# → ``Unknown``."""
    if "#" not in header:
        return "Unknown"
    return header.split("#", 1)[1].strip().split()[0]


def is_compatible(a: str, b: str) -> bool:
    """Strict path-prefix compatibility on canonical slash-form classes.

    Compatible iff the two classes are equal or one is an ancestor of
    the other in the slash-delimited tree. Siblings (CACTA vs hAT) are
    NOT compatible, even though they share the TIR parent — picking
    one over the other based on a blast hit would be guessing.
    """
    if a == b:
        return True
    pa = a.split("/")
    pb = b.split("/")
    n = min(len(pa), len(pb))
    if pa[:n] != pb[:n]:
        return False
    # Equal up to the shorter; one is an ancestor of the other only if
    # one is exactly that shorter prefix.
    return len(pa) == n or len(pb) == n


def fasta_seqid_class_index(path: Path) -> dict[str, str]:
    """``seqid → classification`` from FASTA headers.

    The seqid is the entire token before whitespace, including any
    ``#class`` suffix — mmseqs2 and blast both use that whole string
    as the sequence identifier, so we must too.
    """
    out: dict[str, str] = {}
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                hdr = line[1:].strip()
                seqid = hdr.split()[0]
                out[seqid] = parse_class_from_header(hdr)
    return out


def run(cmd: list[str], log_prefix: str) -> None:
    print(f"[{log_prefix}] {' '.join(cmd)}", file=sys.stderr)
    subprocess.run(cmd, check=True)


_AUDIT_COLS = [
    "rep_id", "rep_class", "status", "reason",
    "conflicting_subject_id", "conflicting_subject_class",
]


def write_empty_outputs(out_fasta: Path, out_tsv: Path, reason: str) -> None:
    out_fasta.write_text("")
    with out_tsv.open("w") as fh:
        csv.writer(fh, delimiter="\t").writerow(_AUDIT_COLS)
    print(f"empty fallback library written: {reason}", file=sys.stderr)


def extract_seqs_by_id(in_fasta: Path, keep_ids: set[str], out_fasta: Path) -> int:
    """Copy FASTA records whose seqid (full pre-whitespace token,
    including any ``#class`` suffix) is in ``keep_ids``."""
    n = 0
    keep = False
    with in_fasta.open() as fin, out_fasta.open("w") as fout:
        for line in fin:
            if line.startswith(">"):
                seqid = line[1:].strip().split()[0]
                keep = seqid in keep_ids
                if keep:
                    n += 1
            if keep:
                fout.write(line)
    return n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--enabled", action="store_true",
                    help="If absent, write empty outputs and exit 0.")
    ap.add_argument("--fallback-fasta", required=True, type=Path)
    ap.add_argument("--ltr-library", required=True, type=Path)
    ap.add_argument("--tir-primary-library", required=True, type=Path)
    ap.add_argument("--line-library", required=True, type=Path)
    ap.add_argument("--custom-library", type=Path, default=None)
    ap.add_argument("--workdir", required=True, type=Path)
    ap.add_argument("--min-multiplicity", required=True, type=int)
    ap.add_argument("--blast-evalue", default="1e-19")
    ap.add_argument("--blast-max-targets", type=int, default=10)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--output-fasta", required=True, type=Path)
    ap.add_argument("--output-dropped-tsv", required=True, type=Path)
    args = ap.parse_args()

    # ── early exits ────────────────────────────────────────────────
    if not args.enabled:
        write_empty_outputs(args.output_fasta, args.output_dropped_tsv,
                            "feature disabled (include_dante_tir_fallback_in_library=false)")
        return 0
    if (not args.fallback_fasta.exists()
            or args.fallback_fasta.stat().st_size == 0):
        write_empty_outputs(args.output_fasta, args.output_dropped_tsv,
                            "no fallback survivors after primary-overlap filter")
        return 0

    args.workdir.mkdir(parents=True, exist_ok=True)

    # ── 1. mmseqs re-cluster ───────────────────────────────────────
    cluster_dir = args.workdir / "cluster"
    cluster_dir.mkdir(exist_ok=True)
    run([
        "mmseqs", "easy-cluster",
        str(args.fallback_fasta),
        str(cluster_dir / "cluster"),
        str(cluster_dir / "tmp"),
        "--threads", str(args.threads),
    ], log_prefix="mmseqs")

    cluster_tsv = cluster_dir / "cluster_cluster.tsv"
    rep_fasta = cluster_dir / "cluster_rep_seq.fasta"
    if (not cluster_tsv.exists() or not rep_fasta.exists()
            or rep_fasta.stat().st_size == 0):
        write_empty_outputs(args.output_fasta, args.output_dropped_tsv,
                            "mmseqs produced no clusters")
        return 0

    # ── 2. Multiplicity floor ──────────────────────────────────────
    cluster_sizes: dict[str, int] = {}
    with cluster_tsv.open() as fh:
        for line in fh:
            rep, _ = line.rstrip("\n").split("\t")
            cluster_sizes[rep] = cluster_sizes.get(rep, 0) + 1

    keep_reps = {rep for rep, n in cluster_sizes.items()
                 if n >= args.min_multiplicity}
    print(f"clusters total: {len(cluster_sizes)}, "
          f"with size >= {args.min_multiplicity}: {len(keep_reps)}",
          file=sys.stderr)
    if not keep_reps:
        write_empty_outputs(args.output_fasta, args.output_dropped_tsv,
                            f"no clusters with size >= {args.min_multiplicity}")
        return 0

    candidates = args.workdir / "candidate_reps.fasta"
    n_candidates = extract_seqs_by_id(rep_fasta, keep_reps, candidates)
    print(f"wrote {n_candidates} candidate reps -> {candidates}",
          file=sys.stderr)

    # ── 3. Canonicalise + concatenate default libraries ────────────
    sources: list[tuple[Path, str]] = [
        (args.ltr_library, "DANTE_LTR"),
        (args.tir_primary_library, "DANTE_TIR"),
        (args.line_library, "DANTE_LINE"),
    ]
    if (args.custom_library is not None
            and args.custom_library.exists()
            and args.custom_library.stat().st_size > 0):
        sources.append((args.custom_library, "custom_library"))

    default_db = args.workdir / "default_db.fasta"
    subject_classes: dict[str, str] = {}
    with default_db.open("w") as out_fh:
        for path, source in sources:
            if not path.exists() or path.stat().st_size == 0:
                print(f"skipping empty/missing {path}", file=sys.stderr)
                continue
            canon_path = args.workdir / f"{path.stem}.canon.fasta"
            canonicalise_fasta_headers(path, canon_path, source=source)
            for sid, cls in fasta_seqid_class_index(canon_path).items():
                subject_classes[sid] = cls
            with canon_path.open() as fin:
                shutil.copyfileobj(fin, out_fh)

    if default_db.stat().st_size == 0:
        # Degenerate case (e.g. cold start with no DANTE_LTR/LINE/primary
        # libraries built yet). Keep all multiplicity-passing reps —
        # there's nothing to BLAST against, so no class-aware filter
        # runs. Fail-open here is consistent with "default db empty
        # means we cannot prove anything is wrong".
        shutil.copy(candidates, args.output_fasta)
        with args.output_dropped_tsv.open("w") as fh:
            csv.writer(fh, delimiter="\t").writerow(_AUDIT_COLS)
        print("default DB empty; all candidates passed by default",
              file=sys.stderr)
        return 0

    # ── 4. BLAST DB + blastn ───────────────────────────────────────
    run(["makeblastdb", "-in", str(default_db), "-dbtype", "nucl"],
        log_prefix="makeblastdb")

    blast_out = args.workdir / "blast.tsv"
    run([
        "blastn", "-task", "blastn",
        "-query", str(candidates),
        "-db", str(default_db),
        "-outfmt", "6",
        "-evalue", args.blast_evalue,
        "-max_target_seqs", str(args.blast_max_targets),
        "-num_threads", str(args.threads),
        "-out", str(blast_out),
    ], log_prefix="blastn")

    # ── 5. Strict class-aware filter ───────────────────────────────
    rep_classes = fasta_seqid_class_index(candidates)
    incompatible: dict[str, list[tuple[str, str]]] = {}
    if blast_out.exists() and blast_out.stat().st_size > 0:
        with blast_out.open() as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 2:
                    continue
                qid, sid = f[0], f[1]
                qcls = rep_classes.get(qid, "Unknown")
                scls = subject_classes.get(sid, "Unknown")
                if not is_compatible(qcls, scls):
                    incompatible.setdefault(qid, []).append((sid, scls))

    drop_set = set(incompatible)
    keep_ids = set(rep_classes) - drop_set
    print(f"reps total: {len(rep_classes)}, "
          f"dropped by class-aware filter: {len(drop_set)}, "
          f"kept: {len(keep_ids)}", file=sys.stderr)

    # ── 6. Write outputs ───────────────────────────────────────────
    n_written = extract_seqs_by_id(candidates, keep_ids, args.output_fasta)

    with args.output_dropped_tsv.open("w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(_AUDIT_COLS)
        for rep_id in sorted(rep_classes):
            cls = rep_classes[rep_id]
            if rep_id in drop_set:
                for sid, scls in incompatible[rep_id]:
                    w.writerow([rep_id, cls, "dropped",
                                "incompatible_blast_hit", sid, scls])
            else:
                w.writerow([rep_id, cls, "kept", "", "", ""])

    print(f"wrote {n_written} sequences -> {args.output_fasta}",
          file=sys.stderr)
    print(f"audit log -> {args.output_dropped_tsv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
