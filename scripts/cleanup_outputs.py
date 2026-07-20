#!/usr/bin/env python3
"""Post-run cleanup of intermediate files under a CARP output_dir.

Called by run_pipeline.py after a SUCCESSFUL run (rc == 0, non-dry-run) when the
config key ``cleanup_intermediates`` is ``minimal`` (default) or ``maximal``.
Deletes per-tool scratch that no downstream Snakemake rule consumes and the
output manifest does not require, freeing space without breaking any consumer.

Guarantees (belt-and-suspenders — the delete lists below are already curated to
unconsumed scratch):
  * NEVER deletes a path listed in the manifest (``manifest.py:OUTPUTS``), a
    symlink target of one, a CI/count-checked file, or the run metadata.
  * NEVER deletes a directory that contains any kept path.
  * Confined to ``output_dir``; a missing candidate is skipped; a delete failure
    warns and continues — the pipeline outcome is never affected.

Modes:
  minimal  — clearly-unconsumed scratch (safe; the default)
  maximal  — minimal + big TideCluster trees + tool workdirs + mmseqs tmp
  none     — no-op

Standalone CLI: ``cleanup_outputs.py <output_dir> [--mode minimal|maximal|none]
[--dry-run]``.
"""
import argparse
import os
import shutil
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
from manifest import OUTPUTS  # noqa: E402  — single source of the keep-list

VALID_MODES = ("minimal", "maximal", "none")

# Not in the manifest but must survive cleanup: files CI count-checks
# (assert_fixture_outputs.sh) + the run metadata. Relative to the output root.
_EXTRA_KEEP = (
    "carp_manifest.json",
    ".classifications_validated",
    "DANTE_TIR/DANTE_TIR_final.fasta",
    "DANTE_TIR/DANTE_TIR_combined.gff3",
    "DANTE_LINE/DANTE_LINE.gff3",
    "DANTE_LTR/LTR_RTs_library.fasta",           # symlink; its target is kept too
)

# minimal: verified unconsumed by any Snakefile rule input (per the cleanup
# audit). Entries may be literal paths or globs.
_MINIMAL = (
    "DANTE/DANTE_filtered.gff3.tmp.gff3",         # within-rule scratch
    "DANTE_TIR/DANTE_TIR.RData",                  # debug workspace (~35-53 MB)
    "DANTE_TIR_FALLBACK/TPase_5prime_alignment.tsv",
    "DANTE_TIR_FALLBACK/TPase_3prime_alignment.tsv",
    "RepeatMasker/genome_cleaned.fasta",          # staged genome copy (~150-210 MB)
    "TideCluster/genome_cleaned.fasta",           # staged genome copy
    "DANTE_LTR/library/mmseqs2/mmseqs_all_seqs.fasta",
    "DANTE_LTR/library/mmseqs2/partitioned_s900_w1000.fasta",
    "DANTE_LTR/library/TE*.fasta",                # TE_DL/DLP/DLT/DLTP/DLplus/all
    "DANTE_LTR/LTR_RTs_library.fasta.reformatted",
    "DANTE_LTR/LTR_RTs_library.fasta.blast.csv",
    "DANTE_LTR/LTR_RTs_library.fasta.filtered_ids",
)

# maximal: additionally purge the big TideCluster trees (kite/tarean can be
# multiple GB; TideCluster_report.html does NOT deep-link them) and the tool
# workdirs / mmseqs tmp. All verified unconsumed downstream.
_MAXIMAL = (
    "TideCluster/*/TideCluster_tarean",
    "TideCluster/*/TideCluster_kite",
    "TideCluster/*/TideCluster_consensus",
    "Libraries/workdir",
    "Libraries/containment_workdir",
    "RepeatMasker/workdir",
    "DANTE_TIR/mmseqs_combined",
    "DANTE_TIR/fallback_library_workdir",
    "DANTE_TIR/mmseqs2",
    "DANTE_LINE/mmseqs",
)


def _human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024:
            return f"{n:.1f}{unit}"
        n /= 1024
    return f"{n:.1f}PB"


def _size_bytes(path: Path) -> int:
    if path.is_symlink():
        return 0
    if path.is_file():
        try:
            return path.stat().st_size
        except OSError:
            return 0
    total = 0
    for p in path.rglob("*"):
        if p.is_file() and not p.is_symlink():
            try:
                total += p.stat().st_size
            except OSError:
                pass
    return total


def build_keep_set(root: Path) -> set:
    """Absolute paths (each entry and its symlink-resolved target) never to
    delete: manifest outputs + CI/count-checked files + run metadata."""
    keep = set()

    def add(rel: str):
        p = root / rel
        keep.add(os.path.normpath(str(p)))
        try:
            keep.add(os.path.normpath(os.path.realpath(str(p))))
        except OSError:
            pass

    for rel in OUTPUTS.values():
        add(rel.rstrip("/"))
    for rel in _EXTRA_KEEP:
        add(rel)
    # Top-level deliverables are often symlinks into per-tool subdirs (container
    # layout); protect every top-level symlink's real target so cleanup can never
    # dangle a deliverable regardless of the symlink-vs-real-file layout.
    try:
        for entry in root.iterdir():
            if entry.is_symlink():
                add(os.path.relpath(str(entry), str(root)))
    except OSError:
        pass
    return keep


def _is_protected(path: Path, keep: set) -> bool:
    ap = os.path.normpath(str(path))
    if ap in keep:
        return True
    try:
        if os.path.normpath(os.path.realpath(str(path))) in keep:
            return True
    except OSError:
        pass
    # never delete a directory that IS an ancestor of a kept path
    prefix = ap + os.sep
    return any(k.startswith(prefix) for k in keep)


def cleanup(output_dir, mode="minimal", dry_run=False, log=print) -> int:
    """Delete the mode's scratch under output_dir. Returns bytes freed."""
    if mode not in ("minimal", "maximal"):
        return 0                                   # none / unknown -> no-op
    root = Path(output_dir).resolve()
    if not root.is_dir():
        log(f"[cleanup] output_dir {root} not found; nothing to do")
        return 0
    keep = build_keep_set(root)
    targets = list(_MINIMAL) + (list(_MAXIMAL) if mode == "maximal" else [])

    freed = 0
    removed = 0
    for rel in targets:
        for path in sorted(root.glob(rel)):
            try:
                path.relative_to(root)             # confine to output_dir
            except ValueError:
                continue
            if _is_protected(path, keep):
                log(f"  [cleanup] KEEP (protected): {path.relative_to(root)}")
                continue
            size = _size_bytes(path)
            rp = path.relative_to(root)
            if dry_run:
                log(f"  [cleanup] would remove {rp} ({_human(size)})")
                freed += size
                removed += 1
                continue
            try:
                if path.is_symlink() or path.is_file():
                    path.unlink()
                elif path.is_dir():
                    shutil.rmtree(path)
                else:
                    continue
                freed += size
                removed += 1
                log(f"  [cleanup] removed {rp} ({_human(size)})")
            except OSError as e:
                log(f"  [cleanup] WARN could not remove {rp}: {e}")
    verb = "would free" if dry_run else "freed"
    log(f"[cleanup] mode={mode}: {removed} item(s), {verb} {_human(freed)}")
    return freed


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Post-run cleanup of CARP intermediates.")
    ap.add_argument("output_dir")
    ap.add_argument("--mode", choices=VALID_MODES, default="minimal")
    ap.add_argument("--dry-run", action="store_true",
                    help="list what would be removed without deleting")
    a = ap.parse_args(argv)
    cleanup(a.output_dir, a.mode, dry_run=a.dry_run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
