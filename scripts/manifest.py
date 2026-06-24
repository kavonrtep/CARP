#!/usr/bin/env python3
"""CARP output manifest — ``<output_dir>/carp_manifest.json``.

A machine-readable statement of the *output contract* a run conforms to, written
at the output root so downstream consumers (e.g. the JBrowse portal) resolve
files by a stable logical name instead of sniffing the directory layout.
See FR-3 in ``docs/carp_feature_requests.md`` and the ``schema_version``
changelog in ``docs/output_schema.md``.

This is SEPARATE from ``run_provenance.json`` and is NOT a completion gate:
completion is decided by the provenance ``exit_status`` + the platform's
exit-code markers + output verification. The manifest only declares schema +
file map, so a manifest write hiccup can never mask a real run outcome.

Written by the wrapper (``run_pipeline.py``), mirroring provenance, so the
manifest is guaranteed present with the right status on success AND failure:
  * ``init_manifest()``     before snakemake -> ``exit_status: "running"`` + full map
  * ``finalise_manifest()`` after snakemake  -> flips ``exit_status`` to
    ``"completed"`` / ``"failed"`` (before scratch->NFS copy-back and the
    exit-code DONE markers, so it travels with the outputs).

Keep the authoritative write here in the wrapper — not scattered across
Snakemake rules — so manifest and provenance can't drift.
"""
import json
import sys
from pathlib import Path

# Reuse the provenance module's helpers (same dir) so the two metadata writers
# share one implementation of timestamp / git-ref / atomic-write / version.
_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
from record_provenance import (  # noqa: E402
    PIPELINE_VERSION, _utcnow_iso, _git_sha, _atomic_write_json)

MANIFEST_VERSION = 1

# OUTPUT-CONTRACT version. Bump ONLY on a breaking change to the layout or
# contents of consumed outputs — NOT on every CARP release. Changelog:
# docs/output_schema.md. (v2 == the rc6 Repeat_density_* layout.)
OUTPUT_SCHEMA_VERSION = "2"

# Stable logical name -> path relative to the output root. Every entry is
# produced unconditionally by a successful run. Logical names stay constant
# across schema versions where the meaning is unchanged; a renamed/removed
# output is a logical-name change PLUS an OUTPUT_SCHEMA_VERSION bump. Directory
# entries end with "/".
OUTPUTS = {
    # ── core annotation ──────────────────────────────────────────────
    "unified_gff3":                 "Repeat_Annotation_Unified.gff3",
    "cleaned_fasta":                "genome_cleaned.fasta",
    "summary_statistics":           "summary_statistics.csv",
    # ── density BigWigs (Unified-sourced; sparse) ────────────────────
    "density_total_dir":            "Repeat_density/",
    "density_by_class_dir":         "Repeat_density_by_class_bigwig/",
    "tandem_unified_split_dir":     "Tandem_repeats_unified_split_by_family_bigwig/",
    "tandem_tidecluster_split_dir": "Tandem_repeats_TideCluster_split_by_family_bigwig/",
    "tandem_tidecluster_total_10k":  "Tandem_repeats_TideCluster_10k.bw",
    "tandem_tidecluster_total_100k": "Tandem_repeats_TideCluster_100k.bw",
    # ── per-tool / per-class GFF3s ───────────────────────────────────
    "dante_filtered_gff3":          "DANTE_filtered.gff3",
    "dante_ltr_gff3":               "DANTE_LTR.gff3",
    "dante_tir_gff3":               "DANTE_TIR.gff3",
    "mobile_elements_gff3":         "Mobile_elements_RepeatMasker.gff3",
    "tandem_tidecluster_gff3":      "Tandem_repeats_TideCluster.gff3",
    "tandem_repeatmasker_gff3":     "Tandem_repeats_RepeatMasker.gff3",
    # ── masking / gaps ───────────────────────────────────────────────
    "masking_bed":                  "all_repeats_for_masking.bed",
    "gaps_bed":                     "gaps_10plus.bed",
    # ── reports ──────────────────────────────────────────────────────
    "report_main":                  "repeat_annotation_report.html",
    "report_tidecluster":           "TideCluster_report.html",
    "report_dante_ltr":             "DANTE_LTR_report.html",
    "report_benchmark":             "benchmark_report.html",
    # ── cross-reference ──────────────────────────────────────────────
    "provenance":                   "run_provenance.json",
}

_MANIFEST_NAME = "carp_manifest.json"


def _git_ref() -> str:
    sha, dirty = _git_sha()
    return sha + ("-dirty" if dirty else "")


def _base_record() -> dict:
    return {
        "manifest_version": MANIFEST_VERSION,
        "schema_version": OUTPUT_SCHEMA_VERSION,
        "carp_version": PIPELINE_VERSION,
        "carp_git_ref": _git_ref(),
        "produced_at": _utcnow_iso(),
        "exit_status": "running",
        "outputs": dict(OUTPUTS),
    }


def init_manifest(output_dir) -> Path:
    """Write carp_manifest.json with ``exit_status: "running"`` and the full
    static outputs map. Called before snakemake. Returns the manifest path."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / _MANIFEST_NAME
    _atomic_write_json(out, _base_record())
    return out


def finalise_manifest(output_dir, exit_status: str) -> Path:
    """Flip ``exit_status`` to "completed"/"failed" after snakemake. If init
    never ran (output dir wiped), emit a complete stub — the manifest must
    exist on every run."""
    out = Path(output_dir) / _MANIFEST_NAME
    if out.exists():
        record = json.loads(out.read_text())
    else:
        record = _base_record()
        record["note"] = "init phase missing; manifest written at finalise."
    record["exit_status"] = exit_status
    _atomic_write_json(out, record)
    return out


if __name__ == "__main__":
    # CLI: manifest.py <output_dir> [exit_status]
    if len(sys.argv) < 2:
        sys.exit("usage: manifest.py <output_dir> [exit_status]")
    if len(sys.argv) >= 3:
        finalise_manifest(sys.argv[1], sys.argv[2])
    else:
        init_manifest(sys.argv[1])
