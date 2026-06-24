#!/usr/bin/env python3
"""
Structural / contract tests for scripts/manifest.py (carp_manifest.json, FR-3).

Validates the manifest shape, the init->finalise exit_status transition, and
the outputs map invariants — without needing a full pipeline run.

Exit code 0 = all pass; 1 = one or more failures (printed to stderr).

Usage:
    python tests/test_manifest.py
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "scripts"))

import manifest  # noqa: E402
sys.path.insert(0, str(REPO))
from version import __version__ as VERSION  # noqa: E402

failures: list[str] = []


def check(cond: bool, msg: str) -> None:
    if not cond:
        failures.append(msg)


def main() -> int:
    with tempfile.TemporaryDirectory() as td:
        # ── init ─────────────────────────────────────────────────────
        p = manifest.init_manifest(td)
        check(p.name == "carp_manifest.json", "manifest filename")
        m = json.loads(p.read_text())

        for key in ("manifest_version", "schema_version", "carp_version",
                    "carp_git_ref", "produced_at", "exit_status", "outputs"):
            check(key in m, f"missing top-level key: {key}")

        check(m["manifest_version"] == manifest.MANIFEST_VERSION,
              "manifest_version value")
        check(isinstance(m["schema_version"], str),
              "schema_version must be a string")
        check(m["schema_version"] == manifest.OUTPUT_SCHEMA_VERSION,
              "schema_version matches module constant")
        check(m["carp_version"] == VERSION,
              f"carp_version {m['carp_version']!r} != version.py {VERSION!r}")
        check(m["exit_status"] == "running",
              "exit_status is 'running' after init")
        check(m["produced_at"].endswith("Z"), "produced_at is UTC (Z)")

        # ── outputs map invariants ───────────────────────────────────
        outs = m["outputs"]
        check(isinstance(outs, dict) and len(outs) > 0, "outputs non-empty dict")
        for name, rel in outs.items():
            check(isinstance(rel, str) and rel != "", f"empty path for {name}")
            check(not rel.startswith("/"), f"path must be relative: {name}={rel}")
            check(".." not in Path(rel).parts, f"no parent escape: {name}={rel}")
        # the core portal-consumed logical names must be present
        for required in ("unified_gff3", "summary_statistics",
                         "density_total_dir", "density_by_class_dir",
                         "tandem_unified_split_dir", "report_main"):
            check(required in outs, f"core logical name missing: {required}")

        # ── finalise: running -> completed ───────────────────────────
        manifest.finalise_manifest(td, "completed")
        m2 = json.loads(p.read_text())
        check(m2["exit_status"] == "completed", "exit_status flips to completed")
        check(m2["outputs"] == outs, "outputs map unchanged by finalise")

    # ── finalise without init writes a complete stub ─────────────────
    with tempfile.TemporaryDirectory() as td2:
        manifest.finalise_manifest(td2, "failed")
        sp = Path(td2, "carp_manifest.json")
        check(sp.exists(), "finalise writes manifest even without init")
        ms = json.loads(sp.read_text())
        check(ms["exit_status"] == "failed", "stub exit_status = failed")
        check(ms["schema_version"] == manifest.OUTPUT_SCHEMA_VERSION,
              "stub has schema_version")

    if failures:
        print("FAIL: test_manifest", file=sys.stderr)
        for f in failures:
            print("  - " + f, file=sys.stderr)
        return 1
    print("OK: test_manifest — all manifest checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
