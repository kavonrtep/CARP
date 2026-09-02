#!/usr/bin/env python3
"""Assert that a run's carp_manifest.json matches what the run produced.

Closes the manifest drift gap: scripts/manifest.py holds a hand-maintained
``OUTPUTS`` map (logical name -> path). Nothing else checks that those paths
correspond to files a successful run actually produces, so renaming or
removing an output in the Snakefile silently leaves the manifest pointing at a
path that no longer exists — and a downstream consumer (e.g. the JBrowse
portal) then resolves a logical name to a missing file.

This script loads ``<output_dir>/carp_manifest.json`` from a completed run and
asserts every entry in ``outputs`` exists on disk (directory entries — those
ending in "/" — must be directories; everything else must exist, following
symlinks, since most top-level outputs are symlinks). It also requires the
manifest's own ``exit_status`` to be "completed", tying the contract to a
real successful run.

Run against the E2E fixture output in CI (gated to the `medium` fixture, the
release gate that produces the full output set):

    python scripts/assert_manifest_outputs.py tests/fixtures/output_medium

Exit 0 = manifest is in sync with produced outputs; 1 = drift (missing paths
or non-completed status), printed to stderr.

This intentionally does NOT enforce the reverse (a produced output absent from
the manifest): the manifest is a curated subset of *consumed* outputs, so
completeness is not a contract. Only the dangerous direction — manifest claims
a path the run did not produce — is enforced.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

_MANIFEST_NAME = "carp_manifest.json"


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: assert_manifest_outputs.py <output_dir>", file=sys.stderr)
        return 2
    output_dir = Path(sys.argv[1])
    manifest_path = output_dir / _MANIFEST_NAME

    if not manifest_path.is_file():
        print(f"FAIL: manifest missing: {manifest_path}", file=sys.stderr)
        return 1

    record = json.loads(manifest_path.read_text())
    problems: list[str] = []

    status = record.get("exit_status")
    if status != "completed":
        problems.append(
            f"exit_status is {status!r}, expected 'completed' "
            "(manifest must reflect a successful run)")

    outputs = record.get("outputs")
    if not isinstance(outputs, dict) or not outputs:
        problems.append("manifest has no non-empty 'outputs' map")
        outputs = {}

    for name, rel in outputs.items():
        target = output_dir / rel
        is_dir_entry = rel.endswith("/")
        # exists() follows symlinks, so a top-level symlink whose target is
        # missing correctly counts as a failure.
        if is_dir_entry:
            ok = target.is_dir()
        else:
            ok = target.exists()
        if not ok:
            kind = "directory" if is_dir_entry else "file"
            problems.append(f"{name}: missing {kind} -> {rel}")

    if problems:
        print(f"FAIL: carp_manifest.json out of sync with {output_dir}",
              file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        print("\nFix: update OUTPUTS in scripts/manifest.py to match the "
              "produced layout (and bump OUTPUT_SCHEMA_VERSION if the change "
              "is a breaking rename/removal — see docs/output_schema.md).",
              file=sys.stderr)
        return 1

    print(f"OK: carp_manifest.json — all {len(outputs)} declared outputs "
          f"present in {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
