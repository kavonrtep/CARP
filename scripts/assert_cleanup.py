#!/usr/bin/env python3
"""Assert the post-run cleanup (run_pipeline.py) did its job on a finished
output_dir — used in CI after a fixture run to prove cleanup deletes the right
scratch WITHOUT ever touching a deliverable.

Reads the mode actually applied from ``run_provenance.json``
(``config.cleanup_intermediates``) and checks:

  1. every manifest output (``manifest.py:OUTPUTS``) still exists, and every
     top-level symlink still resolves (cleanup must never dangle a deliverable);
  2. no deletable scratch for that mode remains — i.e. cleanup ran to completion
     (a would-be re-clean finds nothing left).

Usage: ``assert_cleanup.py <output_dir>``. Exit 0 = OK, 1 = a violation.
"""
import json
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
from manifest import OUTPUTS  # noqa: E402
from cleanup_outputs import (  # noqa: E402
    _MINIMAL, _MAXIMAL, build_keep_set, _is_protected)


def check_cleanup(root: Path):
    """Return (mode, errors). errors is a list of human-readable strings."""
    root = Path(root).resolve()
    errs = []
    prov_path = root / "run_provenance.json"
    if not prov_path.exists():
        return None, [f"{prov_path} not found (needed to read the applied mode)"]
    prov = json.loads(prov_path.read_text())
    mode = str(prov.get("config", {}).get("cleanup_intermediates", "minimal")
               ).strip().lower()

    # 1. every manifest deliverable survived
    for name, rel in OUTPUTS.items():
        p = root / rel.rstrip("/")
        if not p.exists():
            errs.append(f"manifest output missing after cleanup: {name} ({rel})")
    # every top-level symlink still resolves (no dangling deliverable)
    for entry in sorted(root.iterdir()):
        if entry.is_symlink() and not entry.exists():
            errs.append(f"top-level symlink dangles after cleanup: {entry.name}")

    # 2. no scratch for the applied mode remains → cleanup ran to completion
    if mode in ("minimal", "maximal"):
        keep = build_keep_set(root)
        targets = list(_MINIMAL) + (list(_MAXIMAL) if mode == "maximal" else [])
        remaining = []
        for rel in targets:
            for p in sorted(root.glob(rel)):
                if not _is_protected(p, keep):
                    remaining.append(str(p.relative_to(root)))
        if remaining:
            errs.append(f"cleanup mode={mode} left scratch behind: {remaining}")

    return mode, errs


def main(argv=None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    if not argv:
        print("usage: assert_cleanup.py <output_dir>", file=sys.stderr)
        return 2
    root = Path(argv[0])
    mode, errs = check_cleanup(root)
    if errs:
        print(f"FAIL: cleanup assertion (mode={mode}) in {root}:")
        for e in errs:
            print("  - " + e)
        return 1
    print(f"OK: cleanup (mode={mode}) removed all scratch and kept every "
          f"manifest output in {root}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
