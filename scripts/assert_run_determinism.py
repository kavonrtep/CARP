#!/usr/bin/env python3
"""Compare two CARP run output directories and fail if they are not deterministic.

Determinism enforcement: run the same input through the pipeline twice (ideally
under a different PYTHONHASHSEED and thread count), then run this on the two
output dirs. Every determinism-critical output must be byte-identical; any
difference is a determinism regression and exits non-zero.

Scope = the manifest ``OUTPUTS`` map (the curated, downstream-consumed outputs),
imported from ``scripts/manifest.py`` so this needs no per-run manifest file and
auto-covers any output added there. Known-volatile entries are excluded and GFF3
comment/header lines are normalized:

  * excluded: ``provenance`` (run_provenance.json) and ``report_*`` (HTML reports)
    — they legitimately embed timestamps/host/paths.
  * GFF3 files: lines beginning with ``#`` (e.g. ``##date``) are dropped before
    comparing, so a header date does not read as non-determinism; the feature
    lines must still match exactly.
  * everything else (CSV, BED, BigWig, FASTA, per-family track dirs) is compared
    byte-for-byte.

Usage:
    assert_run_determinism.py <output_dir_a> <output_dir_b>
Exit 0 = deterministic; 1 = drift (with the differing outputs listed); 2 = usage.
"""
import hashlib
import importlib.util
import sys
from pathlib import Path

_SCRIPTS = Path(__file__).resolve().parent

# logical keys in OUTPUTS whose *content* is legitimately volatile (timestamps,
# host, absolute paths) and must not be compared byte-for-byte.
_EXCLUDE_KEYS = {"provenance"}
_EXCLUDE_KEY_PREFIXES = ("report_",)


def _load_outputs_map():
    """The OUTPUTS map (logical name -> path relative to output_dir)."""
    spec = importlib.util.spec_from_file_location("manifest", _SCRIPTS / "manifest.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.OUTPUTS


def _digest(path: Path, gff3: bool) -> str:
    if gff3:
        with path.open("rb") as fh:
            data = b"".join(ln for ln in fh if not ln.lstrip().startswith(b"#"))
        return hashlib.md5(data).hexdigest()
    return hashlib.md5(path.read_bytes()).hexdigest()


def _entry_files(base: Path, rel: str):
    """(relpath, Path) for a manifest entry that may be a file or a directory."""
    target = base / rel
    if rel.endswith("/") or target.is_dir():
        if target.is_dir():
            for f in sorted(target.rglob("*")):
                if f.is_file():
                    yield str(f.relative_to(base)), f
    elif target.exists():
        yield rel.rstrip("/"), target


def compare(out_a: Path, out_b: Path):
    outputs = _load_outputs_map()
    diffs, compared = [], 0
    for key, rel in outputs.items():
        if key in _EXCLUDE_KEYS or key.startswith(_EXCLUDE_KEY_PREFIXES):
            continue
        for relpath, fa in _entry_files(out_a, rel):
            fb = out_b / relpath
            gff3 = relpath.endswith(".gff3")
            if not fb.exists():
                diffs.append(f"{relpath}: present in run A, missing in run B")
                continue
            compared += 1
            if _digest(fa, gff3) != _digest(fb, gff3):
                diffs.append(relpath + (" [gff3 data lines]" if gff3 else ""))
    return compared, diffs


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: assert_run_determinism.py <output_dir_a> <output_dir_b>", file=sys.stderr)
        return 2
    out_a, out_b = Path(sys.argv[1]), Path(sys.argv[2])
    for d in (out_a, out_b):
        if not d.is_dir():
            print(f"FAIL: not a directory: {d}", file=sys.stderr)
            return 2
    compared, diffs = compare(out_a, out_b)
    print(f"determinism: compared {compared} manifest output files "
          f"({out_a} vs {out_b})")
    if diffs:
        print(f"NON-DETERMINISM DETECTED — {len(diffs)} output(s) differ between "
              f"two runs of the same input:", file=sys.stderr)
        for d in diffs[:100]:
            print(f"  - {d}", file=sys.stderr)
        print("\nA determinism-critical output changed run-to-run. Fix the source "
              "(see the Determinism section in CLAUDE.md) — do not merge.", file=sys.stderr)
        return 1
    print("OK: all compared outputs are byte-identical (GFF3 headers normalized) "
          "— the run is deterministic.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
