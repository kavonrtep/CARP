#!/usr/bin/env python3
"""
Runs scripts/classification.py against tests/classification_cases.tsv.

Exit code 0 = all pass; 1 = one or more failures (printed to stderr).

Usage:
    python tests/test_classification.py
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "scripts"))

import classification  # noqa: E402
from classification import UnknownClassification, canonicalise, load_vocabulary  # noqa: E402


def parse_cases(path: Path):
    with path.open() as fh:
        for line_no, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line or line.lstrip().startswith("#"):
                continue
            if line.startswith("source\t"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            src = parts[0] or None
            raw = parts[1]
            expected = parts[2]
            notes = parts[3] if len(parts) > 3 else ""
            yield line_no, src, raw, expected, notes


def main() -> int:
    vocab = load_vocabulary()
    cases_path = HERE / "classification_cases.tsv"
    failures: list[str] = []
    passed = 0

    for line_no, src, raw, expected, notes in parse_cases(cases_path):
        if expected == "!!ERROR":
            try:
                result = canonicalise(raw, source=src, vocab=vocab)
            except (UnknownClassification, ValueError):
                passed += 1
                continue
            failures.append(
                f"L{line_no}: expected error for source={src!r} raw={raw!r} "
                f"but got {result!r}"
            )
            continue

        try:
            result = canonicalise(raw, source=src, vocab=vocab)
        except (UnknownClassification, ValueError) as e:
            failures.append(
                f"L{line_no}: unexpected error for source={src!r} raw={raw!r}: {e}"
            )
            continue

        if result != expected:
            failures.append(
                f"L{line_no}: source={src!r} raw={raw!r}\n"
                f"  expected: {expected!r}\n"
                f"  got:      {result!r}\n"
                f"  ({notes})"
            )
        else:
            passed += 1

    total = passed + len(failures)
    print(f"{passed}/{total} cases passed", file=sys.stderr)
    if failures:
        for f in failures:
            print(f"FAIL: {f}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
