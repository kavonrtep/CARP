#!/usr/bin/env python3
"""Unit test for scripts/assert_library_canonical_order.py — the deterministic
regression guard for the RepeatMasker-library canonical-order invariant. Verifies
it (a) passes on a canonically-sorted library, (b) FAILS on a scrambled library
of the same content, and (c) skips (passes) when the library is absent.
Pure Python, no pipeline — runs in the cheap unit env.
"""
import importlib.util
import sys
import tempfile
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location(
    "assert_library_canonical_order", SCRIPTS / "assert_library_canonical_order.py")
alco = importlib.util.module_from_spec(spec)
spec.loader.exec_module(alco)

# Three records whose canonical (sequence-sorted) order differs from this
# input order, so a scrambled write is detectably non-canonical.
SCRAMBLED = (
    ">seqB#Class_I/LINE\nTTTTGGGGTTTTGGGG\n"
    ">seqA#Class_I/LINE\nAAAACCCCAAAACCCC\n"
    ">seqC#Class_I/LINE\nCCCCAAAACCCCAAAA\n"
)


def _write_library(output_dir: Path, content: str) -> None:
    lib = output_dir / alco.DEFAULT_REL
    lib.parent.mkdir(parents=True, exist_ok=True)
    lib.write_text(content)


def _canonicalize(output_dir: Path) -> None:
    """Rewrite the library into canonical order via the real sorter (as the
    reduce_library_containment rule does)."""
    lib = output_dir / alco.DEFAULT_REL
    tmp = str(lib) + ".sorted"
    alco.sort_fasta_by_sequence(str(lib), tmp, threads=1, tmpdir=str(lib.parent))
    Path(tmp).replace(lib)


def main() -> int:
    with tempfile.TemporaryDirectory() as d:
        od = Path(d)

        # (a) canonical library -> OK (0)
        _write_library(od, SCRAMBLED)
        _canonicalize(od)
        assert alco.main([str(od)]) == 0, "canonical library should pass"

        # (b) scrambled library of the same content -> FAIL (1)
        _write_library(od, SCRAMBLED)
        assert alco.main([str(od)]) == 1, "scrambled library should fail"

        # (c) missing library -> SKIP (0)
        (od / alco.DEFAULT_REL).unlink()
        assert alco.main([str(od)]) == 0, "absent library should skip (pass)"

    print("test_assert_library_canonical_order: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
