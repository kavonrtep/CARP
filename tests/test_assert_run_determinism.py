#!/usr/bin/env python3
"""Unit test for scripts/assert_run_determinism.py — the determinism gate's
comparator. Verifies it (a) passes on identical runs, (b) ignores GFF3 header/date
differences, (c) FAILS on a GFF3 feature-line change, (d) FAILS on a data-file
(CSV) change, and (e) ignores volatile excluded entries (reports, provenance).
Pure Python, no pipeline — runs in the cheap unit env.
"""
import importlib.util
import sys
import tempfile
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location("assert_run_determinism", SCRIPTS / "assert_run_determinism.py")
ard = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ard)

# resolve the real manifest paths so the test writes files where the comparator looks
OUTPUTS = ard._load_outputs_map()
UNIFIED = OUTPUTS["unified_gff3"]
SUMMARY = OUTPUTS["summary_statistics"]
REPORT = OUTPUTS["report_main"]
PROV = OUTPUTS["provenance"]

GFF3_A = "##gff-version 3\n##date 2026-07-01\nchr1\tCARP\trepeat\t1\t100\t.\t+\t.\tID=x\n"
GFF3_DATE_ONLY = "##gff-version 3\n##date 2026-07-28\nchr1\tCARP\trepeat\t1\t100\t.\t+\t.\tID=x\n"
GFF3_FEATURE_CHANGED = "##gff-version 3\n##date 2026-07-01\nchr1\tCARP\trepeat\t1\t200\t.\t+\t.\tID=x\n"


def _mk(d: Path, unified, summary="a,b\n1,2\n", report="<html>run at 10:00</html>", prov='{"run_started":"t"}'):
    (d / UNIFIED).parent.mkdir(parents=True, exist_ok=True)
    (d / UNIFIED).write_text(unified)
    (d / SUMMARY).write_text(summary)
    (d / REPORT).write_text(report)
    (d / PROV).write_text(prov)


def _run(a, b):
    return ard.compare(Path(a), Path(b))


def test_comparator():
    with tempfile.TemporaryDirectory() as ta, tempfile.TemporaryDirectory() as tb:
        a, b = Path(ta), Path(tb)

        # (a) identical -> no diffs
        _mk(a, GFF3_A); _mk(b, GFF3_A)
        _, diffs = _run(a, b)
        assert diffs == [], f"identical runs flagged: {diffs}"

        # (b) only the GFF3 ##date header differs -> still no diffs (normalized)
        _mk(b, GFF3_DATE_ONLY)
        _, diffs = _run(a, b)
        assert diffs == [], f"GFF3 header-date change wrongly flagged: {diffs}"

        # (e) volatile excluded entries differ -> still no diffs
        _mk(b, GFF3_DATE_ONLY, report="<html>run at 23:59</html>", prov='{"run_started":"other"}')
        _, diffs = _run(a, b)
        assert diffs == [], f"excluded report/provenance change wrongly flagged: {diffs}"

        # (c) a GFF3 FEATURE line changes -> must be flagged
        _mk(b, GFF3_FEATURE_CHANGED)
        _, diffs = _run(a, b)
        assert any("Repeat_Annotation_Unified.gff3" in d for d in diffs), \
            f"GFF3 feature-line change NOT detected: {diffs}"

        # (d) a data CSV changes -> must be flagged
        _mk(b, GFF3_A, summary="a,b\n9,9\n")
        _, diffs = _run(a, b)
        assert any("summary_statistics.csv" in d for d in diffs), \
            f"summary_statistics change NOT detected: {diffs}"
    print("  test_comparator: PASS")


if __name__ == "__main__":
    test_comparator()
    print("test_assert_run_determinism: ALL PASS")
