#!/usr/bin/env python3
"""User-supplied library validation (`classification.py validate --mode fasta`).

The `custom_library` / `tandem_repeat_library` config inputs are the only
classification-bearing files the pipeline does NOT produce, so nothing
canonicalises them: `concatenate_libraries` cats the custom library into the
RepeatMasker library verbatim. Before the `validate_classifications` rule
covered them, a curated library carrying the legacy flat class `rDNA_45S/18S`
put 9 non-vocabulary records into `Repeat_Annotation_Unified.gff3` beside the
canonical `rDNA/45S_rDNA/*` ones, splitting rDNA into parallel
`summary_statistics.csv` rows that do not sum.

These tests pin the behaviour the rule relies on:
  1. a library whose classes are all canonical exits 0,
  2. a library with a non-vocabulary class exits non-zero and names it,
  3. the specific legacy rDNA spellings that caused the bug are rejected,
  4. an empty / header-only library does not crash.

Exit 0 on success, 1 with a message on failure.
"""
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parent.parent
CLS = ROOT / "scripts" / "classification.py"

FAILURES = []


def check(label, condition, detail=""):
    if condition:
        print(f"OK: {label}")
    else:
        print(f"FAIL: {label} {detail}")
        FAILURES.append(label)


def run_validate(text):
    """Run `validate --mode fasta` on a temporary FASTA; return (rc, output)."""
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as fh:
        fh.write(text)
        path = fh.name
    try:
        p = subprocess.run(
            [sys.executable, str(CLS), "validate", "--mode", "fasta",
             "--source", "RepeatMasker", path],
            capture_output=True, text=True, cwd=str(ROOT),
        )
        return p.returncode, p.stdout + p.stderr
    finally:
        pathlib.Path(path).unlink(missing_ok=True)


def main():
    # 1. All-canonical library passes.
    good = (
        ">seq1#rDNA/45S_rDNA/18S\nACGT\n"
        ">seq2#rDNA/45S_rDNA/ITS1\nACGT\n"
        ">seq3#rDNA/5S_rDNA/5S\nACGT\n"
        ">seq4#Class_I/LTR/Ty1_copia/Ale\nACGT\n"
        ">seq5#Class_II/Subclass_1/TIR/EnSpm_CACTA\nACGT\n"
    )
    rc, out = run_validate(good)
    check("canonical library exits 0", rc == 0, f"(rc={rc}, out={out.strip()[:200]})")

    # 2. A single bad class fails the whole file, and is named in the output.
    bad = good + ">seq6#Not_A_Real_Class/Nope\nACGT\n"
    rc, out = run_validate(bad)
    check("non-vocabulary class exits non-zero", rc != 0, f"(rc={rc})")
    check("offending class is named in the output",
          "Not_A_Real_Class/Nope" in out, f"(out={out.strip()[:200]})")

    # 3. The exact legacy spellings that caused the bug must be rejected.
    #    These are the classes found in the curated Pisum library; the
    #    canonical forms are rDNA/45S_rDNA/<subunit> and rDNA/5S_rDNA/5S.
    legacy = ["rDNA_45S/18S", "rDNA_45S/25S", "rDNA_45S/5.8S",
              "rDNA_45S/IGS", "rDNA_45S/ITS1", "rDNA_45S/ITS2",
              "rDNA_45S", "rDNA_5S"]
    for cls in legacy:
        rc, out = run_validate(f">s#{cls}\nACGT\n")
        check(f"legacy {cls!r} rejected", rc != 0, f"(rc={rc})")

    # 4. The canonical counterparts must pass, so the fix is actually reachable.
    for cls in ["rDNA/45S_rDNA/18S", "rDNA/45S_rDNA/25S", "rDNA/45S_rDNA/5.8S",
                "rDNA/45S_rDNA/IGS", "rDNA/45S_rDNA/ITS1", "rDNA/45S_rDNA/ITS2",
                "rDNA/45S_rDNA", "rDNA/5S_rDNA"]:
        rc, out = run_validate(f">s#{cls}\nACGT\n")
        check(f"canonical {cls!r} accepted", rc == 0, f"(rc={rc}, out={out.strip()[:160]})")

    # 5. Degenerate inputs must not crash (the rule guards on -s, but a file
    #    that exists with only whitespace should still exit cleanly).
    rc, out = run_validate("")
    check("empty library does not crash", rc in (0, 1), f"(rc={rc})")

    print()
    if FAILURES:
        print(f"FAILED: {len(FAILURES)} check(s): {', '.join(FAILURES)}")
        return 1
    print("test_validate_library_classes: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
