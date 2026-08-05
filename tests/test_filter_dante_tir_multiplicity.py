#!/usr/bin/env python3
"""
Unit tests for scripts/filter_dante_tir_by_multiplicity.py — the Multiplicity
floor applied to the primary DANTE_TIR library in rule make_tir_combined_library.

Regression origin: the GFF3<->FASTA join was unsatisfiable. DANTE_TIR names the
FASTA record by stripping ``Class_II_Subclass_1_TIR_`` off the GFF3 ``ID``, so
``ID=Class_II_Subclass_1_TIR_hAT_2`` is ``>hAT_2#...`` in the FASTA, and the
script compared the two verbatim. It therefore kept 0 records on every input
(measured 0/7233 on a 90 Gb genome, 0/N on all 18 runs available locally) since
the script was introduced, and — because dante_tir_min_multiplicity defaults to
3 — that emptied the TIR library on every default run. The rule then touched an
empty library and exited 0, so it was silent all the way to an annotation that
contained no similarity-called DNA transposons and an LTR library whose
Class_II screen had degraded to a no-op.

The pipeline fixture cannot catch this: output_medium's records are all
Multiplicity=1, so the broken result (0 kept) and the correct result (0 kept)
are identical. Hence these direct tests.

Exit 0 = all pass; 1 = one or more failures (printed to stderr).
Usage: python tests/test_filter_dante_tir_multiplicity.py
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
SCRIPT = REPO / "scripts" / "filter_dante_tir_by_multiplicity.py"
sys.path.insert(0, str(REPO / "scripts"))

import filter_dante_tir_by_multiplicity as f  # noqa: E402

failures: list[str] = []


def check(cond: bool, msg: str) -> None:
    if not cond:
        failures.append(msg)


def check_call(fn, args: tuple, expected, msg: str) -> None:
    """check() for a direct API call, recording a missing/raising helper as a
    failure instead of aborting the run — so a regressed script still reports
    every behavioural failure rather than dying on the first AttributeError."""
    try:
        got = fn() if fn is None else fn(*args)
    except Exception as exc:  # noqa: BLE001 - any failure here is a test failure
        failures.append(f"{msg}: raised {type(exc).__name__}: {exc}")
        return
    check(got == expected, f"{msg}: expected {expected!r}, got {got!r}")


def api(name: str):
    """Resolve a helper on the module under test, or None if it is absent."""
    return getattr(f, name, None)


def gff_row(attrs: str, ftype: str = "sequence_feature") -> str:
    return "\t".join(["chr1", "DANTE_TIR", ftype, "1", "100", ".", "+", ".", attrs])


def run_filter(gff_text: str, fasta_text: str, threshold: int,
               keep_missing: bool = False):
    """Run the script as a subprocess; return (returncode, stderr, kept_names)."""
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        (d / "in.gff3").write_text(gff_text)
        (d / "in.fasta").write_text(fasta_text)
        cmd = [sys.executable, str(SCRIPT), "--gff", str(d / "in.gff3"),
               "--fasta-in", str(d / "in.fasta"), "--fasta-out", str(d / "out.fasta"),
               "--min-multiplicity", str(threshold)]
        if keep_missing:
            cmd.append("--keep-missing-multiplicity")
        p = subprocess.run(cmd, capture_output=True, text=True)
        out = d / "out.fasta"
        kept = []
        if out.exists():
            kept = [ln[1:].strip() for ln in out.read_text().splitlines()
                    if ln.startswith(">")]
        return p.returncode, p.stderr, kept


# ---------------------------------------------------------------- join key ---
# 0.3.0-style: no Name= attribute, so the FASTA name must be recovered by
# stripping the classification prefix off ID=. This is the regression.
GFF_LEGACY = "\n".join([
    "##gff-version 3",
    gff_row("ID=Class_II_Subclass_1_TIR_hAT_2;Classification=Class_II_Subclass_1_TIR_hAT;Multiplicity=5"),
    gff_row("ID=Class_II_Subclass_1_TIR_hAT_4;Classification=Class_II_Subclass_1_TIR_hAT;Multiplicity=1"),
    gff_row("ID=Class_II_Subclass_1_TIR_EnSpm_CACTA_9;Classification=Class_II_Subclass_1_TIR_EnSpm_CACTA;Multiplicity=3"),
]) + "\n"

# Raw DANTE_TIR header form (TIR_hAT) — what DANTE_TIR writes.
FASTA_RAW = "".join([
    ">hAT_2#Class_II/Subclass_1/TIR_hAT\nACGT\n",
    ">hAT_4#Class_II/Subclass_1/TIR_hAT\nACGT\n",
    ">EnSpm_CACTA_9#Class_II/Subclass_1/TIR_EnSpm_CACTA\nACGT\n",
])

# Canonicalised header form (TIR/hAT) — what the rule actually feeds the script,
# after classification.py canonicalise-fasta-headers.
FASTA_CANON = "".join([
    ">hAT_2#Class_II/Subclass_1/TIR/hAT\nACGT\n",
    ">hAT_4#Class_II/Subclass_1/TIR/hAT\nACGT\n",
    ">EnSpm_CACTA_9#Class_II/Subclass_1/TIR/EnSpm_CACTA\nACGT\n",
])

for label, fasta in (("raw", FASTA_RAW), ("canonicalised", FASTA_CANON)):
    rc, err, kept = run_filter(GFF_LEGACY, fasta, threshold=3)
    check(rc == 0, f"legacy/{label}: expected exit 0, got {rc} ({err.strip()})")
    check(sorted(kept) == ["EnSpm_CACTA_9#Class_II/Subclass_1/TIR" +
                           ("_" if label == "raw" else "/") + "EnSpm_CACTA",
                           "hAT_2#Class_II/Subclass_1/TIR" +
                           ("_" if label == "raw" else "/") + "hAT"],
          f"legacy/{label}: wrong records kept at threshold 3: {kept}")
    check("joined to a GFF3 parent    : 3" in err,
          f"legacy/{label}: expected all 3 records to join, stderr was:\n{err}")

# The exact regression: keying on the raw GFF3 ID must never be how this works.
name_of = api("fasta_name_from_attrs")
check(name_of is not None,
      "fasta_name_from_attrs() is missing — the GFF3<->FASTA join key is not derived")
if name_of is not None:
    check_call(name_of, ("ID=Class_II_Subclass_1_TIR_hAT_2;Multiplicity=5",), "hAT_2",
               "prefix strip: ID= must reduce to the bare FASTA name")

    # ----------------------------------------------------- Name= join key ----
    # DANTE_TIR >= 0.3.1 publishes the FASTA name as Name=; it must win over ID=.
    check_call(name_of, ("ID=Class_II_Subclass_1_TIR_hAT_2;Name=hAT_2;Multiplicity=5",),
               "hAT_2", "Name=: should be used when present")
    check_call(name_of, ("ID=Class_II_Subclass_1_TIR_WEIRD_1;Name=custom_name_7",),
               "custom_name_7", "Name=: must take precedence over the ID= prefix-strip")
    check_call(name_of, ("Classification=X;Multiplicity=2",), None,
               "a row with neither ID= nor Name= must yield None")

GFF_031 = "\n".join([
    "##gff-version 3",
    gff_row("ID=Class_II_Subclass_1_TIR_hAT_2;Name=hAT_2;Multiplicity=5"),
    gff_row("ID=Class_II_Subclass_1_TIR_hAT_4;Name=hAT_4;Multiplicity=2"),
]) + "\n"
FASTA_031 = ">hAT_2#Class_II/Subclass_1/TIR/hAT\nACGT\n>hAT_4#Class_II/Subclass_1/TIR/hAT\nACGT\n"
rc, err, kept = run_filter(GFF_031, FASTA_031, threshold=3)
check(rc == 0 and kept == ["hAT_2#Class_II/Subclass_1/TIR/hAT"],
      f"0.3.1 Name= join: expected only hAT_2 kept, got {kept}")

# ------------------------------------------------- threshold is inclusive ----
rc, _, kept = run_filter(GFF_LEGACY, FASTA_RAW, threshold=1)
check(len(kept) == 3, f"threshold=1 must keep everything, got {len(kept)}")
rc, _, kept = run_filter(GFF_LEGACY, FASTA_RAW, threshold=5)
check(len(kept) == 1, f"threshold=5 must keep only Multiplicity=5, got {len(kept)}")
rc, _, kept = run_filter(GFF_LEGACY, FASTA_RAW, threshold=6)
check(kept == [], f"threshold above every value must keep nothing, got {kept}")

# ------------------------------------- missing Multiplicity is its own state --
# DANTE_TIR < 0.3.1 clustered before round4(), so round-4 elements carry no
# Multiplicity. Defaulting them to 1 silently dropped 25% of a real library.
GFF_MISSING = "\n".join([
    "##gff-version 3",
    gff_row("ID=Class_II_Subclass_1_TIR_hAT_2;Multiplicity=5"),
    gff_row("ID=Class_II_Subclass_1_TIR_hAT_4"),          # round-4: no Multiplicity
]) + "\n"
FASTA_MISSING = ">hAT_2#Class_II/Subclass_1/TIR/hAT\nACGT\n>hAT_4#Class_II/Subclass_1/TIR/hAT\nACGT\n"

rc, err, kept = run_filter(GFF_MISSING, FASTA_MISSING, threshold=3)
check(len(kept) == 1, f"missing Multiplicity dropped by default: got {kept}")
check("no Multiplicity attribute  : 1" in err,
      f"missing-Multiplicity count must be reported, stderr was:\n{err}")
check("WARNING" in err,
      f"dropping for a missing attribute must warn, stderr was:\n{err}")

rc, err, kept = run_filter(GFF_MISSING, FASTA_MISSING, threshold=3, keep_missing=True)
check(len(kept) == 2,
      f"--keep-missing-multiplicity must retain the round-4 element: got {kept}")
check("WARNING" not in err, "no warning expected when missing records are kept")

# Values are compared numerically, not lexically (a naive string compare makes
# "10" < "3"). Multiplicity=10 must clear a threshold of 3.
GFF_NUM = "##gff-version 3\n" + gff_row(
    "ID=Class_II_Subclass_1_TIR_hAT_2;Multiplicity=10") + "\n"
rc, _, kept = run_filter(GFF_NUM, ">hAT_2#Class_II/Subclass_1/TIR/hAT\nACGT\n", 3)
check(len(kept) == 1, "Multiplicity must be compared numerically (10 >= 3)")

# ---------------------------------------------------- fail loud, not empty ----
# A join that matches nothing is a bug. It must exit non-zero so the calling
# rule's `set -euo pipefail` stops it, rather than touching an empty library.
GFF_UNJOINABLE = "##gff-version 3\n" + gff_row(
    "ID=Class_II_Subclass_1_TIR_hAT_2;Multiplicity=5") + "\n"
rc, err, _ = run_filter(GFF_UNJOINABLE,
                        ">completely_unrelated_name#Class_II/Subclass_1/TIR/hAT\nACGT\n", 3)
check(rc != 0, "a total join failure must exit non-zero")
check("no FASTA record joined" in err,
      f"total join failure must say so, stderr was:\n{err}")

# An empty input FASTA is a legitimate empty result, not a join failure.
rc, err, kept = run_filter(GFF_UNJOINABLE, "", 3)
check(rc == 0 and kept == [], f"empty FASTA input must exit 0, got rc={rc}")

# Multi-line sequences must be carried through whole, not just the header.
GFF_ML = "##gff-version 3\n" + gff_row(
    "ID=Class_II_Subclass_1_TIR_hAT_2;Multiplicity=5") + "\n"
with tempfile.TemporaryDirectory() as td:
    d = Path(td)
    (d / "g.gff3").write_text(GFF_ML)
    (d / "f.fasta").write_text(">hAT_2#Class_II/Subclass_1/TIR/hAT\nAAAA\nCCCC\nGG\n")
    subprocess.run([sys.executable, str(SCRIPT), "--gff", str(d / "g.gff3"),
                    "--fasta-in", str(d / "f.fasta"), "--fasta-out", str(d / "o.fasta"),
                    "--min-multiplicity", "3"], capture_output=True, text=True)
    body = "".join(ln.strip() for ln in (d / "o.fasta").read_text().splitlines()
                   if not ln.startswith(">"))
    check(body == "AAAACCCCGG", f"multi-line sequence must survive intact, got {body!r}")

# Only sequence_feature rows carry Multiplicity; other feature types (e.g. the
# per-element TIR children) must not be treated as parents.
GFF_TYPES = "\n".join([
    "##gff-version 3",
    gff_row("ID=Class_II_Subclass_1_TIR_hAT_2;Multiplicity=5"),
    gff_row("ID=Class_II_Subclass_1_TIR_hAT_9;Multiplicity=9", ftype="terminal_inverted_repeat"),
]) + "\n"
loader = api("load_multiplicity")
if loader is None:
    failures.append("load_multiplicity() is missing from the module under test")
else:
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "g.gff3"
        p.write_text(GFF_TYPES)
        m = loader(p)
    check(set(m) == {"hAT_2"},
          f"only sequence_feature rows may be parents, got {sorted(m)}")


if failures:
    print(f"FAILED {len(failures)} check(s):", file=sys.stderr)
    for msg in failures:
        print(f"  - {msg}", file=sys.stderr)
    sys.exit(1)
print("test_filter_dante_tir_multiplicity: all checks passed")
