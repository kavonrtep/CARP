#!/usr/bin/env python3
"""
Unit tests for scripts/containment_reduce_library.py (the second-round
containment reduction of the RepeatMasker library).

Exercises the greedy within-class containment logic directly against a
synthetic self-blast table — no blastn required, so it runs in the lightweight
unit CI. Verifies the invariants that make the reduction masking-/class-safe.

Exit 0 = all pass; 1 = one or more failures (printed to stderr).
Usage: python tests/test_containment_reduce.py
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "scripts"))

import containment_reduce_library as cr  # noqa: E402

failures: list[str] = []


def check(cond: bool, msg: str) -> None:
    if not cond:
        failures.append(msg)


# Synthetic library: name#class -> sequence (content irrelevant, only length is)
records = [
    ("elemA#ClassX", "A" * 1000),   # full-length element, longest -> always kept
    ("fragA#ClassX", "A" * 300),    # contained in elemA, same class -> REMOVE
    ("fragB#ClassY", "A" * 300),    # contained in elemA but DIFFERENT class -> KEEP
    ("fragC#ClassX", "A" * 300),    # only 50% covered -> KEEP
    ("fragD#ClassX", "A" * 300),    # full coverage but 70% identity -> KEEP
    ("fragE#ClassX", "A" * 300),    # two HSPs that merge to full coverage -> REMOVE
    ("twinA#ClassX", "A" * 500),    # equal-length pair: no strictly-longer
    ("twinB#ClassX", "A" * 500),    # container -> BOTH kept (no mutual removal)
]

# self-blast (qseqid sseqid pident length qlen slen qstart qend). blast uses the
# whole FASTA header token as the seqid, so ids carry the '#class' suffix.
HITS = """\
fragA#ClassX\telemA#ClassX\t95\t300\t300\t1000\t1\t300
fragB#ClassY\telemA#ClassX\t95\t300\t300\t1000\t1\t300
fragC#ClassX\telemA#ClassX\t95\t150\t300\t1000\t1\t150
fragD#ClassX\telemA#ClassX\t70\t300\t300\t1000\t1\t300
fragE#ClassX\telemA#ClassX\t95\t150\t300\t1000\t1\t150
fragE#ClassX\telemA#ClassX\t95\t150\t300\t1000\t151\t300
twinA#ClassX\ttwinB#ClassX\t100\t500\t500\t500\t1\t500
twinB#ClassX\ttwinA#ClassX\t100\t500\t500\t500\t1\t500
"""


def main() -> int:
    with tempfile.TemporaryDirectory() as td:
        hits_path = Path(td) / "self_hits.tsv"
        hits_path.write_text(HITS)

        kept = cr.containment_reduce(records, str(hits_path),
                                     min_cov=0.90, min_pid=80.0)

        check("elemA#ClassX" in kept, "longest element must be kept")
        check("fragA#ClassX" not in kept,
              "fragment fully contained in a longer same-class element must be removed")
        check("fragB#ClassY" in kept,
              "cross-class containment must NOT remove (classification preserved)")
        check("fragC#ClassX" in kept,
              "fragment below coverage threshold must be kept")
        check("fragD#ClassX" in kept,
              "fragment below identity threshold must be kept")
        check("fragE#ClassX" not in kept,
              "merged HSP coverage above threshold must remove the fragment")
        check("twinA#ClassX" in kept and "twinB#ClassX" in kept,
              "equal-length sequences have no strictly-longer container -> both kept (no mutual removal)")

        # Every removed sequence must have a kept, longer, same-class container
        # (the losslessness-by-containment invariant).
        removed = {n for n, _ in records} - kept
        check(removed == {"fragA#ClassX", "fragE#ClassX"},
              f"unexpected removed set: {sorted(removed)}")

        # Helper sanity: classification parsing + fallback shape.
        check(cr.classification("x#Class_I/LTR/Ty1_copia/Ale") == "Class_I/LTR/Ty1_copia/Ale",
              "classification() must return the part after '#'")
        check(cr.classification("noclass") == "", "classification() of a header without '#' is ''")

        # A higher coverage/identity bar keeps more (monotonicity).
        kept_strict = cr.containment_reduce(records, str(hits_path),
                                            min_cov=0.99, min_pid=99.0)
        check(kept <= kept_strict,
              "stricter thresholds must keep a superset (never remove more)")

    if failures:
        print("FAIL: test_containment_reduce", file=sys.stderr)
        for f in failures:
            print("  - " + f, file=sys.stderr)
        return 1
    print("OK: test_containment_reduce — all containment checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
