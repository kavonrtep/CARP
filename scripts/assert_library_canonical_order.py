#!/usr/bin/env python3
"""Assert the RepeatMasker library is in canonical (sequence-sorted) order.

Why this exists
---------------
RepeatMasker builds its search DB with ``makeblastdb``, which assigns OIDs by
input order, and BLAST tie-breaks equal-scoring HSPs by OID — so a library
presented in a different order masks (very slightly) differently. This shipped
as a real bug: ``reduce_library_size.py`` streams mmseqs cluster representatives
in mmseqs' native, thread-scheduling-dependent order (deliberately, to stay
byte-identical to its R reference), so the containment-reduced library that
RepeatMasker indexes had stable *content* but an unstable *order*, and ~600 bp
of ``Class_I/LINE`` masking flipped run-to-run. The fix canonically sorts that
library in the ``reduce_library_containment`` rule.

This check is the deterministic regression guard for that fix. Unlike the
two-run determinism diff — which only catches the drift if mmseqs happens to
reorder between the two runs (probabilistic on small inputs) — this asserts the
**invariant directly**: the produced library must equal its own canonical sort.
It therefore fails 100% of the time if the sort is removed or a new
library→RepeatMasker path is added without one, on any fixture whose library is
non-trivially ordered (the medium fixture builds a real multi-representative
mmseqs LINE library, so its natural order differs from canonical).

Usage
-----
    assert_library_canonical_order.py <output_dir> [--library REL_PATH]

Exit 0 = canonical (or library absent — treated as skip); exit 1 = not
canonical (regression); exit 2 = usage error.
"""
import argparse
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from canonical_sort_fasta import sort_fasta_by_sequence  # noqa: E402

DEFAULT_REL = "Libraries/combined_library_reduced_containment.fasta"


def is_canonically_sorted(library_path: str) -> bool:
    """Return True iff ``library_path`` equals its own canonical sort byte-for-byte.

    ``sort_fasta_by_sequence`` is idempotent (already-sorted, 60-col-wrapped input
    round-trips to identical bytes), so this is a faithful idempotence test.
    """
    with tempfile.TemporaryDirectory() as tmp:
        resorted = os.path.join(tmp, "resorted.fasta")
        sort_fasta_by_sequence(library_path, resorted, threads=1, tmpdir=tmp)
        with open(library_path, "rb") as a, open(resorted, "rb") as b:
            # Stream-compare to stay bounded on multi-GB libraries.
            while True:
                ca = a.read(1 << 20)
                cb = b.read(1 << 20)
                if ca != cb:
                    return False
                if not ca:
                    return True


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Assert the RepeatMasker library is canonically sequence-sorted.")
    ap.add_argument("output_dir", help="pipeline output directory")
    ap.add_argument("--library", default=DEFAULT_REL,
                    help=f"library path relative to output_dir (default: {DEFAULT_REL})")
    args = ap.parse_args(argv)

    library_path = os.path.join(args.output_dir, args.library)
    if not os.path.exists(library_path):
        sys.stderr.write(
            f"SKIP: RepeatMasker library not found ({library_path}); nothing to check.\n")
        return 0

    if is_canonically_sorted(library_path):
        sys.stderr.write(
            f"OK: {library_path} is in canonical (sequence-sorted) order.\n")
        return 0

    sys.stderr.write(
        f"FAIL: {library_path} is NOT in canonical order.\n"
        "  The RepeatMasker library must be canonically sorted so makeblastdb OID\n"
        "  order — and therefore masking — is reproducible run-to-run. The\n"
        "  reduce_library_containment rule sorts it via canonical_sort_fasta.py;\n"
        "  a failure here means that sort was removed or a new library->RepeatMasker\n"
        "  path bypassed it. See CLAUDE.md 'Determinism (enforced)'.\n")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
