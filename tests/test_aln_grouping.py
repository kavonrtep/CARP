#!/usr/bin/env python3
"""
Unit tests for the bounded, clustering-based grouping that keeps the all-vs-all
flank aligner (scripts/global_local_aln.py) from OOMing on very large families.

Exercises the pure-Python grouping logic directly — no parasail/mmseqs required,
so it runs in the lightweight unit CI (parasail is imported lazily by the
module). Verifies the invariants that make the split safe:

  * below the fire-threshold it is a single group (byte-identical single pass),
  * every group is <= max_group_size (memory bound actually holds),
  * grouping is a deterministic function of the cluster assignment
    (independent of dict / input ordering),
  * no member is lost or duplicated,
  * the mmseqs-unavailable fallback still bounds size and stays deterministic.

Exit 0 = all pass; 1 = one or more failures (printed to stderr).
Usage: python tests/test_aln_grouping.py
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "scripts"))

import global_local_aln as g  # noqa: E402

failures: list[str] = []


def check(cond: bool, msg: str) -> None:
    if not cond:
        failures.append(msg)


def _flat(groups):
    return [m for grp in groups for m in grp]


# --------------------------------------------------------------------------
# _bin_clusters: deterministic binning of an mmseqs cluster assignment.
# One oversized cluster (>= 2x max -> split into ~equal chunks) plus several
# smaller clusters (merged in a fixed order up to max).
# --------------------------------------------------------------------------
MAX = 10
clusters = {
    "repA": [f"a{i}" for i in range(25)],  # 25 >= 2*10 -> split
    "repB": [f"b{i}" for i in range(7)],
    "repC": [f"c{i}" for i in range(6)],
    "repD": [f"d{i}" for i in range(3)],
    "repE": ["e0"],
}
all_members = _flat(clusters.values())

binned = g._bin_clusters(clusters, MAX)
# Determinism under a different dict iteration order.
binned_rev = g._bin_clusters(dict(reversed(list(clusters.items()))), MAX)

check(binned == binned_rev,
      "_bin_clusters is not deterministic across dict order")
check(all(len(grp) <= MAX for grp in binned),
      f"_bin_clusters produced a group exceeding max_group_size={MAX}: "
      f"{[len(x) for x in binned]}")
check(sorted(_flat(binned)) == sorted(all_members),
      "_bin_clusters lost or duplicated members")
check(len(_flat(binned)) == len(set(_flat(binned))),
      "_bin_clusters produced duplicate members")
# The 25-member cluster must have been split (no single group holds 25).
check(max(len(x) for x in binned) <= MAX,
      "_bin_clusters failed to split the oversized cluster")


# --------------------------------------------------------------------------
# group_sequences_for_alignment: below-threshold no-op (the "limit where it
# fires"). Must return exactly one group == the input, so the caller takes the
# unchanged single-pass path. No fasta is read on this branch.
# --------------------------------------------------------------------------
ids = [f"s{i}" for i in range(50)]
check(g.group_sequences_for_alignment("does-not-exist.fasta", ids, 1000) == [ids],
      "below-threshold input must be a single unchanged group")
check(g.group_sequences_for_alignment("does-not-exist.fasta", ids, None) == [ids],
      "max_group_size=None must be a single unchanged group")
check(g.group_sequences_for_alignment("does-not-exist.fasta", ids, 0) == [ids],
      "falsy max_group_size must be a single unchanged group")
# Exactly at the threshold is still one group (split only when strictly larger).
ids_at = [f"s{i}" for i in range(10)]
check(g.group_sequences_for_alignment("does-not-exist.fasta", ids_at, 10) == [ids_at],
      "N == max_group_size must not split")


# --------------------------------------------------------------------------
# group_sequences_for_alignment: mmseqs-unavailable fallback. Pointing at a
# non-existent fasta forces mmseqs to fail; the deterministic id-sorted chunk
# fallback must still bound size, stay lossless, and be order-independent.
# --------------------------------------------------------------------------
big = [f"x{i:04d}" for i in range(2500)]
chunks = g.group_sequences_for_alignment("/no/such/file.fasta", big, 1000, verbose=False)
check(all(len(c) <= 1000 for c in chunks),
      f"fallback chunks exceed max_group_size: {[len(c) for c in chunks]}")
check(sorted(_flat(chunks)) == sorted(big),
      "fallback lost or duplicated members")
chunks_rev = g.group_sequences_for_alignment(
    "/no/such/file.fasta", list(reversed(big)), 1000, verbose=False)
check(chunks == chunks_rev,
      "fallback chunking is not deterministic across input order")


if failures:
    print(f"FAILED ({len(failures)}):", file=sys.stderr)
    for f in failures:
        print(f"  - {f}", file=sys.stderr)
    sys.exit(1)
print("test_aln_grouping.py: all checks passed")
