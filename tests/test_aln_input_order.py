#!/usr/bin/env python3
"""
The all-vs-all flank aligner (scripts/global_local_aln.py) must give the same
result whatever order its input FASTA is in (GitHub issue #7).

The callers (dante_tir_fallback, dante_line) build that FASTA from
`seqkit subseq --bed`, which writes records in a different chromosome order on
every call. mmseqs easy-cluster (grouping) and easy-search (prefilter, capped
by --max-seqs) both depend on input order, so before the fix the candidate
pairs -- and the inferred element boundaries -- changed run to run.

Part 1 (always runs; no mmseqs/parasail needed, so it fits the unit CI):
  * mmseqs is stubbed by ORDER-SENSITIVE fakes (they keep the first records of
    the file they are given). The output must still be identical for shuffled
    inputs, which holds only if the engine hands mmseqs a canonical order.
  * the prefilter passes --max-seqs >= the number of input sequences.
  * process_strand_orientation_with_seqkit writes records in id order.

Part 2 (runs only when mmseqs is on PATH; skipped otherwise):
  * real mmseqs grouping + prefilter (parasail stubbed, so the output is the
    candidate-pair set) on a synthetic family set big enough to hit both the
    grouping (max_group_size) and the old 300-hit prefilter cap, in two
    shuffled orders -> the TSVs must be byte-identical.

Exit 0 = all pass; 1 = one or more failures (printed to stderr).
Usage: python tests/test_aln_input_order.py
"""
from __future__ import annotations

import os
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "scripts"))

import global_local_aln as g  # noqa: E402

failures: list[str] = []


def check(cond: bool, msg: str) -> None:
    if not cond:
        failures.append(msg)


def write_fasta(path, records):
    with open(path, "w") as fh:
        for sid, seq in records:
            fh.write(f">{sid}\n{seq}\n")


def ids_in(path):
    with open(path) as fh:
        return [line[1:].strip() for line in fh if line.startswith(">")]


def mutate(seq, rate, rng):
    return "".join(rng.choice("ACGT") if rng.random() < rate else b for b in seq)


def make_records(rng):
    """Families of related flanks plus unrelated singletons. The big family has
    700 members, so a query in it has ~700 prefilter hits -- well over the old
    300 cap (with fewer, symmetrising the pairs hides the cap's effect)."""
    recs = []
    for fam, n in (("big", 700), ("mid", 120), ("small", 40)):
        core = "".join(rng.choice("ACGT") for _ in range(150))
        for i in range(n):
            recs.append((f"{fam}_{i:04d}", mutate(core, 0.03, rng)))
    for i in range(160):
        recs.append((f"single_{i:04d}", "".join(rng.choice("ACGT") for _ in range(150))))
    return recs


tmp = Path(tempfile.mkdtemp(prefix="test_aln_input_order_"))
try:
    rng = random.Random(7)
    records = make_records(rng)
    shuffled_a = records[:]
    random.Random(1).shuffle(shuffled_a)
    shuffled_b = records[:]
    random.Random(2).shuffle(shuffled_b)
    fa_a, fa_b = tmp / "a.fasta", tmp / "b.fasta"
    write_fasta(fa_a, shuffled_a)
    write_fasta(fa_b, shuffled_b)

    # ----------------------------------------------------------------------
    # Part 1a: engine with order-sensitive mmseqs stubs.
    # ----------------------------------------------------------------------
    real_cluster = g._run_mmseqs_easy_cluster
    real_prefilter = g.run_mmseqs_prefilter
    real_compare = g.run_alignment_comparison

    def fake_cluster(fasta_file, threads=1, verbose=True):
        # Greedy "clustering" by file order: blocks of 50 consecutive records,
        # representative = first in block. Order-dependent, like mmseqs.
        ids = ids_in(fasta_file)
        return {ids[i]: ids[i:i + 50] for i in range(0, len(ids), 50)}

    def fake_prefilter(fasta_file, end, min_identity=0.8, threads=1,
                       verbose=True, kmer=0):
        # Keep each query's first 5 targets in file order (a tiny "cap").
        ids = ids_in(fasta_file)
        pairs = set()
        for q in ids:
            for t in [x for x in ids if x != q][:5]:
                pairs.add(tuple(sorted((q, t))))
        return pairs

    def fake_compare(id1, s1, id2, s2, args):
        return {"query_id": id1, "ref_id": id2, "max_score": 100}

    g._run_mmseqs_easy_cluster = fake_cluster
    g.run_mmseqs_prefilter = fake_prefilter
    g.run_alignment_comparison = fake_compare
    try:
        outs = []
        for fa in (fa_a, fa_b):
            out = tmp / f"{fa.stem}_stub.tsv"
            g.run_all_vs_all_alignment(str(fa), str(out), verbose=False,
                                       max_group_size=200)
            outs.append(out.read_bytes())
        check(outs[0] == outs[1],
              "stubbed engine output depends on input order (grouped path)")
        outs = []
        for fa in (fa_a, fa_b):
            out = tmp / f"{fa.stem}_stub1.tsv"
            g.run_all_vs_all_alignment(str(fa), str(out), verbose=False,
                                       max_group_size=None)
            outs.append(out.read_bytes())
        check(outs[0] == outs[1],
              "stubbed engine output depends on input order (single-pass path)")
    finally:
        g._run_mmseqs_easy_cluster = real_cluster
        g.run_mmseqs_prefilter = real_prefilter
        g.run_alignment_comparison = real_compare

    # ----------------------------------------------------------------------
    # Part 1b: prefilter lifts mmseqs' --max-seqs 300 cap to the input size.
    # ----------------------------------------------------------------------
    seen_cmds = []
    real_run = subprocess.run

    def fake_run(cmd, *a, **kw):
        seen_cmds.append(cmd)
        Path(cmd[4]).write_text("")  # empty search_results.tsv
        return subprocess.CompletedProcess(cmd, 0, "", "")

    subprocess.run = fake_run
    try:
        g.run_mmseqs_prefilter(str(fa_a), "5", verbose=False)
    finally:
        subprocess.run = real_run
    check(len(seen_cmds) == 1, f"expected one mmseqs call, got {len(seen_cmds)}")
    if seen_cmds:
        cmd = seen_cmds[0]
        if "--max-seqs" not in cmd:
            failures.append("prefilter does not pass --max-seqs (mmseqs default 300 caps hits)")
        else:
            val = int(cmd[cmd.index("--max-seqs") + 1])
            check(val >= len(records),
                  f"--max-seqs {val} < number of input sequences {len(records)}")

    # ----------------------------------------------------------------------
    # Part 1c: the flank FASTA writer emits records in id order.
    # ----------------------------------------------------------------------
    import dante_line  # noqa: E402

    class P:  # minimal FeatureGroup stand-in (only group_id/strand are read)
        def __init__(self, gid, strand):
            self.group_id, self.strand = gid, strand

    plus_ids = [f"g{i:03d}" for i in range(20)]
    raw = tmp / "raw.fasta"
    order = plus_ids[:]
    random.Random(3).shuffle(order)
    # seqkit-style headers: "<id>_5prime ..."; all plus strand so no revcom call
    write_fasta(raw, [(f"{gid}_5prime", "ACGT" * 10) for gid in order])
    oriented = tmp / "oriented.fasta"
    dante_line.process_strand_orientation_with_seqkit(
        str(raw), str(oriented), [P(gid, "+") for gid in plus_ids])
    check(ids_in(oriented) == sorted(plus_ids),
          f"process_strand_orientation_with_seqkit output not id-sorted: {ids_in(oriented)[:5]}")

    # ----------------------------------------------------------------------
    # Part 2: REAL mmseqs (grouping + prefilter), parasail stubbed so the
    # output is exactly the candidate-pair set (and the test stays fast).
    # Shuffled inputs -> identical TSV. Fails on the pre-fix code.
    # ----------------------------------------------------------------------
    if shutil.which("mmseqs") is None:
        print("test_aln_input_order.py: mmseqs not available -- "
              "skipping the real-mmseqs check (Part 2)")
    else:
        g.run_alignment_comparison = fake_compare
        try:
            for mgs in (None, 400):
                outs = []
                for fa in (fa_a, fa_b):
                    out = tmp / f"{fa.stem}_real_{mgs}.tsv"
                    g.run_all_vs_all_alignment(str(fa), str(out), end="5",
                                               threads=4, verbose=False,
                                               max_group_size=mgs)
                    outs.append(out.read_bytes() if out.exists() else b"")
                check(outs[0] != b"",
                      f"real-mmseqs run produced no pairs (max_group_size={mgs})")
                check(outs[0] == outs[1],
                      f"real-mmseqs candidate pairs depend on input order "
                      f"(max_group_size={mgs})")
        finally:
            g.run_alignment_comparison = real_compare
finally:
    shutil.rmtree(tmp, ignore_errors=True)

if failures:
    print(f"FAILED ({len(failures)}):", file=sys.stderr)
    for f in failures:
        print(f"  - {f}", file=sys.stderr)
    sys.exit(1)
print("test_aln_input_order.py: all checks passed")
