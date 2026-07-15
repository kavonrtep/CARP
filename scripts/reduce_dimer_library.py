#!/usr/bin/env python3
"""
Reduce a TideCluster dimer library, rotation-invariantly, per TRC group.

TideCluster emits each tandem family's consensus as a *dimer* (two monomer
copies), in several rotational phases. Those phase variants are ~the same
sequence but start at different points in the monomer, so a linear aligner
covers only ~one monomer (~50%) when comparing two differently-phased dimers —
which defeats coverage-based clustering of the dimers directly (the old
`mmseqs easy-linclust` on dimers barely reduced such families).

This version exploits the fact that the library entries are *already doubled*:
it aligns each **monomer** (the first half of a dimer) against the **dimers**.
Because the dimer target contains two monomer copies, a monomer at *any* phase
aligns full-length somewhere inside it, so coverage of the monomer is ~1.0 for
same-repeat pairs regardless of rotation. Clustering on those hits collapses all
phase variants of a family, while genuinely distinct families (e.g. different
monomer lengths) still fail the monomer-coverage threshold and stay separate.

Per-TRC, so every TRC keeps >=1 representative and unrelated families never
merge. Representative = the longest dimer in each cluster (most complete
consensus). Validated lossless for masking on tiny_pea and an 800 Mbp genome
(masked bp unchanged within +/-0.15%, while shrinking the library ~9x vs the
old reduction).

Headers follow the RepeatMasker convention ``>name#classification``
(e.g. ``>TRC_1#TRC_1``); the TRC group is the part after ``#``.
"""
import argparse
import multiprocessing
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

# mmseqs nucleotide k-mer length for the monomer-vs-dimer prefilter.
#
# The mmseqs DEFAULT (15) segfaults / aborts when the query (a monomer = first
# half of a dimer) is short relative to k — and the required margin grows with
# k (k=15 dies even on a 21 bp monomer; k=9 is fine at 15 bp; k=7 is fine down
# to 11 bp). A short tandem monomer is exactly the common case here, so we pin
# a small k. A smaller k only makes the prefilter MORE sensitive (more candidate
# pairs); specificity is still enforced by --min-seq-id / coverage and the
# Python identity+coverage check, so clustering quality is unaffected.
MMSEQS_KMER = 7

# Groups whose shortest monomer is below this are skipped (kept unreduced): even
# k=7 crashes at <=10 bp, and such ultra-short satellites are rare in the default
# TideCluster run and carry negligible reduction value. 12 keeps a 1 bp cushion
# over the observed k=7 safe floor (11 bp).
MIN_MONOMER_FOR_MMSEQS = 12

# Two dimers are merged only if the shorter is at least this fraction of the
# longer. Rotational phase variants of a family share the dimer length, so this
# is permissive for the intended merges while it blocks nested/harmonic periods
# (2x, 3x, ... the monomer) from collapsing onto the longest rep — which would
# silently drop masking of the shorter period (see reduce_trc_group). Rotations
# share the dimer length exactly, so this is strict (only ~5% wobble allowed for
# consensus differences); erring high merely keeps more reps (still lossless).
LEN_RATIO_MIN = 0.95


def parse_fasta(path):
    """Yield (header_without_gt, sequence) pairs from a FASTA file."""
    header = None
    parts = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:]
                parts = []
            else:
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def write_fasta(records, path):
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def trc_id_of(header):
    """Return the classification part (after '#') of a FASTA header."""
    idx = header.find("#")
    return header[idx + 1:] if idx != -1 else header


def reduce_trc_group(records, work_dir, min_seq_id, min_cov, threads):
    """Cluster one TRC group rotation-invariantly (monomer-vs-dimer).

    ``records`` is a list of (header, dimer_seq). Returns the representative
    records. Groups of <=1 are returned unchanged.

    Clustering is GREEDY / star, not single-linkage: a sequence is dropped only
    if its monomer aligns (>=min_seq_id identity, >=min_cov query coverage,
    within the length guard) directly to a KEPT representative's dimer, so the
    rep is similar enough to mask every sequence it replaces. (An earlier
    union-find version chained variants transitively, so the single kept rep
    could be too divergent from the far end of the chain.) This keeps the
    reduction lossless for small/medium families; large satellite families can
    still under-mask because RepeatMasker needs the full consensus diversity to
    tile them — a known, accepted limitation documented in CLAUDE.md.
    """
    n = len(records)
    if n <= 1:
        return records

    mono = os.path.join(work_dir, "mono.fasta")   # query: first-half monomers
    dim = os.path.join(work_dir, "dim.fasta")     # target: full dimers
    with open(mono, "w") as fm, open(dim, "w") as fd:
        for i, (_, seq) in enumerate(records):
            fm.write(f">{i}\n{seq[:len(seq) // 2]}\n")
            fd.write(f">{i}\n{seq}\n")

    out_m8 = os.path.join(work_dir, "hits.m8")
    mmseqs_tmp = os.path.join(work_dir, "mmseqs_tmp")
    cmd = [
        "mmseqs", "easy-search", mono, dim, out_m8, mmseqs_tmp,
        "--search-type", "3",            # nucleotide
        "--min-seq-id", str(min_seq_id),
        "--cov-mode", "2",               # coverage of the QUERY (the monomer)
        "-c", str(min_cov),
        "--alignment-mode", "3", "-a",
        "-k", str(MMSEQS_KMER),          # small k: default (15) crashes on short monomers
        "-e", "1e-5",
        "--threads", str(threads),
        "-v", "1",
        "--format-output", "query,target,pident,qcov",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    sys.stderr.write(result.stderr)
    if result.returncode != 0 or not os.path.exists(out_m8):
        tail = "\n".join(result.stderr.strip().splitlines()[-5:])
        raise RuntimeError(
            f"mmseqs easy-search failed (returncode {result.returncode}, "
            f"signal {-result.returncode if result.returncode < 0 else 'n/a'})"
            + (f"; stderr tail:\n{tail}" if tail else "")
        )

    # represents[t] = set of sequences q that target t can represent: q's
    # monomer aligns into t's dimer at >=identity / >=query-coverage AND the two
    # dimers are of near-equal length. Direction matters — t's dimer contains
    # q's monomer, so RepeatMasker masks q's genomic arrays using t.
    represents = defaultdict(set)
    for i in range(n):
        represents[i].add(i)
    with open(out_m8) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            q, t, pid, qcov = f[0], f[1], f[2], f[3]
            q, t = int(q), int(t)
            if q == t:
                continue
            # Length guard: only merge dimers of near-equal length. Rotational
            # phase variants of a family share the dimer length (= 2x the same
            # consensus monomer); a nested/harmonic period (e.g. a 21 bp monomer
            # vs a 105 bp one in the same TRC) would otherwise pass the
            # query-coverage test — the short query aligns fully inside the long
            # dimer — and collapsing it onto the longer rep LOSES masking of the
            # shorter period. Different length => different period => keep separate.
            lq, lt = len(records[q][1]), len(records[t][1])
            if min(lq, lt) < LEN_RATIO_MIN * max(lq, lt):
                continue
            pid = float(pid)
            pid = pid * 100 if pid <= 1 else pid
            if pid >= min_seq_id * 100 and float(qcov) >= min_cov:
                represents[t].add(q)

    # Greedy set cover, longest dimer first (most complete consensus). A rep
    # covers exactly the sequences that align directly to it; anything left
    # uncovered becomes its own rep — so every dropped sequence is guaranteed to
    # have a kept rep that masks it.
    assigned = set()
    reps = []
    for i in sorted(range(n), key=lambda i: (-len(records[i][1]), i)):
        if i in assigned:
            continue
        reps.append(records[i])
        assigned |= represents[i]
    return reps


def _reduce_one_trc(task):
    """Pool worker: reduce ONE TRC group. Returns (trc_id, reps, status, detail).

    ``status`` is "reduced" | "short" | "fallback". mmseqs runs single-threaded
    here because parallelism is across TRC groups; the result is identical to the
    serial multi-threaded path — the reduction is thread-invariant (mmseqs hits
    feed a set-based greedy cover; verified by comparing -t 1 vs -t 4 output).
    Single-threaded is also the more robust path (the code's own retry dropped to
    1 thread precisely because multi-thread mmseqs races/segfaults). The group's
    scratch dir is removed on exit to bound peak disk.
    """
    trc_id, records, work_dir, min_seq_id, min_cov = task
    os.makedirs(work_dir, exist_ok=True)
    try:
        shortest_mono = min(len(s) // 2 for _, s in records) if records else 0
        if len(records) > 1 and shortest_mono < MIN_MONOMER_FOR_MMSEQS:
            return (trc_id, records, "short", shortest_mono)
        try:
            reps = reduce_trc_group(records, work_dir, min_seq_id, min_cov, 1)
            return (trc_id, reps, "reduced", "")
        except Exception as exc:
            # Reduction is a lossless optimisation, never a correctness
            # requirement — keep the group UNREDUCED rather than fail the run.
            return (trc_id, records, "fallback", str(exc))
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser(
        description="Rotation-invariant per-TRC reduction of a TideCluster dimer library."
    )
    parser.add_argument("-i", "--input", required=True, help="Input dimer library FASTA")
    parser.add_argument("-o", "--output", required=True, help="Output reduced library FASTA")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Threads for mmseqs2")
    parser.add_argument("--min_seq_id", type=float, default=0.8,
                        help="Minimum sequence identity (default: 0.8)")
    parser.add_argument("--min_cov", type=float, default=0.8,
                        help="Minimum monomer (query) coverage (default: 0.8)")
    args = parser.parse_args()

    groups = defaultdict(list)
    for header, seq in parse_fasta(args.input):
        groups[trc_id_of(header)].append((header, seq))

    n_input = sum(len(v) for v in groups.values())
    sys.stderr.write(
        f"Reducing dimer library (rotation-invariant): "
        f"{n_input} sequences in {len(groups)} TRC groups\n"
    )

    output_records = []
    fallback_groups = 0
    fallback_seqs = 0
    short_groups = 0
    trc_ids = sorted(groups)
    # Reduce TRC groups in parallel (each mmseqs single-threaded) instead of
    # serially with a multi-threaded mmseqs each. Every group is a tiny,
    # independent mmseqs easy-search, so on a satellite-rich genome the old serial
    # loop paid thousands of ~1-core mmseqs startups back-to-back. Output is
    # identical: the reduction is thread-invariant, pool.map preserves input
    # (sorted-TRC) order, and each group is independent. spawn keeps workers off
    # the parent's heap; each group's records travel as the pickled task.
    with tempfile.TemporaryDirectory() as tmp_root:
        tasks = [
            (trc_id, groups[trc_id],
             os.path.join(tmp_root, trc_id.replace("/", "_")),
             args.min_seq_id, args.min_cov)
            for trc_id in trc_ids
        ]
        ctx = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(max_workers=max(1, args.threads),
                                 mp_context=ctx) as pool:
            results = list(pool.map(_reduce_one_trc, tasks))

    for trc_id, reps, status, detail in results:
        n_rec = len(groups[trc_id])
        if status == "short":
            sys.stderr.write(
                f"  {trc_id}: shortest monomer ~{detail} bp is below the mmseqs "
                f"prefilter floor ({MIN_MONOMER_FOR_MMSEQS} bp); kept unreduced "
                f"({n_rec} sequences)\n"
            )
            short_groups += 1
        elif status == "fallback":
            sys.stderr.write(
                f"  WARNING: {trc_id}: reduction failed ({detail}); keeping all "
                f"{n_rec} sequences unreduced (lossless — masking unaffected)\n"
            )
            fallback_groups += 1
            fallback_seqs += n_rec
        sys.stderr.write(f"  {trc_id}: {n_rec} → {len(reps)} sequences\n")
        output_records.extend(reps)

    write_fasta(output_records, args.output)
    sys.stderr.write(
        f"Dimer library: {n_input} → {len(output_records)} sequences "
        f"after rotation-invariant per-TRC reduction\n"
    )
    if short_groups:
        sys.stderr.write(
            f"NOTE: {short_groups} TRC group(s) had monomers too short for the "
            f"mmseqs prefilter and were kept unreduced (lossless).\n"
        )
    if fallback_groups:
        sys.stderr.write(
            f"NOTE: {fallback_groups} TRC group(s) ({fallback_seqs} sequences) "
            f"could not be reduced (mmseqs failed) and were kept unreduced. "
            f"This is lossless for masking; the only effect is a larger dimer "
            f"library and slightly slower remasking for those groups.\n"
        )


if __name__ == "__main__":
    main()
