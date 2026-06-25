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
import os
import subprocess
import sys
import tempfile
from collections import defaultdict


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


class _UF:
    """Minimal union-find."""
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


def reduce_trc_group(records, work_dir, min_seq_id, min_cov, threads):
    """Cluster one TRC group rotation-invariantly (monomer-vs-dimer).

    ``records`` is a list of (header, dimer_seq). Returns the representative
    records (longest dimer per cluster). Groups of <=1 are returned unchanged.
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

    uf = _UF(n)
    with open(out_m8) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            # default m8: query target pident ... (qcov not in default cols)
            # we request an explicit format below, so parse accordingly
            q, t, pid, qcov = f[0], f[1], f[2], f[3]
            q, t = int(q), int(t)
            if q == t:
                continue
            pid = float(pid)
            pid = pid * 100 if pid <= 1 else pid
            if pid >= min_seq_id * 100 and float(qcov) >= min_cov:
                uf.union(q, t)

    clusters = defaultdict(list)
    for i in range(n):
        clusters[uf.find(i)].append(i)

    reps = []
    for members in clusters.values():
        best = max(members, key=lambda i: len(records[i][1]))  # longest dimer
        reps.append(records[best])
    return reps


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
    with tempfile.TemporaryDirectory() as tmp_root:
        for trc_id in sorted(groups):
            records = groups[trc_id]
            work_dir = os.path.join(tmp_root, trc_id.replace("/", "_"))
            os.makedirs(work_dir, exist_ok=True)
            try:
                reps = reduce_trc_group(
                    records, work_dir, args.min_seq_id, args.min_cov, args.threads
                )
            except Exception as exc:
                # mmseqs can segfault / error on a pathological group. The
                # reduction is a lossless optimisation, never a correctness
                # requirement, so degrade gracefully instead of failing the
                # whole pipeline. First retry single-threaded (many mmseqs
                # segfaults are thread races / memory pressure with many
                # threads); if that still fails, keep the group UNREDUCED —
                # masking is unaffected, the dimer library for this TRC just
                # stays full-size.
                sys.stderr.write(
                    f"  WARNING: {trc_id}: reduction failed ({exc}); "
                    f"retrying single-threaded\n"
                )
                retry_dir = os.path.join(work_dir, "retry_t1")
                os.makedirs(retry_dir, exist_ok=True)
                try:
                    reps = reduce_trc_group(
                        records, retry_dir, args.min_seq_id, args.min_cov, 1
                    )
                except Exception as exc2:
                    sys.stderr.write(
                        f"  WARNING: {trc_id}: reduction failed again ({exc2}); "
                        f"keeping all {len(records)} sequences unreduced "
                        f"(lossless — masking unaffected)\n"
                    )
                    reps = records
                    fallback_groups += 1
                    fallback_seqs += len(records)
            sys.stderr.write(f"  {trc_id}: {len(records)} → {len(reps)} sequences\n")
            output_records.extend(reps)

    write_fasta(output_records, args.output)
    sys.stderr.write(
        f"Dimer library: {n_input} → {len(output_records)} sequences "
        f"after rotation-invariant per-TRC reduction\n"
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
