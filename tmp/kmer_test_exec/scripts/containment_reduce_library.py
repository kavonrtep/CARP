#!/usr/bin/env python3
"""Second-round *containment* reduction of a CARP repeat library.

The per-class CAP3/mmseqs reduction (`reduce_library_size.py`) collapses
near-duplicate elements within a class but leaves many short fragments that are
fully *contained* in a longer element of the same class. RepeatMasker would mask
those fragments' genomic copies via the longer element anyway, so carrying them
only slows the (dominant) RepeatMasker step.

This pass removes each sequence that is >= --min-coverage covered, at
>= --min-identity percent identity, over its own length by a RETAINED, strictly
longer sequence of the SAME classification. Greedy longest-first, so every
removed sequence is contained in a kept one of the same class:

  * masking is preserved  — the container masks the fragment's genomic copies;
  * classification is preserved — container and fragment share the #class.

Validated masked-bp-lossless with RepeatMasker on the Pisum pangenome at
--min-identity 80 --min-coverage 0.90 (about -22 percent library bp, about
-30 percent RepeatMasker wall-time). Output feeds RepeatMasker.

Self-comparison uses blastn (the RepeatMasker-relevant aligner). Lossless
fallback: if blastn is unavailable or fails, the input is copied through
unchanged (the reduction is an optimisation, never a correctness requirement).
"""
import argparse
import os
import subprocess
import sys
from collections import defaultdict


def iter_fasta(path):
    name = None
    parts = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name = line[1:].rstrip("\n").split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, "".join(parts)


def classification(name):
    return name.split("#", 1)[1] if "#" in name else ""


def run_blast(library, workdir, threads):
    """makeblastdb + self blastn; return path to the outfmt-6 hits, or None."""
    os.makedirs(workdir, exist_ok=True)
    db = os.path.join(workdir, "libdb")
    hits = os.path.join(workdir, "self_hits.tsv")
    try:
        subprocess.run(["makeblastdb", "-in", library, "-dbtype", "nucl",
                        "-out", db], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        with open(hits, "w") as out:
            subprocess.run(
                ["blastn", "-query", library, "-db", db,
                 "-outfmt", "6 qseqid sseqid pident length qlen slen qstart qend",
                 "-evalue", "1e-10", "-num_threads", str(threads),
                 "-max_target_seqs", "100"],
                check=True, stdout=out, stderr=subprocess.DEVNULL)
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        sys.stderr.write(f"WARNING: blastn self-comparison failed ({exc}); "
                         f"copying library through unreduced.\n")
        return None
    return hits


def containment_reduce(records, hits, min_cov, min_pid):
    """Return the kept names. Greedy longest-first; a sequence is dropped only
    if a kept, strictly-longer, same-class sequence covers >=min_cov of it at
    >=min_pid identity."""
    length = {n: len(s) for n, s in records}
    cls = {n: classification(n) for n, _ in records}

    # HSP query-intervals per (query, subject), identity-filtered
    ivs = defaultdict(list)
    with open(hits) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            q, s, pid, _ln, _ql, _sl, qs, qe = f
            if q == s or float(pid) < min_pid:
                continue
            qs, qe = int(qs), int(qe)
            if qe < qs:
                qs, qe = qe, qs
            ivs[(q, s)].append((qs, qe))

    def covered(q, s):
        a = sorted(ivs.get((q, s), ()))
        if not a:
            return 0.0
        tot = cs = ce = 0
        cs = ce = None
        for x, y in a:
            if cs is None:
                cs, ce = x, y
            elif x <= ce + 1:
                ce = max(ce, y)
            else:
                tot += ce - cs + 1
                cs, ce = x, y
        tot += ce - cs + 1
        return tot / length[q]

    # candidate containers per query: strictly longer, same class, cov>=min_cov
    containers = defaultdict(list)
    for (q, s) in ivs:
        if length[s] > length[q] and cls[s] == cls[q] and covered(q, s) >= min_cov:
            containers[q].append(s)

    kept = set()
    for n in sorted(length, key=lambda x: (-length[x], x)):
        if any(c in kept for c in containers.get(n, ())):
            continue
        kept.add(n)
    return kept


def write_fasta(records, keep, path):
    with open(path, "w") as fo:
        for name, seq in records:
            if name in keep:
                fo.write(f">{name}\n")
                for i in range(0, len(seq), 60):
                    fo.write(seq[i:i + 60] + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-i", "--input", required=True, help="Input library FASTA")
    p.add_argument("-o", "--output", required=True, help="Output reduced FASTA")
    p.add_argument("-d", "--workdir", default="containment_workdir",
                   help="Scratch dir for blast db/hits")
    p.add_argument("-t", "--threads", type=int, default=1)
    p.add_argument("--min-identity", type=float, default=80.0,
                   help="Min %% identity of containment HSPs [default %(default)s]")
    p.add_argument("--min-coverage", type=float, default=0.90,
                   help="Min fraction of the shorter sequence covered [default %(default)s]")
    args = p.parse_args()

    records = list(iter_fasta(args.input))
    n0 = len(records)
    hits = run_blast(args.input, args.workdir, args.threads)
    if hits is None:
        keep = {n for n, _ in records}            # fallback: keep everything
    else:
        keep = containment_reduce(records, hits, args.min_coverage, args.min_identity)
    write_fasta(records, keep, args.output)

    bp0 = sum(len(s) for _, s in records)
    bp1 = sum(len(s) for n, s in records if n in keep)
    sys.stderr.write(
        f"containment reduction (COV>={args.min_coverage}, PID>={args.min_identity}): "
        f"{n0} -> {len(keep)} seqs ({100*(n0-len(keep))/max(1,n0):.1f}% fewer), "
        f"{bp0} -> {bp1} bp ({100*(bp0-bp1)/max(1,bp0):.1f}% fewer)\n")


if __name__ == "__main__":
    main()
