#!/usr/bin/env python3
"""Is the mmseqs prefilter, not the support rule, what stops LINE extension?

CONTEXT. 71% of LINE elements get zero extension. The suspected cause was the
support rule (min_group_alignments=5, support_fraction=0.5). It is not: no group
that clears the size gate ends up with a k-th largest of 0. What kills the layer
is that most groups never reach 5 partners, because a pair is only ever aligned
if mmseqs easy-search matches the 30 nt ADJACENT TO THE DOMAIN CORE at >=80%
identity (global_local_aln.run_mmseqs_prefilter, prefilter_identity=0.8).

That 0.8 is HARDCODED in run_all_vs_all_alignment and exposed nowhere -- no CLI
flag, no config key -- yet it gates the whole extension layer.

THIS PROBE asks whether 0.8 is defensible, using families defined from the CORE
(the REXdb protein each domain hit) so the test is independent of the flanks it
is judging. For every same-family and different-family pair it measures the
identity of exactly the 30 nt the prefilter compares.

Usage:
  prefilter_identity_probe.py --gff <DANTE_LINE.gff3> --flank <ENDO_RT_5prime.fasta>
"""
import argparse, bisect, itertools, random, re
from collections import defaultdict


def load_families(gff):
    """element id -> tuple of 'DOMAIN:REXdbID', derived from the core only."""
    els, by_seq, dom = {}, defaultdict(list), defaultdict(list)
    for line in open(gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        d = dict(kv.split("=", 1) for kv in f[8].rstrip(";").split(";") if "=" in kv)
        if f[2] == "LINE_element":
            els[d["ID"]] = f[0]
            by_seq[f[0]].append((int(f[3]), int(f[4]), d["ID"]))
        elif f[2] == "protein_domain":
            m = re.search(r"REXdb_ID(\d+)", d.get("Best_Hit", ""))
            dom[f[0]].append((int(f[3]), int(f[4]), d.get("Name", "?"),
                              m.group(1) if m else "?"))
    fam = {}
    for seq, lst in by_seq.items():
        lst.sort()
        starts = [x[0] for x in lst]
        for ds, de, nm, rx in dom.get(seq, []):
            i = bisect.bisect_right(starts, ds) - 1
            if i >= 0 and lst[i][0] <= ds and de <= lst[i][1]:
                fam.setdefault(lst[i][2], set()).add(f"{nm}:{rx}")
    return {k: tuple(sorted(v)) for k, v in fam.items()}


def read_fasta(path):
    out, cur, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if cur:
                out[cur] = "".join(buf)
            cur = line[1:].strip().split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if cur:
        out[cur] = "".join(buf)
    return out


def identity(a, b):
    return sum(1 for x, y in zip(a, b) if x == y) / len(a)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--flank", required=True)
    ap.add_argument("--region", type=int, default=30,
                    help="nt compared by the prefilter (default 30, its value)")
    ap.add_argument("--core-at", choices=("end", "start"), default="end",
                    help="which end of the flank abuts the core; verified "
                         "empirically as 'end' for the 5prime files")
    ap.add_argument("--decoys", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    fam = load_families(a.gff)
    seqs = read_fasta(a.flank)
    groups = defaultdict(list)
    for eid, f in fam.items():
        if eid in seqs and len(seqs[eid]) >= a.region:
            groups[f].append(eid)
    groups = {k: v for k, v in groups.items() if len(v) >= 2}

    def region(eid):
        s = seqs[eid].upper()
        return s[-a.region:] if a.core_at == "end" else s[:a.region]

    same = []
    for k in sorted(groups):
        for x, y in itertools.combinations(sorted(groups[k]), 2):
            same.append(identity(region(x), region(y)))

    rng = random.Random(a.seed)
    flat = [(k, e) for k in sorted(groups) for e in groups[k]]
    diff = []
    for _ in range(a.decoys):
        (k1, x), (k2, y) = rng.sample(flat, 2)
        if k1 != k2:
            diff.append(identity(region(x), region(y)))

    if not same or not diff:
        print("not enough pairs")
        return

    def q(v, p):
        v = sorted(v)
        return v[min(len(v) - 1, int(p * (len(v) - 1)))]

    print(f"Identity of the {a.region} nt the prefilter compares "
          f"(families from the core's REXdb hit)\n")
    print(f"  SAME family  n={len(same):5d}  median {q(same,.5):.2f}  "
          f"p75 {q(same,.75):.2f}  p90 {q(same,.9):.2f}")
    print(f"  DIFF family  n={len(diff):5d}  median {q(diff,.5):.2f}  "
          f"p75 {q(diff,.75):.2f}  p90 {q(diff,.9):.2f}")
    print()
    for t in (0.6, 0.7, 0.8, 0.9):
        s = 100 * sum(1 for x in same if x >= t) / len(same)
        d = 100 * sum(1 for x in diff if x >= t) / len(diff)
        mark = "  <- shipped" if abs(t - 0.8) < 1e-9 else ""
        print(f"  threshold {t:.1f}: admits {s:5.1f}% of SAME-family pairs, "
              f"{d:5.1f}% of DIFF-family pairs{mark}")


main()
