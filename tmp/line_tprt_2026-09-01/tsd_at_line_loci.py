#!/usr/bin/env python3
"""The multiplier that decides Phase 2: how long are recoverable TSDs at REAL
LINE loci, and how many clear the chance floor?

Chain (Phase 1 + Phase 2 machinery, end to end):
  1. anchor the 3' end on a poly-A homopolymer (>=15 A, <=1 mismatch) within
     2500 bp downstream of the domain core;
  2. take the 20 bp immediately 3' of the tail as the candidate TSD;
  3. search UPSTREAM OF THE CORE, bounded by the LINE 5' extension cap (~2 kb),
     for the matching direct repeat -- the bound is what keeps the chance floor
     at 11 bp instead of 12+;
  4. compare each hit to the per-locus decoy null.

A perfect TSD of length L scores exactly L, so the reported score IS the
effective matched length. What matters is the fraction clearing the floor.
"""
import argparse, bisect, random, sys
from collections import defaultdict
sys.path.insert(0, "tmp/line_tprt_2026-09-01")
from tsd_finder import Faidx, rc, scan, null_scores   # noqa: E402

_C = str.maketrans("ACGTN", "TGCAN")
def comp(s): return s.translate(_C)

def longest_run(seq, ch="A", max_mm=1):
    best = 0; at = -1; n = len(seq); i = 0
    while i < n:
        if seq[i] != ch: i += 1; continue
        mm = 0; j = i; last = i
        while j < n:
            if seq[j] == ch: last = j
            else:
                mm += 1
                if mm > max_mm: break
            j += 1
        if last - i + 1 > best: best, at = last - i + 1, i
        i = last + 1 if last > i else i + 1
    return best, at

def load_line(gff):
    els = []; by_seq = defaultdict(list); doms = defaultdict(list)
    for line in open(gff):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9: continue
        if f[2] == "LINE_element":
            els.append((f[0], int(f[3]), int(f[4]), f[6]))
            by_seq[f[0]].append((int(f[3]), int(f[4]), len(els) - 1))
        elif f[2] == "protein_domain":
            doms[f[0]].append((int(f[3]), int(f[4])))
    core = {}
    for seq, lst in by_seq.items():
        lst.sort(); st = [x[0] for x in lst]
        for ds, de in doms.get(seq, []):
            i = bisect.bisect_right(st, ds) - 1
            if i >= 0 and lst[i][0] <= ds and de <= lst[i][1]:
                k = lst[i][2]; c = core.get(k)
                core[k] = (min(c[0], ds), max(c[1], de)) if c else (ds, de)
    return [(els[k][0], els[k][3], core[k][0], core[k][1]) for k in sorted(core)]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True); ap.add_argument("--gff", required=True)
    ap.add_argument("--tail-search", type=int, default=2500)
    ap.add_argument("--min-run", type=int, default=15)
    ap.add_argument("--up-window", type=int, default=2100)
    ap.add_argument("--floor", type=int, default=10, help="null p99 for this window")
    ap.add_argument("--decoys", type=int, default=6)
    ap.add_argument("--limit", type=int, default=400); ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--label", default="")
    a = ap.parse_args()

    fa = Faidx(a.fasta); rng = random.Random(a.seed)
    loci = load_line(a.gff); rng.shuffle(loci); loci = loci[:a.limit]
    n = tails = 0; scores = []; beat_all = 0; above_floor = 0
    for seqname, strand, cs, ce in loci:
        pad5 = a.up_window + 60; pad3 = a.tail_search + 80
        lo = (cs - 1 - pad5) if strand == "+" else (cs - 1 - pad3)
        hi = (ce + pad3) if strand == "+" else (ce + pad5)
        reg = fa.fetch(seqname, max(0, lo), hi)
        if len(reg) < (hi - lo) // 2 or reg.count("N") > len(reg) // 10: continue
        C5o = (cs - 1) - max(0, lo); C3o = (ce - 1) - max(0, lo)
        if strand == "-":
            L = len(reg); reg = rc(reg); C5, C3 = L - 1 - C3o, L - 1 - C5o
        else:
            C5, C3 = C5o, C3o
        if not (0 <= C5 < C3 < len(reg)): continue
        n += 1
        dn = reg[C3 + 1:C3 + 1 + a.tail_search]
        run, at = longest_run(dn, "A")
        if run < a.min_run: continue
        tails += 1
        tail_end = C3 + 1 + at + run
        query = reg[tail_end:tail_end + 20]
        w0 = max(0, C5 - a.up_window)
        win = reg[w0:C5]
        if len(query) < 8 or len(win) < 100: continue
        sc, p, k, ident = scan(query, win, 8, 20)
        scores.append(sc)
        if sc > a.floor: above_floor += 1
        nulls = null_scores(query, fa, seqname, max(0, lo), a.decoys,
                            len(win), rng, 8, 20)
        if nulls and all(sc > x for x in nulls): beat_all += 1
    print(f"### {a.label}   loci={n}")
    if not n: return
    print(f"  poly-A tail (>={a.min_run} A within {a.tail_search} bp): {tails} ({100*tails/n:.1f}% of loci)")
    if not scores:
        print("  no TSD searches performed"); return
    v = sorted(scores); q = lambda pp: v[min(len(v)-1, int(pp*(len(v)-1)))]
    print(f"  TSD match length at those loci (score == matched length when perfect):")
    print(f"     median {q(.5)}   p75 {q(.75)}   p90 {q(.9)}   max {v[-1]}")
    print(f"  clearing the {a.up_window} bp chance floor (score > {a.floor}): "
          f"{above_floor}/{len(scores)} ({100*above_floor/len(scores):.1f}% of tailed loci)")
    print(f"  beating every decoy: {beat_all}/{len(scores)} ({100*beat_all/len(scores):.1f}%)")
    print(f"  >>> END-TO-END YIELD (tail AND TSD above floor): "
          f"{100*above_floor/n:.1f}% of all LINE loci")

if __name__ == "__main__":
    main()
