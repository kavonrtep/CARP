#!/usr/bin/env python3
"""Joint false-positive rate of poly-A AND TSD, measured rather than multiplied.

The 3' side of a LINE genuinely carries a poly-A tail. The 5' side does not, so
any poly-A run found there is spurious by construction. Running the IDENTICAL
chain on both sides therefore measures:

  3' side = real rate (signal + noise)
  5' side = NOISE ALONE, for both tests combined

Chain, per side: find a poly-A homopolymer -> take the 20 bp beyond it as the
candidate TSD -> search the other flank for a matching direct repeat -> require
it to clear the measured chance floor AND beat every decoy window.
"""
import argparse, bisect, random, sys
from collections import defaultdict
sys.path.insert(0, "tmp/line_tprt_2026-09-01")
from tsd_finder import Faidx, rc, scan, null_scores          # noqa: E402
from tsd_at_line_loci import longest_run, load_line          # noqa: E402

_C = str.maketrans("ACGTN", "TGCAN")
def comp(s): return s.translate(_C)

ap = argparse.ArgumentParser()
ap.add_argument("--fasta", required=True); ap.add_argument("--gff", required=True)
ap.add_argument("--tail-search", type=int, default=2500)
ap.add_argument("--min-run", type=int, default=15)
ap.add_argument("--window", type=int, default=2100)
ap.add_argument("--floor", type=int, default=10)
ap.add_argument("--decoys", type=int, default=6)
ap.add_argument("--limit", type=int, default=400); ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--label", default="")
a = ap.parse_args()

fa = Faidx(a.fasta); rng = random.Random(a.seed)
loci = load_line(a.gff); rng.shuffle(loci); loci = loci[:a.limit]

res = {"3prime": dict(n=0, tail=0, floor=0, beat=0),
       "5prime": dict(n=0, tail=0, floor=0, beat=0)}

for seqname, strand, cs, ce in loci:
    pad = a.window + a.tail_search + 200
    lo, hi = cs - 1 - pad, ce + pad
    reg = fa.fetch(seqname, max(0, lo), hi)
    if len(reg) < (hi - lo) // 2 or reg.count("N") > len(reg) // 10:
        continue
    C5o, C3o = (cs - 1) - max(0, lo), (ce - 1) - max(0, lo)
    if strand == "-":
        L = len(reg); reg = rc(reg); C5, C3 = L - 1 - C3o, L - 1 - C5o
    else:
        C5, C3 = C5o, C3o
    if not (0 <= C5 < C3 < len(reg)):
        continue

    for side in ("3prime", "5prime"):
        r = res[side]
        r["n"] += 1
        if side == "3prime":
            out = reg[C3 + 1:C3 + 1 + a.tail_search]          # downstream, real tail
            def anchor_end(at, run): return C3 + 1 + at + run
            win_lo, win_hi = max(0, C5 - a.window), C5        # search upstream
        else:
            # mirror image: walk UPSTREAM from the core, reading outward, so a
            # poly-A here is the same test applied where no tail should exist
            out = reg[max(0, C5 - a.tail_search):C5][::-1]
            def anchor_end(at, run): return C5 - (at + run)
            win_lo, win_hi = C3 + 1, min(len(reg), C3 + 1 + a.window)
        run, at = longest_run(out, "A")
        if run < a.min_run:
            continue
        r["tail"] += 1
        e = anchor_end(at, run)
        query = reg[e:e + 20] if side == "3prime" else reg[max(0, e - 20):e][::-1]
        win = reg[win_lo:win_hi]
        if len(query) < 8 or len(win) < 100:
            continue
        sc, p, k, ident = scan(query, win, 8, 20)
        if sc > a.floor:
            r["floor"] += 1
        nulls = null_scores(query, fa, seqname, max(0, lo), a.decoys, len(win), rng, 8, 20)
        if sc > a.floor and nulls and all(sc > x for x in nulls):
            r["beat"] += 1

print(f"### {a.label}   loci={res['3prime']['n']}")
hdr_a = "3' (real)"; hdr_b = "5' (decoy)"
print(f"  {'':<34}{hdr_a:>12}{hdr_b:>12}{'ratio':>9}")
for lab, k in (("poly-A run found", "tail"),
               ("+ TSD above the chance floor", "floor"),
               ("+ TSD beats every decoy window", "beat")):
    t, d = res["3prime"][k], res["5prime"][k]
    n = res["3prime"]["n"]
    ratio = f"{t/d:.1f}x" if d else ("inf" if t else "-")
    print(f"  {lab:<34}{t:6d} ({100*t/n:4.1f}%){d:6d} ({100*d/n:4.1f}%){ratio:>9}")
