#!/usr/bin/env python3
"""What the TSD-confirmed subset is FOR: per-genome ground truth on element geometry.

A 0.7-5% confirmed rate is not a failure -- it is (a) a measure of how recently
LINEs were active in that genome, and (b) the only set of loci where BOTH
boundaries are observed rather than inferred. That makes it a calibration set:
the true per-side extension measured here is what the inference at the other
95% of loci should be reproducing.

Emits, for every confirmed locus: distance from the domain core 3' end to the
END of the poly-A tail (= true 3' extension) and from the core 5' start out to
the TSD (= true 5' extension).
"""
import argparse, bisect, random, sys
sys.path.insert(0, "tmp/line_tprt_2026-09-01")
from tsd_finder import Faidx, rc, scan, null_scores          # noqa: E402
from tsd_at_line_loci import longest_run, load_line          # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--fasta", required=True); ap.add_argument("--gff", required=True)
ap.add_argument("--tail-search", type=int, default=2500); ap.add_argument("--min-run", type=int, default=15)
ap.add_argument("--up-window", type=int, default=2100); ap.add_argument("--floor", type=int, default=10)
ap.add_argument("--decoys", type=int, default=6); ap.add_argument("--limit", type=int, default=400)
ap.add_argument("--seed", type=int, default=0); ap.add_argument("--label", default="")
ap.add_argument("--tsv", default=None)
a = ap.parse_args()

fa = Faidx(a.fasta); rng = random.Random(a.seed)
loci = load_line(a.gff); rng.shuffle(loci); loci = loci[:a.limit]
rows = []
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
    dn = reg[C3 + 1:C3 + 1 + a.tail_search]
    run, at = longest_run(dn, "A")
    if run < a.min_run: continue
    tail_end = C3 + 1 + at + run
    query = reg[tail_end:tail_end + 20]
    w0 = max(0, C5 - a.up_window); win = reg[w0:C5]
    if len(query) < 8 or len(win) < 100: continue
    sc, p, k, ident = scan(query, win, 8, 20)
    nulls = null_scores(query, fa, seqname, max(0, lo), a.decoys, len(win), rng, 8, 20)
    if sc <= a.floor or not (nulls and all(sc > x for x in nulls)): continue
    ext3 = tail_end - (C3 + 1)          # core 3' end -> end of poly-A tail
    ext5 = C5 - (w0 + p)                # TSD (element start) -> core 5' start
    # core coords in GENOME space, so a locus can be matched to a DANTE_LINE.gff3
    # produced under different parameters (the core does not move, only extensions do)
    rows.append((seqname, strand, sc, k, run, at, ext5, ext3,
                 ext5 + ext3 + (C3 - C5 + 1), cs, ce))

print(f"### {a.label}   confirmed loci = {len(rows)}")
if not rows: sys.exit(0)
def q(v, pp):
    v = sorted(v); return v[min(len(v) - 1, int(pp * (len(v) - 1)))]
e5 = [r[6] for r in rows]; e3 = [r[7] for r in rows]
tot = [r[8] for r in rows]; runs = [r[4] for r in rows]; ks = [r[3] for r in rows]
print(f"  {'':22}{'p10':>7}{'median':>8}{'p90':>7}{'max':>7}")
for nm, v in (("TRUE 5' extension", e5), ("TRUE 3' extension", e3),
              ("full element length", tot), ("poly-A run length", runs),
              ("TSD matched length", ks)):
    print(f"  {nm:<22}{q(v,.1):7d}{q(v,.5):8d}{q(v,.9):7d}{max(v):7d}")
print(f"\n  current shipped caps: 5' = 2000 bp, 3' = 800 bp")
over3 = sum(1 for x in e3 if x > 800); over5 = sum(1 for x in e5 if x > 2000)
print(f"  confirmed loci whose TRUE 3' extension EXCEEDS the 800 bp cap: "
      f"{over3}/{len(e3)} ({100*over3/len(e3):.0f}%)")
print(f"  confirmed loci whose TRUE 5' extension EXCEEDS the 2000 bp cap: "
      f"{over5}/{len(e5)} ({100*over5/len(e5):.0f}%)")
if a.tsv:
    with open(a.tsv, "w") as fh:
        fh.write("seq\tstrand\tscore\ttsd_len\tpolya_len\tpolya_offset\text5\text3\telem_len\tcore_start\tcore_end\n")
        for r in sorted(rows): fh.write("\t".join(map(str, r)) + "\n")
