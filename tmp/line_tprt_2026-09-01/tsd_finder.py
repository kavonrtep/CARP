#!/usr/bin/env python3
"""Phase 0 + 0b: a direct-repeat (TSD) finder with a calibrated null.

WHY THE NULL COMES FIRST. An 8-20 bp direct repeat searched over a multi-kb
window WILL be found by chance, most readily in AT-rich repeat-dense sequence --
exactly where the LINE boundary problem lives. Without a per-locus null,
"TSD-confirmed" is a confident label on noise.

NULL MODEL (decoy windows). The same query is searched in windows drawn from
distant loci on the same sequence. This preserves both the query and realistic
genomic composition/repeat structure -- unlike a shuffle, which destroys the
local structure that makes chance hits likely in the first place, and would
therefore understate the false-positive rate.

VALIDATION (0b) uses ground truth we already have: DANTE_LTR reports TSD= and
DANTE_TIR reports tsd=, so the true insertion start is known. We hide it, anchor
at the element 3' end, search upstream, and ask whether the finder recovers it.
"""
import argparse, bisect, random, sys
from collections import defaultdict
import numpy as np

# ---------------------------------------------------------------- sequence io
class Faidx:
    def __init__(self, fa):
        self.f = open(fa, "rb"); self.idx = {}
        for line in open(fa + ".fai"):
            n, ln, off, lb, lw = line.split("\t")[:5]
            self.idx[n] = (int(ln), int(off), int(lb), int(lw))
    def length(self, name): return self.idx[name][0] if name in self.idx else 0
    def fetch(self, name, a, b):
        if name not in self.idx: return ""
        ln, off, lb, lw = self.idx[name]; a = max(0, a); b = min(ln, b)
        if b <= a: return ""
        s = off + a // lb * lw + a % lb
        e = off + (b - 1) // lb * lw + (b - 1) % lb + 1
        self.f.seek(s)
        return self.f.read(e - s).decode().replace("\n", "").replace("\r", "").upper()

_C = str.maketrans("ACGTN", "TGCAN")
def rc(s): return s.translate(_C)[::-1]

def enc(s):
    a = np.frombuffer(s.encode(), dtype=np.uint8)
    return a

# ---------------------------------------------------------------- the finder
def scan(query, window, min_k=8, max_k=20, mm_penalty=2):
    """Best direct-repeat match of `query`'s prefix anywhere in `window`.

    Score for a prefix of length k with m matches is m - mm_penalty*(k-m), so a
    longer repeat only wins if the extra bases are mostly matches. Returns
    (score, pos, k, identity); pos is the window offset of the match start.
    """
    q = enc(query[:max_k]); w = enc(window)
    n = len(w) - min_k + 1
    if n <= 0 or len(q) < min_k: return (-10**9, -1, 0, 0.0)
    best = np.full(n, -10**9, dtype=np.int32)
    best_k = np.zeros(n, dtype=np.int32)
    m = np.zeros(n, dtype=np.int32)
    for i in range(min(max_k, len(q))):
        seg = w[i:i + n] if i + n <= len(w) else np.pad(w[i:], (0, i + n - len(w)), constant_values=0)
        m += (seg == q[i])
        k = i + 1
        if k < min_k: continue
        sc = m * (1 + mm_penalty) - mm_penalty * k
        upd = sc > best
        best = np.where(upd, sc, best)
        best_k = np.where(upd, k, best_k)
    p = int(np.argmax(best))
    k = int(best_k[p]); s = int(best[p])
    ident = 0.0
    if k:
        mm = (s + mm_penalty * k) / (1 + mm_penalty)
        ident = mm / k
    return (s, p, k, ident)

def null_scores(query, fa, seqname, avoid, n_decoy, wlen, rng, min_k, max_k):
    """Same query, windows from distant loci on the same sequence."""
    L = fa.length(seqname); out = []
    tries = 0
    while len(out) < n_decoy and tries < n_decoy * 8:
        tries += 1
        start = rng.randrange(0, max(1, L - wlen))
        if abs(start - avoid) < 200_000: continue
        w = fa.fetch(seqname, start, start + wlen)
        if len(w) < wlen // 2 or w.count("N") > wlen // 10: continue
        out.append(scan(query, w, min_k, max_k)[0])
    return out

# ---------------------------------------------------------------- ground truth
def attrs(c9):
    d = {}
    for kv in c9.rstrip(";").split(";"):
        if "=" in kv: k, v = kv.split("=", 1); d[k] = v
    return d

def load_truth(gff, kind):
    """Elements with a KNOWN TSD length: (seq, start, end, strand, tsd_len)."""
    out = []
    for line in open(gff):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9: continue
        a = attrs(f[8])
        if kind == "ltr":
            if f[2] != "transposable_element": continue
            t = a.get("TSD", "").split("/")[0]
        else:
            if not (f[2] == "sequence_feature" and f[1] == "DANTE_TIR"): continue
            t = a.get("tsd", "")
        # DANTE writes the SENTINEL "not_found" in this field when it failed to
        # find a TSD. It is 9 characters long, so a naive len() counts it as a
        # 9 bp TSD -- which is exactly a LINE-range length, and would silently
        # fill the validation set with elements that have no TSD at all.
        if not t or not set(t.upper()) <= set("ACGT"): continue
        L = len(t)
        if L < 1: continue
        out.append((f[0], int(f[3]), int(f[4]), f[6], L))
    return out

# ---------------------------------------------------------------- validation
def validate(fa, truth, window, min_k, max_k, n_decoy, seed, limit, tol=1):
    rng = random.Random(seed)
    truth = [t for t in truth if t[4] >= min_k]
    rng.shuffle(truth)
    truth = truth[:limit]
    hit = miss = 0; ranks = []; obs_scores = []; null_all = []
    correct_scores = []; wrong_scores = []
    geom_ok = geom_n = 0; skipped_long = 0
    for seqname, s, e, strand, L in truth:
        # The search runs UPSTREAM FROM THE ANCHOR (the element 3' end), so the
        # window must span the whole element or the true 5' TSD copy is simply
        # not inside it. Silently keeping those loci scores the finder 0% for a
        # reason that has nothing to do with the finder.
        if (e - s + 1) + L > window - 100:
            skipped_long += 1; continue
        pad = window + L + 60
        a = max(0, s - 1 - pad); b = e + pad
        reg = fa.fetch(seqname, a, b)
        if len(reg) < (b - a) // 2 or reg.count("N") > len(reg) // 10: continue
        S0 = (s - 1) - a; E0 = (e - 1) - a      # 0-based local, genome orientation
        if strand == "-":
            Lr = len(reg); reg = rc(reg)
            S, E = Lr - 1 - E0, Lr - 1 - S0     # use len(reg): the fetch may clamp
        else:
            S, E = S0, E0
        if not (0 <= S < E < len(reg)): continue
        query = reg[E + 1:E + 1 + max_k]
        if len(query) < min_k: continue
        wstart = max(0, E - window)
        win = reg[wstart:E]
        if len(win) < min_k: continue
        # geometry self-check, independent of the finder: the two TSD copies
        # must actually match each other in the genome.
        c5 = reg[S - L:S]; c3 = reg[E + 1:E + 1 + L]
        if len(c5) == L and len(c3) == L:
            geom_n += 1
            geom_ok += sum(x == y for x, y in zip(c5, c3)) >= 0.8 * L
        sc, p, k, ident = scan(query, win, min_k, max_k)
        true_pos = (S - L) - wstart
        ok = abs(p - true_pos) <= tol
        hit += ok; miss += (not ok)
        obs_scores.append(sc)
        (correct_scores if ok else wrong_scores).append(sc)
        nulls = null_scores(query, fa, seqname, a, n_decoy, window, rng, min_k, max_k)
        null_all.extend(nulls)
        if nulls:
            ranks.append(sum(1 for x in nulls if x >= sc) / len(nulls))
    return dict(n=hit + miss, hit=hit, miss=miss, ranks=ranks,
                obs=obs_scores, null=null_all, geom_ok=geom_ok, geom_n=geom_n,
                skipped_long=skipped_long,
                correct=correct_scores, wrong=wrong_scores)

def pct(v, p):
    v = sorted(v); return v[min(len(v) - 1, int(p * (len(v) - 1)))] if v else 0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True); ap.add_argument("--gff", required=True)
    ap.add_argument("--kind", choices=("ltr", "tir"), required=True)
    ap.add_argument("--window", type=int, default=2000)
    ap.add_argument("--min-k", type=int, default=8); ap.add_argument("--max-k", type=int, default=20)
    ap.add_argument("--decoys", type=int, default=12)
    ap.add_argument("--limit", type=int, default=200); ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--label", default="")
    a = ap.parse_args()

    fa = Faidx(a.fasta)
    truth = load_truth(a.gff, a.kind)
    r = validate(fa, truth, a.window, a.min_k, a.max_k, a.decoys, a.seed, a.limit)
    if not r["n"]:
        print(f"{a.label}: no usable elements with TSD >= {a.min_k}"); return
    rec = 100 * r["hit"] / r["n"]
    print(f"### {a.label or a.gff}   window={a.window}  k={a.min_k}-{a.max_k}  n={r['n']}"
          f"  (skipped {r['skipped_long']} elements longer than the window)")
    g = 100 * r["geom_ok"] / r["geom_n"] if r["geom_n"] else 0
    print(f"  GEOMETRY self-check: the two TSD copies match in {r['geom_ok']}/{r['geom_n']} ({g:.1f}%)"
          f"   <- must be high, else the coordinates are wrong, not the finder")
    print(f"  RECALL  finder lands on the true TSD: {r['hit']}/{r['n']} ({rec:.1f}%)")
    print(f"  observed score  median {pct(r['obs'],.5)}   correct {pct(r['correct'],.5)}   wrong {pct(r['wrong'],.5)}")
    print(f"  NULL (decoy windows, n={len(r['null'])})  median {pct(r['null'],.5)}  "
          f"p95 {pct(r['null'],.95)}  p99 {pct(r['null'],.99)}  max {max(r['null']) if r['null'] else 0}")
    thr = pct(r["null"], .99)
    tp = sum(1 for s in r["correct"] if s > thr); fp = sum(1 for s in r["wrong"] if s > thr)
    print(f"  at the null p99 threshold (score > {thr}):")
    print(f"     kept {tp+fp}/{r['n']} loci ({100*(tp+fp)/r['n']:.1f}%)   "
          f"PRECISION {100*tp/(tp+fp) if tp+fp else 0:.1f}%   recall-of-true {100*tp/max(1,r['hit']):.1f}%")
    if r["ranks"]:
        conf = sum(1 for x in r["ranks"] if x == 0)
        print(f"  loci beating EVERY decoy: {conf}/{len(r['ranks'])} ({100*conf/len(r['ranks']):.1f}%)")

if __name__ == "__main__":
    main()
