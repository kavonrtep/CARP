#!/usr/bin/env python3
"""Phase 0 deliverable: the chance-score floor as a function of search-window size.

A perfect TSD of length L scores exactly L under the finder's scoring
(3*matches - 2*k = k when perfect). So a TSD is separable from chance only if
L exceeds the chance floor for the window it is searched in. This measures that
floor directly -- no ground truth needed -- by scoring random genomic queries
against random genomic windows of each size.

The practical consequence: the domain anchor is what makes this tractable. It
bounds the 5' search to roughly the per-superfamily extension cap rather than a
whole element length, and the floor falls steeply as the window shrinks.
"""
import random, sys
import numpy as np
sys.path.insert(0, "tmp/line_tprt_2026-09-01")
from tsd_finder import Faidx, scan   # noqa: E402

def curve(fa, seqnames, windows, n=250, min_k=8, max_k=20, seed=0):
    rng = random.Random(seed)
    print(f"  {'window':>8}{'n':>6}{'median':>8}{'p90':>6}{'p95':>6}{'p99':>6}{'max':>6}   "
          f"{'min TSD len separable at p99':>30}")
    for W in windows:
        scores = []
        while len(scores) < n:
            s1 = rng.choice(seqnames); L1 = fa.length(s1)
            s2 = rng.choice(seqnames); L2 = fa.length(s2)
            if L1 < 10_000 or L2 < W + 10_000: continue
            q = fa.fetch(s1, rng.randrange(0, L1 - max_k - 1), 0) or ""
            qp = rng.randrange(0, L1 - max_k - 1)
            q = fa.fetch(s1, qp, qp + max_k)
            wp = rng.randrange(0, L2 - W)
            w = fa.fetch(s2, wp, wp + W)
            if len(q) < max_k or len(w) < W: continue
            if q.count("N") or w.count("N") > W // 20: continue
            scores.append(scan(q, w, min_k, max_k)[0])
        v = sorted(scores); q_ = lambda p: v[min(len(v) - 1, int(p * (len(v) - 1)))]
        print(f"  {W:8d}{len(v):6d}{q_(.5):8d}{q_(.9):6d}{q_(.95):6d}{q_(.99):6d}{v[-1]:6d}"
              f"{q_(.99) + 1:30d}")

import sys as _sys
_args = _sys.argv[1:]
if not _args:
    _args = ["/nfsroot/projects/darwin/runs/tmp3/GCA_982474395.1/carp_output/genome_cleaned.fasta::Lycopus"]
for _spec in _args:
    _fa, _, _lab = _spec.partition("::")
    fa = Faidx(_fa)
    names = [n for n in fa.idx if fa.length(n) > 2_000_000][:40]
    if not names:
        names = sorted(fa.idx, key=lambda n: -fa.length(n))[:40]
    print(f"### {_lab or _fa}  ({len(names)} sequences)")
    curve(fa, names, [2000, 4000, 6000, 8000])
