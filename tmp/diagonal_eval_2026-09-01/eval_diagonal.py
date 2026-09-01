#!/usr/bin/env python3
"""
Evaluate anchored-diagonal drift as a boundary signal for flank extension.

GROUND TRUTH.  DANTE_LTR elements whose Rank contains 'L' have both LTRs, so
their boundary is exact.  Their protein_domain children give the core.  The
true extension = core edge -> element end.  We hand the algorithm only the
core edge and ask it to rediscover that distance from sequence alone.

TWO ARMS, both read off the SAME parasail traceback, production parameters
(match 2, mismatch -2, gap_open 12, gap_extend 3, one end anchored):
  L_score  the shipping method: query length at max cumulative score
  drift    where the alignment leaves the diagonal (gap burst / net offset)

BLIND-SPOT TEST.  For every pair we record whether an overshoot past the true
boundary came with a drift signal.  Overshoots with NO drift are the case Petr
raised: the element inserts at a fixed position inside another repeat, so the
neighbour is colinear too and no gap is needed to resynchronise.

Flank orientation (both members must be read on the element's own strand,
outward from the core, so index 0 is the anchor):
  strand +   f5 = reverse(g_up)        f3 = g_dn
  strand -   f5 = complement(g_dn)     f3 = revcomp(g_up)
where g_up/g_dn are the W bases immediately 5'/3' of the core IN THE GENOME.
"""
import argparse, os, random, re, sys, importlib.util
from collections import defaultdict

def load_gla(repo):
    spec = importlib.util.spec_from_file_location(
        "gla", os.path.join(repo, "scripts", "global_local_aln.py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

class Faidx:
    def __init__(self, fasta):
        self.f = open(fasta, "rb"); self.idx = {}
        for line in open(fasta + ".fai"):
            n, ln, off, lb, lw = line.split("\t")[:5]
            self.idx[n] = (int(ln), int(off), int(lb), int(lw))
    def fetch(self, name, start, end):
        if name not in self.idx: return ""
        ln, off, lb, lw = self.idx[name]
        start, end = max(0, start), min(ln, end)
        if end <= start: return ""
        s = off + start // lb * lw + start % lb
        e = off + (end - 1) // lb * lw + (end - 1) % lb + 1
        self.f.seek(s)
        return self.f.read(e - s).decode().replace("\n", "").replace("\r", "").upper()

_C = str.maketrans("ACGTN", "TGCAN")
def comp(s): return s.translate(_C)
def rc(s):   return s.translate(_C)[::-1]

def rexdb(a):
    m = re.search(r'REXdb_ID(\d+)', a.get('Best_Hit', ''))
    return m.group(1) if m else '?'


def attrs(c9):
    d = {}
    for kv in c9.rstrip(";").split(";"):
        if "=" in kv: k, v = kv.split("=", 1); d[k] = v
    return d

def parse_ltr(gff):
    elems, doms = {}, defaultdict(list)
    for line in open(gff):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9: continue
        a = attrs(f[8])
        if f[2] == "transposable_element":
            elems[a.get("ID")] = dict(seq=f[0], start=int(f[3]), end=int(f[4]),
                                      strand=f[6], rank=a.get("Rank", ""),
                                      cls=a.get("Final_Classification", "?"))
        elif f[2] == "protein_domain" and "Parent" in a:
            doms[a["Parent"]].append((int(f[3]), int(f[4]), a.get("Name", "?"), rexdb(a)))
    return elems, doms


TIR_BY_CLUSTER = [True]

def parse_tir(gff):
    """DANTE_TIR primary elements: boundary is exact (TIRs on both ends + TSD).

    The element is a `sequence_feature` carrying tir_seq5/tir_seq3/tsd; its core
    is the `dante protein_domain` (TPase) lying INSIDE it.  There is no Parent
    attribute, so pair by containment rather than trusting file adjacency.
    """
    import bisect
    elems, raw_doms = {}, defaultdict(list)
    by_seq = defaultdict(list)
    for line in open(gff):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9: continue
        a = attrs(f[8])
        if f[2] == "sequence_feature" and f[1] == "DANTE_TIR":
            if not ("tir_seq5" in a and "tir_seq3" in a and a.get("tsd")): continue
            eid = a.get("ID")
            elems[eid] = dict(seq=f[0], start=int(f[3]), end=int(f[4]), strand=f[6],
                              rank="L",           # exact boundary, same gate as LTR
                              cls=a.get("Classification", "?") + ("|cl" + a.get("Cluster_ID", "0") if TIR_BY_CLUSTER[0] else ""))
            by_seq[f[0]].append((int(f[3]), int(f[4]), eid))
        elif f[2] == "protein_domain":
            raw_doms[f[0]].append((int(f[3]), int(f[4]), a.get("Name", "?"), rexdb(a)))
    doms = defaultdict(list)
    for seq, els in by_seq.items():
        els.sort()
        starts = [e[0] for e in els]
        for ds, de, nm, rx in raw_doms.get(seq, []):
            i = bisect.bisect_right(starts, ds) - 1
            if i >= 0 and els[i][0] <= ds and de <= els[i][1]:
                doms[els[i][2]].append((ds, de, nm, rx))
    return elems, doms

def drift_profile(res):
    """Anchor is at column 0 (end='5': trailing gaps are free).

    offset = query residues consumed - ref residues consumed.  Constant on the
    diagonal; each gap column shifts it by one.  A burst means the two copies
    needed different amounts of sequence to stay aligned -- they have stopped
    being the same element.  Returns (bursts, first_offset_crossing, aln_qlen).
    """
    tb = res.traceback
    q, r = tb.query, tb.ref
    n = len(q); j = n - 1
    while j >= 0 and (q[j] == '-' or r[j] == '-'): j -= 1
    hi = j + 1
    qpos = off = 0
    bursts, first = [], {}
    run_len, run_side, run_q = 0, None, 0
    for i in range(hi):
        side = 'q' if q[i] == '-' else ('r' if r[i] == '-' else None)
        if side is None:
            if run_len: bursts.append((run_q, run_len)); run_len, run_side = 0, None
            qpos += 1
        else:
            if side != run_side:
                if run_len: bursts.append((run_q, run_len))
                run_len, run_side, run_q = 0, side, qpos
            run_len += 1
            if side == 'r': qpos += 1; off += 1
            else: off -= 1
        for t in (10, 20, 50, 100):
            if t not in first and abs(off) >= t: first[t] = qpos
    if run_len: bursts.append((run_q, run_len))
    return bursts, first, qpos


def identity_profile(res, true_bp, acc, lo=-2500, hi=2500, w=100):
    """Accumulate per-column identity binned by distance RELATIVE to the true
    element boundary.  Inside the element (<0) a conserved family should stay
    high; if the interval between core and boundary is variable, identity falls
    BEFORE 0 -- which is the claim under test."""
    tb = res.traceback
    q, r = tb.query, tb.ref
    n = len(q); j = n - 1
    while j >= 0 and (q[j] == '-' or r[j] == '-'): j -= 1
    qpos = 0
    for i in range(j + 1):
        d = qpos - true_bp
        if lo <= d < hi:
            b = (d // w) * w
            a = acc.setdefault(b, [0, 0])
            a[1] += 1
            if q[i] == r[i] and q[i] != '-': a[0] += 1
        if r[i] != '-' or q[i] != '-':
            if q[i] != '-': qpos += 1

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True); ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--repo", default=os.getcwd())
    ap.add_argument("--mode", choices=("ltr","tir"), default="ltr")
    ap.add_argument("--tir-group", choices=("cluster","class"), default="cluster")
    ap.add_argument("--group-by", choices=("default","rexdb"), default="default")
    ap.add_argument("--profile", default=None)
    ap.add_argument("--window", type=int, default=4000)
    ap.add_argument("--max-groups", type=int, default=50)
    ap.add_argument("--max-per-group", type=int, default=25)
    ap.add_argument("--max-pairs", type=int, default=100)
    ap.add_argument("--min-true", type=int, default=150)
    a = ap.parse_args()

    TIR_BY_CLUSTER[0] = (a.tir_group == 'cluster')
    gla = load_gla(a.repo); fa = Faidx(a.fasta)
    elems, doms = (parse_ltr if a.mode == 'ltr' else parse_tir)(a.gff)

    recs = []
    for eid, e in sorted(elems.items()):
        if "L" not in e["rank"]: continue
        d = doms.get(eid)
        if not d: continue
        cs, ce = min(x[0] for x in d), max(x[1] for x in d)
        if not (e["start"] <= cs and ce <= e["end"]): continue
        g_up = fa.fetch(e["seq"], cs - 1 - a.window, cs - 1)
        g_dn = fa.fetch(e["seq"], ce, ce + a.window)
        up, dn = cs - e["start"], e["end"] - ce
        if e["strand"] == "-":
            f5, f3, t5, t3 = comp(g_dn), rc(g_up), dn, up
        else:
            f5, f3, t5, t3 = g_up[::-1], g_dn, up, dn
        recs.append(dict(id=eid, cls=e["cls"],
                         names=tuple(sorted(set(x[2] for x in d))),
                         rex=tuple(sorted(set(f"{x[2]}:{x[3]}" for x in d))),
                         f5=f5, f3=f3, t5=t5, t3=t3))
    sys.stderr.write(f"elements with exact boundaries: {len(recs)}\n")

    groups = defaultdict(list)
    for r in recs:
        gkey = r["rex"] if a.group_by == "rexdb" else (r["names"] if a.mode == "ltr" else ())
        groups[(r["cls"], gkey)].append(r)
    groups = {k: v for k, v in groups.items() if len(v) >= 4}
    order = sorted(groups, key=lambda k: (-len(groups[k]), k))[:a.max_groups]
    sys.stderr.write(f"groups (lineage + identical domain set, n>=4): {len(order)}\n")

    rng = random.Random(0)
    out = open(a.out, "w")
    out.write("\t".join(["group","side","q","r","true_q","true_r","true_min",
                         "L_score","L_score_r","aln_q","nb20","b20","nb50","b50",
                         "off10","off20","off50","off100"]) + "\n")
    func = gla.pick_func("5"); mat = gla.make_matrix(2, -2)
    prof_acc = {}
    nal = 0
    for gk in order:
        mem = sorted(groups[gk], key=lambda r: r["id"])
        if len(mem) > a.max_per_group:
            mem = sorted(rng.sample(mem, a.max_per_group), key=lambda r: r["id"])
        for side, fk, tk in (("5", "f5", "t5"), ("3", "f3", "t3")):
            pool = [m for m in mem if a.min_true <= m[tk] <= a.window
                    and len(m[fk]) >= a.min_true and "N" * 50 not in m[fk]]
            if len(pool) < 2: continue
            pairs = [(i, j) for i in range(len(pool)) for j in range(i + 1, len(pool))]
            if len(pairs) > a.max_pairs: pairs = sorted(rng.sample(pairs, a.max_pairs))
            for i, j in pairs:
                A, B = pool[i], pool[j]
                res = func(A[fk], B[fk], 12, 3, mat)
                sc = gla.per_column_scores_alt(res, match=2, mismatch=-2,
                                               gap_open=12, gap_extend=3, cumulative=False)
                rd = gla.extract_optimal_alignment(res, sc, "5")
                bursts, foff, alnq = drift_profile(res)
                if a.profile: identity_profile(res, min(A[tk], B[tk]), prof_acc)
                b20 = [p for p, L in bursts if L >= 20]
                b50 = [p for p, L in bursts if L >= 50]
                out.write("\t".join(map(str, [
                    f"{gk[0]}|{'+'.join(gk[1])}", side, A["id"], B["id"],
                    A[tk], B[tk], min(A[tk], B[tk]), rd["degapped_query_len"], rd["degapped_ref_len"], alnq,
                    len(b20), b20[0] if b20 else -1,
                    len(b50), b50[0] if b50 else -1,
                    foff.get(10, -1), foff.get(20, -1),
                    foff.get(50, -1), foff.get(100, -1)])) + "\n")
                nal += 1
    out.close(); sys.stderr.write(f"alignments: {nal}\n")
    if a.profile:
        with open(a.profile, "w") as ph:
            ph.write("rel_bin\tmatches\tcolumns\tidentity\n")
            for b in sorted(prof_acc):
                m, t = prof_acc[b]
                ph.write(f"{b}\t{m}\t{t}\t{m/t if t else 0:.4f}\n")

main()
