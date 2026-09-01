#!/usr/bin/env python3
"""Go/no-go for a TPRT-signature boundary step: is a poly-A tail findable?

Every later step (TSD search anchored at the tail, 3'-anchored MSA) depends on
locating the 3' end first. This measures that alone, with the element's own 5'
side as the NEGATIVE CONTROL -- a real poly-A must be enriched downstream of the
core and NOT upstream. Reported against local base composition, because an
absolute A-fraction is meaningless in an AT-rich genome.
"""
import sys, bisect
from collections import defaultdict

class Faidx:
    def __init__(self, fa):
        self.f=open(fa,"rb"); self.idx={}
        for line in open(fa+".fai"):
            n,ln,off,lb,lw=line.split("\t")[:5]
            self.idx[n]=(int(ln),int(off),int(lb),int(lw))
    def fetch(self,name,a,b):
        if name not in self.idx: return ""
        ln,off,lb,lw=self.idx[name]; a=max(0,a); b=min(ln,b)
        if b<=a: return ""
        s=off+a//lb*lw+a%lb; e=off+(b-1)//lb*lw+(b-1)%lb+1
        self.f.seek(s); return self.f.read(e-s).decode().replace("\n","").replace("\r","").upper()

_C=str.maketrans("ACGTN","TGCAN")
def comp(s): return s.translate(_C)
def rc(s): return s.translate(_C)[::-1]

def attrs(c9):
    d={}
    for kv in c9.rstrip(";").split(";"):
        if "=" in kv: k,v=kv.split("=",1); d[k]=v
    return d

def best_a_window(seq, w=20):
    """(max A-fraction, offset) over sliding windows."""
    if len(seq) < w: return 0.0, -1
    cnt = seq[:w].count("A"); best, at = cnt, 0
    for i in range(1, len(seq)-w+1):
        cnt += (seq[i+w-1]=="A") - (seq[i-1]=="A")
        if cnt > best: best, at = cnt, i
    return best/w, at

def main(gff, fasta, window=600, thr=0.70):
    fa=Faidx(fasta)
    els=[]; by_seq=defaultdict(list); doms=defaultdict(list)
    for line in open(gff):
        if line.startswith("#"): continue
        f=line.rstrip("\n").split("\t")
        if len(f)<9: continue
        if f[2]=="LINE_element":
            els.append((f[0],int(f[3]),int(f[4]),f[6])); by_seq[f[0]].append((int(f[3]),int(f[4]),len(els)-1))
        elif f[2]=="protein_domain":
            doms[f[0]].append((int(f[3]),int(f[4])))
    core={}
    for seq,lst in by_seq.items():
        lst.sort(); st=[x[0] for x in lst]
        for ds,de in doms.get(seq,[]):
            i=bisect.bisect_right(st,ds)-1
            if i>=0 and lst[i][0]<=ds and de<=lst[i][1]:
                k=lst[i][2]; c=core.get(k)
                core[k]=(min(c[0],ds),max(c[1],de)) if c else (ds,de)
    hits3=[]; hits5=[]; dists=[]
    for k,(seq,s,e,strand) in enumerate(els):
        if k not in core: continue
        cs,ce=core[k]
        g_up=fa.fetch(seq,cs-1-window,cs-1); g_dn=fa.fetch(seq,ce,ce+window)
        if strand=="-": f3,f5 = rc(g_up), comp(g_dn)
        else:           f3,f5 = g_dn, g_up[::-1]
        if len(f3)<50 or len(f5)<50: continue
        a3,off3=best_a_window(f3); a5,_=best_a_window(f5)
        hits3.append(a3); hits5.append(a5)
        if a3>=thr: dists.append(off3)
    n=len(hits3)
    if not n: print("no usable elements"); return
    p3=100*sum(1 for x in hits3 if x>=thr)/n
    p5=100*sum(1 for x in hits5 if x>=thr)/n
    print(f"  elements with a core: {n}")
    print(f"  3' side (expected)  : {p3:5.1f}%  have a 20bp window >={int(thr*100)}% A")
    print(f"  5' side (CONTROL)   : {p5:5.1f}%  <- background; the gap is the real signal")
    print(f"  enrichment          : {p3/p5:.2f}x" if p5 else "  enrichment: inf (no 5' background)")
    if dists:
        dists.sort()
        print(f"  distance core->tail : median {dists[len(dists)//2]} bp, p90 {dists[int(.9*(len(dists)-1))]} bp")

main(sys.argv[1], sys.argv[2])
