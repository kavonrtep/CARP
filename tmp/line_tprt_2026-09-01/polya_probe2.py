#!/usr/bin/env python3
"""Poly-A anchor detectability, corrected.

v1 was wrong twice: it searched only 600 bp (the RT domain end is NOT the
element end -- the ORF2 C-terminal domain plus the 3'UTR lie between them), and
it scored "70% A in a 20 bp window", which coding AT-richness satisfies. This
version looks for a genuine homopolymer run over a realistic range, with the
element's own 5' side as the negative control.
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

def longest_run(seq, ch, max_mm=1):
    """Longest run of `ch` tolerating up to max_mm mismatches, and where it starts."""
    best=0; at=-1
    n=len(seq); i=0
    while i < n:
        if seq[i]!=ch: i+=1; continue
        mm=0; j=i; last=i
        while j < n:
            if seq[j]==ch: last=j
            else:
                mm+=1
                if mm>max_mm: break
            j+=1
        L=last-i+1
        if L>best: best, at = L, i
        i=last+1 if last>i else i+1
    return best, at

def run(gff, fasta, window, minrun):
    fa=Faidx(fasta)
    els=[]; by_seq=defaultdict(list); doms=defaultdict(list)
    for line in open(gff):
        if line.startswith("#"): continue
        f=line.rstrip("\n").split("\t")
        if len(f)<9: continue
        if f[2]=="LINE_element":
            els.append((f[0],int(f[3]),int(f[4]),f[6])); by_seq[f[0]].append((int(f[3]),int(f[4]),len(els)-1))
        elif f[2]=="protein_domain": doms[f[0]].append((int(f[3]),int(f[4])))
    core={}
    for seq,lst in by_seq.items():
        lst.sort(); st=[x[0] for x in lst]
        for ds,de in doms.get(seq,[]):
            i=bisect.bisect_right(st,ds)-1
            if i>=0 and lst[i][0]<=ds and de<=lst[i][1]:
                k=lst[i][2]; c=core.get(k); core[k]=(min(c[0],ds),max(c[1],de)) if c else (ds,de)
    n=h3=h5=0; dists=[]
    for k,(seq,s,e,strand) in enumerate(els):
        if k not in core: continue
        cs,ce=core[k]
        g_up=fa.fetch(seq,cs-1-window,cs-1); g_dn=fa.fetch(seq,ce,ce+window)
        if strand=="-": f3,f5 = rc(g_up), comp(g_dn)
        else:           f3,f5 = g_dn, g_up[::-1]
        if len(f3)<window//2 or len(f5)<window//2: continue
        n+=1
        r3,at3=longest_run(f3,"A"); r5,_=longest_run(f5,"A")
        if r3>=minrun: h3+=1; dists.append(at3)
        if r5>=minrun: h5+=1
    return n,h3,h5,dists

for lab,acc,root in (("boechera","GCA_018361405.1","/nfsroot/projects/darwin/runs/tmp3"),
                     ("lycopus","GCA_982474395.1","/nfsroot/projects/darwin/runs/tmp3")):
    D=f"{root}/{acc}/carp_output"
    print(f"### {lab}")
    print(f"   {'window':>7}{'minrun':>8}{'n':>6}{'3prime':>9}{'5prime ctl':>12}{'enrich':>8}{'median dist':>13}")
    for w in (1000, 2500):
        for mr in (10, 15):
            n,h3,h5,d = run(f"{D}/DANTE_LINE/DANTE_LINE.gff3", f"{D}/genome_cleaned.fasta", w, mr)
            if not n: continue
            d.sort(); md = d[len(d)//2] if d else -1
            en = (h3/h5) if h5 else float('inf')
            print(f"   {w:7d}{mr:8d}{n:6d}{100*h3/n:8.1f}%{100*h5/n:11.1f}%{en:8.2f}{md:13d}")
