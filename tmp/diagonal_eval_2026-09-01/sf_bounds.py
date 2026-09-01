#!/usr/bin/env python3
"""Per-superfamily distribution of the TRUE core->boundary distance.

Measured from DANTE_TIR primary successes (exact TIR+TSD boundary), so any cap
below the high quantiles here truncates REAL element.  This is the evidence for
per-superfamily extension bounds instead of one shared number.
"""
import sys, bisect
from collections import defaultdict

def attrs(c9):
    d={}
    for kv in c9.rstrip(";").split(";"):
        if "=" in kv: k,v=kv.split("=",1); d[k]=v
    return d

per=defaultdict(list); elen=defaultdict(list)
for gff in sys.argv[1:]:
    els=[]; by_seq=defaultdict(list); doms=defaultdict(list)
    for line in open(gff):
        if line.startswith("#"): continue
        f=line.rstrip("\n").split("\t")
        if len(f)<9: continue
        a=attrs(f[8])
        if f[2]=="sequence_feature" and f[1]=="DANTE_TIR":
            if not ("tir_seq5" in a and "tir_seq3" in a and a.get("tsd")): continue
            els.append((f[0],int(f[3]),int(f[4]),a.get("Classification","?")))
            by_seq[f[0]].append((int(f[3]),int(f[4]),len(els)-1))
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
    for k,(seq,s,e,cls) in enumerate(els):
        if k not in core: continue
        cs,ce=core[k]
        sf=cls.replace("Class_II_Subclass_1_TIR_","")
        per[sf]+= [cs-s, e-ce]
        elen[sf].append(e-s+1)

def q(v,p):
    v=sorted(v); return v[min(len(v)-1,int(p*(len(v)-1)))]
print(f"  {'superfamily':<16}{'sides':>7}{'median':>8}{'p75':>7}{'p90':>7}{'p95':>7}{'p99':>7}{'max':>8}   {'elem p95':>9}")
for sf in sorted(per, key=lambda s: -len(per[s])):
    v=[x for x in per[sf] if x>=0]
    if len(v)<20: continue
    print(f"  {sf:<16}{len(v):7d}{q(v,.5):8d}{q(v,.75):7d}{q(v,.90):7d}{q(v,.95):7d}{q(v,.99):7d}{max(v):8d}   {q(elen[sf],.95):9d}")
print("\n  (current shipping caps: dante_line --max-extension 1500 per side;")
print("   dante_tir_fallback has NO per-side cap, only the vocabulary whole-element bound)")
