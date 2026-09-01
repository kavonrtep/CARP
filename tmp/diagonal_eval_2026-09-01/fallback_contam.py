#!/usr/bin/env python3
"""Contamination of DANTE_TIR_fallback EXTENSIONS on the real population.

The fallback runs exactly where DANTE_TIR failed, so no exact boundary exists.
But contamination is measurable without one: how much of the inferred flank
lands on a structural annotation of an INCOMPATIBLE class (a DANTE_LTR element,
a TideCluster satellite)?  The TPase core is the built-in control -- it should
be clean, and any core signal is the background rate.
"""
import sys, bisect
from collections import defaultdict

def attrs(c9):
    d={}
    for kv in c9.rstrip(";").split(";"):
        if "=" in kv: k,v=kv.split("=",1); d[k]=v
    return d

def load_intervals(path, want_source=None, want_type=None, cls_filter=None):
    iv=defaultdict(list)
    try: fh=open(path)
    except FileNotFoundError: return iv
    for line in fh:
        if line.startswith("#"): continue
        f=line.rstrip("\n").split("\t")
        if len(f)<9: continue
        if want_source and f[1]!=want_source: continue
        if want_type and f[2]!=want_type: continue
        if cls_filter:
            a=attrs(f[8])
            c=(a.get("Final_Classification") or a.get("Classification") or a.get("Name") or "").replace("|","/")
            if not cls_filter(c): continue
        iv[f[0]].append((int(f[3]),int(f[4])))
    return iv

def merge(iv):
    out={}
    for s,v in iv.items():
        v.sort(); m=[]
        for a,b in v:
            if m and a<=m[-1][1]+1: m[-1][1]=max(m[-1][1],b)
            else: m.append([a,b])
        out[s]=(([x[0] for x in m]),([x[1] for x in m]))
    return out

def ov_bp(merged, seq, a, b):
    if seq not in merged: return 0
    st,en=merged[seq]
    i=bisect.bisect_right(st,b)
    tot=0
    j=i-1
    while j>=0 and en[j]>=a:
        tot+=min(b,en[j])-max(a,st[j])+1
        j-=1
        if j>=0 and en[j]<a-100000: break
    return max(0,tot)

def fallback_parts(gff):
    """yield (seq, kind, start, end) for core / 5' ext / 3' ext of each fallback element"""
    els=[]; doms=defaultdict(list); by_seq=defaultdict(list)
    for line in open(gff):
        if line.startswith("#"): continue
        f=line.rstrip("\n").split("\t")
        if len(f)<9 or f[1]!="DANTE_TIR_fallback": continue
        a=attrs(f[8])
        if f[2]=="sequence_feature":
            els.append((f[0],int(f[3]),int(f[4]),a.get("ID")))
            by_seq[f[0]].append((int(f[3]),int(f[4]),len(els)-1))
        elif f[2]=="protein_domain":
            doms[f[0]].append((int(f[3]),int(f[4])))
    # pair domain -> containing element
    idx=defaultdict(lambda:(None,None))
    core={}
    for seq,lst in by_seq.items():
        lst.sort(); starts=[x[0] for x in lst]
        for ds,de in doms.get(seq,[]):
            i=bisect.bisect_right(starts,ds)-1
            if i>=0 and lst[i][0]<=ds and de<=lst[i][1]:
                k=lst[i][2]
                cur=core.get(k)
                core[k]=(min(cur[0],ds),max(cur[1],de)) if cur else (ds,de)
    for k,(seq,s,e,_) in enumerate(els):
        c=core.get(k)
        if not c: continue
        cs,ce=c
        eid=els[k][3]
        yield seq,"core",cs,ce,eid
        if cs>s:  yield seq,"ext",s,cs-1,eid
        if e>ce:  yield seq,"ext",ce+1,e,eid

for arg in sys.argv[1:]:
    lab,_,root=arg.partition("=")
    inc=defaultdict(list)
    for s,v in load_intervals(f"{root}/DANTE_LTR/DANTE_LTR.gff3","dante_ltr","transposable_element").items():
        inc[s]+=v
    for tc in ("TideCluster/default/TideCluster_clustering.gff3",
               "TideCluster/short_monomer/TideCluster_clustering.gff3"):
        for s,v in load_intervals(f"{root}/{tc}").items(): inc[s]+=v
    # RepeatMasker Class_I is far broader coverage and is INDEPENDENT of the
    # fallback (fallback consensi are not in the library by default).
    for s,v in load_intervals(f"{root}/RepeatMasker/RM_on_combined_library.gff3",
                              cls_filter=lambda c: c.startswith("Class_I/")).items():
        inc[s]+=v
    M=merge(inc)
    ref_bp=sum(e-s0+1 for s0,e in ((st,en) for k in M for st,en in zip(*M[k])))
    print(f"\n### {lab}\n   reference (incompatible structure) covers {ref_bp:,} bp")
    agg=defaultdict(lambda:[0,0])
    per_el=defaultdict(lambda:[0,0])          # element -> [ext_bp, contam_bp]
    for seq,kind,a,b,eid in fallback_parts(f"{root}/DANTE_TIR/DANTE_TIR_combined.gff3"):
        o=ov_bp(M,seq,a,b)
        agg[kind][0]+=b-a+1; agg[kind][1]+=o
        if kind=="ext":
            per_el[eid][0]+=b-a+1; per_el[eid][1]+=o
    for kind in ("core","ext"):
        tot,o=agg[kind]
        if tot: print(f"   {kind:5s} {tot:12,d} bp   on incompatible structure {o:11,d} bp  ({100*o/tot:5.2f}%)")
    ext_els=[v for v in per_el.values() if v[0]>0]
    if ext_els:
        for thr in (1,50,200):
            n=sum(1 for e,c in ext_els if c>=thr)
            print(f"   elements whose extension carries >={thr:3d} bp contaminant: "
                  f"{n:5d}/{len(ext_els)} ({100*n/len(ext_els):5.2f}%)")
