#!/usr/bin/env python3
"""Score extension as a RELIABLE CORE, not as a complete element.

Purity  = fraction of emitted extension bp that lies INSIDE the true element.
          Under-extension costs nothing here; only overshoot does.
Contam  = bp emitted past the true boundary -- the quantity that poisons RM.
Shared-extent test: within a family, how variable is the true core->end
distance?  If it is variable (Petr: the TPase->TIR interval often is), then no
consensus can span it and the algorithm SHOULD stop near the family minimum.
"""
import csv, math, sys, statistics as st
from collections import defaultdict

def load(p):
    rows=[r for r in csv.DictReader(open(p),delimiter='\t')]
    for r in rows:
        for k in ("true_q","true_r","true_min","L_score","L_score_r"): r[k]=int(r[k])
    return rows

def sides(rows, sf=0.5, cap=0):
    per=defaultdict(list); truth={}
    for r in rows:
        for sid,L,t in ((r["q"],r["L_score"],r["true_q"]),(r["r"],r["L_score_r"],r["true_r"])):
            key=(r["group"],r["side"],sid); per[key].append(L); truth[key]=t
    out=[]
    for key,L in per.items():
        n=len(L); k=max(3, math.ceil(sf*n)) if sf>0 else 3
        if n<k or n<5: continue
        v=sorted(L,reverse=True)[k-1]
        if cap: v=min(v,cap)
        out.append((key,v,truth[key]))
    return out

for arg in sys.argv[1:]:
    lab,_,path=arg.partition("=")
    rows=load(path)
    res=sides(rows)
    if not res: print(f"{lab}: empty"); continue
    emit=sum(v for _,v,_ in res); inside=sum(min(v,t) for _,v,t in res)
    contam=sum(max(0,v-t) for _,v,t in res)
    dirty=[1 for _,v,t in res if v-t>50]
    print(f"\n########## {lab} ##########")
    print(f"  sides {len(res)}   emitted extension {emit:,} bp")
    print(f"  PURITY  {100*inside/emit:6.2f}%  of emitted bp lies inside the true element")
    print(f"  CONTAM  {contam:,} bp ({100*contam/emit:.2f}% of emitted)   "
          f"sides overshooting >50bp: {len(dirty)} ({100*len(dirty)/len(res):.1f}%)")
    # --- does the family's true extension even HAVE a shared extent? --------
    fam=defaultdict(list)
    for (g,side,sid),v,t in res: fam[(g,side)].append((v,t))
    sp=[]; pos=[]
    for k,vals in fam.items():
        if len(vals)<4: continue
        ts=sorted(t for _,t in vals)
        med=ts[len(ts)//2]; lo=ts[int(.1*(len(ts)-1))]
        if med>0: sp.append((ts[int(.9*(len(ts)-1))]-lo)/med)
        L=st.median([v for v,_ in vals])
        # where does the inferred length sit in the family's true-extent distribution?
        pos.append(100*sum(1 for t in ts if t<=L)/len(ts))
    if sp:
        print(f"  WITHIN-FAMILY SPREAD of true core->end distance: median (p90-p10)/median = {st.median(sp):.2f}")
        print(f"  inferred length sits at the {st.median(pos):.0f}th percentile of its family's true extents")
