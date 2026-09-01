#!/usr/bin/env python3
"""Score boundary estimators against the ground truth produced by eval_diagonal.py."""
import csv, math, sys
from collections import defaultdict

def load(p):
    rows=[r for r in csv.DictReader(open(p),delimiter='\t')]
    for r in rows:
        for k in ("true_q","true_r","true_min","L_score","L_score_r","aln_q",
                  "nb20","b20","nb50","b50","off10","off20","off50","off100"):
            r[k]=int(r[k])
    return rows

def seq_level(rows, sf=0.5, cap=0, cutkey=None):
    per=defaultdict(list); truth={}
    for r in rows:
        cut = r[cutkey] if (cutkey and r[cutkey]>=0) else None
        for sid,L,t in ((r["q"],r["L_score"],r["true_q"]),(r["r"],r["L_score_r"],r["true_r"])):
            if cut is not None: L=min(L,cut)
            key=(r["group"],r["side"],sid); per[key].append(L); truth[key]=t
    out=[]
    for key,L in per.items():
        n=len(L); k=max(3, math.ceil(sf*n)) if sf>0 else 3
        if n<k or n<5: continue
        v=sorted(L,reverse=True)[k-1]
        if cap: v=min(v,cap)
        out.append((v,truth[key]))
    return out

def stats(res):
    if not res: return None
    e=sorted(x-t for x,t in res); q=lambda p: e[int(p*(len(e)-1))]
    return dict(n=len(e), med=q(.5),
                ov3=100*sum(1 for x in e if x>300)/len(e),
                ov1k=100*sum(1 for x in e if x>1000)/len(e),
                un3=100*sum(1 for x in e if x<-300)/len(e),
                rec=100*sum(min(x,t) for x,t in res)/sum(t for _,t in res))

def report(label, path):
    rows=load(path)
    if not rows: print(f"  {label}: empty"); return
    print(f"\n########## {label}   ({len(rows)} pairs) ##########")
    print(f"  {'policy':<26}{'n':>5}{'med err':>9}{'over300':>9}{'over1k':>8}{'under300':>10}{'recovered':>11}")
    for lab,sf,cap in (("shipping sf.5 + cap1500",0.5,1500),("sf.5, no cap",0.5,0),
                       ("sf0 (old) + cap1500",0.0,1500),("sf0, no cap",0.0,0)):
        s=stats(seq_level(rows,sf,cap))
        if s: print(f"  {lab:<26}{s['n']:5d}{s['med']:+9d}{s['ov3']:8.1f}%{s['ov1k']:7.1f}%{s['un3']:9.1f}%{s['rec']:10.1f}%")
    over=[r for r in rows if r["L_score"]-r["true_min"]>300]
    print(f"  overshoot >300bp: {len(over)}/{len(rows)} ({100*len(over)/len(rows):.1f}%)", end="")
    if over:
        a=[r for r in over if 0<=r['off20']<=r['true_min']+200]
        b=[r for r in over if 0<=r['b20'] <=r['true_min']+200]
        print(f"   drift catches: offset>=20 {100*len(a)/len(over):.0f}%  burst>=20 {100*len(b)/len(over):.0f}%")
    else: print()
    b0=stats(seq_level(rows,0.5,1500)); b1=stats(seq_level(rows,0.5,1500,"off20"))
    if b0 and b1:
        print(f"  cost of drift truncation: recovered {b0['rec']:.1f}% -> {b1['rec']:.1f}%   "
              f"over300 {b0['ov3']:.1f}% -> {b1['ov3']:.1f}%")

for arg in sys.argv[1:]:
    lab,_,path = arg.partition("=")
    report(lab, path)
