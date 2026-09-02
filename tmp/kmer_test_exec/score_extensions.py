#!/usr/bin/env python3
"""Score a DANTE_LINE.gff3 against confirmed-boundary ground truth.

PURITY  = emitted extension bp lying INSIDE the true element / emitted bp.
          Primary. An over-extended consensus mislabels whatever its tail
          matches genome-wide; a truncated one only loses a little masking.
COVERAGE= emitted bp inside the true element / true extension bp.
          Secondary: how much real element was recovered.

Loci are matched between the ground truth and a GFF3 produced under different
parameters by the DOMAIN CORE coordinates, which do not move when the extension
parameters change.
"""
import argparse, bisect, csv
from collections import defaultdict


def load_emitted(gff):
    """(seq, core_start, core_end) -> (ext5, ext3)"""
    els, by_seq, dom = {}, defaultdict(list), defaultdict(list)
    for line in open(gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        d = dict(kv.split("=", 1) for kv in f[8].rstrip(";").split(";") if "=" in kv)
        if f[2] == "LINE_element":
            els[d["ID"]] = (f[0], int(d.get("Extension_5prime", 0)),
                            int(d.get("Extension_3prime", 0)))
            by_seq[f[0]].append((int(f[3]), int(f[4]), d["ID"]))
        elif f[2] == "protein_domain":
            dom[f[0]].append((int(f[3]), int(f[4])))
    core = {}
    for seq, lst in by_seq.items():
        lst.sort()
        starts = [x[0] for x in lst]
        for ds, de in dom.get(seq, []):
            i = bisect.bisect_right(starts, ds) - 1
            if i >= 0 and lst[i][0] <= ds and de <= lst[i][1]:
                k = lst[i][2]
                c = core.get(k)
                core[k] = (min(c[0], ds), max(c[1], de)) if c else (ds, de)
    out = {}
    for eid, (seq, e5, e3) in els.items():
        if eid in core:
            out[(seq, core[eid][0], core[eid][1])] = (e5, e3)
    return out


def score(gt_tsv, gff):
    emitted = load_emitted(gff)
    n = matched = 0
    emit_bp = inside_bp = true_bp = 0
    for r in csv.DictReader(open(gt_tsv), delimiter="\t"):
        # Pre-registered filter: a saturated 20 bp match is not a TSD (published
        # range 9-18 nt), and an element over 7500 bp is not a plant LINE (4-7 kb).
        if int(r["score"]) >= 20 or int(r["elem_len"]) > 7500:
            continue
        n += 1
        key = (r["seq"], int(r["core_start"]), int(r["core_end"]))
        if key not in emitted:
            continue
        matched += 1
        e5, e3 = emitted[key]
        t5, t3 = int(r["ext5"]), int(r["ext3"])
        emit_bp += e5 + e3
        inside_bp += min(e5, t5) + min(e3, t3)
        true_bp += t5 + t3
    return dict(n=n, matched=matched, emit=emit_bp, inside=inside_bp, true=true_bp,
                purity=(100 * inside_bp / emit_bp) if emit_bp else None,
                coverage=(100 * inside_bp / true_bp) if true_bp else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--gff", required=True, nargs="+",
                    help="one or more LABEL=path/to/DANTE_LINE.gff3")
    ap.add_argument("--label", default="")
    a = ap.parse_args()
    print(f"### {a.label}")
    print(f"  {'threshold':<12}{'loci':>6}{'emitted bp':>12}{'purity':>9}{'coverage':>10}")
    for spec in a.gff:
        lab, _, path = spec.partition("=")
        s = score(a.truth, path)
        p = f"{s['purity']:.1f}%" if s["purity"] is not None else "n/a"
        c = f"{s['coverage']:.1f}%" if s["coverage"] is not None else "n/a"
        print(f"  {lab:<12}{s['matched']:6d}{s['emit']:12,d}{p:>9}{c:>10}")


main()
