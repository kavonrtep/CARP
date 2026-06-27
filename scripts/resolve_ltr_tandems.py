#!/usr/bin/env python3
"""Collapse tandem LTR-RT (LTR_RT_TR) structures in DANTE_LTR output.

DANTE_LTR emits *overlapping* intact `transposable_element` features when
same-lineage LTR-RTs are arranged head-to-tail sharing an LTR
(`LTR-INT-LTR-INT-LTR`; Macko-Podgórni et al., Mobile DNA 2025 — "tandem LTR-RT"
/ LTR_RT_TR). Left as-is, the shared LTRs are double-annotated in the unified
annotation. This script detects such tandems and writes a **small side file**
containing only the derived **container** features — one per array
(`structure=LTR_RT_TR`, `copy_number=N`, `members=<member element IDs>`, spanning
the array). It does **not** modify DANTE_LTR.gff3: the input is read-only, so the
LTR library / masking track / report keep seeing every individual element exactly
as DANTE_LTR emitted them. make_unified_annotation.R reads both files — it splits
DANTE_LTR.gff3's elements into members (IDs listed here) vs standalone, then
emits one Level-1 container with the member copies nested as Level-2 children.

Detection: a maximal chain of consecutive **same-lineage, same-strand** complete
elements where each adjacent pair's terminal LTRs coincide — element A's
rightmost `long_terminal_repeat` reciprocally overlaps element B's leftmost LTR
by >= --min-ltr-recip (default 0.5). The shared-LTR test distinguishes a tandem
(boundary LTR shared) from a nested insertion (deep overlap, LTRs don't align)
and self-limits array length, so there is no arbitrary kb locus cap. Generic
cross-tool / non-tandem tier-1 overlap resolution is done downstream in
make_unified_annotation.R, not here.

Usage:  resolve_ltr_tandems.py -i DANTE_LTR.gff3 -o DANTE_LTR_tandems.gff3
"""
import argparse
import sys


def parse_attrs(col9):
    out = {}
    for f in col9.rstrip(";").split(";"):
        f = f.strip()
        if f and "=" in f:
            k, v = f.split("=", 1)
            out[k] = v
    return out


def load(path):
    header, recs = [], []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                header.append(ln.rstrip("\n"))
                continue
            if not ln.strip():
                continue
            c = ln.rstrip("\n").split("\t")
            if len(c) != 9:
                continue
            recs.append({"seqid": c[0], "source": c[1], "type": c[2],
                         "start": int(c[3]), "end": int(c[4]), "score": c[5],
                         "strand": c[6], "phase": c[7], "attrs": parse_attrs(c[8]),
                         "col9": c[8]})
    return header, recs


def recip_overlap(a, b):
    """Reciprocal overlap fraction (of the shorter interval) of (start,end) a,b."""
    ov = min(a[1], b[1]) - max(a[0], b[0]) + 1
    if ov <= 0:
        return 0.0
    return ov / min(a[1] - a[0] + 1, b[1] - b[0] + 1)


def shares_ltr(a, b, ltrs, min_recip):
    """True if upstream element a and downstream element b share a boundary LTR:
    a's rightmost LTR reciprocally overlaps b's leftmost LTR."""
    la, lb = ltrs.get(a["attrs"].get("ID")), ltrs.get(b["attrs"].get("ID"))
    if not la or not lb:
        return False
    a_right = max(la, key=lambda l: l[0])     # rightmost LTR of a (its 3' boundary)
    b_left = min(lb, key=lambda l: l[0])      # leftmost LTR of b (its 5' boundary)
    return recip_overlap(a_right, b_left) >= min_recip


def detect_tandems(complete, ltrs, min_recip):
    """Group sorted complete elements into maximal same-lineage shared-LTR chains.
    Returns list of chains (each a list of >=2 element records)."""
    chains, cur = [], []
    for el in complete:
        prev = cur[-1] if cur else None
        if (prev is not None
                and el["attrs"].get("Final_Classification") == prev["attrs"].get("Final_Classification")
                and el["strand"] == prev["strand"]
                and shares_ltr(prev, el, ltrs, min_recip)):
            cur.append(el)
        else:
            if len(cur) >= 2:
                chains.append(cur)
            cur = [el]
    if len(cur) >= 2:
        chains.append(cur)
    return chains


def render(seqid, source, ftype, start, end, score, strand, phase, col9):
    return "\t".join([seqid, source, ftype, str(start), str(end), score,
                      strand, phase, col9])


def main(argv=None):
    ap = argparse.ArgumentParser(description="Collapse tandem LTR-RT structures.")
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--min-ltr-recip", type=float, default=0.5,
                    help="min reciprocal overlap of shared boundary LTRs [0.5]")
    args = ap.parse_args(argv)

    header, recs = load(args.input)
    tes = [r for r in recs if r["type"] == "transposable_element"]
    # LTR children grouped by their Parent element ID
    ltrs = {}
    for r in recs:
        if r["type"] == "long_terminal_repeat":
            ltrs.setdefault(r["attrs"].get("Parent"), []).append((r["start"], r["end"]))

    complete = sorted([t for t in tes
                       if not str(t["attrs"].get("ID", "")).startswith("TE_partial")],
                      key=lambda t: (t["seqid"], t["start"]))
    # chain only within a seqid
    chains, by_seq = [], {}
    for t in complete:
        by_seq.setdefault(t["seqid"], []).append(t)
    for lst in by_seq.values():
        chains.extend(detect_tandems(lst, ltrs, args.min_ltr_recip))

    # One container record per chain. This file holds ONLY the derived containers
    # (the array span + the member element IDs); DANTE_LTR.gff3 is left untouched,
    # so the LTR library / masking / report keep seeing every individual element.
    # make_unified_annotation.R reads both: it splits DANTE_LTR.gff3's elements
    # into members (IDs listed here) vs standalone, and nests members under the
    # container.
    n_members = 0
    containers = []
    for n, ch in enumerate(chains, 1):
        lin = ch[0]["attrs"].get("Final_Classification", "")
        s, e = min(m["start"] for m in ch), max(m["end"] for m in ch)
        mids = [m["attrs"]["ID"] for m in ch]
        n_members += len(mids)
        col9 = (f"ID=LTR_RT_TR_{n:05d};Name={lin};Final_Classification={lin};"
                f"structure=LTR_RT_TR;copy_number={len(ch)};members={','.join(mids)}")
        containers.append((ch[0]["seqid"], s, e, ch[0]["strand"], ch[0]["source"], col9))
    containers.sort(key=lambda c: (c[0], c[1]))

    with open(args.output, "w") as fh:
        fh.write("##gff-version 3\n")
        for seqid, s, e, strand, source, col9 in containers:
            fh.write(render(seqid, source, "transposable_element",
                            s, e, ".", strand, ".", col9) + "\n")

    print(f"resolve_ltr_tandems: {len(chains)} tandem LTR-RT (LTR_RT_TR) "
          f"array(s) covering {n_members} member copies; wrote containers-only "
          f"{args.output} (DANTE_LTR.gff3 untouched).", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
