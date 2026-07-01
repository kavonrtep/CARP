#!/usr/bin/env python3
"""Superfamily-aware array-length filter for TideCluster tc_reannotate.

`tc_reannotate.py` masks the genome with the per-TRC dimer library and then keeps
only regions where a *single TRC* forms an array >= its monomer length. That
filter is **strictly per-TRC** (TideCluster's `merge_overlapping_gff3_intervals`
merges only equal `Cluster_ID`; `filter_intervals` requires equal `Name`;
`split_intervals` drops multi-TRC overlaps as ambiguous). When one real tandem
array is tiled by two or more TRCs of the **same superfamily** (common for short
satellites with near-identical monomers), each per-TRC piece is sub-threshold and
the whole array is discarded — so a TE-derived tandem array is missed and the
underlying TE (e.g. a Tekay LTR-RT) prevails in the unified annotation. Enabling
`-culling_limit` masks the problem by collapsing each locus to one TRC, which
makes the annotation depend on culling.

This filter fixes it deterministically: it qualifies arrays at the **superfamily**
level (so the length test sees the true array across sibling TRCs, independent of
culling), then emits the recovered array geometry labelled with the **dominant
TRC** per segment — preserving the invariant that every satellite feature's
`Name` is a bare `TRC_<n>`. Input is the raw per-TRC hit GFF3 that
`tc_reannotate.py --debug` leaves in its temp dir (`rm.gff3`); output replaces the
per-TRC-filtered `RM_on_TideCluster_Library.gff3` and is byte-format-compatible
with it (`<seqid>\tRepeatMasker\trepeat\t<start>\t<end>\t.\t<strand>\t.\tName=TRC_<n>`).

TRCs absent from the superfamily table are treated as their own singleton group,
so the result is never *less* than the per-TRC filter would keep.
"""
import argparse
import collections


def parse_raw_gff3(path):
    """Yield (seqid, start0, end, strand, trc) from a tc_reannotate rm.gff3."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            name = None
            for kv in f[8].split(";"):
                if kv.startswith("Name="):
                    name = kv[5:].strip()
                    break
            if name is None:
                continue
            yield f[0], int(f[3]) - 1, int(f[4]), f[6], name


def load_superfamily_map(path):
    """Parse TideCluster_trc_superfamilies.csv -> {TRC: 'SF<n>'} (tolerant of
    quoting / column order)."""
    sfmap = {}
    with open(path) as fh:
        for line in fh:
            parts = [p.strip().strip('"') for p in line.replace(",", "\t").split("\t")]
            trc = next((p for p in parts if p.startswith("TRC_")), None)
            sfid = next((p for p in parts if p.isdigit()), None)
            if trc and sfid is not None:
                sfmap[trc] = "SF" + sfid
    return sfmap


def load_monomer_lengths(dimer_fasta):
    """{TRC: monomer_len}. The reduced dimer library holds doubled monomers, so
    monomer = max dimer length for that TRC / 2 (matches tc_reannotate's own
    per-TRC `seq_lengths`)."""
    dimer = {}
    name = None
    length = 0
    def flush(n, l):
        if n is None:
            return
        key = n.split("#")[0]
        dimer[key] = max(dimer.get(key, 0), l)
    with open(dimer_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                flush(name, length)
                name = line[1:].strip()
                length = 0
            else:
                length += len(line.strip())
    flush(name, length)
    return {trc: d / 2.0 for trc, d in dimer.items()}


def merge_plain(intervals):
    """Union-merge sorted (start, end) intervals."""
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def merge_named(intervals):
    """Union-merge (start, end, name); track per-name covered bp for dominant
    labelling. Returns [(start, end, Counter{name: bp})]."""
    out = []
    for s, e, nm in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
            out[-1][2][nm] += e - s
        else:
            out.append([s, e, collections.Counter({nm: e - s})])
    return [(s, e, cov) for s, e, cov in out]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rm-gff3", required=True,
                    help="raw per-TRC hits from tc_reannotate --debug (temp_dir/rm.gff3)")
    ap.add_argument("--superfamily-map", required=True,
                    help="TideCluster_trc_superfamilies.csv")
    ap.add_argument("--dimer-lib", required=True,
                    help="reduced dimer library FASTA (for monomer lengths)")
    ap.add_argument("--output", required=True, help="output GFF3")
    ap.add_argument("--min-units", type=float, default=2.0,
                    help="array must span >= this many monomers of the superfamily "
                         "(default 2, matching tc_reannotate's dimer-length cutoff)")
    args = ap.parse_args()

    sfmap = load_superfamily_map(args.superfamily_map)
    monomer = load_monomer_lengths(args.dimer_lib)

    def sf_of(trc):
        return sfmap.get(trc, trc)  # unmapped TRC = its own singleton group

    # group raw hits by (seqid, strand, superfamily)
    groups = collections.defaultdict(list)
    for seqid, s, e, strand, trc in parse_raw_gff3(args.rm_gff3):
        groups[(seqid, strand, sf_of(trc))].append((s, e, trc))

    # per-superfamily monomer scale (smallest member -> most permissive threshold,
    # matching "2 monomers of the shortest member" as in the per-TRC filter)
    sf_monomer = collections.defaultdict(lambda: float("inf"))
    for trc, m in monomer.items():
        sf_monomer[sf_of(trc)] = min(sf_monomer[sf_of(trc)], m)

    rows = []
    for (seqid, strand, sf), hits in groups.items():
        mono = sf_monomer.get(sf)
        if mono is None or mono == float("inf"):
            mono = min((monomer.get(t, 30.0) for _, _, t in hits), default=30.0)
        thr = args.min_units * mono
        # qualify: union of all sibling-TRC hits, extend to bridge intra-array
        # gaps (<= a dimer, capped at 10% of length, as tc does), keep arrays
        # whose extended span reaches the min-monomer threshold.
        plain = merge_plain([(s, e) for s, e, _ in hits])
        extended = []
        for s, e in plain:
            pad = int(min(2.0 * mono, (e - s) * 0.1))
            extended.append((max(0, s - pad), e + pad))
        qualified = [(s, e) for s, e in merge_plain(extended) if (e - s) >= thr]
        if not qualified:
            continue
        # recover geometry: emit the union of real hits inside qualified arrays,
        # each segment labelled by the TRC covering the most of it.
        for s, e, cov in merge_named(hits):
            if any(a <= s < b or a < e <= b or (s <= a and e >= b) for a, b in qualified):
                dom = max(cov.items(), key=lambda kv: kv[1])[0]
                rows.append((seqid, s, e, strand, dom))

    rows.sort()
    with open(args.output, "w") as out:
        out.write("##gff-version 3\n")
        for seqid, s, e, strand, trc in rows:
            out.write(f"{seqid}\tRepeatMasker\trepeat\t{s + 1}\t{e}\t.\t{strand}\t.\tName={trc}\n")


if __name__ == "__main__":
    main()
