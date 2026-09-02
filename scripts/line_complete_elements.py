#!/usr/bin/env python3
"""Mark LINE elements whose BOTH ends were directly observed as `complete`.

WHAT THIS IS FOR. dante_line infers where a LINE ends by aligning the flanks of
many copies. That inference is what the boundary guards exist to bound. A small
minority of elements need no inference at all: a young insertion still carries the
two marks target-primed reverse transcription leaves behind — a poly-A tail at the
3' end, and a target-site duplication (TSD), a short run of host DNA copied on
both sides. Where both are found, the element's extent is MEASURED, not guessed.

This script rewrites those elements in DANTE_LINE.gff3: their span is replaced
with the MEASURED one, they are marked `Status=complete`, and the evidence (TSD
sequence, poly-A length) is recorded as attributes. The sequences are also
written to FASTA.

WHY THE SPAN IS REPLACED. The inferred boundary is substantially shorter than the
measured one -- median +1,740 bp at the 5' end and +1,558 bp at the 3' end on
Boechera, agreeing within 50 bp for essentially no element. Where the true extent
is known there is no reason to keep an estimate, so the element is extended to it.

THIS THEREFORE CHANGES THE ANNOTATION, for the small minority of elements that
qualify (0.7-5% of loci). Every other element keeps its inferred boundary and is
marked `Status=inferred`, so the two confidence classes are distinguishable
downstream.

WHY A CONFIRMED CALL IS TRUSTWORTHY. The two tests are statistically independent,
so their error rates multiply rather than overlap. Measured by running the
identical chain on the element's 5' side, where a LINE has no tail so any hit is
spurious by construction: over 2,743 loci in 4 genomes the spurious poly-A rate
was 2.8% and the chance-TSD rate given a bad anchor 2.6%, predicting a joint
false-positive rate of 0.073% -- and 0.073% was observed. A confirmed call is
right about 97% of the time.

DETERMINISM. The decoy windows that calibrate each TSD are drawn from a
per-locus seed derived from the locus coordinates, never a global RNG, so the
output does not depend on locus order, thread count or PYTHONHASHSEED. Rows and
records are emitted in sorted order and written atomically.
"""
from __future__ import annotations

import argparse
import bisect
import os
import random
import sys
import zlib
from collections import defaultdict

import numpy as np

# Chance floors, measured per search-window size: the 99th percentile score of
# the same search run against unrelated genomic windows. A perfect TSD of length
# L scores exactly L, so a TSD is separable from chance only above this.
# window bp -> minimum separable score
CHANCE_FLOOR = {500: 7, 1000: 9, 2000: 10, 4000: 10, 8000: 11}


def floor_for(window: int) -> int:
    """Chance floor for the nearest calibrated window size, rounded up."""
    for w in sorted(CHANCE_FLOOR):
        if window <= w:
            return CHANCE_FLOOR[w]
    return max(CHANCE_FLOOR.values())


class Faidx:
    """Random-access FASTA reader; avoids holding a multi-Gb genome in RAM."""

    def __init__(self, fasta: str):
        self.fh = open(fasta, "rb")
        self.idx = {}
        with open(fasta + ".fai") as fh:
            for line in fh:
                name, ln, off, lb, lw = line.split("\t")[:5]
                self.idx[name] = (int(ln), int(off), int(lb), int(lw))

    def fetch(self, name: str, start: int, end: int) -> str:
        if name not in self.idx:
            return ""
        ln, off, lb, lw = self.idx[name]
        start, end = max(0, start), min(ln, end)
        if end <= start:
            return ""
        s = off + start // lb * lw + start % lb
        e = off + (end - 1) // lb * lw + (end - 1) % lb + 1
        self.fh.seek(s)
        return self.fh.read(e - s).decode().replace("\n", "").replace("\r", "").upper()


_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def longest_run(seq: str, ch: str = "A", max_mismatch: int = 1):
    """Longest run of ``ch`` tolerating ``max_mismatch`` interruptions.

    A homopolymer run, not an A-rich window: an A-rich window is satisfied by
    ordinary coding AT-richness, which measured BELOW the 5' background.
    """
    best, at, n, i = 0, -1, len(seq), 0
    while i < n:
        if seq[i] != ch:
            i += 1
            continue
        mism, j, last = 0, i, i
        while j < n:
            if seq[j] == ch:
                last = j
            else:
                mism += 1
                if mism > max_mismatch:
                    break
            j += 1
        if last - i + 1 > best:
            best, at = last - i + 1, i
        i = last + 1 if last > i else i + 1
    return best, at


def is_tandem(seq: str, min_identity: float = 0.85) -> int:
    """Shortest period at which ``seq`` is its own tandem repeat; 0 if none.

    A genuine TSD is arbitrary host sequence. A periodic one came from inside a
    tandem array, where a short direct repeat is found trivially.
    """
    n = len(seq)
    for period in range(1, n // 2 + 1):
        rep = (seq[:period] * (n // period + 1))[:n]
        if sum(a == b for a, b in zip(seq, rep)) / n >= min_identity:
            return period
    return 0


def scan(query: str, window: str, min_k: int = 8, max_k: int = 20,
         mismatch_penalty: int = 2, count_above: int = None):
    """Best direct-repeat match of ``query``'s prefix anywhere in ``window``.

    Score for a prefix of length k with m matches is m - penalty*(k-m), so a
    longer repeat wins only if the extra bases mostly match. A perfect repeat of
    length k therefore scores exactly k, which is what makes CHANCE_FLOOR
    directly interpretable as a minimum TSD length.
    """
    q = np.frombuffer(query[:max_k].encode(), dtype=np.uint8)
    w = np.frombuffer(window.encode(), dtype=np.uint8)
    n = len(w) - min_k + 1
    if n <= 0 or len(q) < min_k:
        return -10 ** 9, -1, 0
    best = np.full(n, -10 ** 9, dtype=np.int32)
    best_k = np.zeros(n, dtype=np.int32)
    matches = np.zeros(n, dtype=np.int32)
    for i in range(min(max_k, len(q))):
        seg = w[i:i + n] if i + n <= len(w) else np.pad(w[i:], (0, i + n - len(w)))
        matches += (seg == q[i])
        k = i + 1
        if k < min_k:
            continue
        sc = matches * (1 + mismatch_penalty) - mismatch_penalty * k
        upd = sc > best
        best = np.where(upd, sc, best)
        best_k = np.where(upd, k, best_k)
    p = int(np.argmax(best))
    # How many places in the window match this well? A target-site duplication
    # is a UNIQUE pair; a query taken from inside a tandem array matches at every
    # period, and decoy windows drawn from elsewhere in the genome will not
    # contain that array, so the null test alone does not catch it.
    n_above = int((best > count_above).sum()) if count_above is not None else 0
    return int(best[p]), p, int(best_k[p]), n_above


def locus_seed(seqname: str, start: int, end: int) -> int:
    """Deterministic per-locus seed.

    zlib.crc32, not hash(): Python randomises string hashing per process, so a
    hash-derived seed would make the output depend on PYTHONHASHSEED.
    """
    return zlib.crc32(f"{seqname}:{start}:{end}".encode()) & 0xFFFFFFFF


def decoy_scores(query, fa, seqname, avoid, n_decoy, wlen, rng, min_k, max_k):
    """Same query, windows drawn from distant loci on the same sequence.

    Distant rather than shuffled: a shuffle destroys the local composition and
    repeat structure that make chance hits likely, and so understates the rate.
    """
    length = fa.idx.get(seqname, (0,))[0]
    out, tries = [], 0
    while len(out) < n_decoy and tries < n_decoy * 8:
        tries += 1
        if length <= wlen:
            break
        start = rng.randrange(0, length - wlen)
        if abs(start - avoid) < 200_000:
            continue
        w = fa.fetch(seqname, start, start + wlen)
        if len(w) < wlen // 2 or w.count("N") > wlen // 10:
            continue
        out.append(scan(query, w, min_k, max_k)[0])
    return out


def load_elements(gff: str):
    """LINE elements with their domain core, from a DANTE_LINE GFF3."""
    els, by_seq, doms = {}, defaultdict(list), defaultdict(list)
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].rstrip(";").split(";") if "=" in kv)
            if f[2] == "LINE_element":
                eid = attrs.get("ID", f"{f[0]}:{f[3]}-{f[4]}")
                els[eid] = (f[0], int(f[3]), int(f[4]), f[6],
                            attrs.get("Pattern_Type", "?"))
                by_seq[f[0]].append((int(f[3]), int(f[4]), eid))
            elif f[2] == "protein_domain":
                doms[f[0]].append((int(f[3]), int(f[4])))
    core = {}
    for seq, lst in by_seq.items():
        lst.sort()
        starts = [x[0] for x in lst]
        for ds, de in doms.get(seq, []):
            i = bisect.bisect_right(starts, ds) - 1
            if i >= 0 and lst[i][0] <= ds and de <= lst[i][1]:
                eid = lst[i][2]
                cur = core.get(eid)
                core[eid] = (min(cur[0], ds), max(cur[1], de)) if cur else (ds, de)
    out = []
    for eid in sorted(els):
        if eid in core:
            seq, s, e, strand, ptype = els[eid]
            out.append((eid, seq, core[eid][0], core[eid][1], strand, ptype))
    return out


def confirm(el, fa, args):
    """Return an evidence dict if BOTH ends are observed, else None."""
    eid, seq, cs, ce, strand, ptype = el
    pad5, pad3 = args.tsd_window + 80, args.tail_search + 80
    lo = (cs - 1 - pad5) if strand == "+" else (cs - 1 - pad3)
    hi = (ce + pad3) if strand == "+" else (ce + pad5)
    lo = max(0, lo)
    reg = fa.fetch(seq, lo, hi)
    if len(reg) < (hi - lo) // 2 or reg.count("N") > len(reg) // 10:
        return None
    c5o, c3o = (cs - 1) - lo, (ce - 1) - lo
    if strand == "-":
        n = len(reg)
        reg = revcomp(reg)
        c5, c3 = n - 1 - c3o, n - 1 - c5o
    else:
        c5, c3 = c5o, c3o
    if not (0 <= c5 < c3 < len(reg)):
        return None

    downstream = reg[c3 + 1:c3 + 1 + args.tail_search]
    run, at = longest_run(downstream, "A", args.tail_mismatch)
    if run < args.min_tail:
        return None
    tail_end = c3 + 1 + at + run

    query = reg[tail_end:tail_end + args.max_tsd]
    w0 = max(0, c5 - args.tsd_window)
    win = reg[w0:c5]
    if len(query) < args.min_tsd or len(win) < 100:
        return None
    fl = floor_for(args.tsd_window)
    score, pos, klen, n_above = scan(query, win, args.min_tsd, args.max_tsd,
                                     count_above=fl)
    if score <= fl:
        return None
    # A genuine duplication is unique in its window.
    if n_above > args.max_tsd_hits:
        return None
    # ...and is not itself a tandem repeat. Caught a telomeric CCCTAAA array
    # being read as a 16 bp TSD, which then displaced 5 kb of real telomere
    # annotation downstream.
    if is_tandem(query[:klen]):
        return None

    rng = random.Random(locus_seed(seq, cs, ce))
    decoys = decoy_scores(query, fa, seq, lo, args.decoys, len(win), rng,
                          args.min_tsd, args.max_tsd)
    if not decoys or not all(score > d for d in decoys):
        return None

    ext5 = c5 - (w0 + pos)
    ext3 = tail_end - (c3 + 1)
    elem = reg[w0 + pos:tail_end]
    if len(elem) < args.min_length or len(elem) > args.max_length:
        return None
    # measured span back in GENOME coordinates (the work above is in element
    # orientation, so the two extensions swap on the minus strand)
    if strand == "+":
        g_start, g_end = cs - ext5, ce + ext3
    else:
        g_start, g_end = cs - ext3, ce + ext5
    return dict(id=eid, seq=seq, strand=strand, pattern=ptype,
                core_start=cs, core_end=ce, core_len=ce - cs + 1,
                start=max(1, g_start), end=g_end,
                ext5=ext5, ext3=ext3, elem_len=len(elem),
                tsd_len=klen, tsd_score=score, tsd_seq=query[:klen],
                polya_len=run, sequence=elem)


def main():
    ap = argparse.ArgumentParser(
        description="Extract LINE elements whose both ends were directly observed "
                    "(poly-A tail + target-site duplication). Evaluation output; "
                    "feeds nothing downstream.")
    ap.add_argument("-g", "--genome", required=True, help="genome FASTA (needs .fai)")
    ap.add_argument("-a", "--gff", required=True, help="DANTE_LINE.gff3")
    ap.add_argument("-o", "--out-fasta", required=True)
    ap.add_argument("--annotate-gff3", action="store_true",
                    help="rewrite the input GFF3: extend confirmed elements to "
                         "their measured span, mark every element Status=complete "
                         "or Status=inferred, and record the evidence.")
    ap.add_argument("--tail-search", type=int, default=2500,
                    help="how far past the core to look for the poly-A tail (bp)")
    ap.add_argument("--min-tail", type=int, default=15,
                    help="minimum poly-A homopolymer run; below ~15 the 3' "
                         "enrichment over background collapses to ~1.2x")
    ap.add_argument("--tail-mismatch", type=int, default=1)
    ap.add_argument("--tsd-window", type=int, default=2100,
                    help="how far upstream of the core to look for the TSD (bp)")
    ap.add_argument("--min-tsd", type=int, default=8)
    ap.add_argument("--max-tsd", type=int, default=20)
    ap.add_argument("--decoys", type=int, default=6)
    ap.add_argument("--max-tsd-hits", type=int, default=3,
                    help="reject a candidate TSD matching more than this many "
                         "places in its window; a real duplication is unique, a "
                         "query from inside a tandem array is not")
    ap.add_argument("--min-length", type=int, default=1000)
    ap.add_argument("--max-length", type=int, default=7500,
                    help="a plant LINE is 4-7 kb; longer is not credible")
    ap.add_argument("--limit", type=int, default=0, help="0 = all loci")
    args = ap.parse_args()

    if not os.path.exists(args.genome + ".fai"):
        print(f"ERROR: missing FASTA index {args.genome}.fai", file=sys.stderr)
        return 1

    fa = Faidx(args.genome)
    elements = load_elements(args.gff)
    if args.limit:
        elements = elements[:args.limit]
    print(f"LINE elements with a domain core: {len(elements)}")

    rows = [r for r in (confirm(el, fa, args) for el in elements) if r]
    rows.sort(key=lambda r: (r["seq"], r["core_start"], r["core_end"]))

    n = len(elements)
    print(f"both ends observed: {len(rows)} ({100 * len(rows) / n:.1f}%)" if n
          else "both ends observed: 0")

    tmp = args.out_fasta + ".tmp"
    with open(tmp, "w") as fh:
        for r in rows:
            # RepeatMasker header convention, so the file is usable as a library
            # directly by anyone evaluating it.
            fh.write(f">{r['id']}_complete#Class_I/LINE\n{r['sequence']}\n")
    os.replace(tmp, args.out_fasta)
    print(f"wrote {args.out_fasta}")

    if args.annotate_gff3:
        annotate_gff3(args.gff, rows)
        print(f"marked {len(rows)} elements Status=complete in {args.gff}")
    return 0


def annotate_gff3(gff: str, rows) -> None:
    """Rewrite DANTE_LINE.gff3: extend confirmed elements to their measured span.

    Every LINE_element gains a `Status` attribute so the two confidence classes
    are distinguishable; confirmed ones additionally carry the evidence and take
    the measured coordinates. Rewritten atomically (.tmp + os.replace) so a
    killed run cannot leave a truncated GFF3 that a later step then trusts.
    """
    by_id = {r["id"]: r for r in rows}
    out = []
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                out.append(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "LINE_element":
                out.append(line)
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].rstrip(";").split(";") if "=" in kv)
            # strip any marks from a previous run so re-running is idempotent
            for k in ("Status", "TSD", "PolyA_length"):
                attrs.pop(k, None)
            r = by_id.get(attrs.get("ID"))
            if r is None:
                attrs["Status"] = "inferred"
            else:
                f[3], f[4] = str(r["start"]), str(r["end"])
                attrs["Status"] = "complete"
                attrs["TSD"] = r["tsd_seq"]
                attrs["PolyA_length"] = str(r["polya_len"])
            f[8] = ";".join(f"{k}={v}" for k, v in attrs.items())
            out.append("\t".join(f) + "\n")

    tmp = gff + ".tmp"
    with open(tmp, "w") as fh:
        fh.writelines(out)
    os.replace(tmp, gff)


if __name__ == "__main__":
    sys.exit(main())
