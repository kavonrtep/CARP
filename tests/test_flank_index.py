#!/usr/bin/env python3
"""Equivalence test for the FeatureIndex flank-clipping optimization.

dante_line.py / dante_tir_fallback.py used to rescan the whole feature list per
pattern/anchor to clip flanking regions -> O(patterns x features). On a large
genome mask_features (raw TideHunter) can hold millions of records, so that
quadratic scan dominates wall-time. FeatureIndex answers the same clip queries
by binary search.

This test proves the indexed result is IDENTICAL to the original brute-force
loops across many random inputs: both strands, boundary clamping, with/without
mask, features on multiple seqnames, and features straddling the region.

Run: python3 tests/test_flank_index.py
"""
import os
import random
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from dante_line import (  # noqa: E402
    FeatureGroup,
    FeatureIndex,
    GFF3Feature,
    create_prime_bed_files,
    get_flanking_regions,
)
import dante_tir_fallback as dtf  # noqa: E402


def mkfeat(seqname, start, end, strand="+"):
    return GFF3Feature(seqname=seqname, source="t", feature="domain", start=start,
                       end=end, score=".", strand=strand, phase=".", attributes="ID=x")


# --- reference brute-force implementations (verbatim original logic) ---------
def ref_get_flanking(pattern, all_features, flank_size, mask_features=None):
    region_start, region_end = pattern.get_region_bounds()
    flank_start = region_start - flank_size
    flank_end = region_end + flank_size
    seq_features = [f for f in all_features if f.seqname == pattern.seqname]
    for feature in seq_features:
        if (feature.end < region_start and feature.end > flank_start and
                feature not in pattern.features):
            flank_start = feature.end + 1
    for feature in seq_features:
        if (feature.start > region_end and feature.start < flank_end and
                feature not in pattern.features):
            flank_end = feature.start - 1
    if mask_features:
        mask_seq = [f for f in mask_features if f.seqname == pattern.seqname]
        for feature in mask_seq:
            if feature.end < region_start and feature.end > flank_start:
                flank_start = feature.end + 1
        for feature in mask_seq:
            if feature.start > region_end and feature.start < flank_end:
                flank_end = feature.start - 1
    return max(1, flank_start), flank_end


def ref_prime_line(patterns, all_features, flank_size, mask_features, seq_lengths):
    """Verbatim original create_prime_bed_files (dante_line) -> list of BED lines."""
    out5, out3 = [], []
    for pattern in patterns:
        seq_length = seq_lengths.get(pattern.seqname) if seq_lengths else None
        if pattern.strand == "+":
            endo, last = pattern.features[0], pattern.features[-1]
            p5s = max(1, endo.start - flank_size); p5e = endo.start - 1
            p3s = last.end + 1; p3e = last.end + flank_size
            sf = [f for f in all_features if f.seqname == pattern.seqname]
            for f in sf:
                if f.end < endo.start and f.end > p5s and f not in pattern.features:
                    p5s = f.end + 1
            for f in sf:
                if f.start > last.end and f.start < p3e and f not in pattern.features:
                    p3e = f.start - 1
            if mask_features:
                ms = [f for f in mask_features if f.seqname == pattern.seqname]
                for f in ms:
                    if f.end < endo.start and f.end > p5s:
                        p5s = f.end + 1
                for f in ms:
                    if f.start > last.end and f.start < p3e:
                        p3e = f.start - 1
        else:
            first, endo = pattern.features[0], pattern.features[-1]
            p5s = endo.end + 1; p5e = endo.end + flank_size
            p3s = max(1, first.start - flank_size); p3e = first.start - 1
            sf = [f for f in all_features if f.seqname == pattern.seqname]
            for f in sf:
                if f.start > endo.end and f.start < p5e and f not in pattern.features:
                    p5e = f.start - 1
            for f in sf:
                if f.end < first.start and f.end > p3s and f not in pattern.features:
                    p3s = f.end + 1
            if mask_features:
                ms = [f for f in mask_features if f.seqname == pattern.seqname]
                for f in ms:
                    if f.start > endo.end and f.start < p5e:
                        p5e = f.start - 1
                for f in ms:
                    if f.end < first.start and f.end > p3s:
                        p3s = f.end + 1
        if seq_length:
            p5e = min(p5e, seq_length); p3e = min(p3e, seq_length)
        if p5e > p5s:
            out5.append(f"{pattern.seqname}\t{p5s}\t{p5e}\t{pattern.group_id}_5prime")
        if p3e > p3s:
            out3.append(f"{pattern.seqname}\t{p3s}\t{p3e}\t{pattern.group_id}_3prime")
    return out5, out3


def ref_prime_fallback(anchors, all_features, flank_size, mask_features, seq_lengths):
    """Verbatim original create_prime_bed_files (dante_tir_fallback) -> BED lines."""
    out5, out3 = [], []
    for anchor in anchors:
        feature = anchor.feature
        seq_length = seq_lengths.get(anchor.seqname) if seq_lengths else None
        if anchor.strand == "+":
            p5s = max(1, feature.start - flank_size); p5e = feature.start - 1
            p3s = feature.end + 1; p3e = feature.end + flank_size
            sf = [f for f in all_features if f.seqname == anchor.seqname]
            for o in sf:
                if o is feature:
                    continue
                if o.end < feature.start and o.end > p5s:
                    p5s = o.end + 1
            for o in sf:
                if o is feature:
                    continue
                if o.start > feature.end and o.start < p3e:
                    p3e = o.start - 1
            if mask_features:
                ms = [f for f in mask_features if f.seqname == anchor.seqname]
                for o in ms:
                    if o.end < feature.start and o.end > p5s:
                        p5s = o.end + 1
                for o in ms:
                    if o.start > feature.end and o.start < p3e:
                        p3e = o.start - 1
        else:
            p5s = feature.end + 1; p5e = feature.end + flank_size
            p3s = max(1, feature.start - flank_size); p3e = feature.start - 1
            sf = [f for f in all_features if f.seqname == anchor.seqname]
            for o in sf:
                if o is feature:
                    continue
                if o.start > feature.end and o.start < p5e:
                    p5e = o.start - 1
            for o in sf:
                if o is feature:
                    continue
                if o.end < feature.start and o.end > p3s:
                    p3s = o.end + 1
            if mask_features:
                ms = [f for f in mask_features if f.seqname == anchor.seqname]
                for o in ms:
                    if o.start > feature.end and o.start < p5e:
                        p5e = o.start - 1
                for o in ms:
                    if o.end < feature.start and o.end > p3s:
                        p3s = o.end + 1
        if seq_length:
            p5e = min(p5e, seq_length); p3e = min(p3e, seq_length)
        if p5e > p5s:
            out5.append(f"{anchor.seqname}\t{p5s}\t{p5e}\t{anchor.group_id}_5prime")
        if p3e > p3s:
            out3.append(f"{anchor.seqname}\t{p3s}\t{p3e}\t{anchor.group_id}_3prime")
    return out5, out3


def rand_features(rng, seqnames, n):
    feats = []
    for _ in range(n):
        sn = rng.choice(seqnames)
        s = rng.randint(1, 20000)
        e = s + rng.randint(0, 800)   # varied lengths incl. long straddlers
        feats.append(mkfeat(sn, s, e))
    return feats


def test_primitive(rng):
    """FeatureIndex primitive vs brute force, standalone."""
    seqnames = ["chr1", "chr2", "c3"]
    for _ in range(300):
        feats = rand_features(rng, seqnames, rng.randint(0, 40))
        idx = FeatureIndex(feats)
        sn = rng.choice(seqnames)
        ref = rng.randint(1, 20000)
        floor = ref - rng.randint(0, 5000)
        ceil = ref + rng.randint(0, 5000)
        # brute force nearest_end_below
        fs = floor
        for f in feats:
            if f.seqname == sn and f.end < ref and f.end > fs:
                fs = f.end + 1
        assert idx.nearest_end_below(sn, ref, floor) == fs
        # brute force nearest_start_above
        ce = ceil
        for f in feats:
            if f.seqname == sn and f.start > ref and f.start < ce:
                ce = f.start - 1
        assert idx.nearest_start_above(sn, ref, ceil) == ce
    print("  primitive: FeatureIndex == brute force (300 trials)")


def test_get_flanking(rng):
    seqnames = ["chr1", "chr2", "c3"]
    for _ in range(200):
        all_f = rand_features(rng, seqnames, rng.randint(0, 60))
        mask_f = rand_features(rng, seqnames, rng.randint(0, 60)) if rng.random() < 0.7 else None
        sn = rng.choice(seqnames)
        # build a 2- or 3-feature pattern on sn
        s = rng.randint(2000, 15000)
        f1 = mkfeat(sn, s, s + rng.randint(50, 300))
        f2 = mkfeat(sn, f1.end + rng.randint(1, 400), f1.end + rng.randint(401, 800))
        pat = FeatureGroup("g1", sn, rng.choice("+-"), [f1, f2], "ENDO-RT")
        flank = rng.choice([2000, 5000, 10000])
        all_idx = FeatureIndex(all_f)
        mask_idx = FeatureIndex(mask_f) if mask_f else None
        got = get_flanking_regions(pat, all_idx, flank, mask_idx)
        exp = ref_get_flanking(pat, all_f, flank, mask_f)
        assert got == exp, f"get_flanking_regions {got} != {exp}"
    print("  get_flanking_regions: indexed == brute force (200 trials, +/- mask)")


def _make_pattern(rng, sn, strand):
    s = rng.randint(2000, 15000)
    f1 = mkfeat(sn, s, s + rng.randint(50, 300), strand)
    f2 = mkfeat(sn, f1.end + rng.randint(1, 400), f1.end + rng.randint(401, 800), strand)
    return FeatureGroup(f"g{rng.randint(0,9999)}", sn, strand, [f1, f2], "ENDO-RT")


def test_prime_line(rng):
    seqnames = ["chr1", "chr2", "c3"]
    for _ in range(80):
        all_f = rand_features(rng, seqnames, rng.randint(0, 60))
        mask_f = rand_features(rng, seqnames, rng.randint(0, 60)) if rng.random() < 0.7 else None
        seq_lengths = {sn: 30000 for sn in seqnames} if rng.random() < 0.5 else None
        patterns = [_make_pattern(rng, rng.choice(seqnames), rng.choice("+-"))
                    for _ in range(rng.randint(1, 8))]
        flank = rng.choice([2000, 5000, 10000])
        with tempfile.TemporaryDirectory() as td:
            b5 = os.path.join(td, "5.bed"); b3 = os.path.join(td, "3.bed")
            create_prime_bed_files(patterns, all_f, flank, b5, b3, mask_f, seq_lengths)
            got5 = open(b5).read().splitlines(); got3 = open(b3).read().splitlines()
        exp5, exp3 = ref_prime_line(patterns, all_f, flank, mask_f, seq_lengths)
        assert got5 == exp5, "prime line 5' BED mismatch"
        assert got3 == exp3, "prime line 3' BED mismatch"
    print("  create_prime_bed_files (dante_line): indexed == brute force (80 trials)")


def test_prime_fallback(rng):
    seqnames = ["chr1", "chr2", "c3"]
    for _ in range(80):
        all_f = rand_features(rng, seqnames, rng.randint(0, 60))
        mask_f = rand_features(rng, seqnames, rng.randint(0, 60)) if rng.random() < 0.7 else None
        seq_lengths = {sn: 30000 for sn in seqnames} if rng.random() < 0.5 else None
        anchors = []
        for _ in range(rng.randint(1, 8)):
            sn = rng.choice(seqnames); strand = rng.choice("+-")
            s = rng.randint(2000, 15000)
            feat = mkfeat(sn, s, s + rng.randint(50, 400), strand)
            anchors.append(dtf.TIRAnchor(group_id=f"a{rng.randint(0,9999)}", seqname=sn,
                                         strand=strand, subtype="hAT", feature=feat))
        flank = rng.choice([2000, 5000, 10000])
        with tempfile.TemporaryDirectory() as td:
            b5 = os.path.join(td, "5.bed"); b3 = os.path.join(td, "3.bed")
            dtf.create_prime_bed_files(anchors, all_f, flank, b5, b3, mask_f, seq_lengths)
            got5 = open(b5).read().splitlines(); got3 = open(b3).read().splitlines()
        exp5, exp3 = ref_prime_fallback(anchors, all_f, flank, mask_f, seq_lengths)
        assert got5 == exp5, "prime fallback 5' BED mismatch"
        assert got3 == exp3, "prime fallback 3' BED mismatch"
    print("  create_prime_bed_files (dante_tir_fallback): indexed == brute force (80 trials)")


if __name__ == "__main__":
    rng = random.Random(1234)
    test_primitive(rng)
    test_get_flanking(rng)
    test_prime_line(rng)
    test_prime_fallback(rng)
    print("test_flank_index: PASSED")
